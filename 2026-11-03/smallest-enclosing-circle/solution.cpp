// 最小圓覆蓋（Smallest Enclosing Circle / Welzl Move-to-Front）
//
// 題意：
//     給定平面上 n 個整數點（|x|, |y| <= 10^4），求覆蓋所有點的最小圓，
//     輸出圓心與半徑（四捨五入到小數點後 2 位，「遠離零」規則；
//     -0.00 正規化為 0.00）。
//     輸入：第一行 n（1 <= n <= 5000），接下來 n 行每行 x y。
//     輸出：一行「cx cy r」，單一空格分隔，各 2 位小數。
//
// 思路：
//     最小覆蓋圓唯一，由至多 3 個點決定（1 點：半徑 0；2 點：直徑圓；
//     3 點：銳角 / 直角三角形的外接圓）。
//     Welzl move-to-front：點隨機打亂後逐個加入；新點不在當前圓內時，
//     它必在最終圓的邊界上 → 帶單約束對前綴再做一遍增量（又遇外點則
//     帶雙約束再掃一遍，依叉積分側取能包住三點的最小外接圓）。
//     期望 O(n)。
//
// 數值細節：
//     * 洗牌用自製 LCG，Python 與 C++ 產生相同點序，輸出逐位元組一致。
//     * 外接圓計算先把座標平移到點集範圍中心，提升精度。
//     * 共線三點無外接圓（行列式為 0），回傳無效由上層跳過。
//
// 應用場景：
//     地理圍欄 / 配送半徑規劃、蜂窩基站最小覆蓋半徑、
//     製造業定位誤差（最小包絡圓）、電腦圖學包圍球預計算。

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <string>
#include <vector>
#include <algorithm>

#ifdef _WIN32
#include <io.h>
#include <fcntl.h>
#endif

typedef long long ll;

const int NMAX = 5000;
const ll XYMAX = 10000;
const double EPS = 1e-14;

struct Point {
    ll x, y;
};

struct Circle {
    double cx, cy, r;
    bool valid;
};

static bool is_in_circle(const Circle &c, double x, double y) {
    if (!c.valid) return false;
    double dx = x - c.cx;
    double dy = y - c.cy;
    return dx * dx + dy * dy <= c.r * c.r * (1.0 + EPS);
}

// 以 a、b 為直徑的圓
static Circle make_diameter(double ax, double ay, double bx, double by) {
    Circle c;
    c.cx = (ax + bx) / 2.0;
    c.cy = (ay + by) / 2.0;
    double dx = ax - bx;
    double dy = ay - by;
    c.r = std::sqrt(dx * dx + dy * dy) / 2.0;
    c.valid = true;
    return c;
}

// 三點外接圓；共線回傳 invalid
static Circle make_circumcircle(double ax, double ay, double bx, double by,
                                double cx, double cy) {
    // 平移到座標範圍中心，提升數值精度
    double ox = (std::min(ax, std::min(bx, cx)) + std::max(ax, std::max(bx, cx))) / 2.0;
    double oy = (std::min(ay, std::min(by, cy)) + std::max(ay, std::max(by, cy))) / 2.0;
    double ax2 = ax - ox, ay2 = ay - oy;
    double bx2 = bx - ox, by2 = by - oy;
    double cx2 = cx - ox, cy2 = cy - oy;
    double d = (ax2 * (by2 - cy2) + bx2 * (cy2 - ay2) + cx2 * (ay2 - by2)) * 2.0;
    Circle c;
    c.valid = false;
    if (d == 0.0) return c;
    double aa = ax2 * ax2 + ay2 * ay2;
    double bb = bx2 * bx2 + by2 * by2;
    double cc = cx2 * cx2 + cy2 * cy2;
    c.cx = ox + (aa * (by2 - cy2) + bb * (cy2 - ay2) + cc * (ay2 - by2)) / d;
    c.cy = oy + (aa * (cx2 - bx2) + bb * (ax2 - cx2) + cc * (bx2 - ax2)) / d;
    double r1 = std::sqrt((c.cx - ax) * (c.cx - ax) + (c.cy - ay) * (c.cy - ay));
    double r2 = std::sqrt((c.cx - bx) * (c.cx - bx) + (c.cy - by) * (c.cy - by));
    double r3 = std::sqrt((c.cx - cx) * (c.cx - cx) + (c.cy - cy) * (c.cy - cy));
    c.r = std::max(r1, std::max(r2, r3));
    c.valid = true;
    return c;
}

// 求覆蓋 points 且邊界經過 p、q 的最小圓
static Circle make_circle_two_points(const std::vector<Point> &points,
                                     double px, double py, double qx, double qy) {
    Circle circ = make_diameter(px, py, qx, qy);
    Circle left, right;
    left.valid = false;
    right.valid = false;
    double left_side = 0.0;   // left 圓心相對弦 p→q 的叉積（>0，越遠越大）
    double right_side = 0.0;  // right 圓心相對弦 p→q 的叉積（<0，越遠越小）
    for (size_t i = 0; i < points.size(); ++i) {
        double sx = (double)points[i].x;
        double sy = (double)points[i].y;
        if (is_in_circle(circ, sx, sy)) continue;
        double cross = (qx - px) * (sy - py) - (qy - py) * (sx - px);
        Circle c = make_circumcircle(px, py, qx, qy, sx, sy);
        if (!c.valid) continue;
        // 候選圓心相對弦 p→q 的叉積（帶符號距離 × |pq|）
        double c_side = (qx - px) * (c.cy - py) - (qy - py) * (c.cx - px);
        // 同側候選中取「圓心離弦最遠」者，才能包住同側所有點
        if (cross > 0.0 && (!left.valid || c_side > left_side)) {
            left = c;
            left_side = c_side;
        } else if (cross < 0.0 && (!right.valid || c_side < right_side)) {
            right = c;
            right_side = c_side;
        }
    }
    if (!left.valid && !right.valid) return circ;
    if (!left.valid) return right;
    if (!right.valid) return left;
    return (left.r <= right.r) ? left : right;
}

// 求覆蓋 points 且邊界經過 p 的最小圓
static Circle make_circle_one_point(const std::vector<Point> &points,
                                    double px, double py) {
    Circle c;
    c.cx = px;
    c.cy = py;
    c.r = 0.0;
    c.valid = true;
    for (size_t i = 0; i < points.size(); ++i) {
        double sx = (double)points[i].x;
        double sy = (double)points[i].y;
        if (!is_in_circle(c, sx, sy)) {
            if (c.r == 0.0) {
                c = make_diameter(px, py, sx, sy);
            } else {
                std::vector<Point> pref(points.begin(), points.begin() + i + 1);
                c = make_circle_two_points(pref, px, py, sx, sy);
            }
        }
    }
    return c;
}

// 自製 LCG 的 Fisher-Yates 洗牌（與 Python 版序列一致）
static std::vector<Point> lcg_shuffle(const std::vector<Point> &points) {
    std::vector<Point> a = points;
    unsigned int state = 20261103u;
    for (int i = (int)a.size() - 1; i > 0; --i) {
        state = (unsigned int)((ll)state * 1103515245LL + 12345LL) & 0x7FFFFFFFu;
        int j = (int)(state % (unsigned)(i + 1));
        std::swap(a[i], a[j]);
    }
    return a;
}

// Welzl move-to-front
static Circle smallest_enclosing_circle(const std::vector<Point> &points) {
    std::vector<Point> shuffled = lcg_shuffle(points);
    Circle c;
    c.valid = false;
    for (size_t i = 0; i < shuffled.size(); ++i) {
        double px = (double)shuffled[i].x;
        double py = (double)shuffled[i].y;
        if (!c.valid || !is_in_circle(c, px, py)) {
            std::vector<Point> pref(shuffled.begin(), shuffled.begin() + i + 1);
            c = make_circle_one_point(pref, px, py);
        }
    }
    return c;
}

// ---- 暴力對拍（O(n^4)，僅小規模）----
static bool covers_all(const std::vector<Point> &pts,
                       double cx, double cy, double r) {
    double tol = std::max(r, 1.0) * 1e-7;
    for (size_t i = 0; i < pts.size(); ++i) {
        double dx = (double)pts[i].x - cx;
        double dy = (double)pts[i].y - cy;
        if (dx * dx + dy * dy > r * r + tol) return false;
    }
    return true;
}

static double brute_force(const std::vector<Point> &pts) {
    int n = (int)pts.size();
    double best = 1e100;
    for (int i = 0; i < n; ++i) {
        for (int j = i; j < n; ++j) {
            double ax = (double)pts[i].x, ay = (double)pts[i].y;
            double bx = (double)pts[j].x, by = (double)pts[j].y;
            double cx = (ax + bx) / 2.0, cy = (ay + by) / 2.0;
            double dx = ax - bx, dy = ay - by;
            double r = std::sqrt(dx * dx + dy * dy) / 2.0;
            if (covers_all(pts, cx, cy, r)) best = std::min(best, r);
        }
    }
    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            for (int k = j + 1; k < n; ++k) {
                Circle c = make_circumcircle((double)pts[i].x, (double)pts[i].y,
                                             (double)pts[j].x, (double)pts[j].y,
                                             (double)pts[k].x, (double)pts[k].y);
                if (!c.valid) continue;
                if (covers_all(pts, c.cx, c.cy, c.r)) best = std::min(best, c.r);
            }
        }
    }
    return best;
}

// ---- 輸出格式化（與 Python 完全一致）----
static double round2(double x) {
    double v;
    if (x >= 0.0) v = (double)(long long)(x * 100.0 + 0.5) / 100.0;
    else          v = (double)(long long)(x * 100.0 - 0.5) / 100.0;
    if (v == 0.0) return 0.0;
    return v;
}

static std::string fmt(double x) {
    char buf[64];
    sprintf(buf, "%.2f", round2(x));
    return std::string(buf);
}

static void emit(const std::string &s) { fwrite(s.data(), 1, s.size(), stdout); }

static void run_tests() {
    // ---- 固定用例 ----
    {
        std::vector<Point> pts; pts.push_back({3, 4});
        Circle c = smallest_enclosing_circle(pts);
        if (c.cx != 3.0 || c.cy != 4.0 || c.r != 0.0) { emit("FAIL 1\n"); return; }
    }
    {
        std::vector<Point> pts; pts.push_back({0, 0}); pts.push_back({4, 0});
        Circle c = smallest_enclosing_circle(pts);
        if (c.cx != 2.0 || c.cy != 0.0 || c.r != 2.0) { emit("FAIL 2\n"); return; }
    }
    {
        // 直角三角形：圓心 (2, 1.5)、半徑 2.5
        std::vector<Point> pts; pts.push_back({0, 0}); pts.push_back({4, 0}); pts.push_back({0, 3});
        Circle c = smallest_enclosing_circle(pts);
        if (c.cx != 2.0 || c.cy != 1.5 || c.r != 2.5) { emit("FAIL 3\n"); return; }
    }
    {
        std::vector<Point> pts; pts.push_back({0, 0}); pts.push_back({2, 0});
        pts.push_back({2, 2}); pts.push_back({0, 2});
        Circle c = smallest_enclosing_circle(pts);
        if (std::abs(c.cx - 1.0) > 1e-9 || std::abs(c.cy - 1.0) > 1e-9 ||
            std::abs(c.r - std::sqrt(2.0)) > 1e-9) { emit("FAIL 4\n"); return; }
    }
    {
        std::vector<Point> pts; pts.push_back({-1, 0}); pts.push_back({0, 0}); pts.push_back({5, 0});
        Circle c = smallest_enclosing_circle(pts);
        if (c.cx != 2.0 || c.cy != 0.0 || c.r != 3.0) { emit("FAIL 5\n"); return; }
    }
    {
        std::vector<Point> pts; pts.push_back({2, 2}); pts.push_back({2, 2}); pts.push_back({2, 2});
        Circle c = smallest_enclosing_circle(pts);
        if (c.cx != 2.0 || c.cy != 2.0 || c.r != 0.0) { emit("FAIL 6\n"); return; }
    }

    // ---- 與暴力對拍 ----
    unsigned int lcg_state = 20261103u;
    auto lcg_next = [&]() -> unsigned int {
        lcg_state = (unsigned int)((ll)lcg_state * 1103515245LL + 12345LL) & 0x7FFFFFFFu;
        return lcg_state;
    };
    for (int it = 0; it < 250; ++it) {
        int n = 1 + (int)(lcg_next() % 14u);
        std::vector<Point> pts;
        for (int i = 0; i < n; ++i)
            pts.push_back({(ll)(lcg_next() % 41u) - 20, (ll)(lcg_next() % 41u) - 20});
        Circle c = smallest_enclosing_circle(pts);
        // 包容性
        for (size_t i = 0; i < pts.size(); ++i) {
            double dx = (double)pts[i].x - c.cx;
            double dy = (double)pts[i].y - c.cy;
            if (dx * dx + dy * dy > c.r * c.r * (1.0 + 1e-9)) {
                emit("FAIL random contain\n"); return;
            }
        }
        double want = brute_force(pts);
        if (std::abs(c.r - want) > std::max(want, 1.0) * 1e-6) {
            emit("FAIL random radius\n"); return;
        }
    }
    // 大座標範圍
    for (int it = 0; it < 80; ++it) {
        int n = 1 + (int)(lcg_next() % 10u);
        std::vector<Point> pts;
        for (int i = 0; i < n; ++i)
            pts.push_back({(ll)(lcg_next() % 20001u) - 10000, (ll)(lcg_next() % 20001u) - 10000});
        Circle c = smallest_enclosing_circle(pts);
        for (size_t i = 0; i < pts.size(); ++i) {
            double dx = (double)pts[i].x - c.cx;
            double dy = (double)pts[i].y - c.cy;
            if (dx * dx + dy * dy > c.r * c.r * (1.0 + 1e-9)) {
                emit("FAIL big contain\n"); return;
            }
        }
        double want = brute_force(pts);
        if (std::abs(c.r - want) > std::max(want, 1.0) * 1e-6) {
            emit("FAIL big radius\n"); return;
        }
    }
}

// ---- IO 模式 ----
static bool parse_ll(const char *s, ll &out) {
    if (!s || !*s) return false;
    int i = 0;
    if (s[0] == '+' || s[0] == '-') i = 1;
    if (!s[i]) return false;
    for (int j = i; s[j]; ++j)
        if (s[j] < '0' || s[j] > '9') return false;
    if (strlen(s) - i > 19) return false;
    long long v = 0;
    for (int j = i; s[j]; ++j) v = v * 10 + (s[j] - '0');
    if (v > 1000000000000000000LL) return false;
    out = (s[0] == '-') ? -v : v;
    return true;
}

int main() {
#ifdef _WIN32
    _setmode(_fileno(stdout), _O_BINARY);
#endif
    std::string raw;
    {
        char buf[1 << 16];
        size_t n;
        while ((n = fread(buf, 1, sizeof(buf), stdin)) > 0) raw.append(buf, n);
    }
    std::vector<std::string> toks;
    {
        std::string cur;
        for (size_t i = 0; i <= raw.size(); ++i) {
            char ch = (i < raw.size()) ? raw[i] : ' ';
            if (ch == ' ' || ch == '\n' || ch == '\r' || ch == '\t') {
                if (!cur.empty()) { toks.push_back(cur); cur.clear(); }
            } else {
                cur += ch;
            }
        }
    }
    if (toks.empty()) {
        run_tests();
        emit("all tests passed\n");
        return 0;
    }

    size_t pos = 0;
    ll n_ll;
    if (pos >= toks.size() || !parse_ll(toks[pos++].c_str(), n_ll) ||
        n_ll < 1 || n_ll > NMAX) { emit("error\n"); return 0; }
    int n = (int)n_ll;

    std::vector<Point> pts;
    for (int i = 0; i < n; ++i) {
        ll x, y;
        if (pos >= toks.size() || !parse_ll(toks[pos++].c_str(), x) ||
            pos >= toks.size() || !parse_ll(toks[pos++].c_str(), y)) {
            emit("error\n"); return 0;
        }
        if (x < -XYMAX || x > XYMAX || y < -XYMAX || y > XYMAX) {
            emit("error\n"); return 0;
        }
        Point p; p.x = x; p.y = y;
        pts.push_back(p);
    }

    Circle c = smallest_enclosing_circle(pts);
    emit(fmt(c.cx) + " " + fmt(c.cy) + " " + fmt(c.r) + "\n");
    return 0;
}
