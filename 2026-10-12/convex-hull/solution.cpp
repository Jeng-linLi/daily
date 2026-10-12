// 計算幾何：凸包（Andrew / Jarvis / Graham）+ 面積 + 旋轉卡殼直徑 + 點在凸包內
//
// 與 solution.py 完全同構：同樣的演算法、同樣的確定性規則、同樣的輸入輸出格式。
// 全程整數運算（叉積 / 面積兩倍 / 距離平方），不使用任何浮點，故兩邊輸出逐字節一致。
// 編譯：g++ -std=c++17 -O2 -Wall solution.cpp -o solution
// 無 stdin 輸入時運行內置斷言測試並輸出 `all tests passed`。

#include <algorithm>
#include <cassert>
#include <cctype>
#include <cstdlib>
#include <iostream>
#include <random>
#include <set>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

using namespace std;

typedef long long ll;
typedef pair<ll, ll> Point;

static const ll COORD_MAX = 1000000LL;   // 座標截斷上限，保證叉積不溢位
static const int N_MAX = 2000;           // 點數上限

// ---------------------------------------------------------------- 基本工具

// (a-o) × (b-o)：>0 表示 b 在 o→a 的左側（逆時針），=0 共線。
static ll cross(Point o, Point a, Point b) {
    return (a.first - o.first) * (b.second - o.second) - (a.second - o.second) * (b.first - o.first);
}

static ll absl(ll x) { return x < 0 ? -x : x; }

// 兩點距離的平方（整數，不開根號）。
static ll dist2(Point a, Point b) {
    ll dx = a.first - b.first, dy = a.second - b.second;
    return dx * dx + dy * dy;
}

// 去重後按 (x, y) 排序 —— 這是三個演算法共用的起手式。
static vector<Point> dedupSorted(const vector<Point> &points) {
    vector<Point> v = points;
    sort(v.begin(), v.end());
    v.erase(unique(v.begin(), v.end()), v.end());
    return v;
}

// ---------------------------------------------------------------- 三種凸包構造

// Andrew's Monotone Chain，O(n log n)，逆時針、共線點剔除。
static vector<Point> andrewHull(const vector<Point> &points) {
    vector<Point> pts = dedupSorted(points);
    if ((int)pts.size() <= 2) return pts;
    vector<Point> lower, upper;
    for (Point p : pts) {
        while (lower.size() >= 2 && cross(lower[lower.size() - 2], lower.back(), p) <= 0)
            lower.pop_back();
        lower.push_back(p);
    }
    for (int i = (int)pts.size() - 1; i >= 0; i--) {
        Point p = pts[i];
        while (upper.size() >= 2 && cross(upper[upper.size() - 2], upper.back(), p) <= 0)
            upper.pop_back();
        upper.push_back(p);
    }
    vector<Point> hull;
    for (int i = 0; i + 1 < (int)lower.size(); i++) hull.push_back(lower[i]);
    for (int i = 0; i + 1 < (int)upper.size(); i++) hull.push_back(upper[i]);
    return hull;
}

// Jarvis march（禮物包裝），O(n·h)，每步選「最順時針」的點，共線取較遠者。
static vector<Point> jarvisHull(const vector<Point> &points) {
    vector<Point> pts = dedupSorted(points);
    int n = (int)pts.size();
    if (n <= 2) return pts;
    vector<Point> hull;
    int cur = 0;                               // 字典序最小點必定是凸包頂點
    while (true) {
        hull.push_back(pts[cur]);
        int nxt = -1;
        for (int j = 0; j < n; j++) {
            if (j == cur) continue;
            if (nxt == -1) { nxt = j; continue; }
            ll c = cross(pts[cur], pts[nxt], pts[j]);
            if (c < 0) nxt = j;                // j 更順時針 → 換成 j
            else if (c == 0 && dist2(pts[cur], pts[j]) > dist2(pts[cur], pts[nxt]))
                nxt = j;                       // 共線 → 取較遠的，跳過中間點
        }
        if (nxt == 0) break;
        cur = nxt;
        if ((int)hull.size() > n) break;       // 保險絲，正常不會觸發
    }
    return hull;
}

// Graham scan，O(n log n)，極角排序只用叉積 + 距離，不使用 atan2。
static vector<Point> grahamHull(const vector<Point> &points) {
    vector<Point> pts = dedupSorted(points);
    if ((int)pts.size() <= 2) return pts;
    Point p0 = pts[0];
    vector<Point> rest(pts.begin() + 1, pts.end());
    sort(rest.begin(), rest.end(), [&](Point a, Point b) {
        ll c = cross(p0, a, b);
        if (c > 0) return true;                // a 的極角較小，排在前面
        if (c < 0) return false;
        ll da = dist2(p0, a), db = dist2(p0, b);
        return da < db;
    });
    vector<Point> stk;
    stk.push_back(p0);
    stk.push_back(rest[0]);
    for (int i = 1; i < (int)rest.size(); i++) {
        Point p = rest[i];
        while (stk.size() >= 2 && cross(stk[stk.size() - 2], stk.back(), p) <= 0) stk.pop_back();
        stk.push_back(p);
    }
    return stk;
}

// ---------------------------------------------------------------- 凸包上的計算

// 鞋帶公式求面積的**兩倍**（整數，避免浮點）。
static ll hullArea2(const vector<Point> &hull) {
    int h = (int)hull.size();
    if (h < 3) return 0;
    ll s = 0;
    for (int i = 0; i < h; i++) {
        Point a = hull[i], b = hull[(i + 1) % h];
        s += a.first * b.second - b.first * a.second;
    }
    return absl(s);
}

// 旋轉卡殼求最遠兩點距離的平方，O(h)。
static ll diameter2(const vector<Point> &hull) {
    int h = (int)hull.size();
    if (h < 2) return 0;
    if (h == 2) return dist2(hull[0], hull[1]);
    ll best = 0;
    int j = 1;
    for (int i = 0; i < h; i++) {
        int ni = (i + 1) % h;
        int steps = 0;
        while (steps < 2 * h) {                // 保險絲：j 全程最多前進 2h 次
            int nj = (j + 1) % h;
            if (absl(cross(hull[i], hull[ni], hull[nj])) > absl(cross(hull[i], hull[ni], hull[j]))) {
                j = nj;
                steps++;
            } else {
                break;
            }
        }
        best = max(best, dist2(hull[i], hull[j]));
        best = max(best, dist2(hull[ni], hull[j]));
    }
    return best;
}

// O(h²) 暴力求最遠點對，僅供對拍。
static ll diameter2Bruteforce(const vector<Point> &hull) {
    int h = (int)hull.size();
    ll best = 0;
    for (int i = 0; i < h; i++)
        for (int j = i + 1; j < h; j++) best = max(best, dist2(hull[i], hull[j]));
    return best;
}

// q 是否在凸包內（含邊界），O(log h)。
static bool pointInConvex(const vector<Point> &hull, Point q) {
    int h = (int)hull.size();
    if (h == 0) return false;
    if (h == 1) return q == hull[0];
    if (h == 2) {
        if (cross(hull[0], hull[1], q) != 0) return false;
        return min(hull[0].first, hull[1].first) <= q.first &&
               q.first <= max(hull[0].first, hull[1].first) &&
               min(hull[0].second, hull[1].second) <= q.second &&
               q.second <= max(hull[0].second, hull[1].second);
    }
    // 楔形外直接排除：q 必須在 p0→p1 的左側、且在 p0→p_{h-1} 的右側
    if (cross(hull[0], hull[1], q) < 0) return false;
    if (cross(hull[0], hull[h - 1], q) > 0) return false;
    // 二分：找最大的 i（1 ≤ i ≤ h-2）使「q 在 p0→p_i 的左側（含共線）」
    int lo = 1, hi = h - 1;
    while (hi - lo > 1) {
        int mid = (lo + hi) / 2;
        if (cross(hull[0], hull[mid], q) >= 0) lo = mid;
        else hi = mid;
    }
    return cross(hull[lo], hull[lo + 1], q) >= 0;
}

// O(h) 暴力判定：凸多邊形內部 ⟺ 在所有邊的左側（含共線）。
static bool pointInConvexBruteforce(const vector<Point> &hull, Point q) {
    int h = (int)hull.size();
    if (h == 0) return false;
    if (h == 1) return q == hull[0];
    if (h == 2) return pointInConvex(hull, q);   // 退化情形直接沿用同一套邏輯
    for (int i = 0; i < h; i++)
        if (cross(hull[i], hull[(i + 1) % h], q) < 0) return false;
    return true;
}

// ---------------------------------------------------------------- 輸出

static string run(const vector<Point> &points, const vector<Point> &queries) {
    vector<Point> hull = andrewHull(points);
    ostringstream oss;
    oss << hull.size() << "\n";
    for (int i = 0; i < (int)hull.size(); i++) {
        if (i) oss << " ";
        oss << hull[i].first << "," << hull[i].second;
    }
    oss << "\n" << hullArea2(hull) << "\n" << diameter2(hull) << "\n";
    for (int i = 0; i < (int)queries.size(); i++) {
        if (i) oss << " ";
        oss << (pointInConvex(hull, queries[i]) ? 1 : 0);
    }
    return oss.str();
}

// ---------------------------------------------------------------- 測試

static void checkHull(const vector<Point> &points) {
    vector<Point> a = andrewHull(points);
    vector<Point> j = jarvisHull(points);
    vector<Point> g = grahamHull(points);
    assert(a == j);
    assert(a == g);

    int h = (int)a.size();
    set<Point> uniq(points.begin(), points.end());
    for (Point p : a) assert(uniq.count(p));
    if (h >= 3)
        for (int i = 0; i < h; i++)
            assert(cross(a[i], a[(i + 1) % h], a[(i + 2) % h]) > 0);
    for (Point p : uniq) assert(pointInConvexBruteforce(a, p));
    assert(diameter2(a) == diameter2Bruteforce(a));
}

static void runTests() {
    mt19937 rng(20261012);
    auto ri = [&](int lo, int hi) { return lo + (int)(rng() % (hi - lo + 1)); };

    // ---- 固定用例 ----
    vector<Point> sq = {{0, 0}, {10, 0}, {10, 10}, {0, 10}};
    assert((andrewHull(sq) == vector<Point>{{0, 0}, {10, 0}, {10, 10}, {0, 10}}));
    assert(hullArea2(andrewHull(sq)) == 200);
    assert(diameter2(andrewHull(sq)) == 200);
    checkHull(sq);

    checkHull({{0, 0}, {10, 0}, {10, 10}, {0, 10}, {5, 5}, {5, 0}, {0, 5}});
    assert((andrewHull({{0, 0}, {10, 0}, {10, 10}, {0, 10}, {5, 5}, {5, 0}, {0, 5}}) ==
            vector<Point>{{0, 0}, {10, 0}, {10, 10}, {0, 10}}));

    assert(andrewHull({}).empty());
    assert((andrewHull({{3, 4}}) == vector<Point>{{3, 4}}));
    assert((andrewHull({{0, 0}, {5, 0}}) == vector<Point>{{0, 0}, {5, 0}}));
    assert((andrewHull({{0, 0}, {1, 0}, {2, 0}, {3, 0}}) == vector<Point>{{0, 0}, {3, 0}}));
    assert((andrewHull({{1, 1}, {1, 1}, {1, 1}}) == vector<Point>{{1, 1}}));
    assert(hullArea2({{0, 0}, {3, 0}}) == 0);
    assert(diameter2({{0, 0}, {3, 0}}) == 9);
    checkHull({});
    checkHull({{3, 4}});
    checkHull({{0, 0}, {5, 0}});
    checkHull({{0, 0}, {1, 0}, {2, 0}, {3, 0}});

    checkHull({{-5, -5}, {5, -5}, {5, 5}, {-5, 5}, {0, 0}});
    assert(hullArea2(andrewHull({{-5, -5}, {5, -5}, {5, 5}, {-5, 5}})) == 200);

    assert(run({{0, 0}, {10, 0}, {10, 10}, {0, 10}}, {{5, 5}, {20, 20}, {10, 10}}) ==
           "4\n0,0 10,0 10,10 0,10\n200\n200\n1 0 1");
    assert(run({}, {}) == "0\n\n0\n0\n");

    // ---- 隨機對拍：小座標，共線與重複點大量出現 ----
    for (int iter = 0; iter < 600; iter++) {
        int n = ri(0, 10);
        vector<Point> pts;
        for (int i = 0; i < n; i++) pts.push_back({ri(-4, 4), ri(-4, 4)});
        checkHull(pts);
    }

    // ---- 隨機對拍：中等座標 ----
    for (int iter = 0; iter < 400; iter++) {
        int n = ri(0, 25);
        vector<Point> pts;
        for (int i = 0; i < n; i++) pts.push_back({ri(-50, 50), ri(-50, 50)});
        checkHull(pts);
    }

    // ---- 隨機對拍：去重後的小座標集合 ----
    for (int iter = 0; iter < 200; iter++) {
        int n = ri(3, 12);
        set<Point> s;
        for (int i = 0; i < n; i++) s.insert({ri(-6, 6), ri(-6, 6)});
        vector<Point> pts(s.begin(), s.end());
        checkHull(pts);
    }

    // ---- 點在凸包內：與暴力逐點對拍 ----
    for (int iter = 0; iter < 400; iter++) {
        int n = ri(0, 12);
        vector<Point> pts;
        for (int i = 0; i < n; i++) pts.push_back({ri(-6, 6), ri(-6, 6)});
        vector<Point> hull = andrewHull(pts);
        for (int t = 0; t < 20; t++) {
            Point q = {ri(-9, 9), ri(-9, 9)};
            assert(pointInConvex(hull, q) == pointInConvexBruteforce(hull, q));
        }
    }
}

// stdin 解析：嚴格只吃 `^[+-]?[0-9]+$` 的 token，非數字一律當作輸入結束。
static bool isInt(const string &s) {
    if (s.empty()) return false;
    int i = 0;
    if (s[0] == '+' || s[0] == '-') i = 1;
    if (i >= (int)s.size()) return false;
    for (int j = i; j < (int)s.size(); j++)
        if (!isdigit((unsigned char)s[j])) return false;
    return true;
}

static ll clampCoord(ll v) { return max(-COORD_MAX, min(COORD_MAX, v)); }

static void runIO(const vector<string> &tokens) {
    auto readInt = [&](size_t i) -> pair<bool, ll> {
        if (i >= tokens.size() || !isInt(tokens[i])) return {false, 0};
        return {true, stoll(tokens[i])};
    };

    auto p0 = readInt(0);
    int n = p0.first ? (int)max<ll>(0, min<ll>(N_MAX, p0.second)) : 0;
    vector<Point> pts;
    for (int i = 0; i < n; i++) {
        auto x = readInt(1 + 2 * (size_t)i);
        auto y = readInt(2 + 2 * (size_t)i);
        if (!x.first || !y.first) break;
        pts.push_back({clampCoord(x.second), clampCoord(y.second)});
    }

    size_t base = 1 + 2 * (size_t)n;
    auto q0 = readInt(base);
    int q = q0.first ? (int)max<ll>(0, min<ll>(N_MAX, q0.second)) : 0;
    vector<Point> queries;
    for (int i = 0; i < q; i++) {
        auto x = readInt(base + 1 + 2 * (size_t)i);
        auto y = readInt(base + 2 + 2 * (size_t)i);
        if (!x.first || !y.first) break;
        queries.push_back({clampCoord(x.second), clampCoord(y.second)});
    }

    cout << run(pts, queries) << "\n";
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    vector<string> tokens;
    bool any = false;
    string cur;
    int c;
    while ((c = cin.get()) != EOF) {
        if (isspace(c)) {
            if (!cur.empty()) { tokens.push_back(cur); cur.clear(); }
        } else {
            cur.push_back((char)c);
            any = true;
        }
    }
    if (!cur.empty()) tokens.push_back(cur);

    if (any) {
        runIO(tokens);
    } else {
        runTests();
        cout << "all tests passed" << endl;
    }
    return 0;
}
