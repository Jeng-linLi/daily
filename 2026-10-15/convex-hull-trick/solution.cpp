// 凸包優化（Convex Hull Trick）與李超線段樹（Li Chao Segment Tree）
//
// 題意：
//     給定 k 條直線 y = m·x + b、q 個查詢點 x，以及一個長度 n 的陣列 a。
//     需要回答三件事：
//       1. 對每個 x，求所有直線在 x 上的最小值 —— 用**李超線段樹** O(log C)，
//          並與 O(k) 的暴力逐點驗證。
//       2. 把直線依斜率排序後建**單調凸包**（deque），查詢點排序後 O(1) 均攤回答。
//       3. 斜率優化 DP：把 a 切成若干段，每段成本 = (段內元素和)² + C，求最小總成本。
//          dp[i] = min_{j<i} { dp[j] + (pre[i] - pre[j])² + C }，
//          暴力 O(n²) 與李超線段樹 O(n log C) 互相驗證。
//
// 思路：
//     ### 為什麼需要「凸包優化」
//     一大類 DP 可以寫成
//         dp[i] = min_{j<i} { dp[j] + f(j)·g(i) + h(i) + w(j) }
//     把 j 視成一條直線：斜率 m_j = f(j)、截距 b_j = dp[j] + w(j)，
//     而 i 的轉移就是在 x = g(i) 上取所有直線的最小值。
//     於是把「O(n) 找最小值」換成「一堆直線求 min」，就有了優化的空間。
//
//     ### 李超線段樹：把「直線」當成區間上的懶標記
//     定義域 x ∈ [L, R] 建立在線段樹上，每個節點存一條「在該區間中點最優」的直線。
//     插入新直線 nw 時，跟節點當前的直線 cur 比中點：
//       - 若 nw 在中點更優 → 交換，nw 變成「被擠掉的那一條」往下傳；
//       - 若 nw 在左端點更優 → 它只可能在左半區間勝出，往左遞迴；
//       - 否則若 nw 在右端點更優 → 往右遞迴；
//       - 兩端都不優 → nw 在整個區間都被 cur 壓住，直接丟掉。
//     正確性來自「兩條直線最多交一次」：在一個區間裡，若 nw 在兩端都輸給 cur，
//     則它在區間內部也不可能贏。每個節點只存一條線，插入 O(log C)、查詢 O(log C)。
//     優點是**完全不要求斜率或查詢點單調**，缺點是需要離散（或座標壓縮）的定義域。
//
//     ### 單調凸包（dq-CHT）：斜率與查詢都單調時可以做到 O(1)
//     把直線依**斜率遞減**排列後建下凸包：斜率最大者在 x → -∞ 時最優，
//     斜率最小者在 x → +∞ 時最優，所以查詢點 x 遞增時，最佳直線在 hull 上的
//     位置只會單調往前走，用一根指標就夠，全部查詢 O(q + k)。
//       - 加入：若最後兩條與新線的交點不滿足「交點遞增」，彈掉最後一條；
//       - 查詢：若下一條比當前好，指標就前進（指標不回頭）。
//     判斷「中間那條是不是多餘的」用整數交叉相乘，完全不碰浮點數
//     （斜率遞減：m1 > m2 > m3）：
//         交點順序亂掉 ⟺ (b2-b1)(m2-m3) ≥ (b3-b2)(m1-m2)
//
//     ### 兩者的取捨
//       李超線段樹：任意順序、任意斜率，O(log C)，程式稍長，定義域要離散。
//       單調凸包　：要求斜率單調 + 查詢單調，O(1) 均攤，常數極小。
//       實務上：DP 轉移順序天然單調時優先單調凸包；否則用李超。
//
// 複雜度（k 條線、q 次查詢、定義域寬度 C、DP 長度 n）：
//     李超線段樹插入 / 查詢   O(log C)，空間 O(k log C)
//     單調凸包建構 O(k log k)（排序）/ O(k)（已排序）；全部查詢 O(q + k)
//     斜率優化 DP             O(n log C)（李超）或 O(n)（單調凸包）

#include <algorithm>
#include <cassert>
#include <cctype>
#include <fcntl.h>
#include <io.h>
#include <iostream>
#include <random>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace {

const long long INF = 4000000000000000000LL;  // 整數哨兵（不用 float，避免污染整數輸出）
const long long XLO = -10000000LL;            // 李超線段樹定義域（離散整數座標）
const long long XHI = 10000000LL;

const int MAXK = 200;                    // 直線數上限
const int MAXQ = 200;                    // 查詢數上限
const int MAXN = 150;                    // DP 陣列長度上限
const long long MLIM = 1000000LL;        // 斜率截斷
const long long BLIM = 1000000000LL;     // 截距截斷
const long long XLIM = 1000000LL;        // 查詢點截斷
const long long ALIM = 10000LL;          // DP 陣列元素截斷
const long long COST = 10;               // 每段固定成本 C

typedef std::pair<long long, long long> Line;  // (m, b)

long long ev(const Line& ln, long long x) { return ln.first * x + ln.second; }

// ------------------------------------------------------- 李超線段樹

class LiChao {
public:
    LiChao(long long xlo, long long xhi) : xlo_(xlo), xhi_(xhi) {
        line_.push_back(Line(0, 0));
        has_.push_back(false);
        lc_.push_back(-1);
        rc_.push_back(-1);
    }

    void add(const Line& nw) { addRec(nw, 0, xlo_, xhi_); }

    long long query(long long x) const { return queryRec(x, 0, xlo_, xhi_); }

private:
    int node() {
        line_.push_back(Line(0, 0));
        has_.push_back(false);
        lc_.push_back(-1);
        rc_.push_back(-1);
        return static_cast<int>(line_.size()) - 1;
    }

    void addRec(Line nw, int v, long long l, long long r) {
        if (!has_[v]) {
            line_[v] = nw;
            has_[v] = true;
            return;
        }
        Line cur = line_[v];
        // 用 l + (r - l) / 2 保證兩個語言在負座標下算出同一個中點
        long long mid = l + (r - l) / 2;
        if (ev(nw, mid) < ev(cur, mid)) {
            std::swap(nw, cur);
            line_[v] = cur;
        }
        if (l == r) return;
        if (ev(nw, l) < ev(cur, l)) {
            if (lc_[v] == -1) lc_[v] = node();
            addRec(nw, lc_[v], l, mid);
        } else if (ev(nw, r) < ev(cur, r)) {
            if (rc_[v] == -1) rc_[v] = node();
            addRec(nw, rc_[v], mid + 1, r);
        }
        // 兩端都不優 → nw 在整個區間都被 cur 壓住（兩線只交一次），丟掉
    }

    long long queryRec(long long x, int v, long long l, long long r) const {
        long long res = has_[v] ? ev(line_[v], x) : INF;
        if (l == r) return res;
        long long mid = l + (r - l) / 2;
        if (x <= mid) {
            if (lc_[v] != -1) {
                long long t = queryRec(x, lc_[v], l, mid);
                if (t < res) res = t;
            }
        } else {
            if (rc_[v] != -1) {
                long long t = queryRec(x, rc_[v], mid + 1, r);
                if (t < res) res = t;
            }
        }
        return res;
    }

    long long xlo_, xhi_;
    std::vector<Line> line_;
    std::vector<bool> has_;
    std::vector<int> lc_, rc_;
};

// ------------------------------------------------------- 單調凸包

bool bad(const Line& l1, const Line& l2, const Line& l3) {
    // 斜率嚴格遞減時，判斷 l2 是否被 l1 與 l3 夾成多餘（全整數運算）
    return (l2.second - l1.second) * (l2.first - l3.first) >=
           (l3.second - l2.second) * (l1.first - l2.first);
}

std::vector<Line> chtBuild(std::vector<Line> lines) {
    // 排序鍵 (-m, b)：斜率由大到小，同斜率時截距小的排前面（後面同斜率的直接捨棄）
    std::sort(lines.begin(), lines.end(), [](const Line& a, const Line& b) {
        if (a.first != b.first) return a.first > b.first;
        return a.second < b.second;
    });
    std::vector<Line> hull;
    for (const Line& ln : lines) {
        if (!hull.empty() && hull.back().first == ln.first) continue;
        while (hull.size() >= 2 && bad(hull[hull.size() - 2], hull[hull.size() - 1], ln)) {
            hull.pop_back();
        }
        hull.push_back(ln);
    }
    return hull;
}

std::vector<long long> chtQuery(const std::vector<Line>& hull,
                                const std::vector<long long>& xs) {
    std::vector<long long> res;
    if (hull.empty()) {
        res.assign(xs.size(), INF);
        return res;
    }
    size_t j = 0;
    for (long long x : xs) {
        while (j + 1 < hull.size() && ev(hull[j], x) >= ev(hull[j + 1], x)) ++j;
        res.push_back(ev(hull[j], x));
    }
    return res;
}

// ------------------------------------------------------- 斜率優化 DP

std::vector<long long> dpLiChao(const std::vector<long long>& a) {
    int n = static_cast<int>(a.size());
    std::vector<long long> pre(n + 1, 0);
    for (int i = 0; i < n; ++i) pre[i + 1] = pre[i] + a[i];
    LiChao lc(XLO, XHI);
    std::vector<long long> dp(n + 1, 0);
    // j = 0 的線：斜率 -2·pre[0] = 0，截距 dp[0] + pre[0]² + C = C
    lc.add(Line(0, COST));
    for (int i = 1; i <= n; ++i) {
        dp[i] = pre[i] * pre[i] + lc.query(pre[i]);
        lc.add(Line(-2 * pre[i], dp[i] + pre[i] * pre[i] + COST));
    }
    return std::vector<long long>(dp.begin() + 1, dp.end());
}

std::vector<long long> dpBrute(const std::vector<long long>& a) {
    int n = static_cast<int>(a.size());
    std::vector<long long> pre(n + 1, 0);
    for (int i = 0; i < n; ++i) pre[i + 1] = pre[i] + a[i];
    std::vector<long long> dp(n + 1, INF);
    dp[0] = 0;
    for (int i = 1; i <= n; ++i) {
        long long best = INF;
        for (int j = 0; j < i; ++j) {
            long long d = pre[i] - pre[j];
            long long v = dp[j] + d * d + COST;
            if (v < best) best = v;
        }
        dp[i] = best;
    }
    return std::vector<long long>(dp.begin() + 1, dp.end());
}

// ------------------------------------------------------- 暴力對拍

std::vector<long long> bruteMin(const std::vector<Line>& lines,
                                const std::vector<long long>& xs) {
    std::vector<long long> out;
    for (long long x : xs) {
        long long best = INF;
        for (const Line& ln : lines) {
            long long v = ev(ln, x);
            if (v < best) best = v;
        }
        out.push_back(best);
    }
    return out;
}

// ------------------------------------------------------- 驅動

std::string join(const std::vector<long long>& v) {
    std::ostringstream oss;
    for (size_t i = 0; i < v.size(); ++i) {
        if (i) oss << ' ';
        oss << v[i];
    }
    return oss.str();
}

std::string run(const std::vector<Line>& lines,
                const std::vector<long long>& xs,
                const std::vector<long long>& a) {
    LiChao lc(XLO, XHI);
    for (const Line& ln : lines) lc.add(ln);
    std::vector<long long> a1;
    for (long long x : xs) a1.push_back(lc.query(x));
    std::vector<long long> a2 = bruteMin(lines, xs);

    std::vector<Line> hull = chtBuild(lines);
    std::vector<long long> sx = xs;
    std::sort(sx.begin(), sx.end());
    std::vector<long long> a3 = chtQuery(hull, sx);
    std::vector<long long> a4 = bruteMin(lines, sx);

    std::vector<long long> d1 = dpLiChao(a);
    std::vector<long long> d2 = dpBrute(a);

    int ok = (a1 == a2 && a3 == a4 && d1 == d2) ? 1 : 0;
    std::ostringstream oss;
    oss << join(a1) << '\n'
        << join(a2) << '\n'
        << join(a3) << '\n'
        << join(d1) << '\n'
        << join(d2) << '\n'
        << ok;
    return oss.str();
}

// ------------------------------------------------------- 內置測試

void runTests() {
    // ---- 固定用例 ----
    std::vector<Line> lines{{1, 0}, {-1, 5}, {0, 3}};
    std::vector<long long> xs{-5, 0, 5, 10};
    std::vector<long long> arr{1, 2, 3};
    std::string out = run(lines, xs, arr);
    assert(out.size() >= 2 && out[out.size() - 1] == '1');

    // 只有一條線：最小值永遠是它自己
    assert(run({{2, 7}}, {0, 1, -1}, {}).substr(0, 6) == "7 9 5\n");

    // 空輸入
    assert(run({}, {}, {}) == "\n\n\n\n\n1");

    // 兩條過原點的線 y = x 與 y = -x：最小值恆為 -|x|
    assert((chtQuery(chtBuild({{1, 0}, {-1, 0}}), {-3, -1, 0, 1, 3}) ==
            std::vector<long long>{-3, -1, 0, -1, -3}));

    // DP：李超版必須與暴力版一致
    assert(dpBrute({1, 2, 3}) == dpLiChao({1, 2, 3}));

    // ---- 隨機對拍 ----
    std::mt19937 rng(20261015);
    for (int iter = 0; iter < 300; ++iter) {
        int k = static_cast<int>(rng() % 31);
        std::vector<Line> ls;
        for (int i = 0; i < k; ++i) {
            ls.push_back({static_cast<long long>(rng() % 41) - 20,
                          static_cast<long long>(rng() % 101) - 50});
        }
        int q = static_cast<int>(rng() % 31);
        std::vector<long long> xq;
        for (int i = 0; i < q; ++i) xq.push_back(static_cast<long long>(rng() % 61) - 30);
        int n = static_cast<int>(rng() % 26);
        std::vector<long long> ar;
        for (int i = 0; i < n; ++i) ar.push_back(static_cast<long long>(rng() % 61) - 30);
        std::string s = run(ls, xq, ar);
        assert(s.size() >= 2 && s[s.size() - 1] == '1');
    }

    // ---- 大座標只用李超（不走 O(k·q) 暴力）----
    std::vector<Line> big;
    for (int i = 0; i < MAXK; ++i) {
        long long m = static_cast<long long>(rng() % (2 * MLIM + 1)) - MLIM;
        long long b = static_cast<long long>(rng() % (2 * BLIM + 1)) - BLIM;
        big.push_back({m, b});
    }
    LiChao lc(XLO, XHI);
    for (const Line& ln : big) lc.add(ln);
    std::vector<long long> pts{XLO, -1, 0, 1, XHI};
    std::vector<long long> expect;
    for (long long x : pts) expect.push_back(bruteMin(big, {x})[0]);
    std::vector<long long> got;
    for (long long x : pts) got.push_back(lc.query(x));
    assert(got == expect);
    // 單調凸包在大樣本上與李超一致（查詢點升序）
    std::vector<long long> sx;
    for (long long x = -XLIM; x <= XLIM; x += 20000) sx.push_back(x);
    std::vector<Line> hull = chtBuild(big);
    std::vector<long long> viaHull = chtQuery(hull, sx);
    std::vector<long long> viaLc;
    for (long long x : sx) viaLc.push_back(lc.query(x));
    assert(viaHull == viaLc);
}

// ------------------------------------------------------- IO

bool isInt(const std::string& t) {
    if (t.empty()) return false;
    size_t i = 0;
    if (t[0] == '+' || t[0] == '-') i = 1;
    if (i >= t.size()) return false;
    for (; i < t.size(); ++i) {
        if (t[i] < '0' || t[i] > '9') return false;
    }
    return true;
}

class Reader {
public:
    explicit Reader(const std::string& raw) {
        std::istringstream iss(raw);
        std::string t;
        while (iss >> t) toks_.push_back(t);
    }
    bool readInt(long long& out) {
        // 非數字 token 一律視為「輸入到此為止」，且**不消費**該 token
        // （與 Python 版一致：後續所有讀取都會失敗）。
        if (pos_ >= toks_.size()) return false;
        const std::string& t = toks_[pos_];
        if (!isInt(t)) return false;
        ++pos_;
        bool neg = false;
        size_t i = 0;
        if (t[0] == '+' || t[0] == '-') {
            neg = (t[0] == '-');
            i = 1;
        }
        const long long LIM = 10000000000000000LL;  // 遠大於所有截斷上限
        long long v = 0;
        for (; i < t.size(); ++i) {
            if (v > LIM / 10) { v = LIM; break; }   // 再乘 10 一定超過 LIM
            v = v * 10 + (t[i] - '0');
            if (v > LIM) { v = LIM; break; }
        }
        out = neg ? -v : v;
        return true;
    }

private:
    std::vector<std::string> toks_;
    size_t pos_ = 0;
};

long long clampAbs(long long v, long long lim) {
    if (v < -lim) return -lim;
    if (v > lim) return lim;
    return v;
}

int clampNonNeg(long long v, int hi) {
    if (v < 0) return 0;
    if (v > hi) return hi;
    return static_cast<int>(v);
}

void runIo(const std::string& raw) {
    Reader rd(raw);
    long long t = 0;
    int k = rd.readInt(t) ? clampNonNeg(t, MAXK) : 0;
    std::vector<Line> lines;
    for (int i = 0; i < k; ++i) {
        long long m, b;
        if (!rd.readInt(m)) break;
        if (!rd.readInt(b)) break;
        lines.push_back({clampAbs(m, MLIM), clampAbs(b, BLIM)});
    }
    int q = rd.readInt(t) ? clampNonNeg(t, MAXQ) : 0;
    std::vector<long long> xs;
    for (int i = 0; i < q; ++i) {
        if (!rd.readInt(t)) break;
        xs.push_back(clampAbs(t, XLIM));
    }
    int n = rd.readInt(t) ? clampNonNeg(t, MAXN) : 0;
    std::vector<long long> a;
    for (int i = 0; i < n; ++i) {
        if (!rd.readInt(t)) break;
        a.push_back(clampAbs(t, ALIM));
    }
    std::cout << run(lines, xs, a) << "\n";
}

}  // namespace

int main() {
#ifdef _WIN32
    _setmode(_fileno(stdout), _O_BINARY);
#endif
    std::string raw((std::istreambuf_iterator<char>(std::cin)),
                    std::istreambuf_iterator<char>());
    bool blank = true;
    for (char c : raw) {
        unsigned char u = static_cast<unsigned char>(c);
        if (u > ' ') { blank = false; break; }
    }
    if (!blank) {
        runIo(raw);
    } else {
        runTests();
        std::cout << "all tests passed" << std::endl;
    }
    return 0;
}
