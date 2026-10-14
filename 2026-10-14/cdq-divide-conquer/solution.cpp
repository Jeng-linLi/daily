// CDQ 分治：逆序對計數 / 三維偏序（陌上花開）/ 分治優化 DP
//
// 題意：
//     給三個長度 n 的整數陣列 a、b、c（第 i 個三元組為 (a_i,b_i,c_i)）與權重 w，輸出：
//       1. b 的逆序對個數；
//       2. 對每個 i，f(i) = #{ j ≠ i : a_j≤a_i, b_j≤b_i, c_j≤c_i }；
//       3. 分治優化 DP：dp[i] = w_i + max{ dp[j] : j<i 且 a_j≤a_i }（無合法 j 時 dp[i]=w_i）；
//       4. max(dp)（n = 0 時為 0）。
//     全部用 CDQ 分治完成，並與 O(n²) 暴力對拍。
//
// 思路：
//     ### CDQ 分治在做什麼
//     普通分治是「切兩半 → 各自遞迴 → 合併」。CDQ 分治多一步：
//     **統計「跨左右兩半」的貢獻**──有多少對 (u,v) 滿足 u 在左、v 在右且滿足某條件。
//     每一對 (u,v) 在遞迴樹上恰好在「它們分開的那一層」被統計一次，
//     於是全體配對的答案就是所有跨半貢獻之和。這就是把 O(n²) 降到 O(n log n) 的關鍵。
//
//     ### ① 逆序對（一維 CDQ）
//     求 #{(i,j): i<j 且 b_i > b_j}。左右兩半都已依 b 排序，歸併時
//     右半元素 b_j 讓左半指標停下來時，左半剩下的每一個都比它大 → 一次加一整段。
//
//     ### ② 三維偏序（CDQ + 樹狀陣列）
//     先依 (a,b,c) 字典序排序。已依 a 升序 ⟹ 左半位置 u < 右半位置 v 必有 a_u ≤ a_v，
//     a 這一維自動滿足，只剩 b、c 兩維。CDQ 的不變量是「返回後區段依 b 排序」，
//     於是用兩個指標掃描、左半走過的元素把 c 加到樹狀陣列（c 離散化），
//     右半每個 v 查「c ≤ c_v 的個數」即為左半對它的貢獻。O(n log² n)。
//
//     **相同三元組要補**：CDQ 只數到「排序後位置 u<v」的配對，題意卻要 j≠i 的全部。
//     若相同三元組佔排序後位置 p..p+s-1，則位置 v 的元素還有 s-1-(v-p) 個組內同伴
//     排在它後面沒被數到。注意分組資訊必須以**原始下標**為鍵保存
//     （CDQ 返回時陣列已依 b 排序，位置 v 早已換人）。
//
//     ### ③ 分治優化 DP
//     dp[i] = w_i + max{ dp[j] : j<i 且 a_j≤a_i }，看成「二維偏序上的最大值」。
//     CDQ(l,r) 的順序（最容易寫錯的地方）：
//       1. 遞迴左半 cdq(l,m)：左半 dp 算完，且依 a 排序；
//       2. 跨半更新：左半預處理「dp 前綴最大值」，右半每個 i 二分找「a ≤ a_i 的最右位置」，
//          O(log n) 拿到 max{dp[j]}；右半此時仍是原始下標順序，正好符合 j<i；
//       3. 再遞迴右半 cdq(m,r)；
//       4. 依 a 歸併，維持不變量。
//     第 2 步必須在第 3 步之前。
//     另外「沒有合法 j」與「有個 dp[j] = 0 的合法 j」語義不同
//     （後者可能是負數），所以要用 NEG_INF 哨兵區分，不能讓 w_i 直接參與 max。
//
// 複雜度：逆序對 O(n log n)；三維偏序 O(n log² n)；分治 DP O(n log n)；暴力 O(n²)。
//
// 輸入（stdin，空白分隔整數）：
//     n
//     a_1 ... a_n
//     b_1 ... b_n
//     c_1 ... c_n
//     w_1 ... w_n
//     非數字 token 視為輸入結束；讀不滿的那一列直接丟棄；n 截斷到 [0,3000]。
// 輸出（stdout，4 行）：
//     1: 逆序對個數  2: f(1)..f(n)  3: dp_1..dp_n  4: max(dp)（n=0 時 0）
// 無 stdin 時跑內置斷言測試並輸出 all tests passed。

#include <cassert>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <limits>
#include <random>
#include <sstream>
#include <string>
#include <vector>

#ifdef _WIN32
#include <fcntl.h>
#include <io.h>
#endif

namespace {

constexpr int MAXN = 3000;
constexpr long long VMAX = 1000000000000LL;
constexpr long long NEG_INF = -1000000000000000000LL;   // -10^18

// ---------------------------------------------------------------- 樹狀陣列


class Fenwick {
public:
    explicit Fenwick(int n) : n_(n), t_(static_cast<std::size_t>(n) + 1, 0) {}

    void add(int i, long long delta) {
        while (i <= n_) {
            t_[static_cast<std::size_t>(i)] += delta;
            i += i & (-i);
        }
    }

    long long sum(int i) const {   // 前綴和 [1, i]
        if (i > n_) i = n_;
        long long res = 0;
        while (i > 0) {
            res += t_[static_cast<std::size_t>(i)];
            i -= i & (-i);
        }
        return res;
    }

private:
    int n_;
    std::vector<long long> t_;
};

// ---------------------------------------------------------------- ① 逆序對


long long inversionCount(std::vector<long long> a) {   // 傳值：內部會重排
    const int n = static_cast<int>(a.size());
    std::vector<long long> tmp(static_cast<std::size_t>(n));
    // 用顯式疊代實作的遞迴 lambda 會很囉嗦，這裡開一個檔案內部的函式物件
    struct Rec {
        std::vector<long long>* a;
        std::vector<long long>* tmp;
        long long run(int l, int r) {
            if (r - l <= 1) return 0;
            int m = (l + r) / 2;
            long long cnt = run(l, m) + run(m, r);
            int i = l, j = m, k = l;
            std::vector<long long>& av = *a;
            std::vector<long long>& tv = *tmp;
            while (i < m && j < r) {
                if (av[static_cast<std::size_t>(i)] <= av[static_cast<std::size_t>(j)]) {
                    tv[static_cast<std::size_t>(k++)] = av[static_cast<std::size_t>(i++)];
                } else {
                    cnt += m - i;                 // 左半剩下的全部 > av[j]
                    tv[static_cast<std::size_t>(k++)] = av[static_cast<std::size_t>(j++)];
                }
            }
            while (i < m) tv[static_cast<std::size_t>(k++)] = av[static_cast<std::size_t>(i++)];
            while (j < r) tv[static_cast<std::size_t>(k++)] = av[static_cast<std::size_t>(j++)];
            for (int t = l; t < r; ++t)
                av[static_cast<std::size_t>(t)] = tv[static_cast<std::size_t>(t)];
            return cnt;
        }
    };
    Rec rec{&a, &tmp};
    return rec.run(0, n);
}

// ---------------------------------------------------------------- ② 三維偏序


int lowerBound(const std::vector<long long>& arr, long long x) {
    int lo = 0, hi = static_cast<int>(arr.size());
    while (lo < hi) {
        int mid = (lo + hi) / 2;
        if (arr[static_cast<std::size_t>(mid)] < x) lo = mid + 1;
        else hi = mid;
    }
    return lo;
}


void compress(const std::vector<long long>& vals, std::vector<int>& out, int& m) {
    std::vector<long long> uniq = vals;
    std::sort(uniq.begin(), uniq.end());
    uniq.erase(std::unique(uniq.begin(), uniq.end()), uniq.end());
    m = static_cast<int>(uniq.size());
    out.resize(vals.size());
    for (std::size_t i = 0; i < vals.size(); ++i)
        out[i] = lowerBound(uniq, vals[i]) + 1;
}


std::vector<long long> partialOrder3d(const std::vector<long long>& a,
                                      const std::vector<long long>& b,
                                      const std::vector<long long>& c) {
    const int n = static_cast<int>(a.size());
    if (n == 0) return {};
    std::vector<int> cr;
    int m = 0;
    compress(c, cr, m);

    std::vector<int> order(static_cast<std::size_t>(n));
    for (int i = 0; i < n; ++i) order[static_cast<std::size_t>(i)] = i;
    std::sort(order.begin(), order.end(), [&](int x, int y) {
        if (a[static_cast<std::size_t>(x)] != a[static_cast<std::size_t>(y)])
            return a[static_cast<std::size_t>(x)] < a[static_cast<std::size_t>(y)];
        if (b[static_cast<std::size_t>(x)] != b[static_cast<std::size_t>(y)])
            return b[static_cast<std::size_t>(x)] < b[static_cast<std::size_t>(y)];
        if (c[static_cast<std::size_t>(x)] != c[static_cast<std::size_t>(y)])
            return c[static_cast<std::size_t>(x)] < c[static_cast<std::size_t>(y)];
        return x < y;
    });

    std::vector<long long> sa(static_cast<std::size_t>(n));
    std::vector<long long> sb(static_cast<std::size_t>(n));
    std::vector<int> sc(static_cast<std::size_t>(n));
    for (int v = 0; v < n; ++v) {
        int i = order[static_cast<std::size_t>(v)];
        sa[static_cast<std::size_t>(v)] = a[static_cast<std::size_t>(i)];
        sb[static_cast<std::size_t>(v)] = b[static_cast<std::size_t>(i)];
        sc[static_cast<std::size_t>(v)] = cr[static_cast<std::size_t>(i)];
    }
    std::vector<long long> g(static_cast<std::size_t>(n), 0);

    // 分組資訊以原始下標為鍵（CDQ 會打亂位置）
    std::vector<int> posOrig(static_cast<std::size_t>(n), 0);
    std::vector<int> grpStart(static_cast<std::size_t>(n), 0);
    std::vector<int> grpSize(static_cast<std::size_t>(n), 0);
    for (int v = 0; v < n; ++v) posOrig[static_cast<std::size_t>(order[static_cast<std::size_t>(v)])] = v;
    int p = 0;
    while (p < n) {
        int q = p + 1;
        while (q < n && sa[static_cast<std::size_t>(q)] == sa[static_cast<std::size_t>(p)] &&
               sb[static_cast<std::size_t>(q)] == sb[static_cast<std::size_t>(p)] &&
               sc[static_cast<std::size_t>(q)] == sc[static_cast<std::size_t>(p)])
            ++q;
        for (int v = p; v < q; ++v) {
            grpStart[static_cast<std::size_t>(order[static_cast<std::size_t>(v)])] = p;
            grpSize[static_cast<std::size_t>(order[static_cast<std::size_t>(v)])] = q - p;
        }
        p = q;
    }

    Fenwick bit(m);
    std::vector<long long> tmpA(static_cast<std::size_t>(n));
    std::vector<long long> tmpB(static_cast<std::size_t>(n));
    std::vector<int> tmpC(static_cast<std::size_t>(n));
    std::vector<int> tmpO(static_cast<std::size_t>(n));
    std::vector<long long> tmpG(static_cast<std::size_t>(n));

    struct Rec {
        std::vector<long long>* sa;
        std::vector<long long>* sb;
        std::vector<int>* sc;
        std::vector<int>* order;
        std::vector<long long>* g;
        std::vector<long long>* tmpA;
        std::vector<long long>* tmpB;
        std::vector<int>* tmpC;
        std::vector<int>* tmpO;
        std::vector<long long>* tmpG;
        Fenwick* bit;

        void run(int l, int r) {
            if (r - l <= 1) return;
            int m = (l + r) / 2;
            run(l, m);
            run(m, r);
            // 跨半統計：左右都已依 b 排序
            int i = l;
            for (int v = m; v < r; ++v) {
                while (i < m && (*sb)[static_cast<std::size_t>(i)] <= (*sb)[static_cast<std::size_t>(v)]) {
                    bit->add((*sc)[static_cast<std::size_t>(i)], 1);
                    ++i;
                }
                (*g)[static_cast<std::size_t>(v)] += bit->sum((*sc)[static_cast<std::size_t>(v)]);
            }
            for (int t = l; t < i; ++t) bit->add((*sc)[static_cast<std::size_t>(t)], -1);
            // 依 (b, c) 歸併
            int pp = l, qq = m, k = l;
            while (pp < m && qq < r) {
                bool takeLeft;
                if ((*sb)[static_cast<std::size_t>(pp)] != (*sb)[static_cast<std::size_t>(qq)])
                    takeLeft = (*sb)[static_cast<std::size_t>(pp)] < (*sb)[static_cast<std::size_t>(qq)];
                else
                    takeLeft = (*sc)[static_cast<std::size_t>(pp)] <= (*sc)[static_cast<std::size_t>(qq)];
                int src = takeLeft ? pp++ : qq++;
                (*tmpA)[static_cast<std::size_t>(k)] = (*sa)[static_cast<std::size_t>(src)];
                (*tmpB)[static_cast<std::size_t>(k)] = (*sb)[static_cast<std::size_t>(src)];
                (*tmpC)[static_cast<std::size_t>(k)] = (*sc)[static_cast<std::size_t>(src)];
                (*tmpO)[static_cast<std::size_t>(k)] = (*order)[static_cast<std::size_t>(src)];
                (*tmpG)[static_cast<std::size_t>(k)] = (*g)[static_cast<std::size_t>(src)];
                ++k;
            }
            while (pp < m) {
                (*tmpA)[static_cast<std::size_t>(k)] = (*sa)[static_cast<std::size_t>(pp)];
                (*tmpB)[static_cast<std::size_t>(k)] = (*sb)[static_cast<std::size_t>(pp)];
                (*tmpC)[static_cast<std::size_t>(k)] = (*sc)[static_cast<std::size_t>(pp)];
                (*tmpO)[static_cast<std::size_t>(k)] = (*order)[static_cast<std::size_t>(pp)];
                (*tmpG)[static_cast<std::size_t>(k)] = (*g)[static_cast<std::size_t>(pp)];
                ++pp;
                ++k;
            }
            while (qq < r) {
                (*tmpA)[static_cast<std::size_t>(k)] = (*sa)[static_cast<std::size_t>(qq)];
                (*tmpB)[static_cast<std::size_t>(k)] = (*sb)[static_cast<std::size_t>(qq)];
                (*tmpC)[static_cast<std::size_t>(k)] = (*sc)[static_cast<std::size_t>(qq)];
                (*tmpO)[static_cast<std::size_t>(k)] = (*order)[static_cast<std::size_t>(qq)];
                (*tmpG)[static_cast<std::size_t>(k)] = (*g)[static_cast<std::size_t>(qq)];
                ++qq;
                ++k;
            }
            for (int t = l; t < r; ++t) {
                (*sa)[static_cast<std::size_t>(t)] = (*tmpA)[static_cast<std::size_t>(t)];
                (*sb)[static_cast<std::size_t>(t)] = (*tmpB)[static_cast<std::size_t>(t)];
                (*sc)[static_cast<std::size_t>(t)] = (*tmpC)[static_cast<std::size_t>(t)];
                (*order)[static_cast<std::size_t>(t)] = (*tmpO)[static_cast<std::size_t>(t)];
                (*g)[static_cast<std::size_t>(t)] = (*tmpG)[static_cast<std::size_t>(t)];
            }
        }
    };
    Rec rec{&sa, &sb, &sc, &order, &g, &tmpA, &tmpB, &tmpC, &tmpO, &tmpG, &bit};
    rec.run(0, n);

    std::vector<long long> res(static_cast<std::size_t>(n), 0);
    for (int v = 0; v < n; ++v) {
        int i = order[static_cast<std::size_t>(v)];
        long long extra = static_cast<long long>(grpSize[static_cast<std::size_t>(i)]) - 1 -
                          (posOrig[static_cast<std::size_t>(i)] - grpStart[static_cast<std::size_t>(i)]);
        res[static_cast<std::size_t>(i)] = g[static_cast<std::size_t>(v)] + extra;
    }
    return res;
}

// ---------------------------------------------------------------- ③ 分治優化 DP


std::vector<long long> cdqDp(const std::vector<long long>& a,
                             const std::vector<long long>& w) {
    const int n = static_cast<int>(a.size());
    if (n == 0) return {};
    std::vector<int> idx(static_cast<std::size_t>(n));
    std::vector<long long> aval(static_cast<std::size_t>(n));
    std::vector<long long> acc(static_cast<std::size_t>(n), NEG_INF);
    std::vector<long long> dval(static_cast<std::size_t>(n), 0);
    std::vector<int> tmpI(static_cast<std::size_t>(n));
    std::vector<long long> tmpA(static_cast<std::size_t>(n));
    std::vector<long long> tmpC(static_cast<std::size_t>(n));
    for (int i = 0; i < n; ++i) {
        idx[static_cast<std::size_t>(i)] = i;
        aval[static_cast<std::size_t>(i)] = a[static_cast<std::size_t>(i)];
    }

    struct Rec {
        std::vector<int>* idx;
        std::vector<long long>* aval;
        std::vector<long long>* acc;
        std::vector<long long>* dval;
        std::vector<int>* tmpI;
        std::vector<long long>* tmpA;
        std::vector<long long>* tmpC;
        const std::vector<long long>* w;

        void run(int l, int r) {
            if (r - l <= 1) return;
            int m = (l + r) / 2;
            // 1) 左半先算完（左半依 a 排序，dp 已是最終值）
            run(l, m);
            for (int t = l; t < m; ++t) {
                long long gain = (*acc)[static_cast<std::size_t>(t)];
                long long base = (*w)[static_cast<std::size_t>((*idx)[static_cast<std::size_t>(t)])];
                (*dval)[static_cast<std::size_t>(t)] = (gain == NEG_INF) ? base : base + gain;
            }
            // 2) 跨半更新：前綴最大值 + 二分
            std::vector<long long> pref(static_cast<std::size_t>(m - l));
            long long best = NEG_INF;
            for (int t = l; t < m; ++t) {
                long long v = (*dval)[static_cast<std::size_t>(t)];
                if (v > best) best = v;
                pref[static_cast<std::size_t>(t - l)] = best;
            }
            std::vector<long long> keys(static_cast<std::size_t>(m - l));
            for (int t = l; t < m; ++t)
                keys[static_cast<std::size_t>(t - l)] = (*aval)[static_cast<std::size_t>(t)];
            for (int v = m; v < r; ++v) {
                int lo = 0, hi = m - l;
                while (lo < hi) {
                    int mid = (lo + hi) / 2;
                    if (keys[static_cast<std::size_t>(mid)] <= (*aval)[static_cast<std::size_t>(v)])
                        lo = mid + 1;
                    else
                        hi = mid;
                }
                if (lo > 0 && pref[static_cast<std::size_t>(lo - 1)] > (*acc)[static_cast<std::size_t>(v)])
                    (*acc)[static_cast<std::size_t>(v)] = pref[static_cast<std::size_t>(lo - 1)];
            }
            // 3) 右半繼續算
            run(m, r);
            // 4) 依 (a, 下標) 歸併
            int p = l, q = m, k = l;
            while (p < m && q < r) {
                bool takeLeft;
                if ((*aval)[static_cast<std::size_t>(p)] != (*aval)[static_cast<std::size_t>(q)])
                    takeLeft = (*aval)[static_cast<std::size_t>(p)] < (*aval)[static_cast<std::size_t>(q)];
                else
                    takeLeft = (*idx)[static_cast<std::size_t>(p)] <= (*idx)[static_cast<std::size_t>(q)];
                int src = takeLeft ? p++ : q++;
                (*tmpI)[static_cast<std::size_t>(k)] = (*idx)[static_cast<std::size_t>(src)];
                (*tmpA)[static_cast<std::size_t>(k)] = (*aval)[static_cast<std::size_t>(src)];
                (*tmpC)[static_cast<std::size_t>(k)] = (*acc)[static_cast<std::size_t>(src)];
                ++k;
            }
            while (p < m) {
                (*tmpI)[static_cast<std::size_t>(k)] = (*idx)[static_cast<std::size_t>(p)];
                (*tmpA)[static_cast<std::size_t>(k)] = (*aval)[static_cast<std::size_t>(p)];
                (*tmpC)[static_cast<std::size_t>(k)] = (*acc)[static_cast<std::size_t>(p)];
                ++p;
                ++k;
            }
            while (q < r) {
                (*tmpI)[static_cast<std::size_t>(k)] = (*idx)[static_cast<std::size_t>(q)];
                (*tmpA)[static_cast<std::size_t>(k)] = (*aval)[static_cast<std::size_t>(q)];
                (*tmpC)[static_cast<std::size_t>(k)] = (*acc)[static_cast<std::size_t>(q)];
                ++q;
                ++k;
            }
            for (int t = l; t < r; ++t) {
                (*idx)[static_cast<std::size_t>(t)] = (*tmpI)[static_cast<std::size_t>(t)];
                (*aval)[static_cast<std::size_t>(t)] = (*tmpA)[static_cast<std::size_t>(t)];
                (*acc)[static_cast<std::size_t>(t)] = (*tmpC)[static_cast<std::size_t>(t)];
            }
        }
    };
    Rec rec{&idx, &aval, &acc, &dval, &tmpI, &tmpA, &tmpC, &w};
    rec.run(0, n);

    std::vector<long long> out(static_cast<std::size_t>(n), 0);
    for (int k = 0; k < n; ++k) {
        int i = idx[static_cast<std::size_t>(k)];
        long long gain = acc[static_cast<std::size_t>(k)];
        long long base = w[static_cast<std::size_t>(i)];
        out[static_cast<std::size_t>(i)] = (gain == NEG_INF) ? base : base + gain;
    }
    return out;
}

// ---------------------------------------------------------------- 主流程


std::string joinLongs(const std::vector<long long>& v) {
    std::ostringstream oss;
    for (std::size_t i = 0; i < v.size(); ++i) {
        if (i) oss << ' ';
        oss << v[i];
    }
    return oss.str();
}


std::string run(const std::vector<long long>& a, const std::vector<long long>& b,
                const std::vector<long long>& c, const std::vector<long long>& w) {
    long long inv = inversionCount(b);
    std::vector<long long> f = partialOrder3d(a, b, c);
    std::vector<long long> dp = cdqDp(a, w);
    long long best = dp.empty() ? 0 : *std::max_element(dp.begin(), dp.end());
    std::ostringstream oss;
    oss << inv << "\n";
    oss << joinLongs(f) << "\n";
    oss << joinLongs(dp) << "\n";
    oss << best;
    return oss.str();
}

// ---------------------------------------------------------------- 對拍用暴力


long long bfInversion(const std::vector<long long>& b) {
    long long cnt = 0;
    int n = static_cast<int>(b.size());
    for (int i = 0; i < n; ++i)
        for (int j = i + 1; j < n; ++j)
            if (b[static_cast<std::size_t>(i)] > b[static_cast<std::size_t>(j)]) ++cnt;
    return cnt;
}


std::vector<long long> bfPartialOrder(const std::vector<long long>& a,
                                      const std::vector<long long>& b,
                                      const std::vector<long long>& c) {
    int n = static_cast<int>(a.size());
    std::vector<long long> res(static_cast<std::size_t>(n), 0);
    for (int i = 0; i < n; ++i) {
        long long cnt = 0;
        for (int j = 0; j < n; ++j)
            if (j != i && a[static_cast<std::size_t>(j)] <= a[static_cast<std::size_t>(i)] &&
                b[static_cast<std::size_t>(j)] <= b[static_cast<std::size_t>(i)] &&
                c[static_cast<std::size_t>(j)] <= c[static_cast<std::size_t>(i)])
                ++cnt;
        res[static_cast<std::size_t>(i)] = cnt;
    }
    return res;
}


std::vector<long long> bfDp(const std::vector<long long>& a,
                            const std::vector<long long>& w) {
    int n = static_cast<int>(a.size());
    std::vector<long long> dp(static_cast<std::size_t>(n), 0);
    for (int i = 0; i < n; ++i) {
        bool has = false;
        long long best = 0;
        for (int j = 0; j < i; ++j)
            if (a[static_cast<std::size_t>(j)] <= a[static_cast<std::size_t>(i)] &&
                (!has || dp[static_cast<std::size_t>(j)] > best)) {
                best = dp[static_cast<std::size_t>(j)];
                has = true;
            }
        dp[static_cast<std::size_t>(i)] = has ? w[static_cast<std::size_t>(i)] + best
                                              : w[static_cast<std::size_t>(i)];
    }
    return dp;
}


std::string bfRun(const std::vector<long long>& a, const std::vector<long long>& b,
                  const std::vector<long long>& c, const std::vector<long long>& w) {
    std::vector<long long> dp = bfDp(a, w);
    long long best = dp.empty() ? 0 : *std::max_element(dp.begin(), dp.end());
    std::ostringstream oss;
    oss << bfInversion(b) << "\n";
    oss << joinLongs(bfPartialOrder(a, b, c)) << "\n";
    oss << joinLongs(dp) << "\n";
    oss << best;
    return oss.str();
}


void check(const std::vector<long long>& a, const std::vector<long long>& b,
           const std::vector<long long>& c, const std::vector<long long>& w) {
    std::string got = run(a, b, c, w);
    std::string want = bfRun(a, b, c, w);
    assert(got == want);
}

// ---------------------------------------------------------------- 測試


void runTests() {
    assert(run({}, {}, {}, {}) == "0\n\n\n0");
    assert(run({1}, {5}, {7}, {3}) == "0\n0\n3\n3");
    assert(run({1, 2}, {2, 1}, {1, 2}, {1, 2}) == "1\n0 0\n1 3\n3");
    assert(run({1, 1}, {1, 1}, {1, 1}, {5, 5}) == "0\n1 1\n5 10\n10");

    check({}, {}, {}, {});
    check({1}, {1}, {1}, {1});
    check({1, 2}, {2, 1}, {1, 2}, {1, 2});
    check({1, 1}, {1, 1}, {1, 1}, {5, 5});
    check({3, 1, 2}, {1, 2, 3}, {2, 2, 2}, {4, -1, 7});
    check({1, 2, 3}, {3, 2, 1}, {1, 1, 1}, {1, 1, 1});

    std::mt19937 rng(20261014);

    for (int iter = 0; iter < 320; ++iter) {
        int n = static_cast<int>(rng() % 13);
        std::vector<long long> a, b, c, w;
        for (int i = 0; i < n; ++i) {
            a.push_back(static_cast<long long>(rng() % 4));
            b.push_back(static_cast<long long>(rng() % 4));
            c.push_back(static_cast<long long>(rng() % 4));
            w.push_back(static_cast<long long>(rng() % 11) - 5);
        }
        check(a, b, c, w);
    }

    for (int iter = 0; iter < 120; ++iter) {
        int n = static_cast<int>(rng() % 41);
        std::vector<long long> a, b, c, w;
        for (int i = 0; i < n; ++i) {
            a.push_back(static_cast<long long>(rng() % 101) - 50);
            b.push_back(static_cast<long long>(rng() % 101) - 50);
            c.push_back(static_cast<long long>(rng() % 101) - 50);
            w.push_back(static_cast<long long>(rng() % 61) - 30);
        }
        check(a, b, c, w);
    }

    for (int iter = 0; iter < 160; ++iter) {
        int n = 1 + static_cast<int>(rng() % 30);
        std::vector<long long> a, b, c, w;
        for (int i = 0; i < n; ++i) {
            a.push_back(static_cast<long long>(rng() % 2));
            b.push_back(static_cast<long long>(rng() % 2));
            c.push_back(static_cast<long long>(rng() % 2));
            w.push_back(static_cast<long long>(rng() % 9) - 4);
        }
        check(a, b, c, w);
    }

    for (int n : {1, 2, 3, 17, 64}) {
        std::vector<long long> inc, dec, zeros, ones, w;
        for (int i = 0; i < n; ++i) {
            inc.push_back(i);
            dec.push_back(n - 1 - i);
            zeros.push_back(0);
            ones.push_back(1);
            w.push_back(i + 1);
        }
        check(inc, inc, inc, ones);
        check(inc, dec, zeros, w);
    }

    // 中量級效能（只確認不崩且與暴力一致）
    {
        int n = 3000;
        std::vector<long long> a, b, c, w;
        for (int i = 0; i < n; ++i) {
            a.push_back(static_cast<long long>(rng() % 2000001) - 1000000);
            b.push_back(static_cast<long long>(rng() % 2000001) - 1000000);
            c.push_back(static_cast<long long>(rng() % 2000001) - 1000000);
            w.push_back(static_cast<long long>(rng() % 201) - 100);
        }
        std::string out = run(a, b, c, w);
        assert(!out.empty());
        assert(inversionCount(b) == bfInversion(b));
        assert(cdqDp(a, w) == bfDp(a, w));
    }
}

// ---------------------------------------------------------------- IO


bool isNumToken(const std::string& s) {
    if (s.empty()) return false;
    std::size_t i = 0;
    if (s[0] == '+' || s[0] == '-') i = 1;
    if (i >= s.size()) return false;
    for (std::size_t j = i; j < s.size(); ++j)
        if (s[j] < '0' || s[j] > '9') return false;
    return true;
}


std::vector<std::string> splitWs(const std::string& raw) {
    std::vector<std::string> out;
    std::size_t i = 0;
    while (i < raw.size()) {
        while (i < raw.size() && static_cast<unsigned char>(raw[i]) <= ' ') ++i;
        std::size_t start = i;
        while (i < raw.size() && static_cast<unsigned char>(raw[i]) > ' ') ++i;
        if (i > start) out.push_back(raw.substr(start, i - start));
    }
    return out;
}


class Reader {
public:
    explicit Reader(const std::vector<std::string>& tok) : tok_(tok) {}

    bool readInt(long long& out) {
        if (pos_ >= tok_.size()) return false;
        if (!isNumToken(tok_[pos_])) return false;
        const std::string& s = tok_[pos_];
        bool neg = false;
        std::size_t i = 0;
        if (s[0] == '+' || s[0] == '-') {
            neg = (s[0] == '-');
            i = 1;
        }
        long long v = 0;
        for (; i < s.size(); ++i) {
            if (v > VMAX / 10) { v = VMAX; break; }
            v = v * 10 + (s[i] - '0');
            if (v > VMAX) { v = VMAX; break; }
        }
        out = neg ? -v : v;
        ++pos_;
        return true;
    }

private:
    const std::vector<std::string>& tok_;
    std::size_t pos_ = 0;
};


void runIo(const std::string& raw) {
    std::vector<std::string> tok = splitWs(raw);
    Reader rd(tok);

    long long t;
    int n = rd.readInt(t) ? static_cast<int>(std::max<long long>(0, std::min<long long>(MAXN, t))) : 0;

    std::vector<std::vector<long long>> cols(4);
    for (int k = 0; k < 4; ++k) {
        for (int i = 0; i < n; ++i) {
            if (!rd.readInt(t)) break;
            cols[static_cast<std::size_t>(k)].push_back(t);
        }
    }
    int m = n;
    for (int k = 0; k < 4; ++k)
        m = std::min(m, static_cast<int>(cols[static_cast<std::size_t>(k)].size()));
    for (int k = 0; k < 4; ++k)
        cols[static_cast<std::size_t>(k)].resize(static_cast<std::size_t>(m));

    std::cout << run(cols[0], cols[1], cols[2], cols[3]) << "\n";
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
        if (static_cast<unsigned char>(c) > ' ') { blank = false; break; }
    }
    if (!blank) {
        runIo(raw);
    } else {
        runTests();
        std::cout << "all tests passed" << std::endl;
    }
    return 0;
}
