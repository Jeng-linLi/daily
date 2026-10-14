// 莫隊算法：區間查詢（普通莫隊 + 奇偶優化）與帶修莫隊
//
// 題意：
//     維護整數陣列 a（長度 n）與兩組操作：
//       A. q 個靜態區間查詢 [l,r]：輸出區間內 ① 不同數的個數 ② 區間和 ③ 眾數的出現次數
//       B. m 個帶修改操作：`1 l r` 查詢（在執行到這裡為止的陣列上回答）、`2 p x` 把 a[p] 改成 x
//     A 用普通莫隊（塊排序 + 奇偶優化），B 用帶修莫隊（三維排序），並與 O(n·q) 暴力對拍。
//
// 思路：
//     ### 莫隊的核心：把「區間查詢」變成「指標在陣列上走路」
//     若知道 [L,R] 的答案，且能 O(1) 地在左右兩端加/刪一個元素，就能 O(|Δ|) 移到相鄰區間。
//     本題三個統計量都滿足：
//       sum      ：add 就 += v，remove 就 -= v
//       distinct ：cnt[v] 記出現次數，0→1 時 +1、1→0 時 -1
//       mode     ：再一層「出現次數的出現次數」freqCnt[c] = 有幾個值恰好出現 c 次。
//                  加入時 cnt 增加、超過最大值就更新；刪除時若最大值那一層被清空就往下退
//                  （while maxc > 0 && freqCnt[maxc] == 0: --maxc）。經典 O(1) 維護法。
//
//     ### 排序策略決定總移動量
//     塊長 B ≈ √n，排序鍵 (l 塊, r)，且奇數塊的 r 用**降序**（奇偶優化）：
//       同塊內 r 單調移動 → 每塊 O(n)，共 O(n^1.5)；換塊 r 最多跳 n，共 O(√n) 次；
//       l 在塊內每次最多走 B → 每次查詢 O(√n)。
//     奇偶優化讓 r 在塊與塊之間來回掃而不是每次跳回起點，常數減半。
//
//     ### 帶修莫隊：多一維「時間」
//     每個查詢看成三維點 (l 塊, r 塊, t)，t = 這個查詢之前發生過幾次修改。
//     塊長取 B ≈ n^(2/3)，排序鍵 (l 塊, r 塊, t)。除了動 L、R，還要把時間推向 t：
//       往前：套用修改（若被改位置在當前區間內，先 remove 舊值、add 新值）
//       往後：撤銷修改（同理）
//     複雜度 O(n^(5/3))。
//
//     ### 為什麼兩個語言能逐位元組一致
//     塊長用整數平方根 / 整數立方根（二分）算，不用浮點；排序鍵是純整數並以查詢編號
//     做最終 tie-break，排序結果唯一；而答案與處理順序無關。
//
// 複雜度：普通莫隊 O((n+q)·√n)；帶修莫隊 O(n^(5/3))；暴力 O(n·q)；空間 O(n+q)。
//
// 輸入（stdin，空白分隔整數）：
//     n
//     a_1 ... a_n
//     q                  靜態查詢個數
//     l_1 r_1            q 組，1-indexed 閉區間
//     m                  帶修改操作數
//     m 組 `op x y`：1 l r 查詢 / 2 p x 修改（1-indexed）
//     非數字 token 視為輸入結束；n、q、m 截斷到 [0,3000]；數值截斷到 ±10^9。
// 輸出（stdout）：
//     1: qA      qA 行：`distinct sum mode`
//     1 行: qB   qB 行：`distinct sum mode`
// 無 stdin 時跑內置斷言測試並輸出 all tests passed。

#include <cassert>
#include <algorithm>
#include <array>
#include <iostream>
#include <iterator>
#include <map>
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
constexpr long long VMAX = 1000000000LL;

using Triple = std::array<long long, 3>;
using Op = std::array<long long, 3>;

// ---------------------------------------------------------------- 整數方根


long long isqrtFloor(long long v) {
    long long lo = 0, hi = 1;
    while (hi * hi <= v) hi *= 2;
    lo = hi / 2;
    while (lo < hi) {
        long long mid = (lo + hi + 1) / 2;
        if (mid * mid <= v) lo = mid;
        else hi = mid - 1;
    }
    return lo;
}


long long icbrtFloor(long long v) {
    long long lo = 0, hi = 1;
    while (hi * hi * hi <= v) hi *= 2;
    lo = hi / 2;
    while (lo < hi) {
        long long mid = (lo + hi + 1) / 2;
        if (mid * mid * mid <= v) lo = mid;
        else hi = mid - 1;
    }
    return lo;
}


int blockSizePlain(int n) { return static_cast<int>(std::max<long long>(1, isqrtFloor(n))); }


int blockSizeModified(int n) {
    if (n <= 1) return 1;
    return static_cast<int>(std::max<long long>(1, icbrtFloor(1LL * n * n)));
}

// ---------------------------------------------------------------- 區間統計量


class Counter {
public:
    Counter(int m, int n)
        : cnt_(static_cast<std::size_t>(m) + 2, 0),
          freq_(static_cast<std::size_t>(n) + 4, 0),
          distinct_(0), total_(0), maxc_(0) {}

    void add(int cid, long long val) {
        int c = cnt_[static_cast<std::size_t>(cid)];
        freq_[static_cast<std::size_t>(c)] -= 1;
        ++c;
        cnt_[static_cast<std::size_t>(cid)] = c;
        freq_[static_cast<std::size_t>(c)] += 1;
        if (c == 1) ++distinct_;
        if (c > maxc_) maxc_ = c;
        total_ += val;
    }

    void remove(int cid, long long val) {
        int c = cnt_[static_cast<std::size_t>(cid)];
        freq_[static_cast<std::size_t>(c)] -= 1;
        --c;
        cnt_[static_cast<std::size_t>(cid)] = c;
        freq_[static_cast<std::size_t>(c)] += 1;
        if (c == 0) --distinct_;
        while (maxc_ > 0 && freq_[static_cast<std::size_t>(maxc_)] == 0) --maxc_;
        total_ -= val;
    }

    Triple triple() const { return {distinct_, total_, static_cast<long long>(maxc_)}; }

private:
    std::vector<int> cnt_;
    std::vector<int> freq_;
    long long distinct_;
    long long total_;
    int maxc_;
};

// ---------------------------------------------------------------- ① 普通莫隊


std::vector<Triple> moQueries(const std::vector<long long>& a,
                              const std::vector<std::pair<int, int>>& queries,
                              const std::vector<int>& cid, int m) {
    const int n = static_cast<int>(a.size());
    const int q = static_cast<int>(queries.size());
    std::vector<Triple> ans(static_cast<std::size_t>(q), {0, 0, 0});
    if (n == 0 || q == 0) return ans;

    const int B = blockSizePlain(n);
    std::vector<int> order(static_cast<std::size_t>(q));
    for (int i = 0; i < q; ++i) order[static_cast<std::size_t>(i)] = i;
    std::sort(order.begin(), order.end(), [&](int x, int y) {
        int bx = queries[static_cast<std::size_t>(x)].first / B;
        int by = queries[static_cast<std::size_t>(y)].first / B;
        if (bx != by) return bx < by;
        long long rx = queries[static_cast<std::size_t>(x)].second;
        long long ry = queries[static_cast<std::size_t>(y)].second;
        if (bx % 2 == 1) {           // 奇偶優化：奇數塊的 r 用降序
            if (rx != ry) return rx > ry;
        } else {
            if (rx != ry) return rx < ry;
        }
        return x < y;
    });

    Counter ctr(m, n);
    int L = 0, R = -1;
    for (int idx : order) {
        int l = queries[static_cast<std::size_t>(idx)].first;
        int r = queries[static_cast<std::size_t>(idx)].second;
        while (L > l) { --L; ctr.add(cid[static_cast<std::size_t>(L)], a[static_cast<std::size_t>(L)]); }
        while (R < r) { ++R; ctr.add(cid[static_cast<std::size_t>(R)], a[static_cast<std::size_t>(R)]); }
        while (L < l) { ctr.remove(cid[static_cast<std::size_t>(L)], a[static_cast<std::size_t>(L)]); ++L; }
        while (R > r) { ctr.remove(cid[static_cast<std::size_t>(R)], a[static_cast<std::size_t>(R)]); --R; }
        ans[static_cast<std::size_t>(idx)] = ctr.triple();
    }
    return ans;
}

// ---------------------------------------------------------------- ② 帶修莫隊


struct Mod {
    int p;
    long long nv;
    int nid;
    long long ov;
    int oid;
};


std::vector<Triple> moWithUpdates(const std::vector<long long>& a,
                                  const std::vector<Op>& ops,
                                  const std::map<long long, int>& idof, int m) {
    const int n = static_cast<int>(a.size());
    std::vector<long long> curVal = a;
    std::vector<int> curId(static_cast<std::size_t>(n));
    for (int i = 0; i < n; ++i)
        curId[static_cast<std::size_t>(i)] = idof.at(a[static_cast<std::size_t>(i)]);

    std::vector<std::array<int, 3>> qs;     // (l0, r0, t)
    std::vector<Mod> mods;
    int t = 0;
    for (const Op& op : ops) {
        long long code = op[0], x = op[1], y = op[2];
        if (code == 1) {
            if (n == 0) {
                qs.push_back({0, -1, t});
            } else {
                long long lo = std::max<long long>(0, std::min<long long>(x - 1, n - 1));
                long long hi = std::max<long long>(0, std::min<long long>(y - 1, n - 1));
                if (lo > hi) std::swap(lo, hi);
                qs.push_back({static_cast<int>(lo), static_cast<int>(hi), t});
            }
        } else if (code == 2) {
            if (x >= 1 && x <= n) {
                int p = static_cast<int>(x - 1);
                mods.push_back({p, y, idof.at(y), curVal[static_cast<std::size_t>(p)],
                                curId[static_cast<std::size_t>(p)]});
                curVal[static_cast<std::size_t>(p)] = y;
                curId[static_cast<std::size_t>(p)] = idof.at(y);
                ++t;
            }
        }
    }
    if (n == 0) return std::vector<Triple>(qs.size(), {0, 0, 0});
    if (qs.empty()) return {};

    // 正規化時模擬過一次修改，這裡復原到初始狀態
    curVal = a;
    for (int i = 0; i < n; ++i)
        curId[static_cast<std::size_t>(i)] = idof.at(a[static_cast<std::size_t>(i)]);

    const int B = blockSizeModified(n);
    std::vector<int> order(qs.size());
    for (std::size_t i = 0; i < qs.size(); ++i) order[i] = static_cast<int>(i);
    std::sort(order.begin(), order.end(), [&](int x, int y) {
        int bx = qs[static_cast<std::size_t>(x)][0] / B;
        int by = qs[static_cast<std::size_t>(y)][0] / B;
        if (bx != by) return bx < by;
        int rx = qs[static_cast<std::size_t>(x)][1] / B;
        int ry = qs[static_cast<std::size_t>(y)][1] / B;
        if (rx != ry) return rx < ry;
        int tx = qs[static_cast<std::size_t>(x)][2];
        int ty = qs[static_cast<std::size_t>(y)][2];
        if (tx != ty) return tx < ty;
        return x < y;
    });

    Counter ctr(m, n);
    int L = 0, R = -1, T = 0;
    std::vector<Triple> res(qs.size(), {0, 0, 0});

    auto changeTime = [&](int target) {
        while (T < target) {
            const Mod& md = mods[static_cast<std::size_t>(T)];
            if (L <= md.p && md.p <= R) {
                ctr.remove(md.oid, md.ov);
                ctr.add(md.nid, md.nv);
            }
            curVal[static_cast<std::size_t>(md.p)] = md.nv;
            curId[static_cast<std::size_t>(md.p)] = md.nid;
            ++T;
        }
        while (T > target) {
            --T;
            const Mod& md = mods[static_cast<std::size_t>(T)];
            if (L <= md.p && md.p <= R) {
                ctr.remove(md.nid, md.nv);
                ctr.add(md.oid, md.ov);
            }
            curVal[static_cast<std::size_t>(md.p)] = md.ov;
            curId[static_cast<std::size_t>(md.p)] = md.oid;
        }
    };

    for (int idx : order) {
        int l = qs[static_cast<std::size_t>(idx)][0];
        int r = qs[static_cast<std::size_t>(idx)][1];
        int tt = qs[static_cast<std::size_t>(idx)][2];
        changeTime(tt);
        while (L > l) { --L; ctr.add(curId[static_cast<std::size_t>(L)], curVal[static_cast<std::size_t>(L)]); }
        while (R < r) { ++R; ctr.add(curId[static_cast<std::size_t>(R)], curVal[static_cast<std::size_t>(R)]); }
        while (L < l) { ctr.remove(curId[static_cast<std::size_t>(L)], curVal[static_cast<std::size_t>(L)]); ++L; }
        while (R > r) { ctr.remove(curId[static_cast<std::size_t>(R)], curVal[static_cast<std::size_t>(R)]); --R; }
        res[static_cast<std::size_t>(idx)] = ctr.triple();
    }
    return res;
}

// ---------------------------------------------------------------- 主流程


std::string run(const std::vector<long long>& a,
                const std::vector<std::pair<int, int>>& queries,
                const std::vector<Op>& ops) {
    std::vector<long long> pool = a;
    for (const Op& op : ops)
        if (op[0] == 2) pool.push_back(op[2]);
    std::vector<long long> uniq = pool;
    std::sort(uniq.begin(), uniq.end());
    uniq.erase(std::unique(uniq.begin(), uniq.end()), uniq.end());
    std::map<long long, int> idof;
    for (std::size_t i = 0; i < uniq.size(); ++i) idof[uniq[i]] = static_cast<int>(i);
    const int m = static_cast<int>(uniq.size());
    std::vector<int> cid(a.size());
    for (std::size_t i = 0; i < a.size(); ++i) cid[i] = idof.at(a[i]);

    std::vector<Triple> ra = moQueries(a, queries, cid, m);
    std::vector<Triple> rb = moWithUpdates(a, ops, idof, m);

    std::ostringstream oss;
    oss << ra.size() << "\n";
    for (const Triple& t : ra) oss << t[0] << " " << t[1] << " " << t[2] << "\n";
    oss << rb.size();
    for (const Triple& t : rb) oss << "\n" << t[0] << " " << t[1] << " " << t[2];
    return oss.str();
}

// ---------------------------------------------------------------- 對拍用暴力


Triple bfTriple(const std::vector<long long>& vals) {
    if (vals.empty()) return {0, 0, 0};
    std::map<long long, long long> cnt;
    long long sum = 0;
    for (long long v : vals) {
        cnt[v] += 1;
        sum += v;
    }
    long long mx = 0;
    for (const auto& kv : cnt) mx = std::max(mx, kv.second);
    return {static_cast<long long>(cnt.size()), sum, mx};
}


std::string bfRun(const std::vector<long long>& a,
                  const std::vector<std::pair<int, int>>& queries,
                  const std::vector<Op>& ops) {
    std::vector<Triple> ra;
    for (const auto& qr : queries) {
        // n = 0 時 a.begin() + 1 是越界的（Python 的 a[0:1] 會安全地給空切片）
        std::vector<long long> seg;
        if (!a.empty()) seg.assign(a.begin() + qr.first, a.begin() + qr.second + 1);
        ra.push_back(bfTriple(seg));
    }

    std::vector<long long> cur = a;
    std::vector<Triple> rb;
    for (const Op& op : ops) {
        long long code = op[0], x = op[1], y = op[2];
        if (code == 1) {
            if (cur.empty()) {
                rb.push_back({0, 0, 0});
            } else {
                long long lo = std::max<long long>(0, std::min<long long>(x - 1, (long long)cur.size() - 1));
                long long hi = std::max<long long>(0, std::min<long long>(y - 1, (long long)cur.size() - 1));
                if (lo > hi) std::swap(lo, hi);
                std::vector<long long> seg(cur.begin() + lo, cur.begin() + hi + 1);
                rb.push_back(bfTriple(seg));
            }
        } else if (code == 2) {
            if (x >= 1 && x <= (long long)cur.size()) cur[static_cast<std::size_t>(x - 1)] = y;
        }
    }

    std::ostringstream oss;
    oss << ra.size() << "\n";
    for (const Triple& t : ra) oss << t[0] << " " << t[1] << " " << t[2] << "\n";
    oss << rb.size();
    for (const Triple& t : rb) oss << "\n" << t[0] << " " << t[1] << " " << t[2];
    return oss.str();
}


void check(const std::vector<long long>& a,
           const std::vector<std::pair<int, int>>& queries,
           const std::vector<Op>& ops) {
    assert(run(a, queries, ops) == bfRun(a, queries, ops));
}

// ---------------------------------------------------------------- 測試


void runTests() {
    assert(run({}, {}, {}) == "0\n0");
    assert(run({1, 2, 3}, {}, {}) == "0\n0");
    assert(run({1, 1, 2}, {{0, 2}}, {}) == "1\n2 4 2\n0");
    assert(run({5}, {{0, 0}}, {{1, 1, 1}}) == "1\n1 5 1\n1\n1 5 1");
    assert(run({1, 2}, {{0, 1}}, {{2, 1, 9}, {1, 1, 2}}) == "1\n2 3 1\n1\n2 11 1");
    assert(run({1, 2, 3, 4}, {{0, 3}},
               {{1, 1, 4}, {2, 2, 9}, {1, 1, 4}, {2, 4, 9}, {1, 1, 4}}) ==
           "1\n4 10 1\n3\n4 10 1\n4 17 1\n3 22 2");

    check({}, {}, {});
    check({1}, {}, {});
    check({1}, {{0, 0}}, {});
    check({1, 2, 3}, {{0, 2}, {1, 1}, {0, 0}}, {});
    check({1, 1, 1, 2, 2}, {{0, 4}, {1, 3}, {2, 2}, {0, 1}}, {});
    check({3, 1, 4, 1, 5, 9, 2, 6}, {{0, 7}, {2, 5}, {0, 0}, {3, 3}, {1, 6}}, {});
    check({1, 2, 3}, {}, {{1, 1, 3}});
    check({1, 2, 3}, {}, {{2, 1, 5}, {1, 1, 3}});
    check({1, 2, 3, 4}, {{0, 3}},
          {{1, 1, 4}, {2, 2, 9}, {1, 1, 4}, {2, 4, 9}, {1, 1, 4}});

    std::mt19937 rng(20261014);

    for (int iter = 0; iter < 260; ++iter) {
        int n = static_cast<int>(rng() % 13);
        std::vector<long long> a(static_cast<std::size_t>(n));
        for (int i = 0; i < n; ++i) a[static_cast<std::size_t>(i)] = static_cast<long long>(rng() % 5);
        int q = static_cast<int>(rng() % 11);
        std::vector<std::pair<int, int>> queries;
        for (int i = 0; i < q; ++i) {
            long long l = static_cast<long long>(rng() % (n + 5)) - 2;
            long long r = static_cast<long long>(rng() % (n + 5)) - 2;
            if (n > 0) {
                l = std::max<long long>(0, std::min<long long>(l, n - 1));
                r = std::max<long long>(0, std::min<long long>(r, n - 1));
            } else {
                l = r = 0;
            }
            if (l > r) std::swap(l, r);
            queries.push_back({static_cast<int>(l), static_cast<int>(r)});
        }
        check(a, queries, {});
    }

    for (int iter = 0; iter < 260; ++iter) {
        int n = static_cast<int>(rng() % 13);
        std::vector<long long> a(static_cast<std::size_t>(n));
        for (int i = 0; i < n; ++i) a[static_cast<std::size_t>(i)] = static_cast<long long>(rng() % 5);
        int m = static_cast<int>(rng() % 13);
        std::vector<Op> ops;
        for (int i = 0; i < m; ++i) {
            if (rng() % 2 == 0)
                ops.push_back({1, static_cast<long long>(rng() % (n + 3)) - 1,
                               static_cast<long long>(rng() % (n + 3)) - 1});
            else
                ops.push_back({2, static_cast<long long>(rng() % (n + 2)),
                               static_cast<long long>(rng() % 5)});
        }
        check(a, {}, ops);
    }

    for (int iter = 0; iter < 160; ++iter) {
        int n = 1 + static_cast<int>(rng() % 14);
        std::vector<long long> a(static_cast<std::size_t>(n));
        for (int i = 0; i < n; ++i) a[static_cast<std::size_t>(i)] = static_cast<long long>(rng() % 7);
        std::vector<std::pair<int, int>> queries;
        int q = static_cast<int>(rng() % 9);
        for (int i = 0; i < q; ++i) {
            int l = static_cast<int>(rng() % n);
            int r = static_cast<int>(rng() % n);
            if (l > r) std::swap(l, r);
            queries.push_back({l, r});
        }
        std::vector<Op> ops;
        int m = static_cast<int>(rng() % 13);
        for (int i = 0; i < m; ++i) {
            if (rng() % 2 == 0)
                ops.push_back({1, static_cast<long long>(rng() % n) + 1,
                               static_cast<long long>(rng() % n) + 1});
            else
                ops.push_back({2, static_cast<long long>(rng() % n) + 1,
                               static_cast<long long>(rng() % 7)});
        }
        check(a, queries, ops);
    }

    for (int n : {1, 2, 3, 4, 8, 9, 16, 27, 64, 100}) {
        std::vector<long long> a(static_cast<std::size_t>(n));
        for (int i = 0; i < n; ++i) a[static_cast<std::size_t>(i)] = static_cast<long long>(rng() % 6);
        std::vector<std::pair<int, int>> queries;
        if (n > 0) queries = {{0, n - 1}, {n / 2, n - 1}, {0, 0}};
        std::vector<Op> ops;
        if (n > 0) ops = {{1, 1, n}, {2, 1, 7}, {1, 1, n}};
        check(a, queries, ops);
    }

    // 中量級效能
    {
        int n = 1500;
        std::vector<long long> a(static_cast<std::size_t>(n));
        for (int i = 0; i < n; ++i) a[static_cast<std::size_t>(i)] = static_cast<long long>(rng() % 301);
        std::vector<std::pair<int, int>> queries;
        for (int i = 0; i < 400; ++i) {
            int l = static_cast<int>(rng() % n);
            int r = static_cast<int>(rng() % n);
            if (l > r) std::swap(l, r);
            queries.push_back({l, r});
        }
        std::vector<Op> ops;
        for (int i = 0; i < 400; ++i) {
            if (rng() % 2 == 0)
                ops.push_back({1, static_cast<long long>(rng() % n) + 1,
                               static_cast<long long>(rng() % n) + 1});
            else
                ops.push_back({2, static_cast<long long>(rng() % n) + 1,
                               static_cast<long long>(rng() % 301)});
        }
        std::string out = run(a, queries, ops);
        assert(out.substr(0, 4) == "400\n");
        for (int k = 0; k < 5; ++k) {
            int i = static_cast<int>(rng() % queries.size());
            std::vector<long long> seg(a.begin() + queries[static_cast<std::size_t>(i)].first,
                                       a.begin() + queries[static_cast<std::size_t>(i)].second + 1);
            Triple t = bfTriple(seg);
            std::ostringstream pat;
            pat << t[0] << " " << t[1] << " " << t[2];
            assert(out.find(pat.str()) != std::string::npos);
        }
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

    auto clampN = [](long long v) { return static_cast<int>(std::max<long long>(0, std::min<long long>(MAXN, v))); };
    auto clampV = [](long long v) { return std::max<long long>(-VMAX, std::min<long long>(VMAX, v)); };

    long long t;
    int n = rd.readInt(t) ? clampN(t) : 0;
    std::vector<long long> a;
    for (int i = 0; i < n; ++i) {
        if (!rd.readInt(t)) break;
        a.push_back(clampV(t));
    }
    n = static_cast<int>(a.size());

    int q = rd.readInt(t) ? clampN(t) : 0;
    std::vector<std::pair<int, int>> queries;
    for (int i = 0; i < q; ++i) {
        long long l, r;
        if (!rd.readInt(l)) break;
        if (!rd.readInt(r)) break;
        if (n == 0) {
            queries.push_back({0, 0});
            continue;
        }
        long long lo = std::max<long long>(0, std::min<long long>(clampV(l) - 1, n - 1));
        long long hi = std::max<long long>(0, std::min<long long>(clampV(r) - 1, n - 1));
        if (lo > hi) std::swap(lo, hi);
        queries.push_back({static_cast<int>(lo), static_cast<int>(hi)});
    }

    int m = rd.readInt(t) ? clampN(t) : 0;
    std::vector<Op> ops;
    for (int i = 0; i < m; ++i) {
        long long o, x, y;
        if (!rd.readInt(o)) break;
        if (!rd.readInt(x)) break;
        if (!rd.readInt(y)) break;
        ops.push_back({std::max<long long>(0, std::min<long long>(2, clampV(o))),
                       clampV(x), clampV(y)});
    }

    std::cout << run(a, queries, ops) << "\n";
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
