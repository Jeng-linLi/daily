// 隨機 Treap：旋轉版 / FHQ（split-merge）版 / 隱式 Treap 區間翻轉
//
// 題意：
//     維護一棵以 key 為序的平衡二元搜尋樹，支援插入、刪除、查詢存在性、
//     第 k 小、排名（有多少元素嚴格小於 x）、中序遍歷；
//     另外用**隱式 Treap** 維護一個序列，支援區間翻轉。
//     同時用「旋轉版」與「FHQ 版」兩套實作互相對拍，並驗證堆性質與 BST 性質。
//
// 思路：
//     ### Treap = Binary Search Tree + Heap
//     每個節點有兩個權值：key 滿足 BST 性質（左子樹 < 自己 ≤ 右子樹）；
//     prio 滿足堆性質（本題用小根堆，父節點 prio 小於兩個兒子）。
//     prio 隨機獨立均勻時，「key 決定的 BST 形狀」等價於「按 prio 由小到大插入」，
//     也就是一棵隨機 BST，期望深度 O(log n)──隨機性本身就是平衡。
//
//     ### 旋轉版
//     先按普通 BST 規則遞迴插到底，回程時若「兒子 prio 比自己小」就旋轉
//     （左兒子小 → 右旋；右兒子小 → 左旋）。旋轉不改中序序列，故不破壞 BST 性質。
//     刪除則是把目標節點旋轉到葉子再丟棄：只有一個兒子就讓兒子頂替；
//     兩個都在就把 prio 較小的兒子旋轉上來，再往另一側繼續刪。
//
//     ### FHQ Treap：只用 split 與 merge
//       merge(a,b)：a 的所有 key ≤ b 的所有 key，prio 小者當根，遞迴合併另一側
//       split(t,key)：切成 (≤ key, > key)
//     於是 insert(x) = split → merge(A,new) → merge(...,B)；
//         erase(x)  = split(t,x) → split(A,x-1) → 丟掉中間那一棵的一個節點。
//
//     ### 隱式 Treap：用「下標」當 key
//     序列第 i 個元素放在中序第 i 個位置的節點上，不存 key，只維護子樹大小 sz。
//       split_size(t,k) 切成「前 k 個」與「剩下的」
//       區間翻轉 [l,r] = 切三段 → 給中段打懶標記 → 合併
//     懶標記 rev 用 XOR 下傳（打兩次等於沒打，正好是翻轉語義）。
//
//     ### 為什麼兩個語言能逐位元組一致
//     Treap 形狀完全由 prio 序列決定。本題用 xorshift64 確定性亂數（種子固定），
//     且「隱式 Treap / FHQ / 旋轉版」各自用獨立亂數流，
//     因此 Python 與 C++ 構造出的樹完全同構。
//
// 複雜度（期望）：插入 / 刪除 / 查詢 / 第 k 小 / 排名 O(log n)；
//                 區間翻轉 O(log n)；中序遍歷 O(n)；空間 O(n)。
//
// 輸入（stdin，空白分隔整數）：
//     n                  初始序列長度
//     v_1 ... v_n        初始序列
//     q                  操作數
//     q 組 `op a b`：
//       1 x 0  insert(x)      2 x 0  erase(x)
//       3 x 0  contains(x)    4 k 0  kth(k)（1-indexed，越界 -1）
//       5 x 0  rank(x)        6 l r  隱式 Treap 翻轉 [l,r]（1-indexed 閉區間）
//     非數字 token 視為輸入結束；n、q 截斷到 [0,5000]；數值截斷到 ±10^12。
// 輸出（stdout，5 行）：
//     1: BST 元素個數      2: 中序遍歷      3: op=3/4/5 的查詢結果
//     4: 隱式 Treap 最終序列 5: 1 = 兩版一致且性質成立
// 無 stdin 時跑內置斷言測試並輸出 all tests passed。

#include <cassert>
#include <algorithm>
#include <array>
#include <cstdint>
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

constexpr int MAXN = 5000;
constexpr long long VMAX = 1000000000000LL;
constexpr std::uint64_t MASK64 = ~0ULL;

constexpr std::uint64_t SEED_SEQ = 0x9E3779B97F4A7C15ULL;
constexpr std::uint64_t SEED_FHQ = 0xC2B2AE3D27D4EB4FULL;
constexpr std::uint64_t SEED_ROT = 0x165667B19E3779F9ULL;

// ---------------------------------------------------------------- 確定性亂數


class Rng {
public:
    explicit Rng(std::uint64_t seed) : s_(seed & MASK64) {}

    std::uint64_t nextU64() {
        std::uint64_t x = s_;
        x ^= (x << 13) & MASK64;
        x ^= x >> 7;
        x ^= (x << 17) & MASK64;
        s_ = x & MASK64;
        return s_;
    }

    long long nextPrio() {
        return static_cast<long long>(((nextU64() >> 33) & 0x7FFFFFFFULL)) + 1;
    }

private:
    std::uint64_t s_;
};

// ---------------------------------------------------------------- Treap


class Treap {
public:
    int newNode(long long key, long long prio) {
        int idx = static_cast<int>(key_.size());
        key_.push_back(key);
        prio_.push_back(prio);
        lc_.push_back(-1);
        rc_.push_back(-1);
        sz_.push_back(1);
        rev_.push_back(0);
        return idx;
    }

    int size(int t) const { return t < 0 ? 0 : sz_[static_cast<std::size_t>(t)]; }

    void pull(int t) {
        sz_[static_cast<std::size_t>(t)] =
            1 + size(lc_[static_cast<std::size_t>(t)]) + size(rc_[static_cast<std::size_t>(t)]);
    }

    void push(int t) {
        if (t < 0 || rev_[static_cast<std::size_t>(t)] == 0) return;
        std::swap(lc_[static_cast<std::size_t>(t)], rc_[static_cast<std::size_t>(t)]);
        if (lc_[static_cast<std::size_t>(t)] >= 0)
            rev_[static_cast<std::size_t>(lc_[static_cast<std::size_t>(t)])] ^= 1;
        if (rc_[static_cast<std::size_t>(t)] >= 0)
            rev_[static_cast<std::size_t>(rc_[static_cast<std::size_t>(t)])] ^= 1;
        rev_[static_cast<std::size_t>(t)] = 0;
    }

    // ---- FHQ ----

    int merge(int a, int b) {
        if (a < 0) return b;
        if (b < 0) return a;
        if (prio_[static_cast<std::size_t>(a)] < prio_[static_cast<std::size_t>(b)]) {
            push(a);
            rc_[static_cast<std::size_t>(a)] = merge(rc_[static_cast<std::size_t>(a)], b);
            pull(a);
            return a;
        }
        push(b);
        lc_[static_cast<std::size_t>(b)] = merge(a, lc_[static_cast<std::size_t>(b)]);
        pull(b);
        return b;
    }

    // 切成 (所有 key ≤ key, 所有 key > key)
    std::pair<int, int> splitKey(int t, long long key) {
        if (t < 0) return {-1, -1};
        push(t);
        if (key_[static_cast<std::size_t>(t)] <= key) {
            auto ab = splitKey(rc_[static_cast<std::size_t>(t)], key);
            rc_[static_cast<std::size_t>(t)] = ab.first;
            pull(t);
            return {t, ab.second};
        }
        auto ab = splitKey(lc_[static_cast<std::size_t>(t)], key);
        lc_[static_cast<std::size_t>(t)] = ab.second;
        pull(t);
        return {ab.first, t};
    }

    // 切成 (中序前 k 個, 剩下的)
    std::pair<int, int> splitSize(int t, int k) {
        if (t < 0) return {-1, -1};
        push(t);
        int ls = size(lc_[static_cast<std::size_t>(t)]);
        if (k <= ls) {
            auto ab = splitSize(lc_[static_cast<std::size_t>(t)], k);
            lc_[static_cast<std::size_t>(t)] = ab.second;
            pull(t);
            return {ab.first, t};
        }
        auto ab = splitSize(rc_[static_cast<std::size_t>(t)], k - ls - 1);
        rc_[static_cast<std::size_t>(t)] = ab.first;
        pull(t);
        return {t, ab.second};
    }

    int fhqInsert(int root, long long key, long long prio) {
        auto ab = splitKey(root, key);
        int node = newNode(key, prio);
        return merge(merge(ab.first, node), ab.second);
    }

    int fhqErase(int root, long long key) {
        auto ab = splitKey(root, key);
        auto cd = splitKey(ab.first, key - 1);
        // d 裡面是「恰好等於 key」的那一批，丟掉其中一個（拿它的根開刀）
        int d = cd.second;
        if (d >= 0) d = merge(lc_[static_cast<std::size_t>(d)], rc_[static_cast<std::size_t>(d)]);
        return merge(merge(cd.first, d), ab.second);
    }

    // ---- 旋轉版 ----

    int rotRight(int t) {
        int l = lc_[static_cast<std::size_t>(t)];
        push(t);
        push(l);
        lc_[static_cast<std::size_t>(t)] = rc_[static_cast<std::size_t>(l)];
        rc_[static_cast<std::size_t>(l)] = t;
        pull(t);
        pull(l);
        return l;
    }

    int rotLeft(int t) {
        int r = rc_[static_cast<std::size_t>(t)];
        push(t);
        push(r);
        rc_[static_cast<std::size_t>(t)] = lc_[static_cast<std::size_t>(r)];
        lc_[static_cast<std::size_t>(r)] = t;
        pull(t);
        pull(r);
        return r;
    }

    int rotInsert(int t, long long key, long long prio) {
        if (t < 0) return newNode(key, prio);
        if (key <= key_[static_cast<std::size_t>(t)]) {
            lc_[static_cast<std::size_t>(t)] = rotInsert(lc_[static_cast<std::size_t>(t)], key, prio);
            pull(t);
            if (prio_[static_cast<std::size_t>(lc_[static_cast<std::size_t>(t)])] <
                prio_[static_cast<std::size_t>(t)])
                t = rotRight(t);
        } else {
            rc_[static_cast<std::size_t>(t)] = rotInsert(rc_[static_cast<std::size_t>(t)], key, prio);
            pull(t);
            if (prio_[static_cast<std::size_t>(rc_[static_cast<std::size_t>(t)])] <
                prio_[static_cast<std::size_t>(t)])
                t = rotLeft(t);
        }
        return t;
    }

    int rotErase(int t, long long key) {
        if (t < 0) return -1;
        if (key < key_[static_cast<std::size_t>(t)]) {
            lc_[static_cast<std::size_t>(t)] = rotErase(lc_[static_cast<std::size_t>(t)], key);
            pull(t);
            return t;
        }
        if (key > key_[static_cast<std::size_t>(t)]) {
            rc_[static_cast<std::size_t>(t)] = rotErase(rc_[static_cast<std::size_t>(t)], key);
            pull(t);
            return t;
        }
        if (lc_[static_cast<std::size_t>(t)] < 0) return rc_[static_cast<std::size_t>(t)];
        if (rc_[static_cast<std::size_t>(t)] < 0) return lc_[static_cast<std::size_t>(t)];
        if (prio_[static_cast<std::size_t>(lc_[static_cast<std::size_t>(t)])] <
            prio_[static_cast<std::size_t>(rc_[static_cast<std::size_t>(t)])]) {
            t = rotRight(t);
            rc_[static_cast<std::size_t>(t)] = rotErase(rc_[static_cast<std::size_t>(t)], key);
            pull(t);
            return t;
        }
        t = rotLeft(t);
        lc_[static_cast<std::size_t>(t)] = rotErase(lc_[static_cast<std::size_t>(t)], key);
        pull(t);
        return t;
    }

    // ---- 查詢 ----

    int contains(int root, long long key) {
        int t = root;
        while (t >= 0) {
            push(t);
            if (key_[static_cast<std::size_t>(t)] == key) return 1;
            t = (key < key_[static_cast<std::size_t>(t)]) ? lc_[static_cast<std::size_t>(t)]
                                                          : rc_[static_cast<std::size_t>(t)];
        }
        return 0;
    }

    long long kth(int root, long long k) {
        if (k < 1 || k > static_cast<long long>(size(root))) return -1;
        int t = root;
        while (t >= 0) {
            push(t);
            int ls = size(lc_[static_cast<std::size_t>(t)]);
            if (k <= ls) {
                t = lc_[static_cast<std::size_t>(t)];
            } else if (k == ls + 1) {
                return key_[static_cast<std::size_t>(t)];
            } else {
                k -= ls + 1;
                t = rc_[static_cast<std::size_t>(t)];
            }
        }
        return -1;
    }

    long long rankOf(int root, long long key) {
        long long res = 0;
        int t = root;
        while (t >= 0) {
            push(t);
            if (key_[static_cast<std::size_t>(t)] < key) {
                res += static_cast<long long>(size(lc_[static_cast<std::size_t>(t)])) + 1;
                t = rc_[static_cast<std::size_t>(t)];
            } else {
                t = lc_[static_cast<std::size_t>(t)];
            }
        }
        return res;
    }

    std::vector<long long> toList(int t) {
        std::vector<long long> out;
        std::vector<int> st;
        int cur = t;
        while (cur >= 0 || !st.empty()) {
            while (cur >= 0) {
                st.push_back(cur);
                push(cur);
                cur = lc_[static_cast<std::size_t>(cur)];
            }
            cur = st.back();
            st.pop_back();
            out.push_back(key_[static_cast<std::size_t>(cur)]);
            cur = rc_[static_cast<std::size_t>(cur)];
        }
        return out;
    }

    // ---- 隱式 Treap ----

    int buildSequence(const std::vector<long long>& vals, Rng& rng) {
        int root = -1;
        for (long long v : vals) root = merge(root, newNode(v, rng.nextPrio()));
        return root;
    }

    int reverseRange(int root, int l, int r) {
        int n = size(root);
        if (n == 0 || l > r || l > n) return root;
        if (r > n) r = n;
        if (l < 1) l = 1;
        auto lr = splitSize(root, r);
        auto md = splitSize(lr.first, l - 1);
        if (md.second >= 0) rev_[static_cast<std::size_t>(md.second)] ^= 1;
        return merge(merge(md.first, md.second), lr.second);
    }

    // ---- 性質檢查 ----

    bool checkBst(int t, bool hasLo, long long lo, bool hasHi, long long hi) {
        if (t < 0) return true;
        push(t);
        long long k = key_[static_cast<std::size_t>(t)];
        if (hasLo && k < lo) return false;
        if (hasHi && k > hi) return false;
        if (sz_[static_cast<std::size_t>(t)] !=
            1 + size(lc_[static_cast<std::size_t>(t)]) + size(rc_[static_cast<std::size_t>(t)]))
            return false;
        return checkBst(lc_[static_cast<std::size_t>(t)], hasLo, lo, true, k) &&
               checkBst(rc_[static_cast<std::size_t>(t)], true, k, hasHi, hi);
    }

    bool checkHeap(int t) {
        if (t < 0) return true;
        int l = lc_[static_cast<std::size_t>(t)];
        int r = rc_[static_cast<std::size_t>(t)];
        if (l >= 0 && prio_[static_cast<std::size_t>(l)] < prio_[static_cast<std::size_t>(t)])
            return false;
        if (r >= 0 && prio_[static_cast<std::size_t>(r)] < prio_[static_cast<std::size_t>(t)])
            return false;
        return checkHeap(l) && checkHeap(r);
    }

private:
    std::vector<long long> key_;
    std::vector<long long> prio_;
    std::vector<int> lc_;
    std::vector<int> rc_;
    std::vector<int> sz_;
    std::vector<int> rev_;
};

using Op = std::array<long long, 3>;

// ---------------------------------------------------------------- 主流程


std::string joinLongs(const std::vector<long long>& v) {
    std::ostringstream oss;
    for (std::size_t i = 0; i < v.size(); ++i) {
        if (i) oss << ' ';
        oss << v[i];
    }
    return oss.str();
}


std::string run(const std::vector<long long>& init, const std::vector<Op>& ops) {
    Treap tr;

    Rng rngSeq(SEED_SEQ);
    int seqRoot = tr.buildSequence(init, rngSeq);

    Rng rngFhq(SEED_FHQ);
    Rng rngRot(SEED_ROT);
    int fhqRoot = -1;
    int rotRoot = -1;

    std::vector<long long> queries;

    for (const Op& op : ops) {
        int code = op[0];
        long long a = op[1];
        long long b = op[2];
        if (code == 1) {
            fhqRoot = tr.fhqInsert(fhqRoot, a, rngFhq.nextPrio());
            rotRoot = tr.rotInsert(rotRoot, a, rngRot.nextPrio());
        } else if (code == 2) {
            fhqRoot = tr.fhqErase(fhqRoot, a);
            rotRoot = tr.rotErase(rotRoot, a);
        } else if (code == 3) {
            queries.push_back(tr.contains(fhqRoot, a));
        } else if (code == 4) {
            queries.push_back(tr.kth(fhqRoot, a));
        } else if (code == 5) {
            queries.push_back(tr.rankOf(fhqRoot, a));
        } else if (code == 6) {
            int n = tr.size(seqRoot);
            if (n > 0) {
                long long x = a < 1 ? 1 : (a > n ? n : a);
                long long y = b < 1 ? 1 : (b > n ? n : b);
                if (x > y) std::swap(x, y);
                seqRoot = tr.reverseRange(seqRoot, static_cast<int>(x), static_cast<int>(y));
            }
        }
    }

    std::vector<long long> inorderFhq = tr.toList(fhqRoot);
    std::vector<long long> inorderRot = tr.toList(rotRoot);
    std::vector<long long> seq = tr.toList(seqRoot);

    // 隱式 Treap 的 key 是序列值、不滿足 BST 性質，因此只檢查它的堆性質
    int ok = (inorderFhq == inorderRot && tr.checkBst(fhqRoot, false, 0, false, 0) &&
              tr.checkHeap(fhqRoot) && tr.checkBst(rotRoot, false, 0, false, 0) &&
              tr.checkHeap(rotRoot) && tr.checkHeap(seqRoot))
                 ? 1
                 : 0;

    std::ostringstream oss;
    oss << tr.size(fhqRoot) << "\n";
    oss << joinLongs(inorderFhq) << "\n";
    oss << joinLongs(queries) << "\n";
    oss << joinLongs(seq) << "\n";
    oss << ok;
    return oss.str();
}

// ---------------------------------------------------------------- 對拍用暴力


std::string bfRun(const std::vector<long long>& init, const std::vector<Op>& ops) {
    std::vector<long long> seq = init;
    std::vector<long long> bst;
    std::vector<long long> queries;

    for (const Op& op : ops) {
        int code = op[0];
        long long a = op[1];
        long long b = op[2];
        if (code == 1) {
            auto it = std::lower_bound(bst.begin(), bst.end(), a);
            bst.insert(it, a);
        } else if (code == 2) {
            auto it = std::lower_bound(bst.begin(), bst.end(), a);
            if (it != bst.end() && *it == a) bst.erase(it);
        } else if (code == 3) {
            queries.push_back(std::binary_search(bst.begin(), bst.end(), a) ? 1 : 0);
        } else if (code == 4) {
            queries.push_back((a >= 1 && a <= static_cast<long long>(bst.size()))
                                  ? bst[static_cast<std::size_t>(a - 1)]
                                  : -1);
        } else if (code == 5) {
            auto it = std::lower_bound(bst.begin(), bst.end(), a);
            queries.push_back(static_cast<long long>(it - bst.begin()));
        } else if (code == 6) {
            int n = static_cast<int>(seq.size());
            if (n > 0) {
                long long x = a < 1 ? 1 : (a > n ? n : a);
                long long y = b < 1 ? 1 : (b > n ? n : b);
                if (x > y) std::swap(x, y);
                std::reverse(seq.begin() + (x - 1), seq.begin() + y);
            }
        }
    }

    std::ostringstream oss;
    oss << bst.size() << "\n";
    oss << joinLongs(bst) << "\n";
    oss << joinLongs(queries) << "\n";
    oss << joinLongs(seq) << "\n";
    oss << 1;
    return oss.str();
}


void check(const std::vector<long long>& init, const std::vector<Op>& ops) {
    std::string got = run(init, ops);
    std::string want = bfRun(init, ops);
    assert(got == want);
}

// ---------------------------------------------------------------- 測試


void runTests() {
    assert(run({1, 2, 3, 4, 5}, {}) == "0\n\n\n1 2 3 4 5\n1");
    assert(run({1, 2, 3, 4, 5}, {{6, 2, 4}}) == "0\n\n\n1 4 3 2 5\n1");
    assert(run({5, 4, 3, 2, 1}, {{1, 3, 0}, {1, 1, 0}, {1, 2, 0},
                                 {3, 2, 0}, {4, 2, 0}, {5, 3, 0}}) ==
           "3\n1 2 3\n1 2 2\n5 4 3 2 1\n1");
    assert(run({1, 2, 3}, {{1, 7, 0}, {1, 7, 0}, {2, 7, 0}, {3, 7, 0}, {4, 1, 0}}) ==
           "1\n7\n1 7\n1 2 3\n1");

    check({}, {});
    check({1}, {});
    check({1, 2, 3}, {{1, 2, 0}});
    check({1, 2, 3}, {{6, 1, 3}});
    check({1, 2, 3}, {{6, 3, 1}});
    check({1, 2, 3, 4, 5}, {{6, 2, 4}, {6, 1, 5}});
    check({1, 2, 3, 4, 5}, {{6, 0, 9}});
    check({10, 20, 30}, {{1, 20, 0}, {1, 20, 0}, {1, 10, 0},
                         {3, 20, 0}, {4, 1, 0}, {4, 3, 0}, {5, 20, 0}});
    check({1}, {{1, 5, 0}, {2, 5, 0}, {2, 5, 0}, {3, 5, 0}});
    check({1, 2}, {{4, 0, 0}, {4, 1, 0}});

    std::mt19937 rng(20261014);
    for (int iter = 0; iter < 220; ++iter) {
        int n = static_cast<int>(rng() % 15);
        std::vector<long long> init(static_cast<std::size_t>(n));
        for (int i = 0; i < n; ++i)
            init[static_cast<std::size_t>(i)] = static_cast<long long>(rng() % 19) - 9;
        int q = static_cast<int>(rng() % 17);
        std::vector<Op> ops;
        for (int i = 0; i < q; ++i) {
            int code = static_cast<int>(rng() % 6) + 1;
            int a = static_cast<int>(rng() % 16) - 3;
            int b = static_cast<int>(rng() % 16) - 3;
            ops.push_back({code, a, b});
        }
        check(init, ops);
    }

    for (int iter = 0; iter < 6; ++iter) {
        int n = 200 + static_cast<int>(rng() % 401);
        std::vector<long long> init(static_cast<std::size_t>(n));
        for (int i = 0; i < n; ++i) init[static_cast<std::size_t>(i)] = i;
        std::vector<Op> ops;
        for (int i = 0; i < 400; ++i)
            ops.push_back({static_cast<int>(rng() % 6) + 1,
                           static_cast<int>(rng() % (n + 51)),
                           static_cast<int>(rng() % (n + 51))});
        check(init, ops);
    }

    // 只翻轉的大序列（不跑 O(q·n) 暴力，只驗證可逆性）
    std::vector<long long> big(800);
    for (int i = 0; i < 800; ++i) big[static_cast<std::size_t>(i)] = i + 1;
    std::string o1 = run(big, {{6, 1, 800}, {6, 100, 200}});
    std::string o2 = run(big, {{6, 100, 200}, {6, 1, 800}});
    std::string base = joinLongs(big);
    assert(o1 != "0\n\n\n" + base + "\n1");
    assert(o2 != "0\n\n\n" + base + "\n1");
    std::string twice = run(big, {{6, 1, 800}, {6, 1, 800}});
    assert(twice == "0\n\n\n" + base + "\n1");
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
        // 手寫解析，避免 stoi 丟出例外（非數字 token 一律視為輸入結束）
        const std::string& s = tok_[pos_];
        bool neg = false;
        std::size_t i = 0;
        if (s[0] == '+' || s[0] == '-') {
            neg = (s[0] == '-');
            i = 1;
        }
        long long v = 0;
        bool big = false;
        for (; i < s.size(); ++i) {
            if (v > (VMAX + 8) / 10) { big = true; break; }
            v = v * 10 + (s[i] - '0');
        }
        if (big) v = VMAX;
        if (v > VMAX) v = VMAX;
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

    auto clamp = [](long long v, long long lo, long long hi) {
        return v < lo ? lo : (v > hi ? hi : v);
    };

    long long t;
    int n = rd.readInt(t) ? static_cast<int>(clamp(t, 0, MAXN)) : 0;
    std::vector<long long> init;
    for (int i = 0; i < n; ++i) {
        if (!rd.readInt(t)) break;
        init.push_back(clamp(t, -VMAX, VMAX));
    }

    int q = rd.readInt(t) ? static_cast<int>(clamp(t, 0, MAXN)) : 0;
    std::vector<Op> ops;
    for (int i = 0; i < q; ++i) {
        long long a, b, c;
        if (!rd.readInt(a)) break;
        if (!rd.readInt(b)) break;
        if (!rd.readInt(c)) break;
        ops.push_back({static_cast<int>(clamp(a, 0, 6)),
                       static_cast<int>(clamp(b, -VMAX, VMAX)),
                       static_cast<int>(clamp(c, -VMAX, VMAX))});
    }

    std::string out = run(init, ops);
    std::cout << out << "\n";
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
