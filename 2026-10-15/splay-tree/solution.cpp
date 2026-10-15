// 伸展樹（Splay Tree）：自調整二元搜尋樹
//
// 題意：
//     維護一個可重複元素的有序集合，支援：
//         insert(x)      插入 x
//         erase(x)       刪除一個 x（不存在則無動作）
//         contains(x)    是否存在
//         kth(k)         第 k 小（1-indexed）
//         rank(x)        嚴格小於 x 的元素個數
//         predecessor(x) 嚴格小於 x 的最大值
//         successor(x)   嚴格大於 x 的最小值
//     所有操作都靠**伸展（splay）**把被訪問的節點一路旋轉到根。
//     另外用一個「有序陣列」參考實作（lower_bound）逐項對拍，驗證 splay 的結果正確。
//
// 思路：
//     ### 核心想法：把剛訪問過的節點搬到根
//     伸展樹不維護任何平衡資訊（沒有 AVL 的高度、沒有紅黑樹的顏色），
//     它只遵守一條規則：每次訪問一個節點，就把它旋轉到根。
//     靠均攤分析（access lemma / 勢能法）可證：m 次操作的總代價是
//     O((n + m) log n)，也就是均攤 O(log n)。
//     代價是單次操作可能退化成 O(n)，好處是常數小、程式短，
//     而且自帶「快取效應」：最近被訪問的元素永遠在根附近，
//     對有局部性（locality）的存取序列特別快——這是它真正的價值。
//
//     ### 三種旋轉情況
//     設 x 是要伸展的節點，f = fa[x]，g = fa[f]：
//       - zig（g 為空）：只轉一次 x。
//       - zig-zig（x 與 f 同側，都是左兒子或都是右兒子）：先轉 f 再轉 x。
//       - zig-zag（x 與 f 異側）：轉兩次 x。
//     zig-zig 與 zig-zag 的差別是伸展樹的精髓：
//     zig-zig 先轉父親會把整條路徑「折半」，把樹壓扁；
//     如果這裡也寫成「轉兩次 x」，均攤複雜度就會退化（競爭比變成 Ω(n)）。
//
//     ### 刪除：伸展後靠前驅頂替
//     先把目標 splay 到根：
//       - 計數 > 1：直接減一（本題用 cnt 支援重複元素）；
//       - 沒有兒子 / 只有一個兒子：讓那個兒子當根；
//       - 兩個兒子都在：取左子樹的最大值（前驅），把它在左子樹裡 splay 到根，
//         再把右子樹接成它的右兒子。這跟 BST 的「中序後繼頂替」是同一套思路，
//         差別在於 splay 版順便把前驅送到根，讓下次存取更快。
//
//     ### 與 Treap / 紅黑樹的比較
//       - Treap 靠隨機優先級保證期望平衡，單次操作嚴格 O(log n)；
//       - 紅黑樹靠顏色規則保證最壞 O(log n)，但常數大、程式長；
//       - Splay 不保證單次操作，但均攤 O(log n) 且有局部性加速，
//         是「序列維護」（區間翻轉 / 區間移動）與 LCT 的基礎結構。
//
// 複雜度（均攤）：插入 / 刪除 / 查詢 / 第 k 小 / 排名 O(log n)；
//                 單次操作最壞 O(n)；空間 O(n)。

#include <algorithm>
#include <array>
#include <cassert>
#include <cctype>
#include <fcntl.h>
#include <io.h>
#include <iostream>
#include <optional>
#include <random>
#include <sstream>
#include <string>
#include <vector>

namespace {

const int MAXN = 3000;           // n 與 q 的截斷上限
const long long VMAX = 1000000000000LL;  // 數值截斷上限
const int NULLID = 0;            // 0 號是空節點哨兵，真實節點從 1 開始編號

// ------------------------------------------------------- 伸展樹

class Splay {
public:
    explicit Splay(int capacity) {
        int cap = capacity + 5;
        ch_.assign(cap, {0, 0});
        fa_.assign(cap, 0);
        key_.assign(cap, 0);
        cnt_.assign(cap, 0);
        sz_.assign(cap, 0);
        root_ = 0;
        tot_ = 0;
    }

    int size() const { return sz_[root_]; }
    int root() const { return root_; }

    void insert(long long key) {
        if (root_ == 0) {
            root_ = newNode(key);
            return;
        }
        int x = root_;
        while (true) {
            if (key == key_[x]) {
                ++cnt_[x];
                pushUp(x);
                splay(x);
                return;
            }
            int d = (key > key_[x]) ? 1 : 0;
            if (ch_[x][d] == 0) {
                int i = newNode(key);
                ch_[x][d] = i;
                fa_[i] = x;
                pushUp(x);
                splay(i);
                return;
            }
            x = ch_[x][d];
        }
    }

    int contains(long long key) {
        int x = root_, last = 0;
        while (x) {
            last = x;
            if (key == key_[x]) {
                splay(x);
                return 1;
            }
            x = ch_[x][(key > key_[x]) ? 1 : 0];
        }
        if (last) splay(last);
        return 0;
    }

    void erase(long long key) {
        int x = root_, last = 0;
        while (x) {
            last = x;
            if (key == key_[x]) break;
            x = ch_[x][(key > key_[x]) ? 1 : 0];
        }
        if (x == 0) {
            if (last) splay(last);
            return;
        }
        splay(x);
        if (cnt_[x] > 1) {
            --cnt_[x];
            pushUp(x);
            return;
        }
        int l = ch_[x][0], r = ch_[x][1];
        if (l == 0) {
            root_ = r;
            fa_[r] = 0;
        } else if (r == 0) {
            root_ = l;
            fa_[l] = 0;
        } else {
            // 取左子樹的最大值（前驅），先在左子樹裡把它 splay 到根
            fa_[l] = 0;
            fa_[r] = 0;
            int y = l;
            while (ch_[y][1]) y = ch_[y][1];
            splay(y);
            ch_[y][1] = r;
            fa_[r] = y;
            pushUp(y);
            root_ = y;
        }
    }

    std::optional<long long> kth(int k) {
        if (k < 1 || k > sz_[root_]) return std::nullopt;
        int x = root_, last = 0;
        while (x) {
            last = x;
            int ls = sz_[ch_[x][0]];
            if (k <= ls) {
                x = ch_[x][0];
            } else if (k > ls + cnt_[x]) {
                k -= ls + cnt_[x];
                x = ch_[x][1];
            } else {
                splay(x);
                return key_[x];
            }
        }
        if (last) splay(last);
        return std::nullopt;
    }

    int rank(long long key) {  // 嚴格小於 key 的元素個數
        int x = root_, last = 0;
        int ans = 0;
        while (x) {
            last = x;
            if (key < key_[x]) {
                x = ch_[x][0];
            } else if (key > key_[x]) {
                ans += sz_[ch_[x][0]] + cnt_[x];
                x = ch_[x][1];
            } else {
                ans += sz_[ch_[x][0]];
                splay(x);
                return ans;
            }
        }
        if (last) splay(last);
        return ans;
    }

    std::optional<long long> predecessor(long long key) {
        int x = root_, last = 0;
        std::optional<long long> ans = std::nullopt;
        while (x) {
            last = x;
            if (key_[x] < key) {
                ans = key_[x];
                x = ch_[x][1];
            } else {
                x = ch_[x][0];
            }
        }
        if (last) splay(last);
        return ans;
    }

    std::optional<long long> successor(long long key) {
        int x = root_, last = 0;
        std::optional<long long> ans = std::nullopt;
        while (x) {
            last = x;
            if (key_[x] > key) {
                ans = key_[x];
                x = ch_[x][0];
            } else {
                x = ch_[x][1];
            }
        }
        if (last) splay(last);
        return ans;
    }

    std::vector<long long> inorder() const {
        std::vector<long long> res;
        std::vector<int> st;
        int x = root_;
        while (x || !st.empty()) {
            while (x) {
                st.push_back(x);
                x = ch_[x][0];
            }
            x = st.back();
            st.pop_back();
            for (int i = 0; i < cnt_[x]; ++i) res.push_back(key_[x]);
            x = ch_[x][1];
        }
        return res;
    }

private:
    int newNode(long long key) {
        ++tot_;
        int i = tot_;
        ch_[i][0] = ch_[i][1] = 0;
        fa_[i] = 0;
        key_[i] = key;
        cnt_[i] = 1;
        sz_[i] = 1;
        return i;
    }

    void pushUp(int x) {
        sz_[x] = sz_[ch_[x][0]] + sz_[ch_[x][1]] + cnt_[x];
    }

    int dir(int f, int x) const { return (ch_[f][1] == x) ? 1 : 0; }

    void rotate(int x) {
        int f = fa_[x], g = fa_[f];
        int d = dir(f, x);
        int b = ch_[x][1 - d];
        ch_[f][d] = b;
        if (b) fa_[b] = f;
        ch_[x][1 - d] = f;
        fa_[f] = x;
        fa_[x] = g;
        if (g) {
            ch_[g][(ch_[g][1] == f) ? 1 : 0] = x;
        } else {
            root_ = x;
        }
        pushUp(f);
        pushUp(x);
    }

    void splay(int x) {  // zig / zig-zig / zig-zag
        while (fa_[x]) {
            int f = fa_[x], g = fa_[f];
            if (g == 0) {
                rotate(x);                       // zig
            } else if (dir(g, f) == dir(f, x)) {
                rotate(f);                       // zig-zig：先轉父親
                rotate(x);
            } else {
                rotate(x);                       // zig-zag：轉兩次自己
                rotate(x);
            }
        }
        root_ = x;
    }

    std::vector<std::array<int, 2>> ch_;
    std::vector<int> fa_;
    std::vector<long long> key_;
    std::vector<int> cnt_;
    std::vector<int> sz_;
    int root_;
    int tot_;
};

// ------------------------------------------------------- 參考實作

class RefSet {
public:
    void insert(long long key) {
        a_.insert(std::lower_bound(a_.begin(), a_.end(), key), key);
    }
    void erase(long long key) {
        auto it = std::lower_bound(a_.begin(), a_.end(), key);
        if (it != a_.end() && *it == key) a_.erase(it);
    }
    int contains(long long key) const {
        auto it = std::lower_bound(a_.begin(), a_.end(), key);
        return (it != a_.end() && *it == key) ? 1 : 0;
    }
    std::optional<long long> kth(int k) const {
        if (k < 1 || k > static_cast<int>(a_.size())) return std::nullopt;
        return a_[k - 1];
    }
    int rank(long long key) const {
        return static_cast<int>(std::lower_bound(a_.begin(), a_.end(), key) - a_.begin());
    }
    std::optional<long long> predecessor(long long key) const {
        auto it = std::lower_bound(a_.begin(), a_.end(), key);
        if (it == a_.begin()) return std::nullopt;
        return *(it - 1);
    }
    std::optional<long long> successor(long long key) const {
        auto it = std::upper_bound(a_.begin(), a_.end(), key);
        if (it == a_.end()) return std::nullopt;
        return *it;
    }
    const std::vector<long long>& data() const { return a_; }

private:
    std::vector<long long> a_;
};

// ------------------------------------------------------- 驅動

std::string fmtOpt(const std::optional<long long>& v) {
    return v.has_value() ? std::to_string(*v) : std::string("-1");
}

std::string run(const std::vector<long long>& init,
                const std::vector<std::pair<int, long long>>& ops) {
    Splay tree(static_cast<int>(init.size() + ops.size()) + 10);
    RefSet ref;
    for (long long v : init) {
        tree.insert(v);
        ref.insert(v);
    }

    std::vector<std::string> out;
    int ok = 1;
    for (const auto& op : ops) {
        int o = op.first;
        long long a = op.second;
        if (o == 1) {
            tree.insert(a);
            ref.insert(a);
        } else if (o == 2) {
            tree.erase(a);
            ref.erase(a);
        } else if (o == 3) {
            int v = tree.contains(a);
            out.push_back(std::to_string(v));
            if (v != ref.contains(a)) ok = 0;
        } else if (o == 4) {
            // k 可能超過 int 能表達的範圍：先夾到 [-2e9, 2e9]，保證與 Python 的
            // 「越界回傳 None」語義一致（|k| 遠大於集合大小）。
            int k = static_cast<int>(
                std::max<long long>(-2000000000LL, std::min<long long>(2000000000LL, a)));
            std::optional<long long> v = tree.kth(k);
            out.push_back(fmtOpt(v));
            if (v != ref.kth(k)) ok = 0;
        } else if (o == 5) {
            int v = tree.rank(a);
            out.push_back(std::to_string(v));
            if (v != ref.rank(a)) ok = 0;
        } else if (o == 6) {
            std::optional<long long> v = tree.predecessor(a);
            out.push_back(fmtOpt(v));
            if (v != ref.predecessor(a)) ok = 0;
        } else if (o == 7) {
            std::optional<long long> v = tree.successor(a);
            out.push_back(fmtOpt(v));
            if (v != ref.successor(a)) ok = 0;
        }
    }

    std::vector<long long> order = tree.inorder();
    if (tree.size() != static_cast<int>(ref.data().size()) || order != ref.data()) ok = 0;
    // BST 性質檢查：中序遍歷必須非遞減（同一 key 的重複元素相鄰）
    for (size_t i = 1; i < order.size(); ++i) {
        if (order[i] < order[i - 1]) { ok = 0; break; }
    }

    std::ostringstream oss;
    oss << tree.size() << '\n';
    for (size_t i = 0; i < order.size(); ++i) {
        if (i) oss << ' ';
        oss << order[i];
    }
    oss << '\n';
    for (size_t i = 0; i < out.size(); ++i) {
        if (i) oss << ' ';
        oss << out[i];
    }
    oss << '\n' << ok;
    return oss.str();
}

// ------------------------------------------------------- 內置測試

void runTests() {
    // ---- 固定用例 ----
    assert(run({}, {}) == "0\n\n\n1");
    assert(run({5, 3, 7}, {}) == "3\n3 5 7\n\n1");
    assert(run({}, {{1, 5}, {1, 3}, {1, 7}}) == "3\n3 5 7\n\n1");
    assert(run({1, 2, 3}, {{3, 2}, {3, 9}, {2, 2}, {3, 2}}) == "2\n1 3\n1 0 0\n1");
    assert(run({1, 2, 3}, {{4, 1}, {4, 3}, {4, 4}, {4, 0}}) == "3\n1 2 3\n1 3 -1 -1\n1");
    assert(run({1, 2, 3}, {{5, 1}, {5, 2}, {5, 3}, {5, 0}, {5, 9}}) ==
           "3\n1 2 3\n0 1 2 0 3\n1");
    assert(run({10, 20, 30}, {{6, 25}, {7, 25}, {6, 10}, {7, 30}}) ==
           "3\n10 20 30\n20 30 -1 -1\n1");
    // 重複元素：insert 三次 5 再刪一次，應該剩下兩個 5
    assert(run({}, {{1, 5}, {1, 5}, {1, 5}, {2, 5}}) == "2\n5 5\n\n1");
    assert(run({}, {{1, 5}, {1, 5}, {1, 5}}) == "3\n5 5 5\n\n1");
    assert(run({2, 2, 1}, {{6, 2}, {7, 2}}) == "3\n1 2 2\n1 -1\n1");

    // ---- 隨機對拍 ----
    std::mt19937 rng(20261015);
    for (int iter = 0; iter < 400; ++iter) {
        int n = static_cast<int>(rng() % 41);
        std::vector<long long> init;
        for (int i = 0; i < n; ++i) init.push_back(static_cast<int>(rng() % 41) - 20);
        int q = static_cast<int>(rng() % 81);
        std::vector<std::pair<int, long long>> ops;
        for (int i = 0; i < q; ++i) {
            ops.push_back({static_cast<int>(rng() % 7) + 1,
                           static_cast<int>(rng() % 51) - 25});
        }
        std::string out = run(init, ops);
        assert(out.size() >= 2 && out[out.size() - 1] == '1');
    }

    // ---- 大序列壓測 ----
    Splay big(2 * MAXN + 20);
    std::vector<long long> vals;
    for (int i = 0; i < MAXN; ++i) {
        vals.push_back(static_cast<long long>(rng() % 2000000001ULL) - 1000000000LL);
    }
    for (long long v : vals) big.insert(v);
    assert(big.size() == MAXN);
    std::vector<long long> sorted = vals;
    std::sort(sorted.begin(), sorted.end());
    assert(big.inorder() == sorted);
    for (long long v : vals) assert(big.contains(v) == 1);
    for (long long v : vals) big.erase(v);
    assert(big.size() == 0);
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
        std::string t = toks_[pos_];
        if (!isInt(t)) return false;
        ++pos_;
        // 手動解析：超過 VMAX 就截斷，避免 stoll 拋例外
        bool neg = false;
        size_t i = 0;
        if (t[0] == '+' || t[0] == '-') {
            neg = (t[0] == '-');
            i = 1;
        }
        long long v = 0;
        for (; i < t.size(); ++i) {
            if (v > VMAX / 10) { v = VMAX; break; }  // 再乘 10 一定超過 VMAX
            v = v * 10 + (t[i] - '0');
            if (v > VMAX) { v = VMAX; break; }
        }
        out = neg ? -v : v;
        return true;
    }

private:
    std::vector<std::string> toks_;
    size_t pos_ = 0;
};

long long clampVal(long long v) {
    if (v < -VMAX) return -VMAX;
    if (v > VMAX) return VMAX;
    return v;
}

int clampInt(long long v, long long lo, long long hi) {
    if (v < lo) return static_cast<int>(lo);
    if (v > hi) return static_cast<int>(hi);
    return static_cast<int>(v);
}

void runIo(const std::string& raw) {
    Reader rd(raw);
    long long t = 0;
    int n = rd.readInt(t) ? clampInt(t, 0, MAXN) : 0;
    std::vector<long long> init;
    for (int i = 0; i < n; ++i) {
        if (!rd.readInt(t)) break;
        init.push_back(clampVal(t));
    }
    int q = rd.readInt(t) ? clampInt(t, 0, MAXN) : 0;
    std::vector<std::pair<int, long long>> ops;
    for (int i = 0; i < q; ++i) {
        long long o, a;
        if (!rd.readInt(o)) break;
        if (!rd.readInt(a)) break;
        ops.push_back({clampInt(o, 0, 7), clampVal(a)});
    }
    std::cout << run(init, ops) << "\n";
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
