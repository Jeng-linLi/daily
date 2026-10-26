// 可持久化線段樹（主席樹）：靜態區間第 k 小 / 前驅 / 後繼
//
// 題意：
//     給定一個靜態陣列 a[1..n]，以及 m 個詢問：
//
//         K l r k    a[l..r] 中第 k 小的值（1 <= k <= r-l+1）
//         L l r x    a[l..r] 中小於 x 的最大值（不存在輸出 none）
//         R l r x    a[l..r] 中大於 x 的最小值（不存在輸出 none）
//
//     限制：1 <= n, m <= 30000，|a_i|, |x| <= 10^9。
//
// 思路：
//     ### 1. 值域離散化 + 前綴版本（主席樹）
//     「區間第 k 小」是經典的靜態問題。先把所有值離散化到值域 1..M，
//     對每個前綴 a[1..i] 維護一棵**權值線段樹** T_i：
//     T_i 在離散值 p 上的計數 = 前 i 個元素中等於 p 的個數。
//     T_i 與 T_{i-1} 只差一次單點 +1，用**函數式（持久化）**寫法
//     每次只新建 O(log M) 個節點、其餘與舊版本共用 —— 這就是「主席樹」。
//
//     ### 2. 區間查詢 = 兩個前綴版本相減
//     區間 [l, r] 的權值樹就是 T_r 「減去」 T_{l-1}。查第 k 小時
//     在兩個版本上同步下降：左子樹的計數差
//         cnt = S[left[T_r]] - S[left[T_{l-1}]]
//     即區間裡落在左半值域的元素個數。cnt >= k 往左走；
//     否則往右走並把 k 換成 k - cnt。每次詢問 O(log M)。
//
//     ### 3. 前驅 / 後繼
//     小於 x 的個數 cnt_lt 可用同樣的「版本相減」在值域
//     [1, bisect_left(x)) 上求和得到；cnt_lt > 0 時前驅就是第 cnt_lt 小。
//     後驅同理：數出 <= x 的個數 cnt_le，若 cnt_le < r-l+1，
//     後驅就是第 cnt_le+1 小。三種詢問最後都歸結到同一個 kth。
//
//     ### 4. 空節點 0 的技巧
//     節點池下標 0 當作空節點：S[0] = L[0] = R[0] = 0。
//     版本 0（空前綴）直接就是空節點，update 從空節點長出第一層，
//     不需要顯式建一棵滿零樹，節點數從 O(M) 降回 O(修改次數·log M)。
//
// 複雜度（M = 離散後的值域大小）：
//     建樹          時間 O(n log M)    空間 O(n log M) 個節點
//     單次詢問      O(log M)
//
// 應用場景：
//     - OLAP / 資料倉儲：靜態事實表上反覆出現的「區間分位數」
//       「區間第 k 大」報表，預處理一次，線上 O(log n) 回答
//     - 遊戲排行榜：靜態戰力榜上查某戰力區間的第 k 名、
//       緊鄰的上/下一名（前驅 / 後繼）
//     - 價格系統：靜態歷史報價序列上查「某時段內第 k 便宜/貴的成交價」
//     - 可持久化陣列：把同樣的「路徑複製」技巧用在**下標域**，
//       即得支援歷史版本單點修改 / 查詢的可持久化陣列（思路直接遷移）
//     - 二維推廣：值域再套一維結構（樹套樹 / 平衡樹）即可回答
//       子矩陣第 k 小

#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <iostream>
#include <optional>
#include <random>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

#ifdef _WIN32
#include <fcntl.h>
#include <io.h>
#endif

namespace {

constexpr int NMAX = 30000;
constexpr int MMAX = 30000;
constexpr long long VMAX = 1000000000LL;

void fail(const std::string& msg) {
    std::cerr << "test failed: " << msg << std::endl;
    std::exit(1);
}

// ---------------------------------------------------------------- 主席樹

class PersistentSegTree {
public:
    explicit PersistentSegTree(const std::vector<long long>& arr) {
        values_.push_back(0);  // 離散值域 1-indexed
        for (long long v : arr) {
            values_.push_back(v);
        }
        std::sort(values_.begin() + 1, values_.end());
        values_.erase(std::unique(values_.begin() + 1, values_.end()),
                      values_.end());
        m_ = static_cast<int>(values_.size()) - 1;
        L_.assign(1, 0);
        R_.assign(1, 0);
        S_.assign(1, 0);
        roots_.assign(1, 0);  // roots_[0] = 空節點
        for (long long v : arr) {
            int p = static_cast<int>(
                        std::lower_bound(values_.begin() + 1, values_.end(), v) -
                        values_.begin());
            roots_.push_back(update(roots_.back(), 1, m_, p));
        }
    }

    long long kth(int l, int r, int k) const {
        int p = kthImpl(roots_[l - 1], roots_[r], 1, m_, k);
        return values_[p];
    }

    std::optional<long long> predecessor(int l, int r, long long x) const {
        // 離散值域 [1, p] 內的值都 < x
        int p = static_cast<int>(
                    std::lower_bound(values_.begin() + 1, values_.end(), x) -
                    values_.begin()) -
                1;
        if (p <= 0) {
            return std::nullopt;
        }
        int cnt = countPrefix(roots_[l - 1], roots_[r], 1, m_, p);
        if (cnt == 0) {
            return std::nullopt;
        }
        return values_[kthImpl(roots_[l - 1], roots_[r], 1, m_, cnt)];
    }

    std::optional<long long> successor(int l, int r, long long x) const {
        // 離散值域 [1, p] 內的值都 <= x
        int p = static_cast<int>(
                    std::upper_bound(values_.begin() + 1, values_.end(), x) -
                    values_.begin()) -
                1;
        int total = r - l + 1;
        int cnt = countPrefix(roots_[l - 1], roots_[r], 1, m_, p);
        if (cnt >= total) {
            return std::nullopt;
        }
        return values_[kthImpl(roots_[l - 1], roots_[r], 1, m_, cnt + 1)];
    }

private:
    int newNode(int u) {
        L_.push_back(L_[u]);
        R_.push_back(R_[u]);
        S_.push_back(S_[u]);
        return static_cast<int>(S_.size()) - 1;
    }

    int update(int u, int lo, int hi, int pos) {
        int v = newNode(u);
        S_[v] += 1;
        if (lo != hi) {
            int mid = (lo + hi) / 2;
            if (pos <= mid) {
                L_[v] = update(L_[u], lo, mid, pos);
            } else {
                R_[v] = update(R_[u], mid + 1, hi, pos);
            }
        }
        return v;
    }

    int kthImpl(int u, int v, int lo, int hi, int k) const {
        if (lo == hi) {
            return lo;
        }
        int mid = (lo + hi) / 2;
        int cnt = S_[L_[v]] - S_[L_[u]];
        if (cnt >= k) {
            return kthImpl(L_[u], L_[v], lo, mid, k);
        }
        return kthImpl(R_[u], R_[v], mid + 1, hi, k - cnt);
    }

    int countPrefix(int u, int v, int lo, int hi, int qr) const {
        if (qr < lo) {
            return 0;
        }
        if (hi <= qr) {
            return S_[v] - S_[u];
        }
        int mid = (lo + hi) / 2;
        return countPrefix(L_[u], L_[v], lo, mid, qr) +
               countPrefix(R_[u], R_[v], mid + 1, hi, qr);
    }

    std::vector<long long> values_;
    int m_ = 0;
    std::vector<int> L_, R_, S_;
    std::vector<int> roots_;
};

using Op = std::tuple<std::string, long long, long long, long long>;

std::vector<std::string> solve(int n, const std::vector<long long>& arr,
                               const std::vector<Op>& ops) {
    PersistentSegTree pst(arr);
    std::vector<std::string> out;
    for (const Op& op : ops) {
        const std::string& t = std::get<0>(op);
        long long l = std::get<1>(op);
        long long r = std::get<2>(op);
        long long x = std::get<3>(op);
        if (t == "K") {
            if (!(1 <= l && l <= r && r <= n && 1 <= x && x <= r - l + 1)) {
                out.push_back("error");
                continue;
            }
            out.push_back(std::to_string(pst.kth(static_cast<int>(l),
                                                 static_cast<int>(r),
                                                 static_cast<int>(x))));
        } else if (t == "L" || t == "R") {
            if (!(1 <= l && l <= r && r <= n) || x > VMAX || x < -VMAX) {
                out.push_back("error");
                continue;
            }
            std::optional<long long> v =
                (t == "L") ? pst.predecessor(static_cast<int>(l),
                                             static_cast<int>(r), x)
                           : pst.successor(static_cast<int>(l),
                                           static_cast<int>(r), x);
            out.push_back(v ? std::to_string(*v) : std::string("none"));
        } else {
            out.push_back("error");
        }
    }
    return out;
}

std::vector<std::string> brute(const std::vector<long long>& arr,
                               const std::vector<Op>& ops) {
    std::vector<std::string> out;
    int n = static_cast<int>(arr.size());
    for (const Op& op : ops) {
        const std::string& t = std::get<0>(op);
        long long l = std::get<1>(op);
        long long r = std::get<2>(op);
        long long x = std::get<3>(op);
        if (!(1 <= l && l <= r && r <= n)) {
            out.push_back("error");
            continue;
        }
        std::vector<long long> seg(arr.begin() + (l - 1), arr.begin() + r);
        std::sort(seg.begin(), seg.end());
        if (t == "K") {
            if (!(1 <= x && x <= static_cast<long long>(seg.size()))) {
                out.push_back("error");
                continue;
            }
            out.push_back(std::to_string(seg[x - 1]));
        } else if (t == "L") {
            if (x > VMAX || x < -VMAX) {
                out.push_back("error");
                continue;
            }
            std::optional<long long> best;
            for (long long v : seg) {
                if (v < x && (!best || v > *best)) {
                    best = v;
                }
            }
            out.push_back(best ? std::to_string(*best) : std::string("none"));
        } else if (t == "R") {
            if (x > VMAX || x < -VMAX) {
                out.push_back("error");
                continue;
            }
            std::optional<long long> best;
            for (long long v : seg) {
                if (v > x && (!best || v < *best)) {
                    best = v;
                }
            }
            out.push_back(best ? std::to_string(*best) : std::string("none"));
        } else {
            out.push_back("error");
        }
    }
    return out;
}

// ---------------------------------------------------------------- 測試

void expectEq(const std::vector<std::string>& got,
              const std::vector<std::string>& want, const std::string& tag) {
    if (got != want) {
        std::string msg = tag + "\n got:";
        for (const std::string& s : got) msg += " " + s;
        msg += "\nwant:";
        for (const std::string& s : want) msg += " " + s;
        fail(msg);
    }
}

void runTests() {
    // ---- 固定用例（全部手工驗算過） ----
    expectEq(solve(1, {5},
                   {{"K", 1, 1, 1}, {"L", 1, 1, 4}, {"L", 1, 1, 6},
                    {"R", 1, 1, 5}, {"R", 1, 1, 4}}),
             {"5", "none", "5", "none", "5"}, "fixed single");
    // a = [3,1,3,2,1]，排序後 [1,1,2,3,3]
    expectEq(solve(5, {3, 1, 3, 2, 1},
                   {{"K", 1, 5, 3}, {"L", 1, 5, 3}, {"R", 1, 5, 3},
                    {"K", 2, 4, 2}, {"L", 2, 4, 3}, {"R", 2, 4, 1}}),
             {"2", "2", "none", "2", "2", "2"}, "fixed dup");
    // 負值：[-5,0,5,-10] 排序 [-10,-5,0,5]
    expectEq(solve(4, {-5, 0, 5, -10},
                   {{"K", 1, 4, 1}, {"K", 1, 4, 4}, {"K", 1, 4, 2},
                    {"L", 1, 4, -5}, {"R", 1, 4, -5}, {"L", 1, 4, -10}}),
             {"-10", "5", "-5", "-10", "0", "none"}, "fixed negative");
    // 全部相同
    expectEq(solve(4, {7, 7, 7, 7},
                   {{"K", 2, 3, 1}, {"L", 2, 3, 8}, {"R", 2, 3, 6},
                    {"L", 2, 3, 7}, {"R", 2, 3, 7}}),
             {"7", "7", "7", "none", "none"}, "fixed all-equal");
    // 非法操作
    expectEq(solve(3, {1, 2, 3},
                   {{"K", 2, 1, 1}, {"K", 1, 3, 4}, {"X", 1, 3, 1}}),
             {"error", "error", "error"}, "fixed invalid");

    // ---- 隨機對拍：主席樹 vs 樸素排序 ----
    std::mt19937 rng(20261026);
    for (int tc = 0; tc < 120; ++tc) {
        int n = static_cast<int>(rng() % 14) + 1;
        std::vector<long long> arr(n);
        for (int i = 0; i < n; ++i) {
            arr[i] = static_cast<long long>(rng() % 51) - 25;
        }
        std::vector<Op> ops;
        int q = static_cast<int>(rng() % 30) + 1;
        for (int i = 0; i < q; ++i) {
            const std::string t = (rng() % 3 == 0) ? "K"
                                : (rng() % 2 == 0) ? "L" : "R";
            int l = static_cast<int>(rng() % n) + 1;
            int r = static_cast<int>(rng() % n) + 1;
            if (l > r) std::swap(l, r);
            long long x;
            if (t == "K") {
                x = static_cast<long long>(rng() % (r - l + 1)) + 1;
            } else {
                x = static_cast<long long>(rng() % 61) - 30;
            }
            ops.emplace_back(t, l, r, x);
        }
        expectEq(solve(n, arr, ops), brute(arr, ops), "random small");
    }

    // ---- 稍大規模：確認持久化節點沒有互相污染 ----
    std::mt19937 big(987654321);
    int n = 800;
    std::vector<long long> arr(n);
    for (int i = 0; i < n; ++i) {
        arr[i] = static_cast<long long>(big() % 2000001) - 1000000;
    }
    std::vector<Op> ops;
    for (int i = 0; i < 800; ++i) {
        const std::string t = (big() % 3 == 0) ? "K" : (big() % 2 == 0) ? "L" : "R";
        int l = static_cast<int>(big() % n) + 1;
        int r = static_cast<int>(big() % n) + 1;
        if (l > r) std::swap(l, r);
        long long x;
        if (t == "K") {
            x = static_cast<long long>(big() % (r - l + 1)) + 1;
        } else {
            x = static_cast<long long>(big() % 2000011) - 1000005;
        }
        ops.emplace_back(t, l, r, x);
    }
    expectEq(solve(n, arr, ops), brute(arr, ops), "random large");
}

// ---------------------------------------------------------------- IO 模式

std::optional<long long> parseInt(const std::string& s) {
    // 與 Python 版一致：可選正負號 + 純數字，空串/非數字/越界（>10^18）皆 None
    if (s.empty()) {
        return std::nullopt;
    }
    size_t i = 0;
    bool neg = false;
    if (s[0] == '+' || s[0] == '-') {
        neg = s[0] == '-';
        i = 1;
    }
    if (i == s.size()) {
        return std::nullopt;
    }
    for (size_t j = i; j < s.size(); ++j) {
        if (s[j] < '0' || s[j] > '9') {
            return std::nullopt;
        }
    }
    long long val = 0;
    for (size_t j = i; j < s.size(); ++j) {
        val = val * 10 + (s[j] - '0');
        if (val > 1000000000000000000LL) {
            return std::nullopt;
        }
    }
    return neg ? -val : val;
}

void runIo(const std::string& raw) {
    std::istringstream in(raw);
    std::vector<std::string> toks;
    std::string t;
    while (in >> t) {
        toks.push_back(t);
    }
    size_t pos = 0;
    auto nxt = [&]() -> const std::string* {
        if (pos >= toks.size()) {
            return nullptr;
        }
        return &toks[pos++];
    };

    auto err = []() {
        std::cout << "error\n";
    };

    const std::string* p = nxt();
    std::optional<long long> n = p ? parseInt(*p) : std::nullopt;
    if (!n || !(*n >= 1 && *n <= NMAX)) {
        err();
        return;
    }
    p = nxt();
    std::optional<long long> m = p ? parseInt(*p) : std::nullopt;
    if (!m || !(*m >= 1 && *m <= MMAX)) {
        err();
        return;
    }
    std::vector<long long> arr;
    bool ok = true;
    for (int i = 0; i < *n; ++i) {
        p = nxt();
        std::optional<long long> v = p ? parseInt(*p) : std::nullopt;
        if (!v || *v > VMAX || *v < -VMAX) {
            ok = false;
            break;
        }
        arr.push_back(*v);
    }
    if (!ok) {
        err();
        return;
    }
    std::vector<std::string> outLines;
    std::vector<Op> ops;
    for (int i = 0; i < *m; ++i) {
        // 每個操作一律消費 4 個 token（類型 + 三個整數），
        // 任一環節不合法都在對應行輸出 error，但 token 仍要前進，
        // 這樣後續操作的解析才不會錯位（兩個語言同一條規則）。
        const std::string* pt = nxt();
        p = nxt();
        std::optional<long long> a = p ? parseInt(*p) : std::nullopt;
        p = nxt();
        std::optional<long long> b = p ? parseInt(*p) : std::nullopt;
        p = nxt();
        std::optional<long long> c = p ? parseInt(*p) : std::nullopt;
        bool validOp = pt && (*pt == "K" || *pt == "L" || *pt == "R");
        if (!validOp || !a || !b || !c) {
            outLines.push_back("error");
            continue;
        }
        ops.emplace_back(*pt, *a, *b, *c);
    }
    // solve 內部會逐條校驗 l/r/k/x 的合法性
    std::vector<std::string> res = solve(static_cast<int>(*n), arr, ops);
    for (const std::string& s : res) {
        outLines.push_back(s);
    }
    for (const std::string& s : outLines) {
        std::cout << s << "\n";
    }
}

}  // namespace

int main() {
#ifdef _WIN32
    _setmode(_fileno(stdout), _O_BINARY);
#endif
    std::string raw((std::istreambuf_iterator<char>(std::cin)),
                    std::istreambuf_iterator<char>());
    bool blank = true;
    for (size_t i = 0; i < raw.size(); ++i) {
        unsigned char u = static_cast<unsigned char>(raw[i]);
        if (u > ' ') {
            blank = false;
            break;
        }
    }
    if (!blank) {
        runIo(raw);
    } else {
        runTests();
        std::cout << "all tests passed" << std::endl;
    }
    return 0;
}
