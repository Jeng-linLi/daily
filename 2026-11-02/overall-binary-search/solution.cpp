// 整體二分（Overall Binary Search / 離線整體二分）：單點修改 + 區間第 k 小
//
// 題意：
//     給定一個長度為 n 的數列 a，初始值為 a[1..n]。接著依序執行 m 個操作：
//         1 p x   —— 把 a[p] 改成 x
//         2 l r k —— 查詢「當下」a[l..r] 中第 k 小的值（1 <= k <= r-l+1）
//     要求把所有查詢的答案依序輸出。
//     限制：1 <= n, m <= 100000，|x| <= 10^9。
//
// 思路：
//     ### 1. 為什麼不能用線段樹套值域樹
//     「動態第 k 小」最直覺的做法是主席樹（可持久化線段樹），但那是**線上**解法：
//     每個版本要 O(log V) 空間，n + m 個版本就是 O((n+m) log V) 個節點。
//     本題是純離線的，所以可以用更省空間的 **整體二分**。
//
//     ### 2. 離線二分的基本套路
//     若我們能回答「查詢範圍 [l,r] 中，值 <= mid 的元素有幾個」，
//     就能判斷第 k 小的值落在 [L,mid] 還是 [mid+1,R]（落在右邊時把 k 減掉左邊的個數）。
//     問題是：每一層二分都要重新回答一次「有幾個」，
//     若每次都重掃全部操作，複雜度會變成 O(m log V * m)，不可行。
//
//     ### 3. 整體二分的核心：所有查詢共用同一棵二分樹
//     把所有操作（含修改）依時間順序排成一個事件序列，遞迴處理值域 [vl, vr]：
//         mid = (vl + vr) / 2
//         從頭到尾掃描事件序列，同時用**樹狀數組**維護「目前值落在 [vl,mid] 的位置數」：
//             * 修改事件：若新值的秩 <= mid，就 BIT[pos] += delta（並把事件留在左半邊）
//                          否則留在右半邊（不動 BIT）
//             * 查詢事件：cnt = BIT.sum(r) - BIT.sum(l-1)
//                          若 k <= cnt  → 事件留在左半邊
//                          否則          → k -= cnt，事件留在右半邊
//         掃完後把這一層加進 BIT 的量**全部回滾**，再各自遞迴。
//     因為事件序列在遞迴中保持時間順序，BIT 在任何時刻都精確反映
//     「該時刻值落在 [vl,mid] 的位置數」，所以每一層只需掃一次序列，
//     總複雜度降為 O((n+m) log V log n)。
//
//     ### 4. 「修改」怎麼變成可以被 BIT 累加的事件
//     BIT 只能做「加一點 / 減一點」，所以把賦值拆成「刪舊值 + 加新值」：
//         a[p]: old → new   等價於   DELETE(p, old)  然後  INSERT(p, new)
//     並且在最前面補上所有初始值的 INSERT(p, a[p])。
//     這樣每個位置的「當前值是否落在 [vl,mid]」永遠是 0 或 1，
//     INSERT/DELETE 把它維護得對稱，兩邊相消就是 0。
//     注意 DELETE 事件必須帶「舊值的秩」，所以要先掃一遍操作序列求出每個位置的前值。
//
//     ### 5. 為什麼答案不會有平手問題
//     第 k 小指的是**值**（不是位置），值相同也只佔一格，
//     所以答案在兩個語言裡都唯一，不需要任何正規化處理。
//
// 複雜度：
//     值域離散化（排序）      O((n+m) log(n+m))
//     整體二分（每層掃一次序列）
//         共 O(log V) 層，每層 O((n+m) log n)   →  O((n+m) log V log n)
//     空間                    O(n + m + V)
//
// 應用場景（Applications）：
//     - **即時排行榜 / 賽事名次**：分數被更新後查詢「第 k 名」
//     - **資料庫的 ORDER BY ... LIMIT k**：批次匯出前 N 名
//     - **版本控制 / 分散式版本號排序**：找出第 k 舊的版本
//     - **監控系統**：取一段時間視窗內的第 k 大延遲（p99 latency）
//     - **資本市場**：組合持股第 k 大的權重（風險集中度監控）

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <optional>
#include <string>
#include <unordered_map>
#include <vector>

#ifdef _WIN32
#include <fcntl.h>
#include <io.h>
#endif

namespace {

constexpr int NMAX = 100000;
constexpr int MMAX = 100000;
constexpr long long VAL_MAX = 1000000000LL;
constexpr long long PARSE_LIMIT = 1000000000000000000LL;

// 樹狀數組（Fenwick Tree）：單點加、前綴和
class Fenwick {
public:
    explicit Fenwick(int n) : n_(n), tree_(n + 1, 0) {}

    void add(int i, int delta) {
        for (; i <= n_; i += i & -i) tree_[i] += delta;
    }

    int sum(int i) const {
        int s = 0;
        for (; i > 0; i -= i & -i) s += tree_[i];
        return s;
    }

private:
    int n_;
    std::vector<int> tree_;
};

constexpr int E_UPDATE = 0;
constexpr int E_QUERY = 1;

// 事件：
//   kind == E_UPDATE：pos=位置，delta=±1，rank=值的秩
//   kind == E_QUERY ：l, r 為查詢區間，k 為第 k 小，id 為查詢編號
struct Event {
    int kind = 0;
    int pos = 0;
    int delta = 0;
    int rank = 0;
    int l = 0;
    int r = 0;
    int k = 0;
    int id = 0;
};

struct Op {
    bool isUpdate = false;
    int a = 0;   // U: 位置 / Q: l
    int b = 0;   // U: 值（long long 以值域存） / Q: r
    int c = 0;   // Q: k
    long long val = 0;  // U: 新值
};

Fenwick* g_bit = nullptr;
std::vector<long long> g_vals;
std::vector<long long>* g_answers = nullptr;

void solve(int vl, int vr, std::vector<Event> evs) {
    if (evs.empty()) return;
    if (vl == vr) {
        long long value = g_vals[vl];
        for (const Event& e : evs) {
            if (e.kind == E_QUERY) (*g_answers)[e.id] = value;
        }
        return;
    }

    int mid = (vl + vr) >> 1;
    std::vector<Event> left, right;
    left.reserve(evs.size());
    right.reserve(evs.size());
    for (const Event& e : evs) {
        if (e.kind == E_UPDATE) {
            if (e.rank <= mid) {
                left.push_back(e);
                g_bit->add(e.pos, e.delta);
            } else {
                right.push_back(e);
            }
        } else {
            int cnt = g_bit->sum(e.r) - g_bit->sum(e.l - 1);
            if (e.k <= cnt) {
                left.push_back(e);
            } else {
                Event ne = e;
                ne.k = e.k - cnt;
                right.push_back(ne);
            }
        }
    }

    // 回滾這一層對 BIT 的所有影響
    for (const Event& e : left) {
        if (e.kind == E_UPDATE) g_bit->add(e.pos, -e.delta);
    }
    std::vector<Event>().swap(evs);  // 釋放父層事件列表，控制記憶體

    solve(vl, mid, std::move(left));
    solve(mid + 1, vr, std::move(right));
}

std::vector<long long> overallBinarySearch(int n, const std::vector<long long>& initial,
                                           const std::vector<Op>& ops) {
    // ---- 1. 值域離散化：初始值 + 所有修改後的值 ----
    std::vector<long long> allVals = initial;
    for (const Op& op : ops) {
        if (op.isUpdate) allVals.push_back(op.val);
    }
    std::sort(allVals.begin(), allVals.end());
    allVals.erase(std::unique(allVals.begin(), allVals.end()), allVals.end());
    std::unordered_map<long long, int> rankOf;
    rankOf.reserve(allVals.size() * 2 + 1);
    for (std::size_t i = 0; i < allVals.size(); ++i) {
        rankOf[allVals[i]] = static_cast<int>(i);
    }
    int vcount = static_cast<int>(allVals.size());

    // ---- 2. 建立事件序列（保持時間順序）----
    std::vector<Event> events;
    events.reserve(static_cast<std::size_t>(n) + ops.size() * 2 + 1);
    std::vector<int> curRank(n + 1, 0);
    for (int i = 1; i <= n; ++i) {
        int r = rankOf[initial[static_cast<std::size_t>(i - 1)]];
        curRank[i] = r;
        Event e;
        e.kind = E_UPDATE;
        e.pos = i;
        e.delta = 1;
        e.rank = r;
        events.push_back(e);  // INSERT(p, a[p])
    }

    int qid = 0;
    for (const Op& op : ops) {
        if (op.isUpdate) {
            int newR = rankOf[op.val];
            int oldR = curRank[op.a];
            if (oldR != newR) {
                Event d;
                d.kind = E_UPDATE;
                d.pos = op.a;
                d.delta = -1;
                d.rank = oldR;
                events.push_back(d);  // DELETE(p, old)
                Event ins;
                ins.kind = E_UPDATE;
                ins.pos = op.a;
                ins.delta = 1;
                ins.rank = newR;
                events.push_back(ins);  // INSERT(p, new)
                curRank[op.a] = newR;
            }
        } else {
            Event q;
            q.kind = E_QUERY;
            q.l = op.a;
            q.r = op.b;
            q.k = op.c;
            q.id = qid++;
            events.push_back(q);
        }
    }

    std::vector<long long> answers(static_cast<std::size_t>(qid), 0);
    if (qid == 0) return answers;

    Fenwick bit(n);
    g_bit = &bit;
    g_vals = allVals;
    g_answers = &answers;
    solve(0, vcount - 1, std::move(events));
    g_bit = nullptr;
    g_answers = nullptr;
    return answers;
}

// 參考實作：直接維護數列 + 排序取第 k 小（O(n log n) 每次查詢，僅供對拍）
std::vector<long long> bruteForce(int n, const std::vector<long long>& initial,
                                  const std::vector<Op>& ops) {
    std::vector<long long> a = initial;
    std::vector<long long> out;
    for (const Op& op : ops) {
        if (op.isUpdate) {
            a[static_cast<std::size_t>(op.a - 1)] = op.val;
        } else {
            std::vector<long long> seg(a.begin() + (op.a - 1), a.begin() + op.b);
            std::sort(seg.begin(), seg.end());
            out.push_back(seg[static_cast<std::size_t>(op.c - 1)]);
        }
    }
    return out;
}

// 簡易 LCG，兩種語言共用同一組常數以確保亂數序列一致
struct Lcg {
    unsigned long long s;
    explicit Lcg(unsigned long long seed) : s(seed) {}
    unsigned long long next() {
        s = (s * 6364136223846793005ULL + 1442695040888963407ULL) &
            0xFFFFFFFFFFFFFFFFULL;
        return s;
    }
    int range(int lo, int hi) {  // [lo, hi]
        return lo + static_cast<int>(
                          next() % static_cast<unsigned long long>(hi - lo + 1));
    }
};

void fail(const char* msg) {
    std::cerr << msg << "\n";
    std::exit(1);
}

void checkEqual(const std::vector<long long>& got, const std::vector<long long>& want,
                const char* msg) {
    if (got != want) fail(msg);
}

void runTests() {
    // ---- 固定用例 ----
    {
        std::vector<Op> ops;
        Op q;
        q.isUpdate = false;
        q.a = 1;
        q.b = 3;
        q.c = 2;
        ops.push_back(q);
        std::vector<long long> init{3, 1, 2};
        std::vector<long long> want{2};
        checkEqual(overallBinarySearch(3, init, ops), want, "case 1 failed");
    }
    {
        std::vector<Op> ops;
        Op q;
        q.a = 1;
        q.b = 4;
        q.c = 3;
        ops.push_back(q);
        std::vector<long long> init{7, 7, 7, 7};
        std::vector<long long> want{7};
        checkEqual(overallBinarySearch(4, init, ops), want, "case 2 failed");
    }
    {
        std::vector<Op> ops;
        Op u;
        u.isUpdate = true;
        u.a = 1;
        u.val = 10;
        ops.push_back(u);
        Op q;
        q.a = 1;
        q.b = 3;
        q.c = 1;
        ops.push_back(q);
        std::vector<long long> init{1, 2, 3};
        std::vector<long long> want{2};
        checkEqual(overallBinarySearch(3, init, ops), want, "case 3 failed");
    }
    {
        std::vector<Op> ops;
        Op u;
        u.isUpdate = true;
        u.a = 3;
        u.val = -1;
        ops.push_back(u);
        Op q;
        q.a = 1;
        q.b = 3;
        q.c = 1;
        ops.push_back(q);
        std::vector<long long> init{5, 6, 7};
        std::vector<long long> want{-1};
        checkEqual(overallBinarySearch(3, init, ops), want, "case 4 failed");
    }
    {
        std::vector<Op> ops;
        Op u;
        u.isUpdate = true;
        u.a = 1;
        u.val = 100;
        ops.push_back(u);
        Op q;
        q.a = 1;
        q.b = 3;
        q.c = 3;
        ops.push_back(q);
        std::vector<long long> init{5, 6, 7};
        std::vector<long long> want{100};
        checkEqual(overallBinarySearch(3, init, ops), want, "case 5 failed");
    }
    {
        std::vector<Op> ops;
        Op q;
        q.a = 1;
        q.b = 5;
        q.c = 1;
        ops.push_back(q);
        std::vector<long long> init{-3, 0, -7, 4, 0};
        std::vector<long long> want{-7};
        checkEqual(overallBinarySearch(5, init, ops), want, "case 6a failed");

        std::vector<Op> ops2;
        Op q2;
        q2.a = 1;
        q2.b = 5;
        q2.c = 5;
        ops2.push_back(q2);
        std::vector<long long> want2{4};
        checkEqual(overallBinarySearch(5, init, ops2), want2, "case 6b failed");
    }
    {
        std::vector<Op> ops;
        Op q;
        q.a = 2;
        q.b = 2;
        q.c = 1;
        ops.push_back(q);
        std::vector<long long> init{9, 1, 4};
        std::vector<long long> want{1};
        checkEqual(overallBinarySearch(3, init, ops), want, "case 7 failed");
    }
    {
        std::vector<Op> ops;
        for (int v : {8, 2, 5}) {
            Op u;
            u.isUpdate = true;
            u.a = 1;
            u.val = v;
            ops.push_back(u);
        }
        Op q;
        q.a = 1;
        q.b = 3;
        q.c = 2;
        ops.push_back(q);
        std::vector<long long> init{1, 2, 3};
        std::vector<long long> want{3};
        checkEqual(overallBinarySearch(3, init, ops), want, "case 8 failed");
    }
    {
        std::vector<Op> ops;
        Op u;
        u.isUpdate = true;
        u.a = 2;
        u.val = 2;
        ops.push_back(u);
        Op q;
        q.a = 1;
        q.b = 3;
        q.c = 1;
        ops.push_back(q);
        std::vector<long long> init{1, 2, 3};
        std::vector<long long> want{1};
        checkEqual(overallBinarySearch(3, init, ops), want, "case 9 failed");
    }
    {
        std::vector<Op> ops;
        Op u1;
        u1.isUpdate = true;
        u1.a = 1;
        u1.val = 5;
        ops.push_back(u1);
        Op u2;
        u2.isUpdate = true;
        u2.a = 2;
        u2.val = 6;
        ops.push_back(u2);
        Op q;
        q.a = 1;
        q.b = 4;
        q.c = 2;
        ops.push_back(q);
        std::vector<long long> init{10, 20, 30, 40};
        std::vector<long long> want{6};
        checkEqual(overallBinarySearch(4, init, ops), want, "case 10 failed");
    }
    {
        std::vector<Op> ops;
        Op q1;
        q1.a = 1;
        q1.b = 3;
        q1.c = 3;
        ops.push_back(q1);
        Op u;
        u.isUpdate = true;
        u.a = 1;
        u.val = 100;
        ops.push_back(u);
        Op q2;
        q2.a = 1;
        q2.b = 3;
        q2.c = 3;
        ops.push_back(q2);
        std::vector<long long> init{1, 2, 3};
        std::vector<long long> want{3, 100};
        checkEqual(overallBinarySearch(3, init, ops), want, "case 11 failed");
    }
    {
        std::vector<Op> ops;
        Op q;
        q.a = 1;
        q.b = 1;
        q.c = 1;
        ops.push_back(q);
        std::vector<long long> init{42};
        std::vector<long long> want{42};
        checkEqual(overallBinarySearch(1, init, ops), want, "case 12 failed");
    }

    // ---- 與參考實作對拍 ----
    {
        Lcg rng(20261102ULL);
        for (int iter = 0; iter < 300; ++iter) {
            int n = rng.range(1, 9);
            std::vector<long long> init;
            for (int i = 0; i < n; ++i) init.push_back(rng.range(-12, 12));
            std::vector<Op> ops;
            int m = rng.range(0, 12);
            for (int i = 0; i < m; ++i) {
                if (rng.range(0, 99) < 45) {
                    Op u;
                    u.isUpdate = true;
                    u.a = rng.range(1, n);
                    u.val = rng.range(-12, 12);
                    ops.push_back(u);
                } else {
                    Op q;
                    q.a = rng.range(1, n);
                    q.b = rng.range(q.a, n);
                    q.c = rng.range(1, q.b - q.a + 1);
                    ops.push_back(q);
                }
            }
            checkEqual(overallBinarySearch(n, init, ops), bruteForce(n, init, ops),
                       "random mismatch");
        }
    }
    {
        Lcg rng(55667788ULL);
        for (int iter = 0; iter < 120; ++iter) {
            int n = rng.range(1, 7);
            std::vector<long long> init;
            for (int i = 0; i < n; ++i) init.push_back(rng.range(-1000, 1000));
            std::vector<Op> ops;
            int m = rng.range(1, 16);
            for (int i = 0; i < m; ++i) {
                if (rng.range(0, 99) < 60) {
                    Op u;
                    u.isUpdate = true;
                    u.a = rng.range(1, n);
                    u.val = rng.range(-1000, 1000);
                    ops.push_back(u);
                } else {
                    Op q;
                    q.a = rng.range(1, n);
                    q.b = rng.range(q.a, n);
                    q.c = rng.range(1, q.b - q.a + 1);
                    ops.push_back(q);
                }
            }
            checkEqual(overallBinarySearch(n, init, ops), bruteForce(n, init, ops),
                       "mid random mismatch");
        }
    }
    {
        Lcg rng(13572468ULL);
        for (int iter = 0; iter < 60; ++iter) {
            int n = rng.range(1, 20);
            std::vector<long long> init;
            for (int i = 0; i < n; ++i) init.push_back(rng.range(-50, 50));
            std::vector<Op> ops;
            int m = rng.range(1, 20);
            for (int i = 0; i < m; ++i) {
                Op q;
                q.a = rng.range(1, n);
                q.b = rng.range(q.a, n);
                q.c = rng.range(1, q.b - q.a + 1);
                ops.push_back(q);
            }
            checkEqual(overallBinarySearch(n, init, ops), bruteForce(n, init, ops),
                       "static random mismatch");
        }
    }
}

// ---------------------------------------------------------------- IO 模式

std::optional<long long> parseInt(const std::string& tok) {
    if (tok.empty()) return std::nullopt;
    std::size_t i = 0;
    if (tok[0] == '+' || tok[0] == '-') i = 1;
    if (i == tok.size()) return std::nullopt;
    for (std::size_t j = i; j < tok.size(); ++j) {
        if (tok[j] < '0' || tok[j] > '9') return std::nullopt;
    }
    __int128 val = 0;
    bool neg = (tok[0] == '-');
    for (std::size_t j = i; j < tok.size(); ++j) {
        val = val * 10 + (tok[j] - '0');
        if (val > (static_cast<__int128>(1) << 100)) return std::nullopt;
    }
    if (neg) val = -val;
    if (val > static_cast<__int128>(PARSE_LIMIT)) return std::nullopt;
    if (val < -static_cast<__int128>(PARSE_LIMIT)) return std::nullopt;
    return static_cast<long long>(val);
}

void runIo(const std::string& raw) {
    std::vector<std::string> toks;
    {
        std::size_t i = 0;
        while (i < raw.size()) {
            while (i < raw.size() && static_cast<unsigned char>(raw[i]) <= ' ') ++i;
            if (i >= raw.size()) break;
            std::size_t j = i;
            while (j < raw.size() && static_cast<unsigned char>(raw[j]) > ' ') ++j;
            toks.push_back(raw.substr(i, j - i));
            i = j;
        }
    }
    std::size_t pos = 0;
    auto nxt = [&toks, &pos]() -> const std::string* {
        if (pos >= toks.size()) return nullptr;
        return &toks[pos++];
    };
    auto bad = []() { std::cout << "error\n"; };

    const std::string* pt = nxt();
    std::optional<long long> nOpt = std::nullopt;
    if (pt) nOpt = parseInt(*pt);
    if (!nOpt || *nOpt < 1 || *nOpt > NMAX) {
        bad();
        return;
    }
    const std::string* pt2 = nxt();
    std::optional<long long> mOpt = std::nullopt;
    if (pt2) mOpt = parseInt(*pt2);
    if (!mOpt || *mOpt < 0 || *mOpt > MMAX) {
        bad();
        return;
    }
    int n = static_cast<int>(*nOpt);
    int m = static_cast<int>(*mOpt);

    std::vector<long long> initial;
    initial.reserve(n);
    for (int i = 0; i < n; ++i) {
        const std::string* t = nxt();
        std::optional<long long> v = std::nullopt;
        if (t) v = parseInt(*t);
        if (!v || *v < -VAL_MAX || *v > VAL_MAX) {
            bad();
            return;
        }
        initial.push_back(*v);
    }

    std::vector<Op> ops;
    ops.reserve(m);
    for (int i = 0; i < m; ++i) {
        const std::string* kind = nxt();
        if (!kind) {
            bad();
            return;
        }
        if (*kind == "1") {
            const std::string* tp = nxt();
            const std::string* tx = nxt();
            std::optional<long long> p = std::nullopt;
            std::optional<long long> x = std::nullopt;
            if (tp) p = parseInt(*tp);
            if (tx) x = parseInt(*tx);
            if (!p || !x || *p < 1 || *p > n || *x < -VAL_MAX || *x > VAL_MAX) {
                bad();
                return;
            }
            Op u;
            u.isUpdate = true;
            u.a = static_cast<int>(*p);
            u.val = *x;
            ops.push_back(u);
        } else if (*kind == "2") {
            const std::string* tl = nxt();
            const std::string* tr = nxt();
            const std::string* tk = nxt();
            std::optional<long long> l = std::nullopt;
            std::optional<long long> r = std::nullopt;
            std::optional<long long> k = std::nullopt;
            if (tl) l = parseInt(*tl);
            if (tr) r = parseInt(*tr);
            if (tk) k = parseInt(*tk);
            if (!l || !r || !k) {
                bad();
                return;
            }
            if (*l < 1 || *l > *r || *r > n || *k < 1 || *k > (*r - *l + 1)) {
                bad();
                return;
            }
            Op q;
            q.isUpdate = false;
            q.a = static_cast<int>(*l);
            q.b = static_cast<int>(*r);
            q.c = static_cast<int>(*k);
            ops.push_back(q);
        } else {
            bad();
            return;
        }
    }

    std::vector<long long> answers = overallBinarySearch(n, initial, ops);
    for (long long a : answers) std::cout << a << "\n";
}

}  // namespace

int main() {
#ifdef _WIN32
    _setmode(_fileno(stdout), _O_BINARY);
#endif
    std::string raw((std::istreambuf_iterator<char>(std::cin)),
                    std::istreambuf_iterator<char>());
    bool blank = true;
    for (std::size_t i = 0; i < raw.size(); ++i) {
        if (static_cast<unsigned char>(raw[i]) > ' ') {
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
