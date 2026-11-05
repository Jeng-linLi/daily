// 線段樹分治 + 可撤銷並查集：離線動態連通性
//
// 題意：
//     n 個點，T 個時刻（1..T）。m 條邊，第 i 條邊 (u, v) 在時刻區間
//     [l_i, r_i] 內存在。q 次詢問 (t, a, b)：時刻 t 時 a 與 b 是否連通。
//
//     輸入：第一行 n T m q；接下來 m 行每行 u v l r（1-indexed）；
//           接下來 q 行每行 t a b。
//     輸出：每個詢問一行 yes 或 no。
//
// 思路：
//     可撤銷並查集不能路徑壓縮（撤銷會破壞壓縮結果），只能按大小合併。
//     把「時間軸」建成線段樹：每條邊的存活區間 [l, r] 拆成 O(log T) 個
//     節點掛上去，然後 DFS 整顆樹——進入節點時把掛在其上的邊並入並查集
//     （記錄棧深），葉子處恰好得到該時刻的邊集、回答詢問，離開時按棧
//     回滾。每條邊只被並入/撤銷 O(log T) 次。
//
// 複雜度：
//     時間 O(m log T log n + q log n)，空間 O(m log T + q)。
//
// 應用場景（Applications）：
//     - 社交網路好友關係時斷時續，查某時刻兩人是否同圈子
//     - 網路鏈路按計畫時段啟停的故障影響分析
//     - 同框架可換資料結構解「線段樹分治 + 線性基」等問題
//
// 內建測試：固定用例 + 200 組隨機資料與「逐時刻 BFS」暴力對拍。

#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <optional>
#include <string>
#include <vector>

#ifdef _WIN32
#include <fcntl.h>
#include <io.h>
#endif

const int NMAX = 300;
const int TMAX = 300;
const long long VAL_LIM = 1000000000000000000LL;

// ----------------------------------------------------- 可撤銷並查集

struct RollbackDSU {
    std::vector<int> fa, sz;
    std::vector<std::pair<int, int>> stk;  // (被改的子節點, 其原父親)

    explicit RollbackDSU(int n) : fa(n), sz(n, 1) {
        for (int i = 0; i < n; ++i) fa[i] = i;
    }

    int find(int x) const {
        while (fa[x] != x) x = fa[x];
        return x;
    }

    bool unite(int a, int b) {
        int ra = find(a), rb = find(b);
        if (ra == rb) return false;
        if (sz[ra] < sz[rb]) std::swap(ra, rb);
        stk.push_back({rb, fa[rb]});
        fa[rb] = ra;
        sz[ra] += sz[rb];
        return true;
    }

    size_t snapshot() const { return stk.size(); }

    void rollback(size_t mark) {
        while (stk.size() > mark) {
            auto [node, oldFa] = stk.back();
            stk.pop_back();
            sz[fa[node]] -= sz[node];
            fa[node] = oldFa;
        }
    }
};

// ----------------------------------------------------- 核心算法

struct DynamicConnectivity {
    int n, T;
    std::vector<std::vector<std::pair<int, int>>> tree;  // 每個線段樹節點掛的邊
    std::vector<std::vector<std::array<int, 3>>> qAt;    // 每個葉子掛的詢問 (idx, a, b)
    std::vector<std::string> out;

    DynamicConnectivity(int n_, int T_) : n(n_), T(T_) {}

    void addInterval(int node, int nl, int nr, int l, int r, std::pair<int, int> e) {
        if (r < nl || nr < l) return;
        if (l <= nl && nr <= r) {
            tree[node].push_back(e);
            return;
        }
        int mid = (nl + nr) / 2;
        addInterval(2 * node, nl, mid, l, r, e);
        addInterval(2 * node + 1, mid + 1, nr, l, r, e);
    }

    void dfs(int node, int nl, int nr, RollbackDSU& dsu) {
        size_t mark = dsu.snapshot();
        for (auto [u, v] : tree[node]) dsu.unite(u, v);
        if (nl == nr) {
            for (auto [idx, a, b] : qAt[node]) {
                out[idx] = (dsu.find(a) == dsu.find(b)) ? "yes" : "no";
            }
        } else {
            int mid = (nl + nr) / 2;
            dfs(2 * node, nl, mid, dsu);
            dfs(2 * node + 1, mid + 1, nr, dsu);
        }
        dsu.rollback(mark);
    }

    // edges: (u, v, l, r) 0-indexed、時刻閉區間 [l, r]
    // queries: (t, a, b)；回傳每個詢問 "yes"/"no"
    std::vector<std::string> solve(
        const std::vector<std::array<int, 4>>& edges,
        const std::vector<std::array<int, 3>>& queries) {
        size_t size = 1;
        while (size < (size_t)T) size *= 2;
        tree.assign(2 * size, {});
        qAt.assign(2 * size, {});
        for (auto [u, v, l, r] : edges) {
            if (l <= r) addInterval(1, 0, (int)size - 1, l, r, {u, v});
        }
        for (size_t idx = 0; idx < queries.size(); ++idx) {
            auto [t, a, b] = queries[idx];
            // 時刻 t 對應的葉子節點編號是 size + t（堆式編號）
            qAt[size + t].push_back({(int)idx, a, b});
        }
        out.assign(queries.size(), "no");
        RollbackDSU dsu(n);
        if (T > 0) dfs(1, 0, (int)size - 1, dsu);
        return out;
    }
};

// ----------------------------------------------------- 暴力對拍

// 逐時刻重建活躍邊集 + DFS 求連通
std::string bruteQuery(int n, int t, int a, int b,
                       const std::vector<std::array<int, 4>>& edges) {
    std::vector<std::vector<int>> adj(n);
    for (auto [u, v, l, r] : edges) {
        if (l <= t && t <= r) {
            adj[u].push_back(v);
            adj[v].push_back(u);
        }
    }
    std::vector<bool> seen(n, false);
    std::vector<int> stk{a};
    seen[a] = true;
    while (!stk.empty()) {
        int x = stk.back();
        stk.pop_back();
        for (int y : adj[x]) {
            if (!seen[y]) {
                seen[y] = true;
                stk.push_back(y);
            }
        }
    }
    return seen[b] ? "yes" : "no";
}

// ----------------------------------------------------- 內建測試

static uint32_t g_rng = 24681357u;
int rnd(int mod) {
    g_rng = (uint32_t)((g_rng * 1103515245ULL + 12345ULL) % 2147483648ULL);
    return (int)(g_rng % (uint32_t)mod);
}

void fail(const std::string& msg) {
    std::cerr << msg << std::endl;
    std::exit(1);
}

void runTests() {
    // 固定用例 1：邊 (0,1) 全程在，邊 (1,2) 只在時刻 2-3
    {
        DynamicConnectivity dc(3, 4);
        std::vector<std::array<int, 4>> edges = {{0, 1, 0, 3}, {1, 2, 2, 3}};
        std::vector<std::array<int, 3>> qs = {{0, 0, 2}, {2, 0, 2}, {3, 1, 2}, {1, 1, 2}};
        auto got = dc.solve(edges, qs);
        std::vector<std::string> want = {"no", "yes", "yes", "no"};
        if (got != want) fail("case1");
    }

    // 固定用例 2：直連邊 [0,1]，繞路鏈 [2,4] —— 全程連通
    {
        DynamicConnectivity dc(3, 5);
        std::vector<std::array<int, 4>> edges = {{0, 1, 0, 1}, {0, 2, 2, 4}, {2, 1, 2, 4}};
        std::vector<std::array<int, 3>> qs = {{0, 0, 1}, {1, 0, 1}, {2, 0, 1}, {4, 0, 1}, {3, 1, 2}};
        auto got = dc.solve(edges, qs);
        std::vector<std::string> want = {"yes", "yes", "yes", "yes", "yes"};
        if (got != want) fail("case2");
    }

    // 固定用例 3：零詢問、零邊
    {
        DynamicConnectivity dc(2, 1);
        auto got = dc.solve({}, {});
        if (!got.empty()) fail("case3");
    }

    // 隨機對拍：n <= 8, T <= 8
    for (int it = 0; it < 200; ++it) {
        int n = rnd(8) + 1;
        int T = rnd(8) + 1;
        int m = rnd(13);
        std::vector<std::array<int, 4>> edges;
        for (int i = 0; i < m; ++i) {
            int u = rnd(n), v = rnd(n);
            if (u == v) continue;
            int l = rnd(T);
            int r = l + rnd(T - l);
            edges.push_back({u, v, l, r});
        }
        int q = rnd(11);
        std::vector<std::array<int, 3>> queries;
        for (int i = 0; i < q; ++i) queries.push_back({rnd(T), rnd(n), rnd(n)});
        DynamicConnectivity dc(n, T);
        auto got = dc.solve(edges, queries);
        for (size_t i = 0; i < queries.size(); ++i) {
            auto [t, a, b] = queries[i];
            std::string want = bruteQuery(n, t, a, b, edges);
            if (got[i] != want) fail("random mismatch at it=" + std::to_string(it));
        }
    }
}

// ----------------------------------------------------- IO 模式

std::optional<long long> parseLL(const std::string& tok) {
    if (tok.empty()) return std::nullopt;
    size_t i = 0;
    if (tok[0] == '+' || tok[0] == '-') i = 1;
    if (i == tok.size()) return std::nullopt;
    for (size_t k = i; k < tok.size(); ++k) {
        if (tok[k] < '0' || tok[k] > '9') return std::nullopt;
    }
    long long val = 0;
    for (size_t k = i; k < tok.size(); ++k) {
        val = val * 10 + (tok[k] - '0');
        if (val > VAL_LIM) return std::nullopt;
    }
    if (tok[0] == '-') val = -val;
    return val;
}

std::vector<std::string> splitTokens(const std::string& raw) {
    std::vector<std::string> toks;
    std::string cur;
    for (char c : raw) {
        if (c == ' ' || c == '\n' || c == '\r' || c == '\t') {
            if (!cur.empty()) toks.push_back(cur);
            cur.clear();
        } else {
            cur.push_back(c);
        }
    }
    if (!cur.empty()) toks.push_back(cur);
    return toks;
}

void runIO(const std::string& raw) {
    std::vector<std::string> toks = splitTokens(raw);
    size_t pos = 0;
    auto nxt = [&]() -> std::string { return pos < toks.size() ? toks[pos++] : std::string(); };
    auto bad = []() { std::cout << "error\n"; };
    auto getInRange = [&](long long lo, long long hi) -> std::optional<long long> {
        auto v = parseLL(nxt());
        if (!v || *v < lo || *v > hi) return std::nullopt;
        return *v;
    };

    auto nv = getInRange(1, NMAX);
    if (!nv) { bad(); return; }
    int n = (int)*nv;
    auto tv = getInRange(1, TMAX);
    if (!tv) { bad(); return; }
    int T = (int)*tv;
    auto mv = getInRange(0, (long long)n * NMAX);
    if (!mv) { bad(); return; }
    int m = (int)*mv;
    auto qv = getInRange(0, 100000);
    if (!qv) { bad(); return; }
    int q = (int)*qv;

    std::vector<std::array<int, 4>> edges;
    for (int i = 0; i < m; ++i) {
        auto u = getInRange(1, n);
        auto v = getInRange(1, n);
        auto l = getInRange(1, T);
        auto r = getInRange(1, T);
        if (!u || !v || !l || !r || *u == *v) { bad(); return; }
        edges.push_back({(int)*u - 1, (int)*v - 1, (int)*l - 1, (int)*r - 1});
    }
    std::vector<std::array<int, 3>> queries;
    for (int i = 0; i < q; ++i) {
        auto t = getInRange(1, T);
        auto a = getInRange(1, n);
        auto b = getInRange(1, n);
        if (!t || !a || !b) { bad(); return; }
        queries.push_back({(int)*t - 1, (int)*a - 1, (int)*b - 1});
    }
    DynamicConnectivity dc(n, T);
    for (const auto& line : dc.solve(edges, queries)) {
        std::cout << line << "\n";
    }
}

int main() {
#ifdef _WIN32
    _setmode(_fileno(stdout), _O_BINARY);
#endif
    std::string raw, line;
    while (std::getline(std::cin, line)) {
        raw += line;
        raw += '\n';
    }
    bool blank = true;
    for (char c : raw) {
        if (c != ' ' && c != '\n' && c != '\r' && c != '\t') {
            blank = false;
            break;
        }
    }
    if (blank) {
        runTests();
        std::cout << "all tests passed" << std::endl;
    } else {
        runIO(raw);
    }
    return 0;
}
