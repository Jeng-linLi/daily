// Prim 最小生成樹（樸素 O(V^2) + 堆優化 O(E log V)，並與 Kruskal 交叉驗證）
//
// 題意：
//     給定一張無向帶權圖（可能不連通），求 MST 總權重、是否連通、MST 具體邊集，
//     不連通時額外給出連通分量個數與最小生成森林總權重。
//     三種獨立實現互相對拍：prim_naive（鄰接矩陣，稠密圖最優）、
//     prim_heap（鄰接表 + 二元堆懶刪除，稀疏圖最優）、kruskal（並查集，第三方基準）。
//
// 思路：
//     ### 貪心原理（cut property）
//     維護已選點集 S。跨越 S 與 V−S 的邊中，權重最小者必屬於某棵 MST。
//     Prim 每輪取這條最小割邊並把新點併入 S，共 V−1 輪。
//     形式與 Dijkstra 很像，差別在 dist 的意義：Dijkstra 是「到起點的路徑總長」，
//     Prim 是「到當前樹的單邊權」。
//     ### 兩個版本怎麼選
//       樸素版每輪 O(V) 掃描 + O(V) 更新 → 總 O(V^2)，與邊數無關，稠密圖最划算。
//       堆版每條邊最多進堆一次 → O(E log V)，稀疏圖明顯更快。
//     ### 為什麼兩版結果能逐字節一致
//       選點：先比 dist，相同則比點編號（小者優先）；
//       鬆弛：只有嚴格 w < dist[v] 才更新 parent（先來先佔）。
//       堆版以 (w, v) 為鍵彈出並跳過已訪問點，等價於「每次取 (dist, 編號) 最小的未訪問點」。
//     ### 不連通圖
//     Prim 從 0 號點出發只覆蓋一個連通分量；若最終選中的點數 < V 即為不連通。
//     最小生成森林則對每個未訪問的點各跑一次 Prim 並加總。
//
// 輸入格式（stdin，全部以空白分隔）：
//     n m
//     m 行：u v w      （無向邊，0-indexed，w 可為負）
// 輸出格式（stdout）：
//     第 1 行：是否連通（1 / 0）
//     第 2 行：prim_naive 的 MST 總權重（不連通為 -1）
//     第 3 行：prim_heap  的 MST 總權重（不連通為 -1）
//     第 4 行：kruskal    的 MST 總權重（不連通為 -1）
//     第 5 行：最小生成森林總權重
//     第 6 行：連通分量個數
//     第 7 行：MST 邊數
//     第 8 行起：每條 MST 邊一行 `u v w`（u < v，按 (u, v) 排序）
// 輸入被截斷時有多少邊讀多少邊。無 stdin 輸入時運行內置斷言測試並輸出 `all tests passed`。

#include <algorithm>
#include <cassert>
#include <iostream>
#include <limits>
#include <optional>
#include <queue>
#include <random>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

using namespace std;

using Edge = tuple<int, int, long long>;     // (u, v, w)
static const long long INF = (long long)1e18;

// ---------------------------------------------------------------- 並查集
class UnionFind {
public:
    vector<int> parent, rankv;
    int count;
    explicit UnionFind(int n) : parent(n), rankv(n, 0), count(n) {
        for (int i = 0; i < n; ++i) parent[(size_t)i] = i;
    }
    int find(int x) {
        while (parent[(size_t)x] != x) {
            parent[(size_t)x] = parent[(size_t)parent[(size_t)x]];
            x = parent[(size_t)x];
        }
        return x;
    }
    bool unionSet(int a, int b) {
        int ra = find(a), rb = find(b);
        if (ra == rb) return false;
        if (rankv[(size_t)ra] < rankv[(size_t)rb]) swap(ra, rb);
        parent[(size_t)rb] = ra;
        if (rankv[(size_t)ra] == rankv[(size_t)rb]) rankv[(size_t)ra]++;
        count--;
        return true;
    }
};

// ---------------------------------------------------------------- 工具
// 過濾自環與越界邊；重邊保留（Prim / Kruskal 都能正確處理）
static vector<Edge> normalizeEdges(int n, const vector<Edge>& edges) {
    vector<Edge> out;
    for (const Edge& e : edges) {
        int u = get<0>(e), v = get<1>(e);
        if (u == v) continue;                          // 自環不可能是生成樹的邊
        if (0 <= u && u < n && 0 <= v && v < n) out.push_back(e);
    }
    return out;
}

// 統一成 u < v 並按 (u, v, w) 排序，讓輸出確定
static vector<Edge> sortMstEdges(const vector<Edge>& mst) {
    vector<Edge> norm;
    for (const Edge& e : mst) {
        int u = get<0>(e), v = get<1>(e);
        norm.emplace_back(min(u, v), max(u, v), get<2>(e));
    }
    sort(norm.begin(), norm.end());
    return norm;
}

// ---------------------------------------------------------------- Prim：樸素 O(V^2)
static tuple<bool, long long, vector<Edge>> primNaive(int n, const vector<Edge>& edges, int root = 0) {
    if (n <= 0) return {true, 0, {}};
    vector<vector<long long>> adj((size_t)n, vector<long long>((size_t)n, INF));
    for (const Edge& e : edges) {
        int u = get<0>(e), v = get<1>(e);
        long long w = get<2>(e);
        if (w < adj[(size_t)u][(size_t)v]) adj[(size_t)u][(size_t)v] = adj[(size_t)v][(size_t)u] = w;
    }

    vector<char> used((size_t)n, 0);
    vector<long long> dist((size_t)n, INF);
    vector<int> parent((size_t)n, -1);
    dist[(size_t)root] = 0;
    long long total = 0;
    vector<Edge> mst;
    int picked = 0;

    for (int round = 0; round < n; ++round) {
        int best = -1;
        for (int v = 0; v < n; ++v) {
            if (used[(size_t)v] || dist[(size_t)v] == INF) continue;
            // 平手取編號小者
            if (best < 0 || dist[(size_t)v] < dist[(size_t)best] ||
                (dist[(size_t)v] == dist[(size_t)best] && v < best))
                best = v;
        }
        if (best < 0) return {false, -1, {}};          // 剩下的點都不可達
        used[(size_t)best] = 1;
        picked++;
        total += dist[(size_t)best];
        if (parent[(size_t)best] >= 0) mst.emplace_back(parent[(size_t)best], best, dist[(size_t)best]);
        for (int v = 0; v < n; ++v) {
            if (!used[(size_t)v] && adj[(size_t)best][(size_t)v] < dist[(size_t)v]) {
                dist[(size_t)v] = adj[(size_t)best][(size_t)v];
                parent[(size_t)v] = best;
            }
        }
    }
    return {picked == n, total, sortMstEdges(mst)};
}

// ---------------------------------------------------------------- Prim：堆優化 O(E log V)
static tuple<bool, long long, vector<Edge>> primHeap(int n, const vector<Edge>& edges, int root = 0) {
    if (n <= 0) return {true, 0, {}};
    vector<vector<pair<int, long long>>> adj((size_t)n);
    for (const Edge& e : edges) {
        int u = get<0>(e), v = get<1>(e);
        long long w = get<2>(e);
        adj[(size_t)u].push_back({v, w});
        adj[(size_t)v].push_back({u, w});
    }

    vector<char> used((size_t)n, 0);
    vector<long long> dist((size_t)n, INF);
    vector<int> parent((size_t)n, -1);
    dist[(size_t)root] = 0;
    // 堆鍵 (w, v)：先比權重，再比點編號
    priority_queue<pair<long long, int>, vector<pair<long long, int>>, greater<pair<long long, int>>> pq;
    pq.push({0, root});
    long long total = 0;
    vector<Edge> mst;
    int picked = 0;

    while (!pq.empty()) {
        auto cur = pq.top();
        pq.pop();
        long long d = cur.first;
        int u = cur.second;
        if (used[(size_t)u]) continue;                 // 過期條目（懶刪除）
        used[(size_t)u] = 1;
        picked++;
        total += d;
        if (parent[(size_t)u] >= 0) mst.emplace_back(parent[(size_t)u], u, d);
        for (const auto& nb : adj[(size_t)u]) {
            int v = nb.first;
            long long w = nb.second;
            if (!used[(size_t)v] && w < dist[(size_t)v]) {
                dist[(size_t)v] = w;
                parent[(size_t)v] = u;
                pq.push({w, v});
            }
        }
    }
    if (picked != n) return {false, -1, {}};
    return {true, total, sortMstEdges(mst)};
}

// ---------------------------------------------------------------- Kruskal（第三方基準）
static tuple<bool, long long, vector<Edge>> kruskal(int n, const vector<Edge>& edges) {
    if (n <= 0) return {true, 0, {}};
    vector<tuple<long long, int, int>> ordered;        // (w, min(u,v), max(u,v))
    for (const Edge& e : edges) {
        int u = get<0>(e), v = get<1>(e);
        ordered.emplace_back(get<2>(e), min(u, v), max(u, v));
    }
    sort(ordered.begin(), ordered.end());
    UnionFind uf(n);
    long long total = 0;
    vector<Edge> mst;
    for (const auto& it : ordered) {
        long long w = get<0>(it);
        int u = get<1>(it), v = get<2>(it);
        if (uf.unionSet(u, v)) {
            total += w;
            mst.emplace_back(u, v, w);
            if ((int)mst.size() == n - 1) break;
        }
    }
    if ((int)mst.size() != n - 1) return {false, -1, {}};
    return {true, total, mst};
}

// ---------------------------------------------------------------- 連通分量 + 最小生成森林
static int connectedComponents(int n, const vector<Edge>& edges) {
    UnionFind uf(n);
    for (const Edge& e : edges) uf.unionSet(get<0>(e), get<1>(e));
    return uf.count;
}

static pair<long long, int> minimumSpanningForest(int n, const vector<Edge>& edges) {
    int comps = connectedComponents(n, edges);
    if (n <= 0) return {0, 0};
    vector<vector<pair<int, long long>>> adj((size_t)n);
    for (const Edge& e : edges) {
        int u = get<0>(e), v = get<1>(e);
        long long w = get<2>(e);
        adj[(size_t)u].push_back({v, w});
        adj[(size_t)v].push_back({u, w});
    }

    vector<char> visited((size_t)n, 0);
    long long total = 0;
    for (int s = 0; s < n; ++s) {
        if (visited[(size_t)s]) continue;
        visited[(size_t)s] = 1;
        priority_queue<pair<long long, int>, vector<pair<long long, int>>, greater<pair<long long, int>>> pq;
        for (const auto& nb : adj[(size_t)s]) pq.push({nb.second, nb.first});
        while (!pq.empty()) {
            auto cur = pq.top();
            pq.pop();
            int u = cur.second;
            if (visited[(size_t)u]) continue;
            visited[(size_t)u] = 1;
            total += cur.first;
            for (const auto& nb : adj[(size_t)u])
                if (!visited[(size_t)nb.first]) pq.push({nb.second, nb.first});
        }
    }
    return {total, comps};
}

// ---------------------------------------------------------------- 小工具
// 只接受 [+-]?digits，與 Python 版的 parse_int 完全一致；其它 token 一律當成「輸入結束」。
static bool tryLL(const string& s, long long& out) {
    if (s.empty()) return false;
    size_t i = 0;
    bool neg = false;
    if (s[0] == '-' || s[0] == '+') {
        neg = (s[0] == '-');
        i = 1;
    }
    if (i >= s.size()) return false;
    long long val = 0;
    for (; i < s.size(); ++i) {
        if (s[i] < '0' || s[i] > '9') return false;
        val = val * 10 + (s[i] - '0');
        if (val > 4000000000000000000LL) return false;      // 溢出保護
    }
    out = neg ? -val : val;
    return true;
}

static mt19937 rngEngine;
static int rndInt(int lo, int hi) { return lo + (int)(rngEngine() % (unsigned)(hi - lo + 1)); }

// ---------------------------------------------------------------- IO 模式
static void runIo(const vector<string>& toks) {
    size_t pos = 0;
    auto nxt = [&]() -> const string* {
        if (pos < toks.size()) return &toks[pos++];
        return nullptr;
    };
    auto nextInt = [&]() -> optional<long long> {
        const string* v = nxt();
        if (!v) return nullopt;
        long long out = 0;
        if (!tryLL(*v, out)) return nullopt;
        return out;
    };
    auto nextIntOr = [&](int defVal) -> int {
        optional<long long> v = nextInt();
        return v ? (int)*v : defVal;
    };

    int n = max(nextIntOr(0), 0);
    int m = max(nextIntOr(0), 0);
    vector<Edge> raw;
    for (int i = 0; i < m; ++i) {
        optional<long long> a = nextInt(), b = nextInt(), c = nextInt();
        if (!a || !b || !c) break;                     // 輸入截斷或出現非數字 token
        raw.emplace_back((int)*a, (int)*b, *c);
    }

    vector<Edge> edges = normalizeEdges(n, raw);
    auto rn = primNaive(n, edges);
    auto rh = primHeap(n, edges);
    auto rk = kruskal(n, edges);
    pair<long long, int> forest = minimumSpanningForest(n, edges);

    cout << ((get<0>(rn) && get<0>(rh) && get<0>(rk)) ? 1 : 0) << '\n';
    cout << get<1>(rn) << '\n';
    cout << get<1>(rh) << '\n';
    cout << get<1>(rk) << '\n';
    cout << forest.first << '\n';
    cout << forest.second << '\n';
    const vector<Edge>& mst = get<2>(rh);
    cout << mst.size() << '\n';
    for (const Edge& e : mst)
        cout << get<0>(e) << ' ' << get<1>(e) << ' ' << get<2>(e) << '\n';
}

// ---------------------------------------------------------------- 測試
static void runTests() {
    // 空圖 / 單點
    assert((primNaive(0, {}) == tuple<bool, long long, vector<Edge>>{true, 0, {}}));
    assert((primHeap(0, {}) == tuple<bool, long long, vector<Edge>>{true, 0, {}}));
    assert((kruskal(0, {}) == tuple<bool, long long, vector<Edge>>{true, 0, {}}));
    assert((primNaive(1, {}) == tuple<bool, long long, vector<Edge>>{true, 0, {}}));
    assert((primHeap(1, {}) == tuple<bool, long long, vector<Edge>>{true, 0, {}}));
    assert((kruskal(1, {}) == tuple<bool, long long, vector<Edge>>{true, 0, {}}));

    // README 示例圖
    {
        vector<Edge> edges = {{0, 1, 2}, {0, 3, 6}, {1, 2, 3}, {1, 3, 8}, {1, 4, 5}, {2, 4, 7}, {3, 4, 9}};
        auto rn = primNaive(5, edges);
        auto rh = primHeap(5, edges);
        auto rk = kruskal(5, edges);
        assert(get<0>(rn) && get<0>(rh) && get<0>(rk));
        assert(get<1>(rn) == 16 && get<1>(rh) == 16 && get<1>(rk) == 16);
        assert(get<2>(rn) == get<2>(rh));
        assert((get<2>(rh) == vector<Edge>{{0, 1, 2}, {0, 3, 6}, {1, 2, 3}, {1, 4, 5}}));
        assert(connectedComponents(5, edges) == 1);
    }

    // 不連通：兩個三角形
    {
        vector<Edge> edges = {{0, 1, 1}, {1, 2, 1}, {0, 2, 5}, {3, 4, 2}, {4, 5, 2}, {3, 5, 9}};
        assert(!get<0>(primNaive(6, edges)));
        assert(get<1>(primNaive(6, edges)) == -1);
        assert(!get<0>(primHeap(6, edges)) && get<1>(primHeap(6, edges)) == -1);
        assert(!get<0>(kruskal(6, edges)) && get<1>(kruskal(6, edges)) == -1);
        pair<long long, int> f = minimumSpanningForest(6, edges);
        assert(f.second == 2);
        assert(f.first == 6);                          // 分量 A: 1+1=2；分量 B: 2+2=4
    }

    // 自環與重邊
    {
        vector<Edge> raw = {{0, 0, -100}, {0, 1, 4}, {0, 1, 2}, {1, 2, 3}, {2, 2, 7}};
        vector<Edge> es = normalizeEdges(3, raw);
        assert(get<1>(primNaive(3, es)) == 5);
        assert(get<1>(primHeap(3, es)) == 5);
        assert(get<1>(kruskal(3, es)) == 5);
        assert((get<2>(primHeap(3, es)) == vector<Edge>{{0, 1, 2}, {1, 2, 3}}));
    }

    // 負權邊：Prim / Kruskal 都照樣成立（最短路才怕負環）
    {
        vector<Edge> edges = {{0, 1, -5}, {1, 2, -7}, {0, 2, 100}, {2, 3, -1}};
        assert(get<1>(primNaive(4, edges)) == -13);
        assert(get<1>(primHeap(4, edges)) == -13);
        assert(get<1>(kruskal(4, edges)) == -13);
        assert(get<2>(primNaive(4, edges)) == get<2>(primHeap(4, edges)));
    }

    // 孤立點
    {
        assert(!get<0>(primNaive(3, {})));
        assert(!get<0>(primHeap(3, {})) && get<1>(primHeap(3, {})) == -1);
        assert(!get<0>(kruskal(3, {})) && get<1>(kruskal(3, {})) == -1);
        pair<long long, int> f = minimumSpanningForest(3, {});
        assert(f.first == 0 && f.second == 3);
    }

    // 鏈狀圖
    {
        vector<Edge> chain;
        for (int i = 0; i < 9; ++i) chain.emplace_back(i, i + 1, i + 1);
        assert(get<1>(primNaive(10, chain)) == 45);
        assert(get<1>(primHeap(10, chain)) == 45);
        assert(get<2>(primNaive(10, chain)) == get<2>(primHeap(10, chain)));
    }

    // 完全圖 K6，邊權 = u+v（大量平手，檢查確定性）
    {
        vector<Edge> full;
        for (int u = 0; u < 6; ++u)
            for (int v = u + 1; v < 6; ++v) full.emplace_back(u, v, u + v);
        auto rn = primNaive(6, full);
        auto rh = primHeap(6, full);
        auto rk = kruskal(6, full);
        assert(get<0>(rn) && get<0>(rh) && get<0>(rk));
        assert(get<1>(rn) == get<1>(rh) && get<1>(rh) == get<1>(rk));
        assert(get<2>(rn) == get<2>(rh));              // 平手規則一致 → 邊集也一致
    }

    // 隨機對拍：三種實現總權重必須一致（MST 權重唯一）
    rngEngine.seed(20261005);
    for (int t = 0; t < 300; ++t) {
        int n = rndInt(1, 9);
        int maxE = n * (n - 1) / 2;
        int m = rndInt(0, maxE);
        vector<pair<int, int>> pool;
        for (int u = 0; u < n; ++u)
            for (int v = u + 1; v < n; ++v) pool.push_back({u, v});
        for (int i = (int)pool.size() - 1; i > 0; --i) swap(pool[(size_t)i], pool[(size_t)rndInt(0, i)]);
        vector<Edge> es;
        for (int i = 0; i < m && i < (int)pool.size(); ++i)
            es.emplace_back(pool[(size_t)i].first, pool[(size_t)i].second, rndInt(-9, 9));

        auto rn = primNaive(n, es);
        auto rh = primHeap(n, es);
        auto rk = kruskal(n, es);
        assert(get<0>(rn) == get<0>(rh) && get<0>(rh) == get<0>(rk));
        assert(get<1>(rn) == get<1>(rh) && get<1>(rh) == get<1>(rk));
        assert(get<2>(rn) == get<2>(rh));
        assert((int)get<2>(rh).size() == (get<0>(rn) ? n - 1 : 0));
        pair<long long, int> f = minimumSpanningForest(n, es);
        if (get<0>(rn)) {
            assert(f.second == 1);
            assert(f.first == get<1>(rn));
        } else {
            assert(f.second > 1);
        }
    }

    // 隨機對拍：最小生成森林 vs 各分量單獨加總
    for (int t = 0; t < 120; ++t) {
        int n = rndInt(1, 10);
        int m = rndInt(0, 12);
        vector<Edge> raw;
        for (int i = 0; i < m; ++i) raw.emplace_back(rndInt(0, n - 1), rndInt(0, n - 1), rndInt(-9, 9));
        vector<Edge> es = normalizeEdges(n, raw);
        pair<long long, int> f = minimumSpanningForest(n, es);
        assert(f.second == connectedComponents(n, es));
        if (f.second == 1) assert(f.first == get<1>(kruskal(n, es)));
    }

    cout << "all tests passed" << '\n';
}

int main() {
    vector<string> toks;
    string t;
    while (cin >> t) toks.push_back(t);
    if (toks.empty()) runTests();
    else runIo(toks);
    return 0;
}
