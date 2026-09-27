// 最小生成树（Kruskal 算法 + 并查集）
// 编译：g++ -std=c++17 -O2 solution.cpp -o solution && ./solution
//
// 思路：贪心 + 并查集。所有边按权值升序排列，依次尝试加入；
//   若两端点当前不连通就选中并合并，否则会成环就丢弃。选够 n-1 条边结束。
//
//   正确性（切分性质）：扫到边 (u, v, w) 且 u、v 未连通时，把 u 所在连通块看成
//   切分的一侧，所有跨切分的边里 (u, v, w) 是当前最小者 —— 更小的边都已考虑过，
//   它们要么不跨这个切分，要么会把点并进来（那样 u、v 就已经连通了）。
//   跨切分的最小边必属于某棵 MST，所以选它不会错。
//
//   排序必须**稳定**（同权边保持输入顺序），这样与 Python 版输出逐字节一致。
//
// 输入（空白分隔）：n m，随后 m 行 `u v w`（顶点 0-based）
// 输出：连通时第一行总权重，随后每行一条选中边 `u v w`；不连通时只输出 disconnected
// 无 stdin 输入时运行内置断言测试。
#include <algorithm>
#include <cassert>
#include <iostream>
#include <limits>
#include <vector>

using namespace std;

struct Edge {
    int u, v;
    long long w;
};

// 并查集：路径压缩 + 按大小合并，同时维护连通块个数
struct DSU {
    vector<int> parent, sz;
    int components;
    explicit DSU(int n) : parent(n), sz(n, 1), components(n) {
        for (int i = 0; i < n; ++i) parent[i] = i;
    }
    int find(int x) {
        while (parent[x] != x) {
            parent[x] = parent[parent[x]];     // 路径压缩（折半）
            x = parent[x];
        }
        return x;
    }
    bool unionSet(int a, int b) {
        int ra = find(a), rb = find(b);
        if (ra == rb) return false;            // 原本就连通，合并无效
        if (sz[ra] < sz[rb]) swap(ra, rb);     // 按大小合并，小树挂到大树下
        parent[rb] = ra;
        sz[ra] += sz[rb];
        --components;
        return true;
    }
};

// 返回 {总权重, 被选中的边, 是否连通}。时间 O(m log m)，空间 O(n + m)
struct KruskalResult {
    long long total;
    vector<Edge> chosen;
    bool connected;
};

KruskalResult kruskal(int n, const vector<Edge>& edges) {
    DSU dsu(n);
    vector<Edge> sorted = edges;
    stable_sort(sorted.begin(), sorted.end(), [](const Edge& a, const Edge& b) {
        return a.w < b.w;                      // 只比权值；同权值靠 stable_sort 保持输入顺序
    });
    KruskalResult res{0, {}, true};
    for (const Edge& e : sorted) {
        if (dsu.unionSet(e.u, e.v)) {
            res.chosen.push_back(e);
            res.total += e.w;
            if (static_cast<int>(res.chosen.size()) == n - 1) break;  // 已是生成树，提前结束
        }
    }
    res.connected = (n <= 1) || (static_cast<int>(res.chosen.size()) == n - 1);
    return res;
}

// 对照用的 Prim（邻接表 + O(n^2) 选最小），用于与 Kruskal 交叉验证总权重
pair<long long, bool> prim(int n, const vector<Edge>& edges) {
    if (n == 0) return {0, true};
    vector<vector<pair<int, long long>>> adj(n);
    for (const Edge& e : edges) {
        if (e.u == e.v) continue;              // 自环对 MST 无意义
        adj[e.u].push_back({e.v, e.w});
        adj[e.v].push_back({e.u, e.w});
    }
    const long long INF = numeric_limits<long long>::max() / 4;
    vector<long long> dist(n, INF);
    vector<bool> used(n, false);
    dist[0] = 0;
    long long total = 0;
    int picked = 0;
    for (int it = 0; it < n; ++it) {
        int best = -1;
        for (int i = 0; i < n; ++i)
            if (!used[i] && (best == -1 || dist[i] < dist[best])) best = i;
        if (best == -1 || dist[best] == INF) return {total, false};  // 够不着，图不连通
        used[best] = true;
        total += dist[best];
        ++picked;
        for (auto& pr : adj[best]) {
            int v = pr.first;
            long long w = pr.second;
            if (!used[v] && w < dist[v]) dist[v] = w;
        }
    }
    return {total, picked == n};
}

// 对照用的指数级枚举：枚举所有边子集，挑权和最小的生成树。仅用于极小规模测试
pair<long long, bool> mstBrute(int n, const vector<Edge>& edges) {
    if (n <= 1) return {0, true};
    int m = static_cast<int>(edges.size());
    const long long INF = numeric_limits<long long>::max() / 4;
    long long best = INF;
    for (int mask = 0; mask < (1 << m); ++mask) {
        if (__builtin_popcount(static_cast<unsigned>(mask)) != n - 1) continue;  // 生成树恰好 n-1 条边
        DSU dsu(n);
        bool ok = true;
        long long total = 0;
        for (int i = 0; i < m; ++i) {
            if (!((mask >> i) & 1)) continue;
            if (dsu.unionSet(edges[i].u, edges[i].v)) {
                total += edges[i].w;
            } else {
                ok = false;                    // 成环，不是树
                break;
            }
        }
        if (ok && dsu.components == 1 && total < best) best = total;
    }
    return {best, best != INF};
}

// 与 Python 版同规模的固定随机序列（LCG），两版各自独立与暴力解对拍
struct LCG {
    unsigned long long s;
    LCG(unsigned long long seed) : s(seed) {}
    int next(int lo, int hi) {
        s = s * 6364136223846793005ULL + 1442695040888963407ULL;
        return lo + static_cast<int>((s >> 33) % static_cast<unsigned long long>(hi - lo + 1));
    }
};

int main() {
    int n, m;
    if (cin >> n >> m) {  // IO 模式
        vector<Edge> edges(m);
        for (int i = 0; i < m; ++i) {
            if (!(cin >> edges[i].u >> edges[i].v >> edges[i].w)) {
                edges[i] = {0, 0, 0};          // 输入被截断时用 0 兜底
            }
        }
        KruskalResult res = kruskal(n, edges);
        if (!res.connected) {
            cout << "disconnected\n";
            return 0;
        }
        cout << res.total << "\n";
        for (const Edge& e : res.chosen) cout << e.u << " " << e.v << " " << e.w << "\n";
        return 0;
    }

    // README 示例：4 点 5 边，MST = (0,1,1) + (1,2,2) + (2,3,3) = 6
    {
        vector<Edge> edges = {{0, 1, 1}, {0, 2, 4}, {1, 2, 2}, {1, 3, 5}, {2, 3, 3}};
        KruskalResult res = kruskal(4, edges);
        assert(res.connected);
        assert(res.total == 6);
        assert(res.chosen.size() == 3);
        assert(res.chosen[0].u == 0 && res.chosen[0].v == 1 && res.chosen[0].w == 1);
        assert(res.chosen[1].u == 1 && res.chosen[1].v == 2 && res.chosen[1].w == 2);
        assert(res.chosen[2].u == 2 && res.chosen[2].v == 3 && res.chosen[2].w == 3);
        assert(prim(4, edges) == make_pair(6LL, true));
        assert(mstBrute(4, edges) == make_pair(6LL, true));
    }

    // 不连通：只有一条边，第三个点孤立
    {
        vector<Edge> edges = {{0, 1, 5}};
        assert(!kruskal(3, edges).connected);
        assert(!prim(3, edges).second);
        assert(!mstBrute(3, edges).second);
    }

    // 退化情形
    {
        vector<Edge> none;
        assert(kruskal(0, none).total == 0);            // 空图视为连通，权重 0
        assert(kruskal(0, none).connected);
        assert(kruskal(1, none).total == 0);            // 单点，不需要边
        assert(kruskal(1, none).connected);
        vector<Edge> one = {{0, 1, 7}};
        assert(kruskal(2, one).total == 7);
        assert(kruskal(2, one).chosen.size() == 1);
        assert(kruskal(2, one).connected);
        assert(!kruskal(2, none).connected);            // 两点无边，不连通
    }

    // 自环与重边：自环必被丢弃，重边只留一条
    {
        vector<Edge> e1 = {{0, 0, 1}, {0, 1, 3}, {0, 1, 3}};
        KruskalResult r1 = kruskal(2, e1);
        assert(r1.total == 3 && r1.chosen.size() == 1);
        vector<Edge> e2 = {{0, 1, 2}, {1, 2, 2}, {0, 2, 2}};
        assert(kruskal(3, e2).total == 4);
    }

    // 负权边同样成立
    {
        vector<Edge> e = {{0, 1, -5}, {1, 2, -1}, {0, 2, 10}};
        KruskalResult r = kruskal(3, e);
        assert(r.total == -6 && r.chosen.size() == 2);
    }

    // 同权边按输入顺序稳定选中（保证两语言输出一致）
    {
        vector<Edge> same = {{2, 3, 1}, {0, 1, 1}, {1, 2, 1}};
        KruskalResult r = kruskal(4, same);
        assert(r.chosen.size() == 3);
        assert(r.chosen[0].u == 2 && r.chosen[0].v == 3);
        assert(r.chosen[1].u == 0 && r.chosen[1].v == 1);
        assert(r.chosen[2].u == 1 && r.chosen[2].v == 2);
    }

    LCG rng(20260927ULL);

    // 随机对拍一：小规模图上 Kruskal、Prim、指数级枚举三者结果一致
    for (int t = 0; t < 150; ++t) {
        int n = rng.next(1, 6);
        int m = rng.next(0, 8);
        vector<pair<int, int>> cand;
        for (int u = 0; u < n; ++u)
            for (int v = u + 1; v < n; ++v) cand.push_back({u, v});
        // 用与 Python 一致的洗牌方式（Fisher-Yates，同样由本 LCG 驱动）
        for (int i = static_cast<int>(cand.size()) - 1; i > 0; --i) {
            int j = rng.next(0, i);
            swap(cand[i], cand[j]);
        }
        vector<Edge> edges;
        for (int i = 0; i < m && i < static_cast<int>(cand.size()); ++i)
            edges.push_back({cand[i].first, cand[i].second, rng.next(1, 20)});

        KruskalResult res = kruskal(n, edges);
        auto pr = prim(n, edges);
        auto br = mstBrute(n, edges);
        assert(res.connected == pr.second);
        assert(res.connected == br.second);
        if (res.connected) {
            assert(res.total == pr.first);
            assert(res.total == br.first);
            assert(static_cast<int>(res.chosen.size()) == n - 1);
            // 选出来的边确实构成生成树：n-1 条且把所有点连成一块
            DSU dsu(n);
            for (const Edge& e : res.chosen) assert(dsu.unionSet(e.u, e.v));
            assert(dsu.components == 1);
        }
    }

    // 随机对拍二：保证连通的随机图（先造链再补随机边），Kruskal 与 Prim 对拍
    for (int t = 0; t < 150; ++t) {
        int n = rng.next(2, 9);
        vector<Edge> edges;
        for (int i = 1; i < n; ++i) {                  // 先连成链，确保一定连通
            int j = rng.next(0, i - 1);
            edges.push_back({j, i, rng.next(1, 30)});
        }
        int extra = rng.next(0, 6);                    // 再补一些随机边
        for (int k = 0; k < extra; ++k) {
            int u = rng.next(0, n - 1), v = rng.next(0, n - 1);
            if (u == v) continue;
            edges.push_back({u, v, rng.next(1, 30)});
        }
        KruskalResult res = kruskal(n, edges);
        auto pr = prim(n, edges);
        assert(res.connected && pr.second);
        assert(res.total == pr.first);
        assert(static_cast<int>(res.chosen.size()) == n - 1);
        DSU dsu(n);
        for (const Edge& e : res.chosen) assert(dsu.unionSet(e.u, e.v));
        assert(dsu.components == 1);
    }

    cout << "all tests passed" << endl;
    return 0;
}
