// Bellman-Ford 单源最短路（支持负权边 + 负环检测 + 路径还原）
// 编译：g++ -std=c++17 -O2 solution.cpp -o solution && ./solution
//
// 思路：Dijkstra 依赖「已出队的点距离不再变小」这个贪心性质，有负权边就不成立。
//   Bellman-Ford 的出发点是一个朴素事实：**一条最短路最多经过 n-1 条边**
//   （再多就一定绕了环；正环/零环可以删掉，负环则根本不存在「最短路」）。
//   于是把所有边整体松弛 n-1 轮即可；再多做一轮还能松弛就说明绕了负环 —— 这既是
//   负环检测，也是 Bellman-Ford 相比 Dijkstra 的核心能力。
//   要点：只松弛 dist[u] 有限的边；某一轮没有更新就提前退出；记录 pre 可还原最短路。
//   SPFA（队列优化）只让「上一轮被更新过的点」继续松弛其出边，负环判据改成
//   「某点入队次数 > n」。
//   注意：负环必须**从 s 可达**才会被检测到；图另一头的负环与 s 无关。
//
// 输入（空白分隔）：n m s / 接着 m 行 u v w
// 输出：第 1 行 1=存在 s 可达的负环、0=不存在；
//       第 2 行（无负环时）dist[0]..dist[n-1]，空格分隔，不可达输出 INF
// 无 stdin 输入时运行内置断言测试。
#include <algorithm>
#include <cassert>
#include <deque>
#include <iostream>
#include <vector>

using namespace std;

const long long INF = (1LL << 60);   // 足够大，且 dist[u] + w 不会溢出

struct Edge {
    int u, v;
    long long w;
};

// Bellman-Ford：返回 has_neg_cycle，距离写入 dist，前驱写入 pre。时间 O(n·m)，空间 O(n)
bool bellmanFord(int n, const vector<Edge>& edges, int s, vector<long long>& dist, vector<int>& pre) {
    dist.assign(n, INF);
    pre.assign(n, -1);
    if (0 <= s && s < n) dist[s] = 0;

    bool hasNeg = false;
    for (int it = 0; it < n; ++it) {
        bool changed = false;
        for (const Edge& e : edges) {
            // 只从已可达的点往外松弛，避免 INF + w 污染结果
            if (dist[e.u] != INF && dist[e.u] + e.w < dist[e.v]) {
                dist[e.v] = dist[e.u] + e.w;
                pre[e.v] = e.u;
                changed = true;
            }
        }
        if (!changed) break;              // 这一轮没人被更新，后面也不可能再变
        if (it == n - 1) hasNeg = true;   // 第 n 轮还能松弛 → 绕了负环
    }
    return hasNeg;
}

// Bellman-Ford 的队列优化版（SPFA）。返回 has_neg_cycle，距离写入 dist
bool spfa(int n, const vector<Edge>& edges, int s, vector<long long>& dist) {
    vector<vector<pair<int, long long>>> adj(n);
    for (const Edge& e : edges) adj[e.u].push_back({e.v, e.w});

    dist.assign(n, INF);
    vector<int> cnt(n, 0);
    vector<char> inq(n, 0);
    deque<int> q;
    if (0 <= s && s < n) {
        dist[s] = 0;
        inq[s] = 1;
        cnt[s] = 1;
        q.push_back(s);
    }

    while (!q.empty()) {
        int u = q.front();
        q.pop_front();
        inq[u] = 0;
        if (dist[u] == INF) continue;
        for (auto [v, w] : adj[u]) {
            if (dist[u] + w < dist[v]) {
                dist[v] = dist[u] + w;
                if (!inq[v]) {
                    inq[v] = 1;
                    ++cnt[v];
                    if (cnt[v] > n) return true;   // 入队超过 n 次 → 有负环
                    q.push_back(v);
                }
            }
        }
    }
    return false;
}

// 沿前驱数组还原 s → t 的最短路；不可达返回空。时间 O(路径长度)
vector<int> buildPath(const vector<int>& pre, int s, int t) {
    if (t < 0 || t >= static_cast<int>(pre.size())) return {};
    vector<int> path;
    int cur = t;
    while (cur != -1) {
        path.push_back(cur);
        if (cur == s) {
            reverse(path.begin(), path.end());
            return path;
        }
        cur = pre[cur];
    }
    return {};
}

// 校验路径：返回 (每条边都存在, 路径总权重)
pair<bool, long long> pathWeight(const vector<Edge>& edges, const vector<int>& path) {
    long long total = 0;
    for (int i = 0; i + 1 < static_cast<int>(path.size()); ++i) {
        int u = path[i], v = path[i + 1];
        bool found = false;
        long long best = 0;
        for (const Edge& e : edges) {
            if (e.u == u && e.v == v && (!found || e.w < best)) {
                best = e.w;
                found = true;
            }
        }
        if (!found) return {false, 0};
        total += best;
    }
    return {true, total};
}

// ---------------- 对照用的其他算法 ----------------

// O(n^2) 版 Dijkstra，仅用于**非负权**图的对拍
vector<long long> dijkstra(int n, const vector<Edge>& edges, int s) {
    vector<vector<pair<int, long long>>> adj(n);
    for (const Edge& e : edges) adj[e.u].push_back({e.v, e.w});
    vector<long long> dist(n, INF);
    vector<char> used(n, 0);
    if (0 <= s && s < n) dist[s] = 0;
    for (int it = 0; it < n; ++it) {
        int u = -1;
        for (int i = 0; i < n; ++i)
            if (!used[i] && dist[i] != INF && (u == -1 || dist[i] < dist[u])) u = i;
        if (u == -1) break;
        used[u] = 1;
        for (auto [v, w] : adj[u])
            if (dist[u] + w < dist[v]) dist[v] = dist[u] + w;
    }
    return dist;
}

// Floyd-Warshall 全源最短路。返回 (图中是否存在任意负环, 距离矩阵)。时间 O(n^3)
pair<bool, vector<vector<long long>>> floydWarshall(int n, const vector<Edge>& edges) {
    vector<vector<long long>> d(n, vector<long long>(n, INF));
    for (int i = 0; i < n; ++i) d[i][i] = 0;
    for (const Edge& e : edges)
        if (0 <= e.u && e.u < n && 0 <= e.v && e.v < n) d[e.u][e.v] = min(d[e.u][e.v], e.w);
    for (int k = 0; k < n; ++k)
        for (int i = 0; i < n; ++i) {
            if (d[i][k] == INF) continue;
            for (int j = 0; j < n; ++j)
                if (d[k][j] != INF && d[i][k] + d[k][j] < d[i][j]) d[i][j] = d[i][k] + d[k][j];
        }
    bool hasNeg = false;
    for (int i = 0; i < n; ++i) if (d[i][i] < 0) hasNeg = true;
    return {hasNeg, d};
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
    int n, m, s;
    if (cin >> n) {  // IO 模式
        if (!(cin >> m)) m = 0;
        if (!(cin >> s)) s = 0;
        vector<Edge> edges;
        for (int i = 0; i < m; ++i) {
            int u = 0, v = 0;
            long long w = 0;
            if (!(cin >> u)) u = 0;
            if (!(cin >> v)) v = 0;
            if (!(cin >> w)) w = 0;
            if (0 <= u && u < n && 0 <= v && v < n) edges.push_back({u, v, w});
        }
        vector<long long> dist;
        vector<int> pre;
        bool hasNeg = bellmanFord(n, edges, s, dist, pre);
        cout << (hasNeg ? 1 : 0) << "\n";
        if (!hasNeg) {
            for (int i = 0; i < n; ++i) {          // n = 0 时输出空行
                if (i) cout << " ";
                if (dist[i] == INF) cout << "INF";
                else cout << dist[i];
            }
            cout << "\n";
        }
        return 0;
    }

    // 示例一：含负权边但无负环
    //   0 -(4)-> 1 -(2)-> 3 -(2)-> 4 -(1)-> 1（正环，不影响）
    //   0 -(2)-> 2 -(-3)-> 1，所以 0→1 走 0→2→1 只需 -1
    {
        vector<Edge> e1 = {{0, 1, 4}, {0, 2, 2}, {2, 1, -3}, {1, 3, 2}, {2, 3, 5}, {3, 4, 2}, {4, 1, 1}};
        vector<long long> dist, distS;
        vector<int> pre;
        bool neg = bellmanFord(5, e1, 0, dist, pre);
        assert(!neg);
        assert((dist == vector<long long>{0, -1, 2, 1, 3}));
        assert(spfa(5, e1, 0, distS) == false);
        assert(distS == dist);
        auto [ok, w] = pathWeight(e1, buildPath(pre, 0, 1));
        assert(ok && w == -1);
        assert((buildPath(pre, 0, 1) == vector<int>{0, 2, 1}));
    }

    // 示例二：负环 1 -> 2 -> 1，权值和 -2，且从 0 可达
    {
        vector<Edge> e2 = {{0, 1, 1}, {1, 2, -3}, {2, 1, 1}};
        vector<long long> dist, distS;
        vector<int> pre;
        assert(bellmanFord(3, e2, 0, dist, pre) == true);
        assert(spfa(3, e2, 0, distS) == true);
    }

    // 负环存在但**从 s 不可达**：不影响 s 的最短路，检测也应当报告「无」
    {
        vector<Edge> e3 = {{0, 1, 1}, {2, 3, -5}, {3, 2, 2}};
        vector<long long> dist, distS;
        vector<int> pre;
        assert(bellmanFord(4, e3, 0, dist, pre) == false);
        assert(dist[0] == 0 && dist[1] == 1);
        assert(dist[2] == INF && dist[3] == INF);
        assert(spfa(4, e3, 0, distS) == false);
    }

    // 自环：负权 → 负环；零权 / 正权 → 不是负环
    {
        vector<long long> d;
        vector<int> p;
        assert(bellmanFord(1, {{0, 0, -1}}, 0, d, p) == true);
        assert(bellmanFord(1, {{0, 0, 0}}, 0, d, p) == false);
        assert(bellmanFord(1, {{0, 0, 5}}, 0, d, p) == false);
    }

    // 边界：单点无边 / 无边图 / 空图
    {
        vector<long long> d;
        vector<int> p;
        vector<Edge> none;
        assert(bellmanFord(1, none, 0, d, p) == false);
        assert(d[0] == 0 && p[0] == -1);
        bellmanFord(3, none, 0, d, p);
        assert(d[0] == 0 && d[1] == INF && d[2] == INF);
        bellmanFord(0, none, 0, d, p);
        assert(d.empty() && p.empty());
    }

    // 重边取最小：两条 0→1，权值 7 和 3
    {
        vector<long long> d;
        vector<int> p;
        assert(bellmanFord(2, {{0, 1, 7}, {0, 1, 3}}, 0, d, p) == false);
        assert((d == vector<long long>{0, 3}));
    }

    // 路径还原：链状图 0→1→2→3
    {
        vector<Edge> e6 = {{0, 1, 2}, {1, 2, 3}, {2, 3, 4}};
        vector<long long> d;
        vector<int> pre;
        bellmanFord(4, e6, 0, d, pre);
        assert((d == vector<long long>{0, 2, 5, 9}));
        assert((buildPath(pre, 0, 3) == vector<int>{0, 1, 2, 3}));
        assert((buildPath(pre, 0, 0) == vector<int>{0}));
        // 不可达的点还原不出路径
        vector<Edge> e7 = {{2, 3, 1}};
        bellmanFord(4, e7, 0, d, pre);
        assert(buildPath(pre, 0, 3).empty());
    }

    LCG rng(20260928ULL);

    // 随机对拍一：**非负权**图，与 O(n^2) Dijkstra 比对
    for (int t = 0; t < 400; ++t) {
        int n2 = rng.next(1, 8);
        int m2 = rng.next(0, n2 * 2);
        vector<Edge> edges;
        for (int i = 0; i < m2; ++i)
            edges.push_back({rng.next(0, n2 - 1), rng.next(0, n2 - 1), rng.next(0, 9)});
        int s2 = rng.next(0, n2 - 1);
        vector<long long> dist;
        vector<int> pre;
        bool neg = bellmanFord(n2, edges, s2, dist, pre);
        assert(!neg);                            // 非负权不可能有负环
        assert(dist == dijkstra(n2, edges, s2));
    }

    // 随机对拍二：**允许负权**，与 SPFA 互相印证，并用 Floyd 与最优性条件兜底
    for (int t = 0; t < 400; ++t) {
        int n2 = rng.next(1, 7);
        int m2 = rng.next(0, n2 * 2);
        vector<Edge> edges;
        for (int i = 0; i < m2; ++i)
            edges.push_back({rng.next(0, n2 - 1), rng.next(0, n2 - 1), rng.next(-6, 9)});
        int s2 = rng.next(0, n2 - 1);

        vector<long long> dist, distS;
        vector<int> pre;
        bool neg = bellmanFord(n2, edges, s2, dist, pre);
        bool negS = spfa(n2, edges, s2, distS);
        assert(neg == negS);                     // 两种判据必须一致

        if (neg) continue;                       // 有负环时距离无意义，只校验判据一致

        assert(dist == distS);                   // 两版距离必须一致
        assert(dist[s2] == 0);

        // 最优性条件：无负环时不存在还能被松弛的边（三角不等式成立）
        for (const Edge& e : edges)
            if (dist[e.u] != INF) assert(dist[e.u] + e.w >= dist[e.v]);

        // 每个可达点的距离都必须对应一条真实存在的、权重相等的路径
        for (int v = 0; v < n2; ++v) {
            if (dist[v] == INF) {
                assert(buildPath(pre, s2, v).empty());
                continue;
            }
            vector<int> p = buildPath(pre, s2, v);
            assert(!p.empty() && p[0] == s2 && p.back() == v);
            auto [ok, w] = pathWeight(edges, p);
            assert(ok && w == dist[v]);
        }

        // 与 Floyd-Warshall 交叉验证（只在整张图都没有负环时才可比）
        auto [hasNegAny, mat] = floydWarshall(n2, edges);
        if (!hasNegAny)
            for (int v = 0; v < n2; ++v) assert(dist[v] == mat[s2][v]);
    }

    cout << "all tests passed" << endl;
    return 0;
}
