// Bellman-Ford 單源最短路（支持負權邊 + 負環檢測 + 路徑還原）
// 編譯：g++ -std=c++17 -O2 solution.cpp -o solution && ./solution
//
// 思路：Dijkstra 依賴「已出隊的點距離不再變小」這個貪心性質，有負權邊就不成立。
//   Bellman-Ford 的出發點是一個樸素事實：**一條最短路最多經過 n-1 條邊**
//   （再多就一定繞了環；正環/零環可以刪掉，負環則根本不存在「最短路」）。
//   於是把所有邊整體鬆弛 n-1 輪即可；再多做一輪還能鬆弛就說明繞了負環 —— 這既是
//   負環檢測，也是 Bellman-Ford 相比 Dijkstra 的核心能力。
//   要點：只鬆弛 dist[u] 有限的邊；某一輪沒有更新就提前退出；記錄 pre 可還原最短路。
//   SPFA（隊列優化）只讓「上一輪被更新過的點」繼續鬆弛其出邊，負環判據改成
//   「某點入隊次數 > n」。
//   注意：負環必須**從 s 可達**才會被檢測到；圖另一頭的負環與 s 無關。
//
// 輸入（空白分隔）：n m s / 接着 m 行 u v w
// 輸出：第 1 行 1=存在 s 可達的負環、0=不存在；
//       第 2 行（無負環時）dist[0]..dist[n-1]，空格分隔，不可達輸出 INF
// 無 stdin 輸入時運行內置斷言測試。
#include <algorithm>
#include <cassert>
#include <deque>
#include <iostream>
#include <vector>

using namespace std;

const long long INF = (1LL << 60);   // 足夠大，且 dist[u] + w 不會溢出

struct Edge {
    int u, v;
    long long w;
};

// Bellman-Ford：返回 has_neg_cycle，距離寫入 dist，前驅寫入 pre。時間 O(n·m)，空間 O(n)
bool bellmanFord(int n, const vector<Edge>& edges, int s, vector<long long>& dist, vector<int>& pre) {
    dist.assign(n, INF);
    pre.assign(n, -1);
    if (0 <= s && s < n) dist[s] = 0;

    bool hasNeg = false;
    for (int it = 0; it < n; ++it) {
        bool changed = false;
        for (const Edge& e : edges) {
            // 只從已可達的點往外鬆弛，避免 INF + w 污染結果
            if (dist[e.u] != INF && dist[e.u] + e.w < dist[e.v]) {
                dist[e.v] = dist[e.u] + e.w;
                pre[e.v] = e.u;
                changed = true;
            }
        }
        if (!changed) break;              // 這一輪沒人被更新，後面也不可能再變
        if (it == n - 1) hasNeg = true;   // 第 n 輪還能鬆弛 → 繞了負環
    }
    return hasNeg;
}

// Bellman-Ford 的隊列優化版（SPFA）。返回 has_neg_cycle，距離寫入 dist
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
                    if (cnt[v] > n) return true;   // 入隊超過 n 次 → 有負環
                    q.push_back(v);
                }
            }
        }
    }
    return false;
}

// 沿前驅數組還原 s → t 的最短路；不可達返回空。時間 O(路徑長度)
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

// 校驗路徑：返回 (每條邊都存在, 路徑總權重)
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

// ---------------- 對照用的其他算法 ----------------

// O(n^2) 版 Dijkstra，僅用於**非負權**圖的對拍
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

// Floyd-Warshall 全源最短路。返回 (圖中是否存在任意負環, 距離矩陣)。時間 O(n^3)
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

// 與 Python 版同規模的固定隨機序列（LCG），兩版各自獨立與暴力解對拍
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
            for (int i = 0; i < n; ++i) {          // n = 0 時輸出空行
                if (i) cout << " ";
                if (dist[i] == INF) cout << "INF";
                else cout << dist[i];
            }
            cout << "\n";
        }
        return 0;
    }

    // 示例一：含負權邊但無負環
    //   0 -(4)-> 1 -(2)-> 3 -(2)-> 4 -(1)-> 1（正環，不影響）
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

    // 示例二：負環 1 -> 2 -> 1，權值和 -2，且從 0 可達
    {
        vector<Edge> e2 = {{0, 1, 1}, {1, 2, -3}, {2, 1, 1}};
        vector<long long> dist, distS;
        vector<int> pre;
        assert(bellmanFord(3, e2, 0, dist, pre) == true);
        assert(spfa(3, e2, 0, distS) == true);
    }

    // 負環存在但**從 s 不可達**：不影響 s 的最短路，檢測也應當報告「無」
    {
        vector<Edge> e3 = {{0, 1, 1}, {2, 3, -5}, {3, 2, 2}};
        vector<long long> dist, distS;
        vector<int> pre;
        assert(bellmanFord(4, e3, 0, dist, pre) == false);
        assert(dist[0] == 0 && dist[1] == 1);
        assert(dist[2] == INF && dist[3] == INF);
        assert(spfa(4, e3, 0, distS) == false);
    }

    // 自環：負權 → 負環；零權 / 正權 → 不是負環
    {
        vector<long long> d;
        vector<int> p;
        assert(bellmanFord(1, {{0, 0, -1}}, 0, d, p) == true);
        assert(bellmanFord(1, {{0, 0, 0}}, 0, d, p) == false);
        assert(bellmanFord(1, {{0, 0, 5}}, 0, d, p) == false);
    }

    // 邊界：單點無邊 / 無邊圖 / 空圖
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

    // 重邊取最小：兩條 0→1，權值 7 和 3
    {
        vector<long long> d;
        vector<int> p;
        assert(bellmanFord(2, {{0, 1, 7}, {0, 1, 3}}, 0, d, p) == false);
        assert((d == vector<long long>{0, 3}));
    }

    // 路徑還原：鏈狀圖 0→1→2→3
    {
        vector<Edge> e6 = {{0, 1, 2}, {1, 2, 3}, {2, 3, 4}};
        vector<long long> d;
        vector<int> pre;
        bellmanFord(4, e6, 0, d, pre);
        assert((d == vector<long long>{0, 2, 5, 9}));
        assert((buildPath(pre, 0, 3) == vector<int>{0, 1, 2, 3}));
        assert((buildPath(pre, 0, 0) == vector<int>{0}));
        // 不可達的點還原不出路徑
        vector<Edge> e7 = {{2, 3, 1}};
        bellmanFord(4, e7, 0, d, pre);
        assert(buildPath(pre, 0, 3).empty());
    }

    LCG rng(20260928ULL);

    // 隨機對拍一：**非負權**圖，與 O(n^2) Dijkstra 比對
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
        assert(!neg);                            // 非負權不可能有負環
        assert(dist == dijkstra(n2, edges, s2));
    }

    // 隨機對拍二：**允許負權**，與 SPFA 互相印證，並用 Floyd 與最優性條件兜底
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
        assert(neg == negS);                     // 兩種判據必須一致

        if (neg) continue;                       // 有負環時距離無意義，只校驗判據一致

        assert(dist == distS);                   // 兩版距離必須一致
        assert(dist[s2] == 0);

        // 最優性條件：無負環時不存在還能被鬆弛的邊（三角不等式成立）
        for (const Edge& e : edges)
            if (dist[e.u] != INF) assert(dist[e.u] + e.w >= dist[e.v]);

        // 每個可達點的距離都必須對應一條真實存在的、權重相等的路徑
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

        // 與 Floyd-Warshall 交叉驗證（只在整張圖都沒有負環時才可比）
        auto [hasNegAny, mat] = floydWarshall(n2, edges);
        if (!hasNegAny)
            for (int v = 0; v < n2; ++v) assert(dist[v] == mat[s2][v]);
    }

    cout << "all tests passed" << endl;
    return 0;
}
