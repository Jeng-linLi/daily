// Floyd-Warshall 全源最短路（含負環檢測、路徑還原、傳遞閉包）
// 編譯：g++ -std=c++17 -O2 solution.cpp -o solution && ./solution
//
// 思路：Floyd-Warshall 本質是集合上的動態規劃：
//     dist[k][i][j] = 從 i 到 j、中間點只允許取自 {0..k-1} 的最短路長度
//   轉移只有兩種選擇：不經過 k（沿用上一層），或經過 k（dist[i][k] + dist[k][j]），取 min。
//   「第 k 層只依賴第 k-1 層」允許把第一維滾動掉、原地更新，代價是三重循環順序
//   必須是 k → i → j，順序寫錯結果就錯。
//
//   負環檢測：跑完後若存在 dist[i][i] < 0，說明從 i 繞一圈回到自己還能更短 —— 這就是負環。
//   Floyd 檢測的是**全圖任意位置**的負環，不像 Bellman-Ford 只檢測單源可達的那些。
//   路徑還原：next[i][j] 記從 i 走向 j 的第一步；被 i→k→j 更新時改寫成 next[i][k]。
//
//   複雜度 O(n^3) 時間、O(n^2) 空間，勝在寫起來極短且一次算出全部點對。
//
// 輸入（空白分隔）：n m / 接着 m 行 u v w
// 輸出：第 1 行負環標記；第 2..n+1 行距離矩陣（不可達打印 INF）；最後一行圖的直徑
// 無 stdin 輸入時運行內置斷言測試。
#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <limits>
#include <random>
#include <sstream>
#include <string>
#include <vector>

using namespace std;

using ll = long long;
const ll INF = (1LL << 60);        // 足夠大，且 dist + w 不會溢出
const string INF_STR = "INF";

struct Edge {
    int u, v;
    ll w;
};

// Floyd-Warshall：返回是否有負環，dist 與 nxt 作爲輸出參數
bool floydWarshall(int n, const vector<Edge>& edges, vector<vector<ll>>& dist, vector<vector<int>>& nxt) {
    dist.assign(n, vector<ll>(n, INF));
    nxt.assign(n, vector<int>(n, -1));
    for (int i = 0; i < n; ++i) {
        dist[i][i] = 0;            // 自己到自己距離 0（負自環會把它改小）
        nxt[i][i] = i;             // 保證「dist 有限 ⟹ next 有效」這條不變量
    }
    for (const Edge& e : edges) {
        if (e.u < 0 || e.u >= n || e.v < 0 || e.v >= n) continue;   // 越界邊忽略
        if (e.w < dist[e.u][e.v]) {
            dist[e.u][e.v] = e.w;
            nxt[e.u][e.v] = e.v;
        }
    }

    // 三重循環順序必須是 k → i → j
    for (int k = 0; k < n; ++k) {
        for (int i = 0; i < n; ++i) {
            if (dist[i][k] == INF) continue;                       // 小剪枝
            ll dik = dist[i][k];
            for (int j = 0; j < n; ++j) {
                if (dist[k][j] == INF) continue;
                ll nd = dik + dist[k][j];
                if (nd < dist[i][j]) {
                    dist[i][j] = nd;
                    nxt[i][j] = nxt[i][k];                         // 第一步 = i 走向 k 的第一步
                }
            }
        }
    }

    for (int i = 0; i < n; ++i)
        if (dist[i][i] < 0) return true;
    return false;
}

// 用 next 矩陣還原 u → v 的一條最短路；不可達（或負環繞圈）返回空
vector<int> restorePath(const vector<vector<int>>& nxt, int u, int v) {
    int n = (int)nxt.size();
    if (u < 0 || u >= n || v < 0 || v >= n) return {};
    if (nxt[u][v] == -1) return {};
    vector<int> path{u};
    int cur = u;
    while (cur != v) {
        cur = nxt[cur][v];
        if (cur == -1) return {};
        path.push_back(cur);
        if ((int)path.size() > n + 1) return {};                   // 負環時防禦死循環
    }
    return path;
}

// 所有有限距離中的最大值；沒有則返回 INF
ll diameter(const vector<vector<ll>>& dist) {
    ll best = INF;
    for (const auto& row : dist)
        for (ll d : row)
            if (d != INF && (best == INF || d > best)) best = d;
    return best;
}

// Warshall 傳遞閉包：把 Floyd 的 min/+ 換成 or/and，同一個三重循環骨架
vector<vector<int>> warshallClosure(int n, const vector<Edge>& edges) {
    vector<vector<int>> reach(n, vector<int>(n, 0));
    for (int i = 0; i < n; ++i) reach[i][i] = 1;
    for (const Edge& e : edges)
        if (e.u >= 0 && e.u < n && e.v >= 0 && e.v < n) reach[e.u][e.v] = 1;
    for (int k = 0; k < n; ++k)
        for (int i = 0; i < n; ++i) {
            if (!reach[i][k]) continue;
            for (int j = 0; j < n; ++j)
                if (reach[k][j]) reach[i][j] = 1;
        }
    return reach;
}

// ---------------- 對照實現（Bellman-Ford，用於對拍） ----------------

pair<vector<ll>, bool> bellmanFordFrom(int n, const vector<Edge>& edges, int s) {
    vector<ll> dist(n, INF);
    dist[s] = 0;
    for (int it = 0; it < n; ++it) {
        bool changed = false;
        for (const Edge& e : edges) {
            if (e.u < 0 || e.u >= n || e.v < 0 || e.v >= n) continue;
            if (dist[e.u] != INF && dist[e.u] + e.w < dist[e.v]) {
                dist[e.v] = dist[e.u] + e.w;
                changed = true;
            }
        }
        if (!changed) break;
        if (it == n - 1) return {dist, true};                      // 第 n 輪還能鬆弛 → 負環
    }
    return {dist, false};
}

// 跑 n 次 Bellman-Ford 得到全源最短路
pair<vector<vector<ll>>, bool> bruteAllPairs(int n, const vector<Edge>& edges) {
    vector<vector<ll>> all;
    bool neg = false;
    for (int s = 0; s < n; ++s) {
        auto [d, hasNeg] = bellmanFordFrom(n, edges, s);
        all.push_back(d);
        neg = neg || hasNeg;
    }
    return {all, neg};
}

// ---------------- IO ----------------

static string fmtDist(ll x) { return (x == INF) ? INF_STR : to_string(x); }

static void runIo(const string& data) {
    istringstream iss(data);
    ll n = 0, m = 0;
    if (!(iss >> n)) n = 0;
    if (!(iss >> m)) m = 0;
    vector<Edge> edges;
    for (ll i = 0; i < m; ++i) {
        ll u = 0, v = 0, w = 0;
        if (!(iss >> u)) u = 0;
        if (!(iss >> v)) v = 0;
        if (!(iss >> w)) w = 0;
        edges.push_back({(int)u, (int)v, w});
    }

    vector<vector<ll>> dist;
    vector<vector<int>> nxt;
    bool hasNeg = floydWarshall((int)n, edges, dist, nxt);
    cout << (hasNeg ? 1 : 0) << '\n';
    for (int i = 0; i < (int)n; ++i) {
        for (int j = 0; j < (int)n; ++j) {
            if (j) cout << ' ';
            cout << fmtDist(dist[i][j]);
        }
        cout << '\n';
    }
    cout << fmtDist(diameter(dist)) << '\n';
}

static void runTests() {
    // README 示例：4 個點、5 條邊
    // 有負權邊但沒有負環
    vector<Edge> edges = {{0, 1, 3}, {0, 2, 7}, {1, 2, -2}, {1, 3, 5}, {2, 3, 1}};
    vector<vector<ll>> dist;
    vector<vector<int>> nxt;
    bool hasNeg = floydWarshall(4, edges, dist, nxt);
    assert(!hasNeg);
    ll expected[4][4] = {
        {0, 3, 1, 2},
        {INF, 0, -2, -1},
        {INF, INF, 0, 1},
        {INF, INF, INF, 0},
    };
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j)
            assert(dist[i][j] == expected[i][j]);
    assert(diameter(dist) == 3);

    // 路徑還原：0 → 3 的最短路是 0 → 1 → 2 → 3
    vector<int> p = restorePath(nxt, 0, 3);
    assert((p == vector<int>{0, 1, 2, 3}));
    ll total = 0;
    for (size_t i = 0; i + 1 < p.size(); ++i)
        for (const Edge& e : edges)
            if (e.u == p[i] && e.v == p[i + 1]) { total += e.w; break; }
    assert(total == dist[0][3]);

    // 傳遞閉包與距離矩陣一致：有限距離 ⟺ 可達
    auto reach = warshallClosure(4, edges);
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j)
            assert((reach[i][j] != 0) == (dist[i][j] != INF));

    // 負環：1 → 2 → 3 → 1 權和 -1
    vector<vector<ll>> d2;
    vector<vector<int>> n2;
    assert(floydWarshall(4, {{0, 1, 1}, {1, 2, 1}, {2, 3, -3}, {3, 1, 1}}, d2, n2));

    // 負自環也算負環
    vector<vector<ll>> ds;
    vector<vector<int>> ns;
    assert(floydWarshall(2, {{0, 0, -1}}, ds, ns));

    // 零環 / 正環不算負環
    vector<vector<ll>> dz, dp;
    vector<vector<int>> nz, np;
    assert(!floydWarshall(3, {{0, 1, 1}, {1, 2, 1}, {2, 0, -2}}, dz, nz));
    assert(!floydWarshall(3, {{0, 1, 1}, {1, 2, 1}, {2, 0, 1}}, dp, np));

    // 空圖 / 單點 / 沒有邊
    vector<vector<ll>> d0;
    vector<vector<int>> n0;
    assert(!floydWarshall(0, {}, d0, n0));
    assert(d0.empty() && diameter(d0) == INF);

    vector<vector<ll>> d1;
    vector<vector<int>> n1;
    assert(!floydWarshall(1, {}, d1, n1));
    assert(d1.size() == 1 && d1[0][0] == 0 && diameter(d1) == 0);
    assert((restorePath(n1, 0, 0) == vector<int>{0}));

    vector<vector<ll>> dEmpty;
    vector<vector<int>> nEmpty;
    assert(!floydWarshall(3, {}, dEmpty, nEmpty));
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            assert(dEmpty[i][j] == (i == j ? 0 : INF));
    assert(diameter(dEmpty) == 0);

    // 不可達：單向邊
    vector<vector<ll>> d3;
    vector<vector<int>> n3;
    assert(!floydWarshall(3, {{0, 1, 5}}, d3, n3));
    assert(d3[0][1] == 5 && d3[1][0] == INF);
    assert(restorePath(n3, 1, 0).empty());
    assert(diameter(d3) == 5);

    // 重複邊取最小的那條
    vector<vector<ll>> d4;
    vector<vector<int>> n4;
    assert(!floydWarshall(2, {{0, 1, 9}, {0, 1, 4}}, d4, n4));
    assert(d4[0][1] == 4);

    mt19937 rng(20260929);

    // 隨機對拍：與跑 n 次 Bellman-Ford 的結果逐格比對
    for (int iter = 0; iter < 400; ++iter) {
        int n = 1 + (int)(rng() % 6);
        int m = (int)(rng() % 11);
        vector<Edge> es;
        for (int i = 0; i < m; ++i)
            es.push_back({(int)(rng() % n), (int)(rng() % n), (ll)(rng() % 19) - 6});

        vector<vector<ll>> distR;
        vector<vector<int>> nxtR;
        bool neg = floydWarshall(n, es, distR, nxtR);
        auto [ref, refNeg] = bruteAllPairs(n, es);
        assert(neg == refNeg);
        if (!neg) {
            for (int i = 0; i < n; ++i)
                for (int j = 0; j < n; ++j)
                    assert(distR[i][j] == ref[i][j]);
            // 三角不等式
            for (int i = 0; i < n; ++i)
                for (int j = 0; j < n; ++j)
                    for (int k = 0; k < n; ++k)
                        if (distR[i][k] != INF && distR[k][j] != INF)
                            assert(distR[i][j] <= distR[i][k] + distR[k][j]);
            // 傳遞閉包一致性
            auto rc = warshallClosure(n, es);
            for (int i = 0; i < n; ++i)
                for (int j = 0; j < n; ++j)
                    assert((rc[i][j] != 0) == (distR[i][j] != INF));
            // 還原出的路徑合法且總長等於 dist
            for (int i = 0; i < n; ++i)
                for (int j = 0; j < n; ++j) {
                    vector<int> path = restorePath(nxtR, i, j);
                    if (distR[i][j] == INF) {
                        assert(path.empty());
                    } else if (i == j) {
                        assert((path == vector<int>{i}));
                    } else {
                        assert(!path.empty() && path.front() == i && path.back() == j);
                        ll cost = 0;
                        bool ok = true;
                        for (size_t t = 0; t + 1 < path.size(); ++t) {
                            ll bestW = INF;
                            for (const Edge& e : es)
                                if (e.u == path[t] && e.v == path[t + 1]) bestW = min(bestW, e.w);
                            if (bestW == INF) { ok = false; break; }
                            cost += bestW;
                        }
                        assert(ok && cost == distR[i][j]);
                    }
                }
        }
    }

    cout << "all tests passed" << '\n';
}

int main() {
    string data, line;
    bool hasInput = false;
    while (getline(cin, line)) {
        data += line;
        data += '\n';
        if (!line.empty()) hasInput = true;
    }
    if (hasInput) runIo(data);
    else runTests();
    return 0;
}
