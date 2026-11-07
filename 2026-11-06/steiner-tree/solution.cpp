// Steiner Tree（斯坦納樹）— 最小連通終端集子圖
//
// 題意：
//   給定一張帶非負邊權的無向圖，以及 k 個「終端節點」(terminals)。
//   求一棵連通子圖（樹），使所有終端都被連通，且總邊權和最小。
//   允許額外包含非終端節點（Steiner 點）以換取更小權重。
//
// 思路（狀壓 DP + Dijkstra）：
//   把 k 個終端編號 0..k-1，定義
//       dp[mask][v] = 連通 mask 所指終端集合、且經過頂點 v 的最小權值。
//   1) 初始化：dp[1<<i][term[i]] = 0。
//   2) 枚舉 mask（由小到大）：
//      a. 子樹合併：在頂點 v 把兩個互斥子集拼起來
//         dp[mask][v] = min_{sub} dp[sub][v] + dp[mask^sub][v]
//      b. 邊鬆弛：以 dp[mask][·] 為初始距離跑 Dijkstra，
//         沿邊 (v,u,w) 更新 dp[mask][u] = min(..., dp[mask][v]+w)
//   3) 答案 = min_v dp[(1<<k)-1][v]；全無窮則 -1。
//
// 時間複雜度：O(3^k · n + 2^k · m · log n)
// 空間複雜度：O(2^k · n)
//
// 輸入格式（stdin 模式）：
//   n m k
//   接下來 m 行：u v w   （0-based，無向邊）
//   最後 1 行：k 個終端頂點
//   輸出：最小總權值，不可連通則 -1。
// 無 stdin 時跑內建斷言並印出 "all tests passed"。

#include <iostream>
#include <vector>
#include <queue>
#include <string>
#include <sstream>
#include <tuple>
#include <algorithm>
#include <numeric>
#include <cassert>
using namespace std;

using ll = long long;
const ll INF = (ll)1e18;

struct Edge { int to; ll w; };

// 回傳連通所有 terminals 的最小總權值；不可連通回傳 INF。
ll steiner_tree(int n, const vector<tuple<int,int,ll>>& edges, const vector<int>& terminals) {
    int k = (int)terminals.size();
    vector<vector<Edge>> adj(n);
    for (auto& e : edges) {
        int u = get<0>(e), v = get<1>(e); ll w = get<2>(e);
        adj[u].push_back({v, w});
        adj[v].push_back({u, w});
    }

    int full = 1 << k;
    vector<vector<ll>> dp(full, vector<ll>(n, INF));

    for (int i = 0; i < k; ++i)
        dp[1 << i][terminals[i]] = 0;

    for (int mask = 1; mask < full; ++mask) {
        // (a) 子樹合併：在頂點 v 把兩個互斥子集拼起來
        for (int sub = (mask - 1) & mask; sub; sub = (sub - 1) & mask) {
            int other = mask ^ sub;
            for (int v = 0; v < n; ++v) {
                if (dp[sub][v] < INF && dp[other][v] < INF) {
                    ll cand = dp[sub][v] + dp[other][v];
                    if (cand < dp[mask][v]) dp[mask][v] = cand;
                }
            }
        }
        // (b) 邊鬆弛：以 dp[mask] 為初始距離跑 Dijkstra
        priority_queue<pair<ll,int>, vector<pair<ll,int>>, greater<pair<ll,int>>> pq;
        for (int v = 0; v < n; ++v)
            if (dp[mask][v] < INF) pq.push({dp[mask][v], v});
        while (!pq.empty()) {
            auto [d, v] = pq.top(); pq.pop();
            if (d > dp[mask][v]) continue;
            for (auto& e : adj[v]) {
                ll nd = d + e.w;
                if (nd < dp[mask][e.to]) {
                    dp[mask][e.to] = nd;
                    pq.push({nd, e.to});
                }
            }
        }
    }

    return *min_element(dp[full - 1].begin(), dp[full - 1].end());
}

// ---- 以下為測試用暴力解（窮舉邊子集），僅供比對 ----
namespace brute {
    struct DSU {
        vector<int> p;
        DSU(int n): p(n) { iota(p.begin(), p.end(), 0); }
        int find(int x){ return p[x]==x ? x : p[x]=find(p[x]); }
        void uni(int a,int b){ p[find(a)]=find(b); }
    };
    ll solve(int n, const vector<tuple<int,int,ll>>& edges, const vector<int>& terminals) {
        int m = (int)edges.size();
        ll best = INF;
        for (int r = 1; r < (1 << m); ++r) {
            int cnt = __builtin_popcount(r);
            // 檢查這些邊是否連通所有終端
            DSU dsu(n);
            ll w = 0;
            for (int i = 0; i < m; ++i) if (r & (1 << i)) {
                int u = get<0>(edges[i]), v = get<1>(edges[i]);
                dsu.uni(u, v); w += get<2>(edges[i]);
            }
            bool ok = true;
            int root = dsu.find(terminals[0]);
            for (int t : terminals) if (dsu.find(t) != root) { ok = false; break; }
            if (ok && w < best) best = w;
            (void)cnt;
        }
        return best;
    }
}

void run_tests() {
    // 測試 1：基礎題（預期 5）
    {
        int n = 5;
        vector<tuple<int,int,ll>> edges = {{0,1,1},{1,2,1},{0,3,4},{3,4,1},{4,2,4}};
        vector<int> term = {2,3,4};
        ll a = steiner_tree(n, edges, term);
        assert(a == 5);
        assert(brute::solve(n, edges, term) == 5);
    }
    // 測試 2：不連通（預期 INF）
    {
        int n = 4;
        vector<tuple<int,int,ll>> edges = {{0,1,1},{2,3,2}};
        vector<int> term = {0,3};
        ll a = steiner_tree(n, edges, term);
        assert(a >= INF);
        assert(brute::solve(n, edges, term) >= INF);
    }
    // 測試 3：終端為全部節點（預期 10）
    {
        int n = 4;
        vector<tuple<int,int,ll>> edges = {{0,1,2},{1,2,3},{2,3,5}};
        vector<int> term = {0,1,2,3};
        ll a = steiner_tree(n, edges, term);
        assert(a == 10);
        assert(brute::solve(n, edges, term) == 10);
    }
    // 測試 4：直接邊連（預期 1）
    {
        int n = 5;
        vector<tuple<int,int,ll>> edges = {{0,1,1},{1,2,1},{0,3,1},{3,4,1},{4,2,1}};
        vector<int> term = {2,4};
        ll a = steiner_tree(n, edges, term);
        assert(a == 1);
        assert(brute::solve(n, edges, term) == 1);
    }
    // 測試 5：繞路用 Steiner 點省錢（預期 4）
    {
        int n = 5;
        vector<tuple<int,int,ll>> edges = {{0,1,1},{1,2,5},{2,3,1},{1,4,1},{4,2,1}};
        vector<int> term = {0,3};
        ll a = steiner_tree(n, edges, term);
        assert(a == 4);
        assert(brute::solve(n, edges, term) == 4);
    }
    cout << "all tests passed" << endl;
}

void solve_io() {
    int n, m, k;
    if (!(cin >> n >> m >> k)) return;
    vector<tuple<int,int,ll>> edges(m);
    for (int i = 0; i < m; ++i) {
        int u, v; ll w;
        cin >> u >> v >> w;
        edges[i] = {u, v, w};
    }
    vector<int> term(k);
    for (int i = 0; i < k; ++i) cin >> term[i];
    ll ans = steiner_tree(n, edges, term);
    cout << (ans >= INF ? -1 : ans) << endl;
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    // 整段讀入 stdin；無輸入則跑內建測試
    stringstream ss;
    ss << cin.rdbuf();
    string all = ss.str();
    if (all.empty()) {
        run_tests();
        return 0;
    }
    istringstream iss(all);
    int n, m, k;
    if (!(iss >> n >> m >> k)) { run_tests(); return 0; }
    vector<tuple<int,int,ll>> edges(m);
    for (int i = 0; i < m; ++i) {
        int u, v; ll w;
        iss >> u >> v >> w;
        edges[i] = {u, v, w};
    }
    vector<int> term(k);
    for (int i = 0; i < k; ++i) iss >> term[i];
    ll ans = steiner_tree(n, edges, term);
    cout << (ans >= INF ? -1 : ans) << endl;
    return 0;
}
