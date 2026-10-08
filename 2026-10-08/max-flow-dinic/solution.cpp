// Dinic 最大流與二分圖匹配（最小割 / 邊不相交路徑 / 匈牙利算法）
//
// 與 solution.py 完全同構：同樣的演算法、同樣的確定性規則、同樣的輸入輸出格式。
// 編譯：g++ -std=c++17 -O2 -Wall solution.cpp -o solution
// 無 stdin 輸入時運行內置斷言測試並輸出 `all tests passed`。

#include <algorithm>
#include <cassert>
#include <cctype>
#include <climits>
#include <cstdint>
#include <deque>
#include <iostream>
#include <random>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

using namespace std;

using ll = long long;

static const ll INF = (ll)4e18;        // 大於任何可能的流量，且不溢出 int64
static const ll CAP_MAX = (ll)1e15;    // 容量上限：保證 Σ 流量仍在 int64 內（與 Python 一致）
static const int MAX_N = 200;          // 防止 IO 模式下圖太大
static const int MAX_M = 5000;

// ---------------------------------------------------------------- 工具

// 嚴格整數規則 ^[+-]?[0-9]+$，與 Python 的 parse_int 完全一致。
static bool tryLL(const string &s, ll &out) {
    if (s.empty()) return false;
    size_t i = 0;
    if (s[0] == '+' || s[0] == '-') i = 1;
    if (i >= s.size()) return false;
    ll val = 0;
    for (; i < s.size(); i++) {
        if (!isdigit((unsigned char)s[i])) return false;
        val = val * 10 + (s[i] - '0');
        if (val > CAP_MAX) val = CAP_MAX;      // 截斷，避免溢位
    }
    out = (s[0] == '-') ? -val : val;
    return true;
}

static ll clampCap(ll v) {
    if (v > CAP_MAX) return CAP_MAX;
    if (v < -CAP_MAX) return -CAP_MAX;
    return v;
}

// ---------------------------------------------------------------- Dinic

struct Edge {
    int to;      // 終點
    ll cap;      // 殘量
    int rev;     // 反向邊在 g[to] 中的下標
};

class Dinic {
public:
    int n;
    vector<vector<Edge>> g;
    vector<tuple<int, int, ll>> edges;         // 原始邊 (u, v, c)

    explicit Dinic(int n_) : n(n_), g(n_) {}

    int addEdge(int u, int v, ll c) {
        if (u < 0 || u >= n || v < 0 || v >= n) return -1;
        if (c < 0) c = 0;
        Edge a{v, c, (int)g[v].size()};
        Edge b{u, 0, (int)g[u].size()};
        g[u].push_back(a);
        g[v].push_back(b);
        edges.emplace_back(u, v, c);
        return (int)edges.size() - 1;
    }

    // 在殘量網絡上分層；level[t] < 0 表示已無增廣路。
    vector<int> bfs(int s, int t) {
        vector<int> level(n, -1);
        level[s] = 0;
        deque<int> q;
        q.push_back(s);
        while (!q.empty()) {
            int u = q.front();
            q.pop_front();
            for (const Edge &e : g[u]) {
                if (e.cap > 0 && level[e.to] < 0) {
                    level[e.to] = level[u] + 1;
                    q.push_back(e.to);
                }
            }
        }
        (void)t;
        return level;
    }

    // 在層次圖上找一條增廣路並推送；返回推送量，0 表示本輪阻塞流已推完。
    ll blocking(int s, int t, vector<int> &level, vector<int> &it) {
        vector<int> stk;
        vector<pair<int, int>> path;           // (u, 邊在 g[u] 中的下標)
        stk.push_back(s);
        while (!stk.empty()) {
            int u = stk.back();
            if (u == t) {
                ll f = INF;
                for (auto &pr : path) {
                    ll c = g[pr.first][pr.second].cap;
                    if (c < f) f = c;
                }
                for (auto &pr : path) {
                    Edge &e = g[pr.first][pr.second];
                    e.cap -= f;
                    g[e.to][e.rev].cap += f;   // 反向邊 +f，保留「反悔」的餘地
                }
                return f;
            }
            bool advanced = false;
            while (it[u] < (int)g[u].size()) {
                int idx = it[u];
                int v = g[u][idx].to;
                if (g[u][idx].cap > 0 && level[v] == level[u] + 1) {
                    path.push_back({u, idx});
                    stk.push_back(v);
                    advanced = true;
                    break;
                }
                it[u]++;
            }
            if (!advanced) {
                level[u] = -1;                 // 該點在層次圖上走不到 t，剪掉
                stk.pop_back();
                if (!path.empty()) {
                    pair<int, int> pr = path.back();
                    path.pop_back();
                    it[pr.first]++;
                }
            }
        }
        return 0;
    }

    ll maxFlow(int s, int t) {
        if (s < 0 || s >= n || t < 0 || t >= n || s == t) return 0;
        ll total = 0;
        while (true) {
            vector<int> level = bfs(s, t);
            if (level[t] < 0) break;
            vector<int> it(n, 0);
            while (true) {
                ll f = blocking(s, t, level, it);
                if (f == 0) break;
                total += f;
            }
        }
        return total;
    }

    // 殘量網絡中從 s 可達的點（升序），即最小割的源側 S。
    vector<int> reachableFrom(int s) {
        vector<int> seen(n, 0);
        if (s >= 0 && s < n) {
            seen[s] = 1;
            deque<int> q;
            q.push_back(s);
            while (!q.empty()) {
                int u = q.front();
                q.pop_front();
                for (const Edge &e : g[u]) {
                    if (e.cap > 0 && !seen[e.to]) {
                        seen[e.to] = 1;
                        q.push_back(e.to);
                    }
                }
            }
        }
        vector<int> out;
        for (int i = 0; i < n; i++)
            if (seen[i]) out.push_back(i);
        return out;
    }

    // 返回 (割容量, 源側點集升序, 割邊[(u,v,c)] 按 (u,v) 升序)。
    tuple<ll, vector<int>, vector<tuple<int, int, ll>>> minCut(int s) {
        vector<int> side = reachableFrom(s);
        vector<int> inSide(n, 0);
        for (int x : side) inSide[x] = 1;
        ll cap = 0;
        vector<tuple<int, int, ll>> cuts;
        for (auto &ed : edges) {
            int u = get<0>(ed), v = get<1>(ed);
            ll c = get<2>(ed);
            if (inSide[u] && !inSide[v]) {
                cap += c;
                cuts.emplace_back(u, v, c);
            }
        }
        sort(cuts.begin(), cuts.end(), [](const tuple<int, int, ll> &a, const tuple<int, int, ll> &b) {
            if (get<0>(a) != get<0>(b)) return get<0>(a) < get<0>(b);
            return get<1>(a) < get<1>(b);
        });
        sort(side.begin(), side.end());
        return make_tuple(cap, side, cuts);
    }
};

static Dinic buildFlow(int n, const vector<tuple<int, int, ll>> &edges) {
    Dinic net(n);
    for (auto &ed : edges) net.addEdge(get<0>(ed), get<1>(ed), get<2>(ed));
    return net;
}

// ---------------------------------------------------------------- 對拍用的其它算法

// Edmonds-Karp：每次用 BFS 找最短增廣路，O(V E²)，作為 Dinic 的獨立對拍。
static ll edmondsKarp(int n, const vector<tuple<int, int, ll>> &edges, int s, int t) {
    if (n <= 0 || s < 0 || s >= n || t < 0 || t >= n || s == t) return 0;
    vector<vector<ll>> cap(n, vector<ll>(n, 0));
    for (auto &ed : edges) {
        int u = get<0>(ed), v = get<1>(ed);
        ll c = get<2>(ed);
        if (u >= 0 && u < n && v >= 0 && v < n && c > 0) cap[u][v] += c;
    }
    ll flow = 0;
    while (true) {
        vector<int> pre(n, -1);
        pre[s] = -2;
        deque<int> q;
        q.push_back(s);
        while (!q.empty() && pre[t] == -1) {
            int u = q.front();
            q.pop_front();
            for (int v = 0; v < n; v++)
                if (pre[v] == -1 && cap[u][v] > 0) {
                    pre[v] = u;
                    q.push_back(v);
                }
        }
        if (pre[t] == -1) break;
        ll f = INF;
        int v = t;
        while (v != s) {
            f = min(f, cap[pre[v]][v]);
            v = pre[v];
        }
        v = t;
        while (v != s) {
            int p = pre[v];
            cap[p][v] -= f;
            cap[v][p] += f;
            v = p;
        }
        flow += f;
    }
    return flow;
}

// 枚舉所有「含 s 不含 t」的點集求最小割容量（n ≤ 14），驗證最大流最小割定理。
static ll minCutBruteforce(int n, const vector<tuple<int, int, ll>> &edges, int s, int t) {
    if (n <= 0 || s < 0 || s >= n || t < 0 || t >= n || s == t) return 0;
    ll best = INF;
    int full = 1 << n;
    for (int mask = 0; mask < full; mask++) {
        if (!((mask >> s) & 1) || ((mask >> t) & 1)) continue;
        ll total = 0;
        for (auto &ed : edges) {
            int u = get<0>(ed), v = get<1>(ed);
            ll c = get<2>(ed);
            if (((mask >> u) & 1) && !((mask >> v) & 1)) total += c;
        }
        if (total < best) best = total;
    }
    return best == INF ? 0 : best;
}

// 無向圖邊不相交路徑數：每條邊拆成兩條容量 1 的有向邊後求最大流。
static ll edgeDisjointPaths(int n, const vector<tuple<int, int, ll>> &edges, int s, int t) {
    if (n <= 0 || s < 0 || s >= n || t < 0 || t >= n || s == t) return 0;
    Dinic net(n);
    for (auto &ed : edges) {
        int u = get<0>(ed), v = get<1>(ed);
        net.addEdge(u, v, 1);
        net.addEdge(v, u, 1);
    }
    return net.maxFlow(s, t);
}

// ---------------------------------------------------------------- 二分圖匹配

// 匈牙利 / Kuhn 算法。adj[u] 為左點 u 可連的右點列表。
// 返回 (最大匹配數, 匹配邊 [(左點, 右點)]，按左點升序)。
static bool tryKuhn(int u, const vector<vector<int>> &adj, vector<int> &matchR, vector<int> &seen) {
    for (int v : adj[u]) {
        if (seen[v]) continue;
        seen[v] = 1;
        if (matchR[v] == -1 || tryKuhn(matchR[v], adj, matchR, seen)) {
            matchR[v] = u;
            return true;
        }
    }
    return false;
}

static pair<int, vector<pair<int, int>>> bipartiteMatching(const vector<vector<int>> &adj, int nl, int nr) {
    if (nl <= 0 || nr <= 0) return {0, {}};
    vector<int> matchR(nr, -1);
    int count = 0;
    for (int u = 0; u < nl; u++) {
        vector<int> seen(nr, 0);
        if (tryKuhn(u, adj, matchR, seen)) count++;
    }
    vector<pair<int, int>> pairs;
    for (int v = 0; v < nr; v++)
        if (matchR[v] != -1) pairs.push_back({matchR[v], v});
    sort(pairs.begin(), pairs.end());
    return {count, pairs};
}

// 把二分圖建模成流網絡求最大匹配：S→左(1)、左→右(1)、右→T(1)。
static int matchingViaDinic(const vector<vector<int>> &adj, int nl, int nr) {
    if (nl <= 0 || nr <= 0) return 0;
    int total = nl + nr;
    int src = total, dst = total + 1;
    Dinic net(total + 2);
    for (int u = 0; u < nl; u++) net.addEdge(src, u, 1);
    for (int v = 0; v < nr; v++) net.addEdge(nl + v, dst, 1);
    for (int u = 0; u < nl; u++)
        for (int v : adj[u]) net.addEdge(u, nl + v, 1);
    return (int)net.maxFlow(src, dst);
}

// 暴力枚舉匹配（nl ≤ 8），作為匈牙利算法的對拍。
static void matchRec(int u, int used, const vector<vector<int>> &adj, int nl, int &best) {
    if (u == nl) {
        int cnt = 0;
        for (int x = used; x; x >>= 1) cnt += (x & 1);
        if (cnt > best) best = cnt;
        return;
    }
    matchRec(u + 1, used, adj, nl, best);
    for (int v : adj[u])
        if (!((used >> v) & 1)) matchRec(u + 1, used | (1 << v), adj, nl, best);
}

static int matchingBruteforce(const vector<vector<int>> &adj, int nl, int nr) {
    (void)nr;
    if (nl <= 0 || nr <= 0) return 0;
    int best = 0;
    matchRec(0, 0, adj, nl, best);
    return best;
}

// 過濾越界與重複的邊，讓兩個語言的行為一致。
static vector<vector<int>> normalizeAdj(const vector<vector<int>> &adj, int nl, int nr) {
    vector<vector<int>> out(nl);
    for (int u = 0; u < nl; u++) {
        const vector<int> &row = (u < (int)adj.size()) ? adj[u] : vector<int>();
        for (int v : row) {
            if (v < 0 || v >= nr) continue;
            bool dup = false;
            for (int x : out[u])
                if (x == v) dup = true;
            if (!dup) out[u].push_back(v);
        }
    }
    return out;
}

// ---------------------------------------------------------------- IO 與測試

static void runIO(const string &raw) {
    istringstream iss(raw);
    vector<string> toks;
    string tk;
    while (iss >> tk) toks.push_back(tk);
    size_t pos = 0;
    auto nxt = [&]() -> ll {
        ll v = 0;
        if (pos < toks.size()) {
            ll parsed = 0;
            if (tryLL(toks[pos], parsed)) v = clampCap(parsed);
        }
        pos++;
        return v;
    };

    ll nv = nxt(), mv = nxt(), sv = nxt(), tv = nxt();
    int n = (nv < 0) ? 0 : (int)min<ll>(nv, MAX_N);
    int m = (mv < 0) ? 0 : (int)min<ll>(mv, MAX_M);
    int s = (int)sv, t = (int)tv;
    vector<tuple<int, int, ll>> edges;
    for (int i = 0; i < m; i++) {
        int u = (int)nxt(), v = (int)nxt();
        ll c = nxt();
        if (c < 0) c = 0;
        edges.emplace_back(u, v, c);
    }

    Dinic net = buildFlow(n, edges);
    ll flow = net.maxFlow(s, t);
    auto cutRes = net.minCut(s);
    ll cutCap = get<0>(cutRes);
    vector<int> side = get<1>(cutRes);
    vector<tuple<int, int, ll>> cuts = get<2>(cutRes);

    ll nlv = nxt(), nrv = nxt(), kv = nxt();
    int nl = (nlv < 0) ? 0 : (int)min<ll>(nlv, MAX_N);
    int nr = (nrv < 0) ? 0 : (int)min<ll>(nrv, MAX_N);
    int k = (kv < 0) ? 0 : (int)min<ll>(kv, MAX_M);
    vector<vector<int>> adj(nl);
    for (int i = 0; i < k; i++) {
        int u = (int)nxt(), v = (int)nxt();
        if (u >= 0 && u < nl && v >= 0 && v < nr) adj[u].push_back(v);
    }
    adj = normalizeAdj(adj, nl, nr);
    auto res = bipartiteMatching(adj, nl, nr);

    ostringstream oss;
    oss << flow << "\n" << cutCap << "\n" << side.size() << "\n";
    for (size_t i = 0; i < side.size(); i++) {
        if (i) oss << " ";
        oss << side[i];
    }
    oss << "\n" << cuts.size() << "\n";
    for (auto &ed : cuts) oss << get<0>(ed) << " " << get<1>(ed) << " " << get<2>(ed) << "\n";
    oss << res.first << "\n";
    for (auto &pr : res.second) oss << pr.first << " " << pr.second << "\n";
    cout << oss.str();
}

static void runTests() {
    // ---- 經典範例：4 點圖，最大流 = 5（對應最小割 {0} / {1,2,3}）----
    vector<tuple<int, int, ll>> e0 = {{0, 1, 3}, {0, 2, 2}, {1, 2, 1}, {1, 3, 2}, {2, 3, 4}};
    {
        Dinic net0 = buildFlow(4, e0);
        assert(net0.maxFlow(0, 3) == 5);
        assert(net0.maxFlow(0, 3) == 0);        // 殘量網絡已滿，再跑一次為 0
    }
    {
        Dinic net0b = buildFlow(4, e0);
        assert(net0b.maxFlow(0, 3) == 5);
        auto cr = net0b.minCut(0);
        assert(get<0>(cr) == 5);
        assert((get<1>(cr) == vector<int>{0}));
        assert(get<2>(cr).size() == 2);
        assert(get<0>(get<2>(cr)[0]) == 0 && get<1>(get<2>(cr)[0]) == 1 && get<2>(get<2>(cr)[0]) == 3);
        assert(get<0>(get<2>(cr)[1]) == 0 && get<1>(get<2>(cr)[1]) == 2 && get<2>(get<2>(cr)[1]) == 2);
    }

    // ---- 二分圖：左 3 右 3，最大匹配 = 3（完美匹配）----
    vector<vector<int>> adj1 = {{0, 1}, {0, 2}, {1}};
    {
        auto r1 = bipartiteMatching(adj1, 3, 3);
        assert(r1.first == 3);
        assert((r1.second == vector<pair<int, int>>{{0, 0}, {1, 2}, {2, 1}}));
    }
    assert(matchingViaDinic(adj1, 3, 3) == 3);
    assert(matchingBruteforce(adj1, 3, 3) == 3);

    // ---- 非完美匹配：左 3 右 2，最多 2 ----
    vector<vector<int>> adj2 = {{0, 1}, {0}, {1}};
    {
        auto r2 = bipartiteMatching(adj2, 3, 2);
        assert(r2.first == 2);
        assert(r2.second.size() == 2);
        assert(r2.second[0].first != r2.second[1].first);
    }
    assert(matchingViaDinic(adj2, 3, 2) == 2);
    assert(matchingBruteforce(adj2, 3, 2) == 2);

    // ---- 空圖 / 退化情形 ----
    {
        Dinic empty0 = buildFlow(0, {});
        assert(empty0.maxFlow(0, 0) == 0);
    }
    assert(bipartiteMatching({}, 0, 0).first == 0);
    assert(bipartiteMatching({{}}, 1, 0).first == 0);
    assert(matchingViaDinic({}, 0, 0) == 0);
    assert(edmondsKarp(0, {}, 0, 0) == 0);
    assert(minCutBruteforce(0, {}, 0, 0) == 0);
    assert(edgeDisjointPaths(3, {{0, 1, 1}}, 0, 2) == 0);

    // ---- 邊不相交路徑：兩條 0→3 的路 ----
    vector<tuple<int, int, ll>> eu = {{0, 1, 1}, {0, 2, 1}, {1, 3, 1}, {2, 3, 1}};
    assert(edgeDisjointPaths(4, eu, 0, 3) == 2);

    // ---- 隨機對拍：Dinic vs Edmonds-Karp vs 暴力最小割 ----
    mt19937 rng(20261008);
    for (int iter = 0; iter < 120; iter++) {
        int n = 2 + (int)(rng() % 7);
        int m = (int)(rng() % 13);
        vector<tuple<int, int, ll>> es;
        for (int i = 0; i < m; i++) {
            int u = (int)(rng() % (unsigned)n);
            int v = (int)(rng() % (unsigned)n);
            int c = (int)(rng() % 10);
            es.emplace_back(u, v, c);
        }
        int s = (int)(rng() % (unsigned)n);
        int t = (int)(rng() % (unsigned)n);
        if (s == t) continue;
        Dinic netA = buildFlow(n, es);
        ll f1 = netA.maxFlow(s, t);
        ll f2 = edmondsKarp(n, es, s, t);
        ll f3 = minCutBruteforce(n, es, s, t);
        assert(f1 == f2 && f2 == f3);
        ll outCap = 0;
        for (auto &ed : es)
            if (get<0>(ed) == s) outCap += get<2>(ed);
        assert(f1 <= outCap);
        Dinic netB = buildFlow(n, es);
        ll fl = netB.maxFlow(s, t);
        auto cr = netB.minCut(s);
        assert(fl == get<0>(cr));
        bool hasS = false, hasT = false;
        for (int x : get<1>(cr)) {
            if (x == s) hasS = true;
            if (x == t) hasT = true;
        }
        assert(hasS && !hasT);
        ll sum = 0;
        for (auto &ed : get<2>(cr)) sum += get<2>(ed);
        assert(sum == get<0>(cr));
    }

    // ---- 隨機對拍：匈牙利 vs Dinic vs 暴力匹配 ----
    for (int iter = 0; iter < 120; iter++) {
        int nl = 1 + (int)(rng() % 7);
        int nr = 1 + (int)(rng() % 7);
        int k = (int)(rng() % 13);
        vector<vector<int>> adj(nl);
        for (int i = 0; i < k; i++) {
            int u = (int)(rng() % (unsigned)nl);
            int v = (int)(rng() % (unsigned)nr);
            adj[u].push_back(v);
        }
        adj = normalizeAdj(adj, nl, nr);
        auto rr = bipartiteMatching(adj, nl, nr);
        int c1 = rr.first;
        int c2 = matchingViaDinic(adj, nl, nr);
        int c3 = matchingBruteforce(adj, nl, nr);
        assert(c1 == c2 && c2 == c3);
        assert((int)rr.second.size() == c1);
        vector<int> seenL(nl, 0), seenR(nr, 0);
        for (auto &pr : rr.second) {
            assert(!seenL[pr.first] && !seenR[pr.second]);
            seenL[pr.first] = seenR[pr.second] = 1;
            bool exists = false;
            for (int v : adj[pr.first])
                if (v == pr.second) exists = true;
            assert(exists);
        }
    }
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    string raw, line;
    bool any = false;
    while (getline(cin, line)) {
        raw += line;
        raw += "\n";
        for (char ch : line)
            if (!isspace((unsigned char)ch)) any = true;
    }
    if (any) {
        runIO(raw);
    } else {
        runTests();
        cout << "all tests passed" << endl;
    }
    return 0;
}
