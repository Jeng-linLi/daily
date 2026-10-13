// 最小費用最大流（Min-Cost Max-Flow, MCMF）
//
// 題意：
//     有向圖上每條邊有容量 cap 與單位費用 cost。
//     在所有「s→t 流量最大」的流之中，找總費用最小的那一條，並輸出每條邊的流量。
//
// 思路：
//     ### 連續最短路增廣（Successive Shortest Path, SSP）
//     Ford-Fulkerson 的推廣：每次不沿「任意」殘量路徑增廣，而是沿「費用最小」的
//     殘量路徑增廣到該路徑飽和。
//
//     > 定理：每次沿當前殘量網路中以費用為權的最短路增廣，得到的每個中間流都是
//     > 「同流量下費用最小」的流；因此結束時即為最小費用最大流。
//
//     ### 殘量網路與反向邊
//     正向邊 (u→v, cap, cost) 配一條反向邊 (v→u, 0, -cost)，增廣時正減反加。
//     反向邊的負費用讓它能夠「退流」，也讓殘量網路出現負權邊。
//
//     ### 兩個版本
//       SPFA（Bellman-Ford 佇列版）  每次 O(VE)，可處理負權邊，最直覺
//       Dijkstra + Johnson 勢能      每次 O(E log V)：先用一次 SPFA 求初始勢能 h，
//                                    再用既約費用 c'(u,v) = cost + h[u] - h[v] ≥ 0 跑 Dijkstra
//     本題兩個版本互相對拍，答案必須一致。
//
//     ### 應用：指派問題
//     s → 每個人(cap 1, cost 0)、人 → 工作(cap 1, cost c[i][j])、工作 → t(cap 1, cost 0)。
//     最大流必為 n，最小費用即最小總費用指派（匈牙利算法是它的 O(n³) 特例）。
//
// 複雜度：SSP+SPFA O(F·V·E)、SSP+Dijkstra O(F·E log V)，F 為最大流值；空間 O(V+E)。
//
// 輸入（stdin，空白分隔）：n m s t，接著 m 條 u v cap cost（點號 0-indexed）。
//     非數字 token 視為輸入結束；n ≤ 100、m ≤ 400、cap ∈ [0,10^6]、cost ∈ [-10^6,10^6]；
//     點號越界的邊視為無效（流量記 0）。
// 輸出（stdout，3 行）：
//     1: flow cost     （SPFA 版最大流值與最小費用）
//     2: m 個整數       （每條輸入邊的流量，按輸入順序）
//     3: flow2 cost2   （Dijkstra + 勢能版）
// 無 stdin 時跑內置斷言測試並輸出 all tests passed。

#include <cassert>
#include <cerrno>
#include <cstdlib>
#include <algorithm>
#include <deque>
#include <functional>
#include <iostream>
#include <iterator>
#include <queue>
#include <random>
#include <string>
#include <tuple>
#include <vector>

#ifdef _WIN32
#include <fcntl.h>
#include <io.h>
#endif

namespace {

constexpr int MAXN = 100;
constexpr int MAXM = 400;
constexpr long long CAP_MAX = 1000000LL;
constexpr long long COST_ABS = 1000000LL;
constexpr long long INF = (1LL << 60);

// ---------------------------------------------------------------- 費用流
struct Edge {
    int to;
    long long cap;
    long long cost;
};

class MinCostFlow {
public:
    explicit MinCostFlow(int n) : n_(n), adj_(n) {}

    int addEdge(int u, int v, long long cap, long long cost) {
        int idx = static_cast<int>(edges_.size());
        edges_.push_back({v, cap, cost});
        edges_.push_back({u, 0, -cost});
        adj_[u].push_back(idx);
        adj_[v].push_back(idx ^ 1);
        return idx;
    }

    const std::vector<Edge>& edges() const { return edges_; }

    // 從 s 出發的最短路（允許負權邊）
    void spfa(int s, std::vector<long long>& dist, std::vector<int>& prevv,
              std::vector<int>& preve) const {
        dist.assign(n_, INF);
        prevv.assign(n_, -1);
        preve.assign(n_, -1);
        std::vector<char> inq(n_, 0);
        dist[s] = 0;
        std::deque<int> dq;
        dq.push_back(s);
        inq[s] = 1;
        while (!dq.empty()) {
            int u = dq.front();
            dq.pop_front();
            inq[u] = 0;
            for (int ei : adj_[u]) {
                const Edge& e = edges_[ei];
                if (e.cap > 0 && dist[u] + e.cost < dist[e.to]) {
                    dist[e.to] = dist[u] + e.cost;
                    prevv[e.to] = u;
                    preve[e.to] = ei;
                    if (!inq[e.to]) {
                        inq[e.to] = 1;
                        dq.push_back(e.to);
                    }
                }
            }
        }
    }

    // Bellman-Ford 負費用環檢測（只看從 s 可達的子圖）
    // 初始圖若含負環，費用可沿負環任意壓低，最小費用無下界、SSP 也不會終止。
    // 本題偵測到此情況一律回傳 (0, 0)，兩個語言用同一套判定。
    bool hasNegativeCycle(int s) const {
        std::vector<long long> dist(n_, INF);
        dist[s] = 0;
        for (int round = 0; round < n_ - 1; ++round) {
            bool changed = false;
            for (int u = 0; u < n_; ++u) {
                if (dist[u] == INF) continue;
                for (int ei : adj_[u]) {
                    const Edge& e = edges_[ei];
                    if (e.cap > 0 && dist[u] + e.cost < dist[e.to]) {
                        dist[e.to] = dist[u] + e.cost;
                        changed = true;
                    }
                }
            }
            if (!changed) break;
        }
        for (int u = 0; u < n_; ++u) {
            if (dist[u] == INF) continue;
            for (int ei : adj_[u]) {
                const Edge& e = edges_[ei];
                if (e.cap > 0 && dist[u] + e.cost < dist[e.to]) return true;
            }
        }
        return false;
    }

    // SPFA 版 SSP；會**修改**殘量網路（Python 版同）
    std::pair<long long, long long> minCostFlowSpfa(int s, int t, long long maxf) {
        if (s == t) return {0, 0};
        if (hasNegativeCycle(s)) return {0, 0};
        long long flow = 0, cost = 0;
        std::vector<long long> dist;
        std::vector<int> prevv, preve;
        while (flow < maxf) {
            spfa(s, dist, prevv, preve);
            if (dist[t] == INF) break;
            long long d = maxf - flow;
            for (int v = t; v != s; v = prevv[v]) d = std::min(d, edges_[preve[v]].cap);
            flow += d;
            cost += d * dist[t];
            for (int v = t; v != s; v = prevv[v]) {
                int ei = preve[v];
                edges_[ei].cap -= d;
                edges_[ei ^ 1].cap += d;
            }
        }
        return {flow, cost};
    }

    // Dijkstra + Johnson 勢能版；會**修改**殘量網路
    std::pair<long long, long long> minCostFlowDijkstra(int s, int t, long long maxf) {
        if (s == t) return {0, 0};
        if (hasNegativeCycle(s)) return {0, 0};
        std::vector<long long> h, dist;
        std::vector<int> prevv, preve;
        spfa(s, h, prevv, preve);                      // 初始勢能（可處理負邊）
        long long flow = 0, cost = 0;
        while (flow < maxf) {
            dist.assign(n_, INF);
            prevv.assign(n_, -1);
            preve.assign(n_, -1);
            dist[s] = 0;
            std::priority_queue<std::pair<long long, int>,
                                std::vector<std::pair<long long, int>>,
                                std::greater<std::pair<long long, int>>> pq;
            pq.push({0, s});
            while (!pq.empty()) {
                long long d = pq.top().first;
                int u = pq.top().second;
                pq.pop();
                if (d != dist[u]) continue;
                for (int ei : adj_[u]) {
                    const Edge& e = edges_[ei];
                    if (e.cap <= 0 || h[u] == INF || h[e.to] == INF) continue;
                    long long nd = d + e.cost + h[u] - h[e.to];     // 既約費用 ≥ 0
                    if (nd < dist[e.to]) {
                        dist[e.to] = nd;
                        prevv[e.to] = u;
                        preve[e.to] = ei;
                        pq.push({nd, e.to});
                    }
                }
            }
            if (dist[t] == INF) break;
            for (int v = 0; v < n_; ++v)
                if (dist[v] < INF) h[v] += dist[v];
            long long d = maxf - flow;
            for (int v = t; v != s; v = prevv[v]) d = std::min(d, edges_[preve[v]].cap);
            flow += d;
            cost += d * (h[t] - h[s]);                 // 真實最短路長 = h[t] - h[s]
            for (int v = t; v != s; v = prevv[v]) {
                int ei = preve[v];
                edges_[ei].cap -= d;
                edges_[ei ^ 1].cap += d;
            }
        }
        return {flow, cost};
    }

private:
    int n_;
    std::vector<Edge> edges_;
    std::vector<std::vector<int>> adj_;
};

// ---------------------------------------------------------------- 最大流（對拍）
class Dinic {
public:
    explicit Dinic(int n) : n_(n), graph_(n) {}

    void addEdge(int u, int v, long long c) {
        to_.push_back(v);
        cap_.push_back(c);
        graph_[u].push_back(static_cast<int>(to_.size()) - 1);
        to_.push_back(u);
        cap_.push_back(0);
        graph_[v].push_back(static_cast<int>(to_.size()) - 1);
    }

    long long maxFlow(int s, int t) {
        if (s == t) return 0;
        long long flow = 0;
        while (true) {
            std::vector<int> level(n_, -1);
            level[s] = 0;
            std::deque<int> q;
            q.push_back(s);
            while (!q.empty()) {
                int u = q.front();
                q.pop_front();
                for (int ei : graph_[u])
                    if (cap_[ei] > 0 && level[to_[ei]] < 0) {
                        level[to_[ei]] = level[u] + 1;
                        q.push_back(to_[ei]);
                    }
            }
            if (level[t] < 0) return flow;
            std::vector<int> it(n_, 0);
            while (true) {
                long long f = dfs(s, t, INF, level, it);
                if (f == 0) break;
                flow += f;
            }
        }
    }

private:
    int n_;
    std::vector<std::vector<int>> graph_;
    std::vector<int> to_;
    std::vector<long long> cap_;

    long long dfs(int u, int t, long long f, const std::vector<int>& level,
                  std::vector<int>& it) {
        if (u == t) return f;
        while (it[u] < static_cast<int>(graph_[u].size())) {
            int ei = graph_[u][it[u]];
            int v = to_[ei];
            if (cap_[ei] > 0 && level[v] == level[u] + 1) {
                long long d = dfs(v, t, std::min(f, cap_[ei]), level, it);
                if (d > 0) {
                    cap_[ei] -= d;
                    cap_[ei ^ 1] += d;
                    return d;
                }
            }
            it[u] += 1;
        }
        return 0;
    }
};

long long edmondsKarp(int n, const std::vector<std::tuple<int, int, long long>>& es,
                      int s, int t) {
    if (s == t) return 0;
    std::vector<std::vector<long long>> cap(n, std::vector<long long>(n, 0));
    for (const auto& e : es) cap[std::get<0>(e)][std::get<1>(e)] += std::get<2>(e);
    long long flow = 0;
    while (true) {
        std::vector<int> prev(n, -1);
        prev[s] = s;
        std::deque<int> q;
        q.push_back(s);
        while (!q.empty() && prev[t] == -1) {
            int u = q.front();
            q.pop_front();
            for (int v = 0; v < n; ++v)
                if (prev[v] == -1 && cap[u][v] > 0) {
                    prev[v] = u;
                    q.push_back(v);
                }
        }
        if (prev[t] == -1) return flow;
        long long d = INF;
        for (int v = t; v != s; v = prev[v]) d = std::min(d, cap[prev[v]][v]);
        for (int v = t; v != s; v = prev[v]) {
            cap[prev[v]][v] -= d;
            cap[v][prev[v]] += d;
        }
        flow += d;
    }
}

// ---------------------------------------------------------------- 暴力對拍
using Edge4 = std::tuple<int, int, long long, long long>;   // (u, v, cap, cost)

// 從 s 沿「容量 > 0」的邊可達的點：費用流只會在這些點上推流，
// 不會跑去 s 到不了的負費用環上繞（暴力對拍必須排除這種循環流）
std::vector<char> reachableFrom(int n, const std::vector<Edge4>& es, int s) {
    std::vector<std::vector<int>> adj(n);
    for (const auto& e : es)
        if (std::get<2>(e) > 0) adj[std::get<0>(e)].push_back(std::get<1>(e));
    std::vector<char> seen(n, 0);
    seen[s] = 1;
    std::deque<int> q;
    q.push_back(s);
    while (!q.empty()) {
        int u = q.front();
        q.pop_front();
        for (int v : adj[u])
            if (!seen[v]) { seen[v] = 1; q.push_back(v); }
    }
    return seen;
}

void bfEnum(std::size_t idx, const std::vector<Edge4>& es, int n, int s, int t,
            long long maxflow, const std::vector<char>& reach,
            std::vector<long long>& cur, long long& best) {
    if (idx == es.size()) {
        std::vector<long long> net(n, 0);
        long long total = 0;
        for (std::size_t i = 0; i < es.size(); ++i) {
            net[std::get<0>(es[i])] -= cur[i];
            net[std::get<1>(es[i])] += cur[i];
            total += cur[i] * std::get<3>(es[i]);
        }
        for (int v = 0; v < n; ++v)
            if (v != s && v != t && net[v] != 0) return;
        if (-net[s] != maxflow) return;
        if (best == INF || total < best) best = total;
        return;
    }
    long long hi = reach[std::get<0>(es[idx])] ? std::get<2>(es[idx]) : 0;
    for (long long f = 0; f <= hi; ++f) {
        cur[idx] = f;
        bfEnum(idx + 1, es, n, s, t, maxflow, reach, cur, best);
    }
    cur[idx] = 0;
}

std::pair<long long, long long> bfMinCostMaxFlow(int n, const std::vector<Edge4>& es,
                                                 int s, int t) {
    if (es.empty()) return {0, 0};
    std::vector<std::tuple<int, int, long long>> caps;
    for (const auto& e : es) caps.emplace_back(std::get<0>(e), std::get<1>(e), std::get<2>(e));
    long long maxflow = edmondsKarp(n, caps, s, t);
    long long best = INF;
    std::vector<long long> cur(es.size(), 0);
    std::vector<char> reach = reachableFrom(n, es, s);
    bfEnum(0, es, n, s, t, maxflow, reach, cur, best);
    return {maxflow, best == INF ? 0 : best};
}

// ---------------------------------------------------------------- 指派問題
std::pair<long long, std::vector<int>> minCostAssignment(
    const std::vector<std::vector<long long>>& cost) {
    int n = static_cast<int>(cost.size());
    if (n == 0) return {0, {}};
    int s = 2 * n, t = 2 * n + 1;
    MinCostFlow g(2 * n + 2);
    std::vector<int> fwd(static_cast<std::size_t>(n) * n, -1);
    for (int i = 0; i < n; ++i) g.addEdge(s, i, 1, 0);
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            fwd[static_cast<std::size_t>(i) * n + j] = g.addEdge(i, n + j, 1, cost[i][j]);
    for (int j = 0; j < n; ++j) g.addEdge(n + j, t, 1, 0);

    std::pair<long long, long long> res = g.minCostFlowSpfa(s, t, INF);
    assert(res.first == n);

    std::vector<int> who(n, -1);
    const std::vector<Edge>& ed = g.edges();
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j) {
            int idx = fwd[static_cast<std::size_t>(i) * n + j];
            if (ed[static_cast<std::size_t>(idx)].cap == 0) who[i] = j;   // 正向邊滿載
        }
    return {res.second, who};
}

void bfAssignmentPerm(const std::vector<std::vector<long long>>& cost, int i,
                      std::vector<int>& perm, std::vector<char>& used, long long& best) {
    int n = static_cast<int>(cost.size());
    if (i == n) {
        long long total = 0;
        for (int k = 0; k < n; ++k) total += cost[k][perm[k]];
        if (total < best) best = total;
        return;
    }
    for (int j = 0; j < n; ++j)
        if (!used[j]) {
            used[j] = 1;
            perm[i] = j;
            bfAssignmentPerm(cost, i + 1, perm, used, best);
            used[j] = 0;
        }
}

long long bfAssignment(const std::vector<std::vector<long long>>& cost) {
    int n = static_cast<int>(cost.size());
    if (n == 0) return 0;
    long long best = INF;
    std::vector<int> perm(n, 0);
    std::vector<char> used(n, 0);
    bfAssignmentPerm(cost, 0, perm, used, best);
    return best;
}

// ---------------------------------------------------------------- 輸出
std::string solve(int n, const std::vector<Edge4>& edges, int s, int t) {
    MinCostFlow g1(n), g2(n);
    for (const auto& e : edges) {
        g1.addEdge(std::get<0>(e), std::get<1>(e), std::get<2>(e), std::get<3>(e));
        g2.addEdge(std::get<0>(e), std::get<1>(e), std::get<2>(e), std::get<3>(e));
    }
    std::pair<long long, long long> r1 = g1.minCostFlowSpfa(s, t, INF);
    std::pair<long long, long long> r2 = g2.minCostFlowDijkstra(s, t, INF);

    std::string out = std::to_string(r1.first) + " " + std::to_string(r1.second) + "\n";
    for (std::size_t i = 0; i < edges.size(); ++i) {
        if (i) out += " ";
        out += std::to_string(std::get<2>(edges[i]) - g1.edges()[2 * i].cap);
    }
    out += "\n";
    out += std::to_string(r2.first) + " " + std::to_string(r2.second);
    return out;
}

void checkGraph(int n, const std::vector<Edge4>& edges, int s, int t) {
    MinCostFlow g1(n), g2(n);
    for (const auto& e : edges) {
        g1.addEdge(std::get<0>(e), std::get<1>(e), std::get<2>(e), std::get<3>(e));
        g2.addEdge(std::get<0>(e), std::get<1>(e), std::get<2>(e), std::get<3>(e));
    }
    if (n > 0 && s != t && g1.hasNegativeCycle(s)) {   // 費用無下界，只驗證兩個版本一致
        assert(g1.minCostFlowSpfa(s, t, INF) == std::make_pair(0LL, 0LL));
        assert(g2.minCostFlowDijkstra(s, t, INF) == std::make_pair(0LL, 0LL));
        return;
    }
    std::pair<long long, long long> r1 = g1.minCostFlowSpfa(s, t, INF);
    std::pair<long long, long long> r2 = g2.minCostFlowDijkstra(s, t, INF);
    assert(r1.first == r2.first);
    assert(r1.second == r2.second);

    Dinic din(n);
    for (const auto& e : edges) din.addEdge(std::get<0>(e), std::get<1>(e), std::get<2>(e));
    assert(r1.first == din.maxFlow(s, t));

    if (s == t) {
        assert(r1.first == 0 && r1.second == 0);
        return;
    }
    std::pair<long long, long long> bf = bfMinCostMaxFlow(n, edges, s, t);
    assert(r1.first == bf.first && r1.second == bf.second);
}

// ---------------------------------------------------------------- 解析
bool isSpaceCh(unsigned char c) {
    return c == ' ' || c == '\t' || c == '\n' || c == '\r' || c == '\v' || c == '\f';
}

std::vector<std::string> splitWs(const std::string& raw) {
    std::vector<std::string> out;
    std::string cur;
    for (unsigned char c : raw) {
        if (isSpaceCh(c)) {
            if (!cur.empty()) { out.push_back(cur); cur.clear(); }
        } else {
            cur.push_back(static_cast<char>(c));
        }
    }
    if (!cur.empty()) out.push_back(cur);
    return out;
}

bool isNumToken(const std::string& tok) {
    if (tok.empty()) return false;
    std::size_t i = (tok[0] == '+' || tok[0] == '-') ? 1 : 0;
    if (i >= tok.size()) return false;
    for (; i < tok.size(); ++i)
        if (tok[i] < '0' || tok[i] > '9') return false;
    return true;
}

long long parseLLClamp(const std::string& tok) {
    errno = 0;
    char* end = nullptr;
    long long v = std::strtoll(tok.c_str(), &end, 10);
    if (end == tok.c_str()) return 0;
    return v;                       // 溢位時 strtoll 已飽和到 LLONG_MAX / LLONG_MIN
}

// ---------------------------------------------------------------- 測試
void runTests() {
    std::mt19937 rng(20261013);

    // ---- 固定用例（與 Python 版逐字節一致）----
    assert(solve(4, {Edge4{0, 1, 3, 1}, Edge4{0, 2, 2, 4}, Edge4{1, 3, 2, 2},
                     Edge4{2, 3, 3, 1}}, 0, 3) == "4 16\n2 2 2 2\n4 16");
    assert(solve(2, {Edge4{0, 1, 5, 3}}, 0, 1) == "5 15\n5\n5 15");
    assert(solve(3, {}, 0, 2) == "0 0\n\n0 0");
    assert(solve(2, {Edge4{1, 0, 4, -2}}, 0, 1) == "0 0\n0\n0 0");

    checkGraph(4, {Edge4{0, 1, 3, 1}, Edge4{0, 2, 2, 4}, Edge4{1, 3, 2, 2},
                   Edge4{2, 3, 3, 1}}, 0, 3);
    checkGraph(2, {Edge4{0, 1, 5, 3}}, 0, 1);
    checkGraph(3, {}, 0, 2);
    checkGraph(2, {Edge4{1, 0, 4, -2}}, 0, 1);

    // ---- 隨機對拍：小圖（費用非負，不會有負環）+ 暴力枚舉 ----
    for (int iter = 0; iter < 400; ++iter) {
        int n = 2 + static_cast<int>(rng() % 4);
        int m = static_cast<int>(rng() % 6);
        std::vector<Edge4> es;
        for (int i = 0; i < m; ++i)
            es.emplace_back(static_cast<int>(rng() % static_cast<unsigned>(n)),
                            static_cast<int>(rng() % static_cast<unsigned>(n)),
                            rng() % 4, static_cast<long long>(rng() % 7));
        int s = static_cast<int>(rng() % static_cast<unsigned>(n));
        int t = static_cast<int>(rng() % static_cast<unsigned>(n));
        checkGraph(n, es, s, t);
    }

    // ---- 隨機對拍：負費用邊（邊一律 u < v，圖無環 ⟹ 保證無負環）----
    const long long pool[5] = {-3, -1, 0, 2, 5};
    for (int iter = 0; iter < 200; ++iter) {
        int n = 3 + static_cast<int>(rng() % 4);
        int m = 1 + static_cast<int>(rng() % 6);
        std::vector<Edge4> es;
        for (int i = 0; i < m; ++i) {
            int u = static_cast<int>(rng() % static_cast<unsigned>(n - 1));
            int v = u + 1 + static_cast<int>(rng() % static_cast<unsigned>(n - 1 - u));
            es.emplace_back(u, v, rng() % 4, pool[rng() % 5]);
        }
        int s = static_cast<int>(rng() % static_cast<unsigned>(n));
        int t = static_cast<int>(rng() % static_cast<unsigned>(n));
        checkGraph(n, es, s, t);
    }

    // ---- 隨機對拍：任意邊 + 任意費用（可能含負環，只驗證兩個版本一致）----
    for (int iter = 0; iter < 150; ++iter) {
        int n = 2 + static_cast<int>(rng() % 4);
        int m = static_cast<int>(rng() % 6);
        std::vector<Edge4> es;
        for (int i = 0; i < m; ++i)
            es.emplace_back(static_cast<int>(rng() % static_cast<unsigned>(n)),
                            static_cast<int>(rng() % static_cast<unsigned>(n)),
                            rng() % 4, static_cast<long long>(rng() % 11) - 4);
        int s = static_cast<int>(rng() % static_cast<unsigned>(n));
        int t = static_cast<int>(rng() % static_cast<unsigned>(n));
        checkGraph(n, es, s, t);
    }

    // ---- 隨機對拍：中圖（比對兩版本 + Dinic + 費用恆等式）----
    for (int iter = 0; iter < 120; ++iter) {
        int n = 2 + static_cast<int>(rng() % 19);
        int m = static_cast<int>(rng() % 41);
        std::vector<Edge4> es;
        for (int i = 0; i < m; ++i)
            es.emplace_back(static_cast<int>(rng() % static_cast<unsigned>(n)),
                            static_cast<int>(rng() % static_cast<unsigned>(n)),
                            rng() % 9, static_cast<long long>(rng() % 21) - 10);
        int s = static_cast<int>(rng() % static_cast<unsigned>(n));
        int t = static_cast<int>(rng() % static_cast<unsigned>(n));
        MinCostFlow g1(n), g2(n);
        for (const auto& e : es) {
            g1.addEdge(std::get<0>(e), std::get<1>(e), std::get<2>(e), std::get<3>(e));
            g2.addEdge(std::get<0>(e), std::get<1>(e), std::get<2>(e), std::get<3>(e));
        }
        if (n > 0 && s != t && g1.hasNegativeCycle(s)) {
            assert(g1.minCostFlowSpfa(s, t, INF) == std::make_pair(0LL, 0LL));
            assert(g2.minCostFlowDijkstra(s, t, INF) == std::make_pair(0LL, 0LL));
            continue;
        }
        std::pair<long long, long long> r1 = g1.minCostFlowSpfa(s, t, INF);
        std::pair<long long, long long> r2 = g2.minCostFlowDijkstra(s, t, INF);
        assert(r1.first == r2.first && r1.second == r2.second);
        Dinic din(n);
        for (const auto& e : es) din.addEdge(std::get<0>(e), std::get<1>(e), std::get<2>(e));
        assert(r1.first == din.maxFlow(s, t));
        long long total = 0;
        for (std::size_t i = 0; i < es.size(); ++i)
            total += (std::get<2>(es[i]) - g1.edges()[2 * i].cap) * std::get<3>(es[i]);
        assert(total == r1.second);
    }

    // ---- 應用：指派問題對拍 ----
    for (int iter = 0; iter < 200; ++iter) {
        int n = 1 + static_cast<int>(rng() % 6);
        std::vector<std::vector<long long>> cost(static_cast<std::size_t>(n),
                                                 std::vector<long long>(n));
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j)
                cost[i][j] = static_cast<long long>(rng() % 26) - 5;
        std::pair<long long, std::vector<int>> got = minCostAssignment(cost);
        assert(got.first == bfAssignment(cost));
        std::vector<int> sortedWho = got.second;
        std::sort(sortedWho.begin(), sortedWho.end());
        for (int i = 0; i < n; ++i) assert(sortedWho[static_cast<std::size_t>(i)] == i);
        long long recomputed = 0;
        for (int i = 0; i < n; ++i)
            recomputed += cost[i][got.second[static_cast<std::size_t>(i)]];
        assert(recomputed == got.first);
    }

    assert(minCostAssignment({{4, 1, 3}, {2, 0, 5}, {3, 2, 2}}).first == 5);
    assert(minCostAssignment({{10}}).first == 10);
    assert(minCostAssignment({}).first == 0);
}

void runIo(const std::string& raw) {
    std::vector<std::string> tok = splitWs(raw);
    auto readInt = [&](std::size_t i, bool& ok) -> long long {
        ok = i < tok.size() && isNumToken(tok[i]);
        return ok ? parseLLClamp(tok[i]) : 0;
    };
    auto clamp = [](long long v, long long lo, long long hi) {
        return v < lo ? lo : (v > hi ? hi : v);
    };
    bool ok;
    long long nRaw = readInt(0, ok);
    int n = ok ? static_cast<int>(clamp(nRaw, 0, MAXN)) : 0;
    long long mRaw = readInt(1, ok);
    int m = ok ? static_cast<int>(clamp(mRaw, 0, MAXM)) : 0;
    long long sRaw = readInt(2, ok);
    int s = ok ? static_cast<int>(clamp(sRaw, 0, MAXN - 1)) : 0;
    long long tRaw = readInt(3, ok);
    int t = ok ? static_cast<int>(clamp(tRaw, 0, MAXN - 1)) : 0;

    std::vector<Edge4> edges;
    int declared = 0;
    for (int i = 0; i < m; ++i) {
        bool o1, o2, o3, o4;
        long long u = readInt(static_cast<std::size_t>(4 + 4 * i), o1);
        long long v = readInt(static_cast<std::size_t>(5 + 4 * i), o2);
        long long cap = readInt(static_cast<std::size_t>(6 + 4 * i), o3);
        long long w = readInt(static_cast<std::size_t>(7 + 4 * i), o4);
        if (!o1 || !o2 || !o3 || !o4) break;
        declared += 1;
        if (u < 0 || u >= n || v < 0 || v >= n) continue;
        edges.emplace_back(static_cast<int>(u), static_cast<int>(v),
                           clamp(cap, 0, CAP_MAX), clamp(w, -COST_ABS, COST_ABS));
    }

    std::string out;
    if (n == 0 || s >= n || t >= n) {
        out = "0 0\n";
        for (int i = 0; i < declared; ++i) {
            if (i) out += " ";
            out += "0";
        }
        out += "\n0 0";
    } else {
        out = solve(n, edges, s, t);
    }
    out += "\n";
    std::cout << out;
}

}  // namespace

int main() {
#ifdef _WIN32
    // Windows 文字模式會把 \n 轉成 \r\n，與 Python 的二進位輸出對不起來
    _setmode(_fileno(stdout), _O_BINARY);
#endif
    std::string raw((std::istreambuf_iterator<char>(std::cin)),
                    std::istreambuf_iterator<char>());
    bool blank = true;
    for (unsigned char c : raw)
        if (!isSpaceCh(c)) { blank = false; break; }
    if (!blank) {
        runIo(raw);
    } else {
        runTests();
        std::cout << "all tests passed" << std::endl;
    }
    return 0;
}
