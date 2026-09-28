// Dijkstra 單源最短路（非負權圖，二叉堆實現）
// 編譯：g++ -std=c++17 -O2 solution.cpp -o solution && ./solution
#include <iostream>
#include <vector>
#include <queue>
#include <tuple>
#include <limits>
#include <cmath>
#include <cassert>

using namespace std;

const double INF = numeric_limits<double>::infinity();

// graph[u] = {(v, w), ...} 爲鄰接表；返回長度 n 的距離數組
vector<double> dijkstra(int n, const vector<vector<pair<int, double>>>& graph, int start) {
    vector<double> dist(n, INF);
    vector<bool> visited(n, false);
    dist[start] = 0.0;

    // 小根堆：greater 使隊首爲距離最小的元素
    // 允許同一節點多次入堆，靠 visited 跳過過期條目（惰性刪除）
    priority_queue<pair<double, int>, vector<pair<double, int>>, greater<pair<double, int>>> pq;
    pq.push({0.0, start});

    while (!pq.empty()) {
        auto [d, u] = pq.top();
        pq.pop();
        if (visited[u]) continue;  // 該節點已定型，堆裏這條是過期記錄
        visited[u] = true;
        for (auto [v, w] : graph[u]) {
            if (visited[v]) continue;
            double nd = d + w;
            if (nd < dist[v]) {  // 鬆弛
                dist[v] = nd;
                pq.push({nd, v});
            }
        }
    }
    return dist;
}

// 便捷接口：傳入無向邊列表 (u, v, w)，返回 start→target 的最短距離；不可達返回 -1
double shortestPath(int n, const vector<tuple<int, int, double>>& edges, int start, int target) {
    vector<vector<pair<int, double>>> graph(n);
    for (auto& e : edges) {
        int u = get<0>(e), v = get<1>(e);
        double w = get<2>(e);
        graph[u].push_back({v, w});
        graph[v].push_back({u, w});
    }
    double d = dijkstra(n, graph, start)[target];
    return isinf(d) ? -1.0 : d;
}

int main() {
    const int n = 5;
    vector<tuple<int, int, double>> edges = {
        {0, 1, 4}, {0, 2, 1},
        {2, 1, 2}, {2, 3, 5},
        {1, 3, 1}, {3, 4, 3},
    };
    // 0→2→1→3→4 = 1+2+1+3 = 7
    assert(fabs(shortestPath(n, edges, 0, 4) - 7.0) < 1e-9);
    // 0→2→1 = 1+2 = 3，比直連的 4 更短
    assert(fabs(shortestPath(n, edges, 0, 1) - 3.0) < 1e-9);
    assert(fabs(shortestPath(n, edges, 4, 4) - 0.0) < 1e-9);
    // 不連通時應返回 -1
    assert(shortestPath(2, {}, 0, 1) == -1.0);

    cout << "all tests passed" << endl;
    return 0;
}
