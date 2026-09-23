// 拓扑排序（Topological Sort，Kahn 算法 + DFS 逆后序）
// 编译：g++ -std=c++17 -O2 solution.cpp -o solution && ./solution
//
// 题意：n 个点 m 条边的有向图（点编号 0..n-1），求一个顺序使每条边 u -> v 都满足
//   u 排在 v 之前；有环则无解。
//
// 思路（两种等价实现）：
//   1. Kahn 算法（BFS / 剥洋葱）：不断拿掉入度为 0 的点，并把它指向的点的入度减一。
//      拿掉了全部 n 个点 -> 拿掉的先后顺序即拓扑序；中途再也找不到入度 0 的点却还有
//      剩余 -> 剩下的点互相卡住，必定在环里。用最小堆挑入度 0 的点，即得字典序最小解。
//   2. DFS 逆后序：一个点的所有后继都访问完之后才把它压栈，最后倒序输出。
//      判环靠三色标记：递归栈上的点是灰色，走到灰色点说明有回边，即有环。
//
//   两者差异：Kahn 是迭代的（无递归深度问题）且天然判环；DFS 版更短但要小心爆栈，
//   且得到的序一般不是字典序最小的。这里两种都用「升序邻接表」驱动，保证结果可复现。
//
// 输入：第一行 n m；接下来 m 行 u v（有向边 u -> v）
// 输出：一行，拓扑序（字典序最小，空格分隔）；若有环输出 -1
// 无 stdin 输入时运行内置断言测试。
#include <algorithm>
#include <cassert>
#include <iostream>
#include <queue>
#include <vector>

using namespace std;

using Edge = pair<int, int>;

// 建邻接表（边的顺序即输入顺序）
vector<vector<int>> buildAdj(int n, const vector<Edge>& edges) {
    vector<vector<int>> adj(n);
    for (const Edge& e : edges) adj[e.first].push_back(e.second);
    return adj;
}

// Kahn 算法 + 最小堆：返回 {是否有解, 拓扑序}，有解时是字典序最小的那个。
// 时间 O((n+m) log n)，空间 O(n + m)
pair<bool, vector<int>> kahnTopologicalSort(int n, const vector<Edge>& edges) {
    vector<vector<int>> adj = buildAdj(n, edges);
    vector<int> indeg(n, 0);
    for (const Edge& e : edges) indeg[e.second]++;

    // 最小堆：每步取编号最小的入度 0 点 -> 结果字典序最小
    priority_queue<int, vector<int>, greater<int>> pq;
    for (int i = 0; i < n; ++i)
        if (indeg[i] == 0) pq.push(i);

    vector<int> order;
    order.reserve(n);
    while (!pq.empty()) {
        int u = pq.top();
        pq.pop();
        order.push_back(u);
        for (int v : adj[u]) {
            if (--indeg[v] == 0) pq.push(v);
        }
    }

    // 拿掉的点不足 n 个 -> 剩下的点都还在环里
    if (static_cast<int>(order.size()) != n) return {false, {}};
    return {true, order};
}

// DFS 逆后序：返回 {是否有解, 拓扑序}。迭代实现，避免深图爆栈。时间 O(n + m)
pair<bool, vector<int>> dfsTopologicalSort(int n, const vector<Edge>& edges) {
    vector<vector<int>> adj = buildAdj(n, edges);
    for (auto& lst : adj) sort(lst.begin(), lst.end());  // 邻接表升序，结果可复现

    const int WHITE = 0, GRAY = 1, BLACK = 2;
    vector<int> color(n, WHITE), post;
    post.reserve(n);

    vector<pair<int, size_t>> stack;  // (当前点, 下一条要走的边下标)
    for (int start = 0; start < n; ++start) {
        if (color[start] != WHITE) continue;
        color[start] = GRAY;
        stack.push_back({start, 0});
        while (!stack.empty()) {
            int u = stack.back().first;
            size_t idx = stack.back().second;
            if (idx < adj[u].size()) {
                stack.back().second = idx + 1;
                int v = adj[u][idx];
                if (color[v] == GRAY) return {false, {}};   // 回边 -> 有环
                if (color[v] == WHITE) {
                    color[v] = GRAY;
                    stack.push_back({v, 0});
                }
            } else {
                color[u] = BLACK;
                post.push_back(u);      // 后继都已完成，本点才算完成
                stack.pop_back();
            }
        }
    }

    reverse(post.begin(), post.end());  // 逆后序即拓扑序
    return {true, post};
}

// 校验：是 n 个点的一个排列，且每条边 u -> v 都满足 u 在 v 之前
bool isValidTopologicalOrder(int n, const vector<Edge>& edges, const vector<int>& order) {
    if (static_cast<int>(order.size()) != n) return false;
    vector<int> pos(n, -1);
    for (int i = 0; i < n; ++i) {
        if (order[i] < 0 || order[i] >= n || pos[order[i]] != -1) return false;
        pos[order[i]] = i;
    }
    for (const Edge& e : edges)
        if (pos[e.first] > pos[e.second]) return false;
    return true;
}

// 对照用的全排列枚举，返回字典序最小的拓扑序。仅用于 n 很小的测试
pair<bool, vector<int>> bruteLexTopologicalOrder(int n, const vector<Edge>& edges) {
    vector<int> perm(n);
    for (int i = 0; i < n; ++i) perm[i] = i;
    sort(perm.begin(), perm.end());
    do {  // next_permutation 按字典序枚举，第一个合法的即字典序最小
        vector<int> pos(n, -1);
        for (int i = 0; i < n; ++i) pos[perm[i]] = i;
        bool ok = true;
        for (const Edge& e : edges) {
            if (pos[e.first] > pos[e.second]) { ok = false; break; }
        }
        if (ok) return {true, perm};
    } while (next_permutation(perm.begin(), perm.end()));
    return {false, {}};
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
        for (int i = 0; i < m; ++i) cin >> edges[i].first >> edges[i].second;
        auto res = kahnTopologicalSort(n, edges);
        if (!res.first) {
            cout << -1 << "\n";
        } else {
            for (size_t i = 0; i < res.second.size(); ++i) {
                if (i) cout << " ";
                cout << res.second[i];
            }
            cout << "\n";
        }
        return 0;
    }

    // README 中的示例：0 1 2 3 5 4
    {
        vector<Edge> edges = {{0, 1}, {0, 2}, {1, 3}, {2, 3}, {2, 4}, {5, 4}};
        auto k = kahnTopologicalSort(6, edges);
        assert(k.first);
        assert(k.second == (vector<int>{0, 1, 2, 3, 5, 4}));
        auto b = bruteLexTopologicalOrder(6, edges);
        assert(b.first && b.second == k.second);
        auto d = dfsTopologicalSort(6, edges);
        assert(d.first && isValidTopologicalOrder(6, edges, d.second));
    }

    // 有环：0 -> 1 -> 2 -> 0
    {
        vector<Edge> cyc = {{0, 1}, {1, 2}, {2, 0}};
        assert(!kahnTopologicalSort(3, cyc).first);
        assert(!dfsTopologicalSort(3, cyc).first);
        assert(!bruteLexTopologicalOrder(3, cyc).first);
    }

    // 自环
    assert(!kahnTopologicalSort(1, {{0, 0}}).first);
    assert(!dfsTopologicalSort(1, {{0, 0}}).first);

    // 空图：任意排列都是拓扑序。Kahn 用最小堆，得到字典序最小的 0 1 2 3；
    // DFS 版按 start = 0,1,2,... 依次完成，逆后序恰好是 3 2 1 0，同样是合法拓扑序
    assert(kahnTopologicalSort(0, {}).first);
    assert(kahnTopologicalSort(0, {}).second.empty());
    assert(kahnTopologicalSort(4, {}).second == (vector<int>{0, 1, 2, 3}));
    assert(dfsTopologicalSort(4, {}).second == (vector<int>{3, 2, 1, 0}));
    assert(isValidTopologicalOrder(4, {}, dfsTopologicalSort(4, {}).second));

    // 单点无边
    assert(kahnTopologicalSort(1, {}).second == (vector<int>{0}));

    // 重边不应导致同一个点被重复输出（入度按边的条数计数，减到 0 才入堆一次）
    {
        vector<Edge> dup = {{0, 1}, {0, 1}};
        assert(kahnTopologicalSort(2, dup).second == (vector<int>{0, 1}));
        assert(dfsTopologicalSort(2, dup).second == (vector<int>{0, 1}));
    }

    // 链：1 -> 2 -> 3，孤立点 0 排最前
    {
        vector<Edge> chain = {{1, 2}, {2, 3}};
        assert(kahnTopologicalSort(4, chain).second == (vector<int>{0, 1, 2, 3}));
    }

    // 环 + 无环部分混合：只要有一个环就整体无解
    {
        vector<Edge> mix = {{0, 1}, {1, 2}, {2, 1}, {0, 3}};
        assert(!kahnTopologicalSort(4, mix).first);
        assert(!dfsTopologicalSort(4, mix).first);
    }

    // 与全排列暴力解随机对拍：校验 Kahn 的字典序最小性、DFS 结果的合法性、
    // 以及两者「有解 / 无解」的判断必须一致
    LCG rng(20260923ULL);
    for (int t = 0; t < 200; ++t) {
        int nn = rng.next(1, 6);
        vector<Edge> edges;
        for (int u = 0; u < nn; ++u)
            for (int v = 0; v < nn; ++v)
                if (u != v && rng.next(0, 99) < 25) edges.push_back({u, v});

        auto k = kahnTopologicalSort(nn, edges);
        auto d = dfsTopologicalSort(nn, edges);
        auto b = bruteLexTopologicalOrder(nn, edges);
        assert(k.first == d.first);
        assert(k.first == b.first);                       // 有解/无解判断一致
        if (k.first) {
            assert(k.second == b.second);                 // Kahn = 字典序最小
            assert(isValidTopologicalOrder(nn, edges, k.second));
            assert(isValidTopologicalOrder(nn, edges, d.second));  // DFS 结果也合法
        }
    }

    // 较大规模的链状图，验证迭代版 DFS 不会爆栈
    {
        const int bigN = 20000;
        vector<Edge> bigEdges;
        bigEdges.reserve(bigN - 1);
        for (int i = 0; i + 1 < bigN; ++i) bigEdges.push_back({i, i + 1});
        auto k = kahnTopologicalSort(bigN, bigEdges);
        auto d = dfsTopologicalSort(bigN, bigEdges);
        assert(k.first && d.first);
        vector<int> expect(bigN);
        for (int i = 0; i < bigN; ++i) expect[i] = i;
        assert(k.second == expect);
        assert(d.second == expect);
    }

    cout << "all tests passed" << endl;
    return 0;
}
