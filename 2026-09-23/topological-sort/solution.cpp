// 拓撲排序（Topological Sort，Kahn 算法 + DFS 逆後序）
// 編譯：g++ -std=c++17 -O2 solution.cpp -o solution && ./solution
//
// 題意：n 個點 m 條邊的有向圖（點編號 0..n-1），求一個順序使每條邊 u -> v 都滿足
//   u 排在 v 之前；有環則無解。
//
// 思路（兩種等價實現）：
//   1. Kahn 算法（BFS / 剝洋蔥）：不斷拿掉入度爲 0 的點，並把它指向的點的入度減一。
//      拿掉了全部 n 個點 -> 拿掉的先後順序即拓撲序；中途再也找不到入度 0 的點卻還有
//      剩餘 -> 剩下的點互相卡住，必定在環裏。用最小堆挑入度 0 的點，即得字典序最小解。
//   2. DFS 逆後序：一個點的所有後繼都訪問完之後才把它壓棧，最後倒序輸出。
//      判環靠三色標記：遞歸棧上的點是灰色，走到灰色點說明有回邊，即有環。
//
//   兩者差異：Kahn 是迭代的（無遞歸深度問題）且天然判環；DFS 版更短但要小心爆棧，
//   且得到的序一般不是字典序最小的。這裡兩種都用「升序鄰接表」驅動，保證結果可復現。
//
// 輸入：第一行 n m；接下來 m 行 u v（有向邊 u -> v）
// 輸出：一行，拓撲序（字典序最小，空格分隔）；若有環輸出 -1
// 無 stdin 輸入時運行內置斷言測試。
#include <algorithm>
#include <cassert>
#include <iostream>
#include <queue>
#include <vector>

using namespace std;

using Edge = pair<int, int>;

// 建鄰接表（邊的順序即輸入順序）
vector<vector<int>> buildAdj(int n, const vector<Edge>& edges) {
    vector<vector<int>> adj(n);
    for (const Edge& e : edges) adj[e.first].push_back(e.second);
    return adj;
}

// Kahn 算法 + 最小堆：返回 {是否有解, 拓撲序}，有解時是字典序最小的那個。
// 時間 O((n+m) log n)，空間 O(n + m)
pair<bool, vector<int>> kahnTopologicalSort(int n, const vector<Edge>& edges) {
    vector<vector<int>> adj = buildAdj(n, edges);
    vector<int> indeg(n, 0);
    for (const Edge& e : edges) indeg[e.second]++;

    // 最小堆：每步取編號最小的入度 0 點 -> 結果字典序最小
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

    // 拿掉的點不足 n 個 -> 剩下的點都還在環裏
    if (static_cast<int>(order.size()) != n) return {false, {}};
    return {true, order};
}

// DFS 逆後序：返回 {是否有解, 拓撲序}。迭代實現，避免深圖爆棧。時間 O(n + m)
pair<bool, vector<int>> dfsTopologicalSort(int n, const vector<Edge>& edges) {
    vector<vector<int>> adj = buildAdj(n, edges);
    for (auto& lst : adj) sort(lst.begin(), lst.end());  // 鄰接表升序，結果可復現

    const int WHITE = 0, GRAY = 1, BLACK = 2;
    vector<int> color(n, WHITE), post;
    post.reserve(n);

    vector<pair<int, size_t>> stack;  // (當前點, 下一條要走的邊下標)
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
                if (color[v] == GRAY) return {false, {}};   // 回邊 -> 有環
                if (color[v] == WHITE) {
                    color[v] = GRAY;
                    stack.push_back({v, 0});
                }
            } else {
                color[u] = BLACK;
                post.push_back(u);      // 後繼都已完成，本點才算完成
                stack.pop_back();
            }
        }
    }

    reverse(post.begin(), post.end());  // 逆後序即拓撲序
    return {true, post};
}

// 校驗：是 n 個點的一個排列，且每條邊 u -> v 都滿足 u 在 v 之前
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

// 對照用的全排列枚舉，返回字典序最小的拓撲序。僅用於 n 很小的測試
pair<bool, vector<int>> bruteLexTopologicalOrder(int n, const vector<Edge>& edges) {
    vector<int> perm(n);
    for (int i = 0; i < n; ++i) perm[i] = i;
    sort(perm.begin(), perm.end());
    do {  // next_permutation 按字典序枚舉，第一個合法的即字典序最小
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

    // 有環：0 -> 1 -> 2 -> 0
    {
        vector<Edge> cyc = {{0, 1}, {1, 2}, {2, 0}};
        assert(!kahnTopologicalSort(3, cyc).first);
        assert(!dfsTopologicalSort(3, cyc).first);
        assert(!bruteLexTopologicalOrder(3, cyc).first);
    }

    // 自環
    assert(!kahnTopologicalSort(1, {{0, 0}}).first);
    assert(!dfsTopologicalSort(1, {{0, 0}}).first);

    // 空圖：任意排列都是拓撲序。Kahn 用最小堆，得到字典序最小的 0 1 2 3；
    // DFS 版按 start = 0,1,2,... 依次完成，逆後序恰好是 3 2 1 0，同樣是合法拓撲序
    assert(kahnTopologicalSort(0, {}).first);
    assert(kahnTopologicalSort(0, {}).second.empty());
    assert(kahnTopologicalSort(4, {}).second == (vector<int>{0, 1, 2, 3}));
    assert(dfsTopologicalSort(4, {}).second == (vector<int>{3, 2, 1, 0}));
    assert(isValidTopologicalOrder(4, {}, dfsTopologicalSort(4, {}).second));

    // 單點無邊
    assert(kahnTopologicalSort(1, {}).second == (vector<int>{0}));

    // 重邊不應導致同一個點被重複輸出（入度按邊的條數計數，減到 0 才入堆一次）
    {
        vector<Edge> dup = {{0, 1}, {0, 1}};
        assert(kahnTopologicalSort(2, dup).second == (vector<int>{0, 1}));
        assert(dfsTopologicalSort(2, dup).second == (vector<int>{0, 1}));
    }

    // 鏈：1 -> 2 -> 3，孤立點 0 排最前
    {
        vector<Edge> chain = {{1, 2}, {2, 3}};
        assert(kahnTopologicalSort(4, chain).second == (vector<int>{0, 1, 2, 3}));
    }

    // 環 + 無環部分混合：只要有一個環就整體無解
    {
        vector<Edge> mix = {{0, 1}, {1, 2}, {2, 1}, {0, 3}};
        assert(!kahnTopologicalSort(4, mix).first);
        assert(!dfsTopologicalSort(4, mix).first);
    }

    // 與全排列暴力解隨機對拍：校驗 Kahn 的字典序最小性、DFS 結果的合法性、
    // 以及兩者「有解 / 無解」的判斷必須一致
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
        assert(k.first == b.first);                       // 有解/無解判斷一致
        if (k.first) {
            assert(k.second == b.second);                 // Kahn = 字典序最小
            assert(isValidTopologicalOrder(nn, edges, k.second));
            assert(isValidTopologicalOrder(nn, edges, d.second));  // DFS 結果也合法
        }
    }

    // 較大規模的鏈狀圖，驗證迭代版 DFS 不會爆棧
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
