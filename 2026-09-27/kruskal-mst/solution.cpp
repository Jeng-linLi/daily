// 最小生成樹（Kruskal 算法 + 併查集）
// 編譯：g++ -std=c++17 -O2 solution.cpp -o solution && ./solution
//
// 思路：貪心 + 併查集。所有邊按權值升序排列，依次嘗試加入；
//   若兩端點當前不連通就選中並合併，否則會成環就丟棄。選夠 n-1 條邊結束。
//
//   正確性（切分性質）：掃到邊 (u, v, w) 且 u、v 未連通時，把 u 所在連通塊看成
//   切分的一側，所有跨切分的邊裏 (u, v, w) 是當前最小者 —— 更小的邊都已考慮過，
//   它們要麼不跨這個切分，要麼會把點並進來（那樣 u、v 就已經連通了）。
//   跨切分的最小邊必屬於某棵 MST，所以選它不會錯。
//
//   排序必須**穩定**（同權邊保持輸入順序），這樣與 Python 版輸出逐字節一致。
//
// 輸入（空白分隔）：n m，隨後 m 行 `u v w`（頂點 0-based）
// 輸出：連通時第一行總權重，隨後每行一條選中邊 `u v w`；不連通時只輸出 disconnected
// 無 stdin 輸入時運行內置斷言測試。
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

// 併查集：路徑壓縮 + 按大小合併，同時維護連通塊個數
struct DSU {
    vector<int> parent, sz;
    int components;
    explicit DSU(int n) : parent(n), sz(n, 1), components(n) {
        for (int i = 0; i < n; ++i) parent[i] = i;
    }
    int find(int x) {
        while (parent[x] != x) {
            parent[x] = parent[parent[x]];     // 路徑壓縮（折半）
            x = parent[x];
        }
        return x;
    }
    bool unionSet(int a, int b) {
        int ra = find(a), rb = find(b);
        if (ra == rb) return false;            // 原本就連通，合併無效
        if (sz[ra] < sz[rb]) swap(ra, rb);     // 按大小合併，小樹掛到大樹下
        parent[rb] = ra;
        sz[ra] += sz[rb];
        --components;
        return true;
    }
};

// 返回 {總權重, 被選中的邊, 是否連通}。時間 O(m log m)，空間 O(n + m)
struct KruskalResult {
    long long total;
    vector<Edge> chosen;
    bool connected;
};

KruskalResult kruskal(int n, const vector<Edge>& edges) {
    DSU dsu(n);
    vector<Edge> sorted = edges;
    stable_sort(sorted.begin(), sorted.end(), [](const Edge& a, const Edge& b) {
        return a.w < b.w;                      // 只比權值；同權值靠 stable_sort 保持輸入順序
    });
    KruskalResult res{0, {}, true};
    for (const Edge& e : sorted) {
        if (dsu.unionSet(e.u, e.v)) {
            res.chosen.push_back(e);
            res.total += e.w;
            if (static_cast<int>(res.chosen.size()) == n - 1) break;  // 已是生成樹，提前結束
        }
    }
    res.connected = (n <= 1) || (static_cast<int>(res.chosen.size()) == n - 1);
    return res;
}

// 對照用的 Prim（鄰接表 + O(n^2) 選最小），用於與 Kruskal 交叉驗證總權重
pair<long long, bool> prim(int n, const vector<Edge>& edges) {
    if (n == 0) return {0, true};
    vector<vector<pair<int, long long>>> adj(n);
    for (const Edge& e : edges) {
        if (e.u == e.v) continue;              // 自環對 MST 無意義
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
        if (best == -1 || dist[best] == INF) return {total, false};  // 夠不着，圖不連通
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

// 對照用的指數級枚舉：枚舉所有邊子集，挑權和最小的生成樹。僅用於極小規模測試
pair<long long, bool> mstBrute(int n, const vector<Edge>& edges) {
    if (n <= 1) return {0, true};
    int m = static_cast<int>(edges.size());
    const long long INF = numeric_limits<long long>::max() / 4;
    long long best = INF;
    for (int mask = 0; mask < (1 << m); ++mask) {
        if (__builtin_popcount(static_cast<unsigned>(mask)) != n - 1) continue;  // 生成樹恰好 n-1 條邊
        DSU dsu(n);
        bool ok = true;
        long long total = 0;
        for (int i = 0; i < m; ++i) {
            if (!((mask >> i) & 1)) continue;
            if (dsu.unionSet(edges[i].u, edges[i].v)) {
                total += edges[i].w;
            } else {
                ok = false;                    // 成環，不是樹
                break;
            }
        }
        if (ok && dsu.components == 1 && total < best) best = total;
    }
    return {best, best != INF};
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
        for (int i = 0; i < m; ++i) {
            if (!(cin >> edges[i].u >> edges[i].v >> edges[i].w)) {
                edges[i] = {0, 0, 0};          // 輸入被截斷時用 0 兜底
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

    // README 示例：4 點 5 邊，MST = (0,1,1) + (1,2,2) + (2,3,3) = 6
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

    // 不連通：只有一條邊，第三個點孤立
    {
        vector<Edge> edges = {{0, 1, 5}};
        assert(!kruskal(3, edges).connected);
        assert(!prim(3, edges).second);
        assert(!mstBrute(3, edges).second);
    }

    // 退化情形
    {
        vector<Edge> none;
        assert(kruskal(0, none).total == 0);            // 空圖視爲連通，權重 0
        assert(kruskal(0, none).connected);
        assert(kruskal(1, none).total == 0);            // 單點，不需要邊
        assert(kruskal(1, none).connected);
        vector<Edge> one = {{0, 1, 7}};
        assert(kruskal(2, one).total == 7);
        assert(kruskal(2, one).chosen.size() == 1);
        assert(kruskal(2, one).connected);
        assert(!kruskal(2, none).connected);            // 兩點無邊，不連通
    }

    // 自環與重邊：自環必被丟棄，重邊只留一條
    {
        vector<Edge> e1 = {{0, 0, 1}, {0, 1, 3}, {0, 1, 3}};
        KruskalResult r1 = kruskal(2, e1);
        assert(r1.total == 3 && r1.chosen.size() == 1);
        vector<Edge> e2 = {{0, 1, 2}, {1, 2, 2}, {0, 2, 2}};
        assert(kruskal(3, e2).total == 4);
    }

    // 負權邊同樣成立
    {
        vector<Edge> e = {{0, 1, -5}, {1, 2, -1}, {0, 2, 10}};
        KruskalResult r = kruskal(3, e);
        assert(r.total == -6 && r.chosen.size() == 2);
    }

    // 同權邊按輸入順序穩定選中（保證兩語言輸出一致）
    {
        vector<Edge> same = {{2, 3, 1}, {0, 1, 1}, {1, 2, 1}};
        KruskalResult r = kruskal(4, same);
        assert(r.chosen.size() == 3);
        assert(r.chosen[0].u == 2 && r.chosen[0].v == 3);
        assert(r.chosen[1].u == 0 && r.chosen[1].v == 1);
        assert(r.chosen[2].u == 1 && r.chosen[2].v == 2);
    }

    LCG rng(20260927ULL);

    // 隨機對拍一：小規模圖上 Kruskal、Prim、指數級枚舉三者結果一致
    for (int t = 0; t < 150; ++t) {
        int n = rng.next(1, 6);
        int m = rng.next(0, 8);
        vector<pair<int, int>> cand;
        for (int u = 0; u < n; ++u)
            for (int v = u + 1; v < n; ++v) cand.push_back({u, v});
        // 用與 Python 一致的洗牌方式（Fisher-Yates，同樣由本 LCG 驅動）
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
            // 選出來的邊確實構成生成樹：n-1 條且把所有點連成一塊
            DSU dsu(n);
            for (const Edge& e : res.chosen) assert(dsu.unionSet(e.u, e.v));
            assert(dsu.components == 1);
        }
    }

    // 隨機對拍二：保證連通的隨機圖（先造鏈再補隨機邊），Kruskal 與 Prim 對拍
    for (int t = 0; t < 150; ++t) {
        int n = rng.next(2, 9);
        vector<Edge> edges;
        for (int i = 1; i < n; ++i) {                  // 先連成鏈，確保一定連通
            int j = rng.next(0, i - 1);
            edges.push_back({j, i, rng.next(1, 30)});
        }
        int extra = rng.next(0, 6);                    // 再補一些隨機邊
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
