// Tarjan 強連通分量：SCC 分解 / 縮點 DAG / 使全圖強連通的最少加邊數
// 編譯：g++ -std=c++17 -O2 -Wall solution.cpp -o solution && ./solution
//
// 思路：
//   強連通分量把有向圖「壓縮」成一棵 DAG：分量內互相可達，分量之間單向連通。
//   Tarjan 做一次 DFS，給每個點記 dfn（訪問時間戳）與 low
//   （沿樹邊往下走、最多再走一條回邊/橫叉邊，能到達的最小 dfn），並把訪問到的點壓棧。
//   當 low[v] == dfn[v] 時，說明 v 的子樹裡沒有任何邊連回 v 的祖先，
//   於是以 v 為根的子樹中「還留在棧裡」的點剛好構成一個 SCC，全部彈出。
//   每個點進棧出棧各一次，複雜度 O(n + m)。
//
//   與 Kosaraju（兩次 DFS + 反圖）相比，Tarjan 只要一次 DFS、不需要反圖，常數更小。
//   這裏把 Tarjan 寫成**迭代版**（顯式棧），避免圖退化成鏈時遞迴深度過深；
//   C++ 與 Python 版用同一套迭代寫法，行為逐字節一致。
//
//   第 5 問的結論：縮點 DAG 上有 src 個入度 0 的點、snk 個出度 0 的點，
//   則最少需要加 max(src, snk) 條邊（k = 1 時為 0）——每次加邊都能把一個源點與一個匯點接上。
//
// 輸入（空白分隔）：n m / u1 v1 / u2 v2 / …
// 輸出：分量數 k / 各點所屬分量編號 / 縮點 DAG 邊數 / 最大分量點數 / 源點數 匯點數 最少加邊數

#include <algorithm>
#include <cassert>
#include <iostream>
#include <numeric>
#include <queue>
#include <random>
#include <set>
#include <sstream>
#include <string>
#include <vector>

using namespace std;

// Tarjan（迭代版）：回傳分量數 k，comp[v] 為點 v（1 … n）所屬的原始分量編號
static int tarjanScc(int n, const vector<vector<int>>& adj, vector<int>& comp) {
    vector<int> dfn(n + 1, 0), low(n + 1, 0);
    vector<bool> onStack(n + 1, false);
    comp.assign(n + 1, -1);
    vector<int> stk;                 // Tarjan 的「待彈出」棧
    int timer = 0, ncomp = 0;

    for (int root = 1; root <= n; ++root) {
        if (dfn[root]) continue;
        // 顯式棧元素 = (當前節點, 下一條要處理的邊的下標)
        vector<pair<int, int>> work;
        work.push_back({root, 0});
        while (!work.empty()) {
            int v = work.back().first;
            // 注意：不要持有 work.back() 的引用——push_back 可能讓 vector 重新配置記憶體
            if (work.back().second == 0) {   // 第一次進入 v
                ++timer;
                dfn[v] = low[v] = timer;
                stk.push_back(v);
                onStack[v] = true;
            }
            bool descended = false;
            while (work.back().second < (int)adj[v].size()) {
                int w = adj[v][work.back().second];
                work.back().second++;
                if (dfn[w] == 0) {    // 樹邊：先下去處理 w
                    work.push_back({w, 0});
                    descended = true;
                    break;
                } else if (onStack[w]) {
                    low[v] = min(low[v], dfn[w]);
                }
            }
            if (descended) continue;
            // v 的所有邊都處理完了
            if (low[v] == dfn[v]) {   // v 是一個 SCC 的根，彈出整個分量
                while (true) {
                    int x = stk.back();
                    stk.pop_back();
                    onStack[x] = false;
                    comp[x] = ncomp;
                    if (x == v) break;
                }
                ++ncomp;
            }
            work.pop_back();
            if (!work.empty()) {      // 把 low[v] 回傳給父節點
                int u = work.back().first;
                low[u] = min(low[u], low[v]);
            }
        }
    }
    return ncomp;
}

// 重新編號：按「分量內最小點號」遞增，讓輸出與 DFS 順序無關。回傳長度 n+1 的陣列（下標 0 未使用）
static vector<int> relabel(const vector<int>& comp, int k, int n) {
    vector<int> first(k, n + 1);
    for (int v = 1; v <= n; ++v) first[comp[v]] = min(first[comp[v]], v);
    vector<int> order(k);
    iota(order.begin(), order.end(), 0);
    sort(order.begin(), order.end(), [&](int a, int b) { return first[a] < first[b]; });
    vector<int> newId(k, -1);
    for (int i = 0; i < k; ++i) newId[order[i]] = i;
    vector<int> res(n + 1, -1);
    for (int v = 1; v <= n; ++v) res[v] = newId[comp[v]];
    return res;
}

// Kosaraju：第二次 DFS 在反圖上做，測試裡當 Tarjan 的對拍基準
static int kosarajuScc(int n, const vector<vector<int>>& adj, vector<int>& comp) {
    vector<vector<int>> radj(n + 1);
    for (int v = 1; v <= n; ++v)
        for (int w : adj[v]) radj[w].push_back(v);

    vector<int> order;
    vector<bool> vis(n + 1, false);
    for (int s = 1; s <= n; ++s) {
        if (vis[s]) continue;
        vector<pair<int, int>> st;
        st.push_back({s, 0});
        vis[s] = true;
        while (!st.empty()) {
            int v = st.back().first;
            if (st.back().second < (int)adj[v].size()) {
                int w = adj[v][st.back().second];
                st.back().second++;
                if (!vis[w]) {
                    vis[w] = true;
                    st.push_back({w, 0});
                }
            } else {
                order.push_back(v);
                st.pop_back();
            }
        }
    }

    comp.assign(n + 1, -1);
    int k = 0;
    for (int idx = (int)order.size() - 1; idx >= 0; --idx) {
        int s = order[idx];
        if (comp[s] != -1) continue;
        vector<int> st{s};
        comp[s] = k;
        while (!st.empty()) {
            int v = st.back();
            st.pop_back();
            for (int w : radj[v])
                if (comp[w] == -1) {
                    comp[w] = k;
                    st.push_back(w);
                }
        }
        ++k;
    }
    return k;
}

// 暴力法：對每個點做 BFS，兩點互相可達則同屬一個分量（O(n·(n+m))，只用於測試）
static vector<int> bruteForceScc(int n, const vector<vector<int>>& adj) {
    vector<vector<bool>> reach(n + 1, vector<bool>(n + 1, false));
    for (int s = 1; s <= n; ++s) {
        reach[s][s] = true;
        queue<int> q;
        q.push(s);
        while (!q.empty()) {
            int v = q.front();
            q.pop();
            for (int w : adj[v])
                if (!reach[s][w]) {
                    reach[s][w] = true;
                    q.push(w);
                }
        }
    }
    vector<int> comp(n + 1, -1);
    int k = 0;
    for (int v = 1; v <= n; ++v) {
        if (comp[v] != -1) continue;
        for (int w = v; w <= n; ++w)
            if (reach[v][w] && reach[w][v]) comp[w] = k;
        ++k;
    }
    return comp;
}

// 縮點後 DAG 的統計：回傳 {邊數, 最大分量大小, 源點數, 匯點數}
static vector<int> condensationInfo(int n, const vector<vector<int>>& adj,
                                    const vector<int>& comp, int k) {
    vector<int> size(k, 0);
    for (int v = 1; v <= n; ++v) size[comp[v]]++;
    set<pair<int, int>> edges;
    vector<int> inDeg(k, 0), outDeg(k, 0);
    for (int v = 1; v <= n; ++v)
        for (int w : adj[v]) {
            int a = comp[v], b = comp[w];
            if (a != b && edges.insert({a, b}).second) {
                outDeg[a]++;
                inDeg[b]++;
            }
        }
    int src = 0, snk = 0;
    for (int i = 0; i < k; ++i) {
        if (inDeg[i] == 0) src++;
        if (outDeg[i] == 0) snk++;
    }
    int biggest = 0;
    for (int i = 0; i < k; ++i) biggest = max(biggest, size[i]);
    return {(int)edges.size(), biggest, src, snk};
}

// 依題目格式算出五行輸出
static vector<string> solve(int n, int m, const vector<pair<int, int>>& edges) {
    vector<vector<int>> adj(n + 1);
    for (auto& e : edges) adj[e.first].push_back(e.second);

    vector<int> comp;
    int k = tarjanScc(n, adj, comp);
    comp = relabel(comp, k, n);                 // 重編號，保證輸出唯一
    vector<int> info = condensationInfo(n, adj, comp, k);
    int add = (k <= 1) ? 0 : max(info[2], info[3]);

    ostringstream line2;
    for (int v = 1; v <= n; ++v) {
        if (v > 1) line2 << ' ';
        line2 << comp[v];
    }
    return {to_string(k), line2.str(), to_string(info[0]), to_string(info[1]),
            to_string(info[2]) + " " + to_string(info[3]) + " " + to_string(add)};
}

// ---------- 輸入輸出 ----------

static void runIo(const string& data) {
    istringstream iss(data);
    vector<int> tk;
    int x;
    while (iss >> x) tk.push_back(x);
    int n = tk.empty() ? 0 : tk[0];
    int m = tk.size() >= 2 ? tk[1] : 0;
    vector<pair<int, int>> edges;
    for (int i = 0; i < m && 2 + 2 * i + 1 < (int)tk.size(); ++i) {
        int u = tk[2 + 2 * i], v = tk[2 + 2 * i + 1];
        if (1 <= u && u <= n && 1 <= v && v <= n) edges.push_back({u, v});
    }
    for (const string& line : solve(n, m, edges)) cout << line << '\n';
}

// ---------- 內置測試 ----------

static vector<vector<int>> buildAdj(int n, const vector<pair<int, int>>& edges) {
    vector<vector<int>> adj(n + 1);
    for (auto& e : edges) adj[e.first].push_back(e.second);
    return adj;
}

static void runTests() {
    vector<pair<int, int>> e1 = {{1, 2}, {2, 3}, {3, 1}, {3, 5}};
    int n1 = 5;
    vector<vector<int>> adj = buildAdj(n1, e1);

    vector<int> comp;
    int k = tarjanScc(n1, adj, comp);
    assert(k == 3);
    vector<int> compK;
    int kk = kosarajuScc(n1, adj, compK);
    assert(kk == k);

    vector<int> bf = bruteForceScc(n1, adj);       // 與暴力法對拍：同分量 ⟺ 互相可達
    for (int v = 1; v <= n1; ++v)
        for (int w = 1; w <= n1; ++w) {
            assert((comp[v] == comp[w]) == (bf[v] == bf[w]));
            assert((compK[v] == compK[w]) == (bf[v] == bf[w]));
        }

    vector<int> lab = relabel(comp, k, n1);
    assert((vector<int>(lab.begin() + 1, lab.end()) == vector<int>{0, 0, 0, 1, 2}));
    vector<int> info = condensationInfo(n1, adj, lab, k);
    assert((info == vector<int>{1, 3, 2, 2}));     // 只有 0→2 一條縮點邊；源點 {0,1}，匯點 {1,2}
    assert(max(info[2], info[3]) == 2);

    assert((solve(0, 0, {}) == vector<string>{"0", "", "0", "0", "0 0 0"}));
    assert((solve(1, 0, {}) == vector<string>{"1", "0", "0", "1", "1 1 0"}));
    assert((solve(3, 3, {{1, 2}, {2, 3}, {3, 1}}) == vector<string>{"1", "0 0 0", "0", "3", "1 1 0"}));
    assert((solve(4, 3, {{1, 2}, {2, 3}, {3, 4}}) == vector<string>{"4", "0 1 2 3", "3", "1", "1 1 1"}));

    mt19937 rng(20260930);
    for (int t = 0; t < 500; ++t) {
        int n = (int)(rng() % 9);
        int m = (int)(rng() % 13);
        vector<pair<int, int>> edges;
        for (int i = 0; i < m && n > 0; ++i)
            edges.push_back({(int)(rng() % n) + 1, (int)(rng() % n) + 1});
        adj = buildAdj(n, edges);

        int k1 = tarjanScc(n, adj, comp);
        int k2 = kosarajuScc(n, adj, compK);
        assert(k1 == k2);

        bf = bruteForceScc(n, adj);
        for (int v = 1; v <= n; ++v)
            for (int w = 1; w <= n; ++w) {
                assert((comp[v] == comp[w]) == (bf[v] == bf[w]));
                assert((compK[v] == compK[w]) == (bf[v] == bf[w]));
            }

        lab = relabel(comp, k1, n);
        vector<int> labK = relabel(compK, k2, n);
        assert(lab == labK);                       // 重編號後兩個算法應給出同一個陣列
        vector<int> first(k1, n + 1);
        for (int v = 1; v <= n; ++v) first[lab[v]] = min(first[lab[v]], v);
        for (int i = 1; i < k1; ++i) assert(first[i - 1] < first[i]);   // 按最小點號遞增

        info = condensationInfo(n, adj, lab, k1);
        int biggest = 0;
        for (int i = 0; i < k1; ++i) {
            int c = 0;
            for (int v = 1; v <= n; ++v)
                if (lab[v] == i) c++;
            biggest = max(biggest, c);
        }
        assert(info[1] == biggest);
        assert(0 <= info[0] && info[0] <= k1 * (k1 - 1));
        if (k1 > 0) assert(info[2] >= 1 && info[3] >= 1);   // DAG 必有源點與匯點
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
