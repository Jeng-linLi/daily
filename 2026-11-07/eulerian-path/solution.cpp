// 歐拉路徑 / 歐拉迴路（Hierholzer 演算法）
//
// 題意：給定一張圖（可帶自環與重邊），求一條經過每一條邊恰好一次的
//       走訪（歐拉路徑）；若起點與終點相同則為歐拉迴路。不存在輸出 -1。
//
// 判定條件：
//   無向圖：度數為奇數的點恰為 0 個（迴路）或 2 個（路徑），
//           且所有度數 > 0 的點在同一連通分量內。
//   有向圖：每點 out-in ∈ {-1, 0, 1}，且 +1 與 -1 的點各至多一個，
//           且所有 out+in > 0 的點弱連通。
//
// IO 模式（stdin 有內容時）：
//   T           （1 = 無向圖，2 = 有向圖）
//   n m
//   m 行：u v    （0-based，允許自環與重邊）
//   存在則輸出一行 m+1 個點，否則輸出 -1。
// 無 stdin 時跑內建斷言（與指數回溯暴力對拍存在性）並輸出 all tests passed。
#include <cstdio>
#include <cstring>
#include <cctype>
#include <sstream>
#include <iostream>
#include <string>
#include <vector>

using namespace std;
typedef long long ll;

static vector<int> solveUndirected(int n, const vector<pair<int,int>>& edges) {
    int m = (int)edges.size();
    if (m == 0) return vector<int>{0};
    vector<int> deg(n, 0);
    vector<vector<pair<int,int> > > adj(n); // (對方點, 弧編號 2i / 2i+1)
    for (int i = 0; i < m; i++) {
        int u = edges[i].first, v = edges[i].second;
        adj[u].push_back(make_pair(v, 2*i));
        adj[v].push_back(make_pair(u, 2*i+1));
        deg[u]++; deg[v]++;
    }
    int firstOdd = -1, oddCnt = 0;
    for (int v = 0; v < n; v++) if (deg[v] & 1) { oddCnt++; if (firstOdd < 0) firstOdd = v; }
    if (oddCnt != 0 && oddCnt != 2) return vector<int>();
    int actFirst = -1, actCnt = 0;
    for (int v = 0; v < n; v++) if (deg[v] > 0) { actCnt++; if (actFirst < 0) actFirst = v; }
    vector<char> seen(n, 0);
    vector<int> st; st.push_back(actFirst); seen[actFirst] = 1;
    int cntSeen = 0;
    while (!st.empty()) {
        int u = st.back(); st.pop_back(); cntSeen++;
        for (size_t j = 0; j < adj[u].size(); j++)
            if (!seen[adj[u][j].first]) { seen[adj[u][j].first] = 1; st.push_back(adj[u][j].first); }
    }
    if (cntSeen != actCnt) return vector<int>();
    int start = (oddCnt == 2) ? firstOdd : actFirst;
    vector<char> used(m, 0);
    vector<size_t> it(n, 0);
    vector<int> stk; stk.push_back(start);
    vector<int> path;
    while (!stk.empty()) {
        int u = stk.back();
        bool advanced = false;
        while (it[u] < adj[u].size()) {
            int v = adj[u][it[u]].first;
            int arc = adj[u][it[u]].second;
            it[u]++;
            if (!used[arc >> 1]) {
                used[arc >> 1] = 1;
                stk.push_back(v);
                advanced = true;
                break;
            }
        }
        if (!advanced) { path.push_back(stk.back()); stk.pop_back(); }
    }
    if ((int)path.size() != m + 1) return vector<int>();
    for (int i = 0; i < (int)path.size() / 2; i++)
        swap(path[i], path[path.size() - 1 - i]);
    return path;
}

static vector<int> solveDirected(int n, const vector<pair<int,int>>& edges) {
    int m = (int)edges.size();
    if (m == 0) return vector<int>{0};
    vector<int> outdeg(n, 0), indeg(n, 0);
    vector<vector<pair<int,int> > > adj(n);
    vector<vector<int> > radj(n);
    for (int i = 0; i < m; i++) {
        int u = edges[i].first, v = edges[i].second;
        adj[u].push_back(make_pair(v, i));
        radj[v].push_back(u);
        outdeg[u]++; indeg[v]++;
    }
    int plus = 0, minus = 0;
    for (int v = 0; v < n; v++) {
        int b = outdeg[v] - indeg[v];
        if (b < -1 || b > 1) return vector<int>();
        if (b == 1) plus++;
        if (b == -1) minus++;
    }
    if (plus != minus) return vector<int>();
    int actFirst = -1, actCnt = 0;
    for (int v = 0; v < n; v++) if (outdeg[v] + indeg[v] > 0) { actCnt++; if (actFirst < 0) actFirst = v; }
    vector<char> seen(n, 0);
    vector<int> st; st.push_back(actFirst); seen[actFirst] = 1;
    int cntSeen = 0;
    while (!st.empty()) {
        int u = st.back(); st.pop_back(); cntSeen++;
        for (size_t j = 0; j < adj[u].size(); j++)
            if (!seen[adj[u][j].first]) { seen[adj[u][j].first] = 1; st.push_back(adj[u][j].first); }
        for (size_t j = 0; j < radj[u].size(); j++)
            if (!seen[radj[u][j]]) { seen[radj[u][j]] = 1; st.push_back(radj[u][j]); }
    }
    if (cntSeen != actCnt) return vector<int>();
    int start = -1;
    for (int v = 0; v < n; v++) if (outdeg[v] - indeg[v] == 1) { start = v; break; }
    if (start < 0) start = actFirst;
    vector<char> used(m, 0);
    vector<size_t> it(n, 0);
    vector<int> stk; stk.push_back(start);
    vector<int> path;
    while (!stk.empty()) {
        int u = stk.back();
        bool advanced = false;
        while (it[u] < adj[u].size()) {
            int v = adj[u][it[u]].first;
            int eid = adj[u][it[u]].second;
            it[u]++;
            if (!used[eid]) { used[eid] = 1; stk.push_back(v); advanced = true; break; }
        }
        if (!advanced) { path.push_back(stk.back()); stk.pop_back(); }
    }
    if ((int)path.size() != m + 1) return vector<int>();
    for (int i = 0; i < (int)path.size() / 2; i++)
        swap(path[i], path[path.size() - 1 - i]);
    return path;
}

// ---------- 暴力對拍（指數回溯，僅用於小規模） ----------

struct BruteState { vector<char> used; };

static bool bruteUnd(const vector<pair<int,int> >& edges, BruteState& bs,
                     int u, int cnt, vector<int>& out) {
    int m = (int)edges.size();
    if (cnt == m) { out.push_back(u); return true; }
    for (int i = 0; i < m; i++) {
        if (bs.used[i]) continue;
        int nxt = -1;
        if (edges[i].first == u) nxt = edges[i].second;
        else if (edges[i].second == u) nxt = edges[i].first;
        if (nxt < 0) continue;
        bs.used[i] = 1;
        if (bruteUnd(edges, bs, nxt, cnt + 1, out)) { out.insert(out.begin(), u); return true; }
        bs.used[i] = 0;
    }
    return false;
}

static vector<int> bruteUndirected(int n, const vector<pair<int,int> >& edges) {
    int m = (int)edges.size();
    if (m == 0) return vector<int>{0};
    vector<int> deg(n, 0);
    for (size_t i = 0; i < edges.size(); i++) { deg[edges[i].first]++; deg[edges[i].second]++; }
    int oddCnt = 0, firstOdd = -1;
    for (int v = 0; v < n; v++) if (deg[v] & 1) { oddCnt++; if (firstOdd < 0) firstOdd = v; }
    if (oddCnt != 0 && oddCnt != 2) return vector<int>();
    BruteState bs; bs.used.assign(m, 0);
    vector<int> starts;
    if (oddCnt == 2) starts.push_back(firstOdd);
    else for (int v = 0; v < n; v++) if (deg[v] > 0) starts.push_back(v);
    vector<int> out;
    for (size_t s = 0; s < starts.size(); s++) {
        out.clear();
        if (bruteUnd(edges, bs, starts[s], 0, out)) return out;
    }
    return vector<int>();
}

static bool bruteDir(const vector<pair<int,int> >& edges, BruteState& bs,
                     int u, int cnt, vector<int>& out) {
    int m = (int)edges.size();
    if (cnt == m) { out.push_back(u); return true; }
    for (int i = 0; i < m; i++) {
        if (!bs.used[i] && edges[i].first == u) {
            bs.used[i] = 1;
            if (bruteDir(edges, bs, edges[i].second, cnt + 1, out)) { out.insert(out.begin(), u); return true; }
            bs.used[i] = 0;
        }
    }
    return false;
}

static vector<int> bruteDirected(int n, const vector<pair<int,int> >& edges) {
    int m = (int)edges.size();
    if (m == 0) return vector<int>{0};
    vector<int> outdeg(n, 0), indeg(n, 0);
    for (size_t i = 0; i < edges.size(); i++) { outdeg[edges[i].first]++; indeg[edges[i].second]++; }
    int plus = 0, minus = 0;
    for (int v = 0; v < n; v++) {
        int b = outdeg[v] - indeg[v];
        if (b < -1 || b > 1) return vector<int>();
        if (b == 1) plus++;
        if (b == -1) minus++;
    }
    if (plus != minus) return vector<int>();
    vector<int> starts;
    for (int v = 0; v < n; v++) if (outdeg[v] - indeg[v] == 1) starts.push_back(v);
    if (starts.empty()) for (int v = 0; v < n; v++) if (outdeg[v] > 0) starts.push_back(v);
    BruteState bs; bs.used.assign(m, 0);
    vector<int> out;
    for (size_t s = 0; s < starts.size(); s++) {
        out.clear();
        if (bruteDir(edges, bs, starts[s], 0, out)) return out;
    }
    return vector<int>();
}

static bool validatePath(int m, const vector<pair<int,int> >& edges,
                         const vector<int>& path, bool directed) {
    if (path.empty()) return true;
    if ((int)path.size() != m + 1) return false;
    vector<pair<pair<int,int>,int> > cnt;
    for (int i = 0; i < m; i++) {
        pair<int,int> key = directed ? edges[i]
            : make_pair(min(edges[i].first, edges[i].second), max(edges[i].first, edges[i].second));
        bool found = false;
        for (size_t j = 0; j < cnt.size(); j++)
            if (cnt[j].first == key) { cnt[j].second++; found = true; break; }
        if (!found) cnt.push_back(make_pair(key, 1));
    }
    for (int i = 0; i < m; i++) {
        int a = path[i], b = path[i + 1];
        pair<int,int> key = directed ? make_pair(a, b) : make_pair(min(a, b), max(a, b));
        bool okStep = false;
        for (size_t j = 0; j < cnt.size(); j++)
            if (cnt[j].first == key && cnt[j].second > 0) { cnt[j].second--; okStep = true; break; }
        if (!okStep) return false;
    }
    for (size_t j = 0; j < cnt.size(); j++) if (cnt[j].second != 0) return false;
    return true;
}

// 固定用例（與 Python run_tests 相同的集合）
static const int FIXED[][3] = {
    {1, 3, 3},   // (t, n, m)
    {1, 3, 2},
    {1, 4, 3},
    {1, 5, 6},
    {1, 4, 4},
    {1, 2, 2},
    {1, 2, 3},
    {2, 3, 3},
    {2, 3, 3},
    {2, 4, 3},
    {2, 3, 3},
};
static const int FIXED_EDGES[][2] = {
    {0,1},{1,2},{2,0},          // 用例 1：三角形迴路
    {0,1},{1,2},                // 用例 2：鏈
    {0,1},{0,2},{0,3},          // 用例 3：三叉星（中心 0，度 3 -> 路徑）
    {0,1},{1,2},{2,0},{2,3},{3,4},{4,2},  // 用例 4：8 字迴路
    {0,1},{0,2},{0,3},{3,1},    // 用例 5：四奇點無解
    {0,0},{0,1},                // 用例 6：自環 + 邊
    {0,1},{0,1},{0,1},          // 用例 7：重邊
    {0,1},{1,2},{2,0},          // 用例 8：有向三角形迴路
    {0,1},{1,2},{0,2},          // 用例 9：有向路徑
    {0,1},{1,2},{3,2},          // 用例 10：bal 檢查失敗
    {0,1},{1,0},{1,2},          // 用例 11：有向路徑
};

static void runAllTests() {
    unsigned long long seed = 20261107ULL;

    // 固定用例：存在性與暴力一致 + 路徑合法
    int nf = sizeof(FIXED) / sizeof(FIXED[0]);
    int ebase = 0;
    for (int c = 0; c < nf; c++) {
        int t = FIXED[c][0], n = FIXED[c][1], m = FIXED[c][2];
        vector<pair<int,int> > es;
        for (int i = 0; i < m; i++)
            es.push_back(make_pair(FIXED_EDGES[ebase + i][0], FIXED_EDGES[ebase + i][1]));
        ebase += m;
        vector<int> got = (t == 1) ? solveUndirected(n, es) : solveDirected(n, es);
        vector<int> want = (t == 1) ? bruteUndirected(n, es) : bruteDirected(n, es);
        if (got.empty() != want.empty()) { printf("FAIL fixed existence case %d\n", c); return; }
        if (!validatePath(m, es, got, t == 2)) { printf("FAIL fixed path case %d\n", c); return; }
    }

    // 隨機用例：n, m <= 8，與回溯暴力對拍存在性
    for (int iter = 0; iter < 300; iter++) {
        seed = (seed * 1103515245ULL + 12345ULL) & 0x7FFFFFFFULL;
        int t = (int)(seed % 2) + 1;
        seed = (seed * 1103515245ULL + 12345ULL) & 0x7FFFFFFFULL;
        int n = (int)(seed % 7) + 1;
        seed = (seed * 1103515245ULL + 12345ULL) & 0x7FFFFFFFULL;
        int m = (int)(seed % 9);
        vector<pair<int,int> > es;
        for (int i = 0; i < m; i++) {
            seed = (seed * 1103515245ULL + 12345ULL) & 0x7FFFFFFFULL;
            int a = (int)(seed % n);
            seed = (seed * 1103515245ULL + 12345ULL) & 0x7FFFFFFFULL;
            int b = (int)(seed % n);
            es.push_back(make_pair(a, b));
        }
        vector<int> got = (t == 1) ? solveUndirected(n, es) : solveDirected(n, es);
        vector<int> want = (t == 1) ? bruteUndirected(n, es) : bruteDirected(n, es);
        if (got.empty() != want.empty()) { printf("FAIL random existence iter %d\n", iter); return; }
        if (!validatePath(m, es, got, t == 2)) { printf("FAIL random path iter %d\n", iter); return; }
    }

    // 大圖冒煙：隨機長度的環（必有歐拉迴路）
    for (int iter = 0; iter < 20; iter++) {
        seed = (seed * 1103515245ULL + 12345ULL) & 0x7FFFFFFFULL;
        int k = (int)(seed % 500) + 2;
        vector<pair<int,int> > cyc;
        for (int i = 0; i < k; i++) cyc.push_back(make_pair(i, (i + 1) % k));
        vector<int> p = solveUndirected(k, cyc);
        if (p.empty() || (int)p.size() != k + 1) { printf("FAIL cycle smoke iter %d\n", iter); return; }
        if (!validatePath(k, cyc, p, false)) { printf("FAIL cycle validate iter %d\n", iter); return; }
    }

    printf("all tests passed\n");
}

int main() {
    // 一次性讀完 stdin
    stringstream buf;
    buf << cin.rdbuf();
    string input = buf.str();
    bool anyTok = false;
    for (size_t i = 0; i < input.size(); i++)
        if (!isspace((unsigned char)input[i])) { anyTok = true; break; }
    if (!anyTok) { runAllTests(); return 0; }
    istringstream it(input);
    int t, n, m;
    if (!(it >> t >> n >> m)) { printf("-1\n"); return 0; }
    vector<pair<int,int> > edges;
    for (int i = 0; i < m; i++) {
        int u, v;
        it >> u >> v;
        edges.push_back(make_pair(u, v));
    }
    vector<int> path = (t == 1) ? solveUndirected(n, edges) : solveDirected(n, edges);
    if (path.empty()) printf("-1\n");
    else {
        ostringstream oss;
        for (size_t i = 0; i < path.size(); i++) {
            if (i) oss << ' ';
            oss << path[i];
        }
        oss << '\n';
        fputs(oss.str().c_str(), stdout);
    }
    return 0;
}
