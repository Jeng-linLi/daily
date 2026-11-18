// Hopcroft-Karp 二分圖最大匹配（BFS 分層 + DFS 找極大不相交增廣路）。
//
// IO 模式（stdin）：
//   第一行 nL nR m，接著 m 行 u v（0-based）。可讀多筆測資直到 EOF。
//   每筆輸出 4 行：最大匹配數 / 匹配對 u-v / 最小點覆蓋大小 / 最小點覆蓋 Lu Rv。
// stdin 無內容時跑內建對拍（vs Kuhn + 2^n 暴力點覆蓋），輸出 all tests passed。
//
// 複雜度 O(E·sqrt(V))，space O(V+E)。

#include <cassert>
#include <cstdio>
#include <string>
#include <vector>
#include <queue>
#include <algorithm>
#include <functional>
#include <sstream>
#include <iostream>

#ifdef _WIN32
#include <fcntl.h>
#include <io.h>
#endif

using namespace std;

struct HK {
    int nL, nR;
    vector<vector<int> > adj;
    vector<int> pairU, pairV, dist;
    static const int INF = 1000000000;

    HK(int nl, int nr) : nL(nl), nR(nr), adj(nl), pairU(nl, -1), pairV(nr, -1), dist(nl) {}

    bool bfs() {
        queue<int> q;
        bool found = false;
        for (int u = 0; u < nL; u++) {
            if (pairU[u] == -1) { dist[u] = 0; q.push(u); }
            else dist[u] = INF;
        }
        while (!q.empty()) {
            int u = q.front(); q.pop();
            for (size_t i = 0; i < adj[u].size(); i++) {
                int v = adj[u][i];
                int pu = pairV[v];
                if (pu == -1) found = true;
                else if (dist[pu] == INF) { dist[pu] = dist[u] + 1; q.push(pu); }
            }
        }
        return found;
    }

    bool dfs(int u) {
        for (size_t i = 0; i < adj[u].size(); i++) {
            int v = adj[u][i];
            int pu = pairV[v];
            if (pu == -1 || (dist[pu] == dist[u] + 1 && dfs(pu))) {
                pairU[u] = v;
                pairV[v] = u;
                return true;
            }
        }
        dist[u] = INF;
        return false;
    }

    int solve() {
        int matching = 0;
        while (bfs())
            for (int u = 0; u < nL; u++)
                if (pairU[u] == -1 && dfs(u)) matching++;
        return matching;
    }
};

// Kőnig：從未匹配左點出發交錯 BFS，覆蓋 = (左部扣掉 Z) 聯集 (右部交 Z)。
static void minVertexCover(const HK &h, vector<int> &left, vector<int> &right) {
    vector<char> seenL(h.nL, 0), seenR(h.nR, 0);
    queue<int> q;
    for (int u = 0; u < h.nL; u++)
        if (h.pairU[u] == -1) { seenL[u] = 1; q.push(u); }
    while (!q.empty()) {
        int u = q.front(); q.pop();
        for (size_t i = 0; i < h.adj[u].size(); i++) {
            int v = h.adj[u][i];
            if (!seenR[v]) {
                seenR[v] = 1;
                int pu = h.pairV[v];
                if (pu != -1 && !seenL[pu]) { seenL[pu] = 1; q.push(pu); }
            }
        }
    }
    for (int u = 0; u < h.nL; u++) if (!seenL[u]) left.push_back(u);
    for (int v = 0; v < h.nR; v++) if (seenR[v]) right.push_back(v);
}

// ----------------------------- 暴力參考 -----------------------------

static int kuhn(int nL, int nR, const vector<vector<int> > &adj) {
    vector<int> pairV(nR, -1);

    function<bool(int, vector<char>&)> tryK = [&](int u, vector<char> &seen) {
        for (size_t i = 0; i < adj[u].size(); i++) {
            int v = adj[u][i];
            if (seen[v]) continue;
            seen[v] = 1;
            if (pairV[v] == -1 || tryK(pairV[v], seen)) { pairV[v] = u; return true; }
        }
        return false;
    };

    int res = 0;
    for (int u = 0; u < nL; u++) {
        vector<char> seen(nR, 0);
        if (tryK(u, seen)) res++;
    }
    return res;
}

static int bruteCover(int nL, int nR, const vector<pair<int,int> > &edges) {
    int n = nL + nR;
    int best = n;
    for (int mask = 0; mask < (1 << n); mask++) {
        int k = __builtin_popcount((unsigned)mask);
        if (k >= best) continue;
        bool ok = true;
        for (size_t i = 0; i < edges.size(); i++) {
            int u = edges[i].first, v = edges[i].second;
            if (!((mask >> u) & 1) && !((mask >> (nL + v)) & 1)) { ok = false; break; }
        }
        if (ok) best = k;
    }
    return best;
}

// ------------------------------- 格式化 -------------------------------

static string formatResult(const HK &h, const vector<int> &left, const vector<int> &right) {
    vector<pair<int,int> > pairs;
    for (int u = 0; u < h.nL; u++)
        if (h.pairU[u] != -1) pairs.push_back({u, h.pairU[u]});
    sort(pairs.begin(), pairs.end());
    ostringstream out;
    out << pairs.size() << '\n';
    for (size_t i = 0; i < pairs.size(); i++) {
        if (i) out << ' ';
        out << pairs[i].first << '-' << pairs[i].second;
    }
    if (pairs.empty()) out << '-';
    out << '\n';
    out << (left.size() + right.size()) << '\n';
    bool first = true;
    for (size_t i = 0; i < left.size(); i++) {
        if (!first) out << ' ';
        first = false;
        out << 'L' << left[i];
    }
    for (size_t i = 0; i < right.size(); i++) {
        if (!first) out << ' ';
        first = false;
        out << 'R' << right[i];
    }
    if (first) out << '-';
    out << '\n';
    return out.str();
}

// ------------------------------- IO 模式 -------------------------------

static unsigned long long g_state = 20261118ULL;
static int rnd(int m) {
    g_state = (g_state * 1103515245ULL + 12345ULL) & 0x7FFFFFFFULL;
    return (int)(g_state % (unsigned)m);
}

static string run_io(const vector<string> &toks) {
    size_t pos = 0;
    string out;
    auto nxt = [&]() -> const string* {
        if (pos >= toks.size()) return nullptr;
        return &toks[pos++];
    };
    while (true) {
        const string *a = nxt(), *b = nxt(), *c = nxt();
        if (!a || !b || !c) break;
        int nL, nR, m;
        try { nL = stoi(*a); nR = stoi(*b); m = stoi(*c); }
        catch (...) { break; }
        if (nL < 0 || nL > 2000 || nR < 0 || nR > 2000 || m < 0 || m > 100000) break;
        HK h(nL, nR);
        bool ok = true;
        for (int i = 0; i < m; i++) {
            const string *x = nxt(), *y = nxt();
            if (!x || !y) { ok = false; break; }
            int u, v;
            try { u = stoi(*x); v = stoi(*y); }
            catch (...) { ok = false; break; }
            if (0 <= u && u < nL && 0 <= v && v < nR) h.adj[u].push_back(v);
        }
        if (!ok) break;
        for (int u = 0; u < nL; u++) {
            sort(h.adj[u].begin(), h.adj[u].end());
            h.adj[u].erase(unique(h.adj[u].begin(), h.adj[u].end()), h.adj[u].end());
        }
        h.solve();
        vector<int> left, right;
        minVertexCover(h, left, right);
        string block = formatResult(h, left, right);
        if (!block.empty() && block[block.size()-1] == '\n') block.erase(block.size()-1);
        out += block + "\n";
    }
    return out;
}

// ------------------------------ 內建測試 ------------------------------

static void run_tests() {
    // 固定情境（輸出已實跑確認）
    {
        HK h(3, 3);
        h.adj[0] = {1, 0}; h.adj[1] = {0}; h.adj[2] = {2};
        assert(h.solve() == 3);
        vector<int> left, right;
        minVertexCover(h, left, right);
        assert(formatResult(h, left, right) == "3\n0-1 1-0 2-2\n3\nL0 L1 L2\n");
    }

    // 隨機對拍
    for (int tc = 0; tc < 500; tc++) {
        int nL = rnd(5) + 1, nR = rnd(5) + 1;
        HK h(nL, nR);
        vector<pair<int,int> > edges;
        for (int u = 0; u < nL; u++)
            for (int v = 0; v < nR; v++)
                if (rnd(3) == 0) { h.adj[u].push_back(v); edges.push_back({u,v}); }
        int match = h.solve();
        assert(match == kuhn(nL, nR, h.adj));
        for (int u = 0; u < nL; u++)
            if (h.pairU[u] != -1) {
                assert(h.pairV[h.pairU[u]] == u);
                assert(find(h.adj[u].begin(), h.adj[u].end(), h.pairU[u]) != h.adj[u].end());
            }
        vector<int> left, right;
        minVertexCover(h, left, right);
        int cover = (int)left.size() + (int)right.size();
        assert(cover == match);
        assert(cover == bruteCover(nL, nR, edges));
        for (size_t i = 0; i < edges.size(); i++) {
            int u = edges[i].first, v = edges[i].second;
            assert(find(left.begin(), left.end(), u) != left.end() ||
                   find(right.begin(), right.end(), v) != right.end());
        }
    }

    // 邊界
    {
        HK h(0, 0);
        assert(h.solve() == 0);
        HK h2(2, 2);
        h2.adj[0] = {0}; h2.adj[1] = {1};
        assert(h2.solve() == 2);
    }
}

int main() {
#ifdef _WIN32
    _setmode(_fileno(stdout), _O_BINARY);
#endif
    ostringstream buf;
    buf << cin.rdbuf();
    string input = buf.str();
    if (!input.empty() && input.find_first_not_of(" \t\r\n") != string::npos) {
        istringstream iss(input);
        vector<string> toks;
        string x;
        while (iss >> x) toks.push_back(x);
        fputs(run_io(toks).c_str(), stdout);
    } else {
        run_tests();
        fputs("all tests passed\n", stdout);
    }
    return 0;
}
