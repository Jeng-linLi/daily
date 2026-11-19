#include <algorithm>
#include <cassert>
#include <cctype>
#include <cstdint>
#include <iostream>
#include <iterator>
#include <queue>
#include <sstream>
#include <string>
#include <tuple>
#include <utility>
#include <vector>
#ifdef _WIN32
#include <fcntl.h>
#include <io.h>
#endif
using namespace std;

const long long INF = (1LL << 60);

// Edmonds-Karp 最小割：回傳 (flow, s 側可達集合)。cap 為 n×n 容量矩陣（合法餘量）
pair<long long, vector<char>> maxflowMincutSide(int n, vector<vector<long long>>& cap, int s, int t) {
    if (s == t) return {0, vector<char>(n, 0)};
    long long flow = 0;
    while (true) {
        vector<int> prev(n, -1);
        prev[s] = -2;
        queue<int> q;
        q.push(s);
        while (!q.empty() && prev[t] == -1) {
            int v = q.front(); q.pop();
            for (int u = 0; u < n; ++u) {
                if (cap[v][u] > 0 && prev[u] == -1) {
                    prev[u] = v;
                    q.push(u);
                    if (u == t) break;
                }
            }
        }
        if (prev[t] == -1) break;
        long long aug = INF;
        for (int v = t; v != s; v = prev[v]) aug = min(aug, cap[prev[v]][v]);
        for (int v = t; v != s; v = prev[v]) {
            cap[prev[v]][v] -= aug;
            cap[v][prev[v]] += aug;
        }
        flow += aug;
    }
    vector<char> seen(n, 0);
    seen[s] = 1;
    queue<int> q; q.push(s);
    while (!q.empty()) {
        int v = q.front(); q.pop();
        for (int u = 0; u < n; ++u)
            if (cap[v][u] > 0 && !seen[u]) { seen[u] = 1; q.push(u); }
    }
    return {flow, seen};
}

// 構造 Gomory-Hu：回傳 (parent, wpar)，wpar[v] = v 與 parent[v] 的邊權
pair<vector<int>, vector<long long>> buildTree(int n, const vector<tuple<int,int,long long>>& edges) {
    vector<vector<long long>> adj(n, vector<long long>(n, 0));
    for (auto [u, v, w] : edges) {
        if (u == v) continue;
        adj[u][v] += w;
        adj[v][u] += w;
    }
    vector<int> parent(n, 0);
    vector<long long> wpar(n, 0);
    if (n > 0) parent[0] = -1;
    for (int s = 1; s < n; ++s) {
        int t = parent[s];
        auto cap = adj;  // 最小割會消耗容量，必須傳副本
        auto [flow, side] = maxflowMincutSide(n, cap, s, t);
        wpar[s] = flow;
        for (int i = s + 1; i < n; ++i)
            if (parent[i] == t && side[i]) parent[i] = s;
    }
    return {parent, wpar};
}

// 樹上 s-t 路徑最小邊權
long long pathMin(const vector<int>& parent, const vector<long long>& wpar, int s, int t) {
    if (s == t) return 0;
    auto depth = [&](int v) {
        int d = 0;
        while (parent[v] != -1) { v = parent[v]; ++d; }
        return d;
    };
    int ds = depth(s), dt = depth(t);
    int a = s, b = t;
    long long best = INF;
    while (ds > dt) { best = min(best, wpar[a]); a = parent[a]; --ds; }
    while (dt > ds) { best = min(best, wpar[b]); b = parent[b]; --dt; }
    while (a != b) {
        best = min(best, wpar[a]); best = min(best, wpar[b]);
        a = parent[a]; b = parent[b];
    }
    return best;
}

bool parseLL(const string& t, long long& out) {
    if (t.empty()) return false;
    size_t i = 0;
    if (t[0] == '+' || t[0] == '-') i = 1;
    if (i >= t.size()) return false;
    for (size_t k = i; k < t.size(); ++k) if (t[k] < '0' || t[k] > '9') return false;
    try { out = stoll(t); } catch (...) { return false; }
    return true;
}

string runIO(const vector<string>& tok) {
    size_t pos = 0;
    auto nxt = [&]() -> string { return pos < tok.size() ? tok[pos++] : string(); };
    long long val;
    string tn = nxt();
    if (tn.empty() || !parseLL(tn, val)) return "";
    int n = (int)min(max(val, 0LL), 60LL);
    string tm = nxt(), tq = nxt();
    long long vm, vq;
    if (tm.empty() || tq.empty() || !parseLL(tm, vm) || !parseLL(tq, vq)) return "";
    int m = (int)min(max(vm, 0LL), 2000LL);
    int q = (int)min(max(vq, 0LL), 2000LL);
    vector<tuple<int,int,long long>> edges;
    for (int i = 0; i < m; ++i) {
        string a = nxt(), b = nxt(), c = nxt();
        long long u, v, w;
        if (a.empty() || b.empty() || c.empty()) break;
        if (!parseLL(a, u) || !parseLL(b, v) || !parseLL(c, w)) break;
        if (0 <= u && u < n && 0 <= v && v < n)
            edges.push_back({(int)u, (int)v, min(max(w, 0LL), 100000000LL)});
    }
    vector<pair<int,int>> queries;
    for (int i = 0; i < q; ++i) {
        string a = nxt(), b = nxt();
        long long s, t;
        if (a.empty() || b.empty()) break;
        if (!parseLL(a, s) || !parseLL(b, t)) break;
        if (0 <= s && s < n && 0 <= t && t < n) queries.push_back({(int)s, (int)t});
    }
    if (n == 0) return string(queries.size(), '\n') == "" ? "" : string(queries.size(), '\n');
    auto [parent, wpar] = buildTree(n, edges);
    string out;
    for (auto [s, t] : queries) out += to_string(pathMin(parent, wpar, s, t)) + "\n";
    return out;
}

uint32_t gstate = 20261119;
uint32_t grnd(uint32_t m) { gstate = gstate * 1103515245u + 12345u; return gstate % m; }

void runTests() {
    {
        vector<tuple<int,int,long long>> e = {{0,1,2},{1,2,3},{0,2,4}};
        auto [p, wp] = buildTree(3, e);
        vector<vector<long long>> adj(3, vector<long long>(3, 0));
        for (auto [u, v, w] : e) { adj[u][v] += w; adj[v][u] += w; }
        for (int s = 0; s < 3; ++s)
            for (int t = 0; t < 3; ++t) {
                auto cp = adj;
                long long want = maxflowMincutSide(3, cp, s, t).first;
                assert(pathMin(p, wp, s, t) == want);
            }
    }
    for (int tc = 0; tc < 120; ++tc) {
        int n = (int)(grnd(7) + 2);
        vector<tuple<int,int,long long>> edges;
        vector<vector<long long>> adj(n, vector<long long>(n, 0));
        for (int u = 0; u < n; ++u)
            for (int v = u + 1; v < n; ++v)
                if (grnd(2) == 0) {
                    long long w = grnd(15) + 1;
                    edges.push_back({u, v, w});
                    adj[u][v] += w; adj[v][u] += w;
                }
        auto [p, wp] = buildTree(n, edges);
        for (int s = 0; s < n; ++s)
            for (int t = 0; t < n; ++t) {
                auto cp = adj;
                long long want = maxflowMincutSide(n, cp, s, t).first;
                assert(pathMin(p, wp, s, t) == want);
            }
    }
}

int main() {
#ifdef _WIN32
    _setmode(_fileno(stdin), _O_BINARY);
    _setmode(_fileno(stdout), _O_BINARY);
#endif
    string data((istreambuf_iterator<char>(cin)), istreambuf_iterator<char>());
    vector<string> tok;
    istringstream iss(data);
    string x;
    while (iss >> x) tok.push_back(x);
    bool has = false;
    for (char c : data) if (!isspace((unsigned char)c)) { has = true; break; }
    if (has) cout << runIO(tok);
    else { runTests(); cout << "all tests passed\n"; }
    return 0;
}
