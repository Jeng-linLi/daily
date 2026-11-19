#include <algorithm>
#include <cassert>
#include <cctype>
#include <cstdint>
#include <iostream>
#include <iterator>
#include <numeric>
#include <sstream>
#include <string>
#include <utility>
#include <vector>
#ifdef _WIN32
#include <fcntl.h>
#include <io.h>
#endif
using namespace std;

// Lengauer–Tarjan 支配樹。n 點有向圖，回傳 idom（idom[s] = -1，不可達為 -1）。
vector<int> dominatorTree(int n, const vector<pair<int,int>>& edges, int s) {
    vector<int> idom(n, -1);
    if (n == 0) return idom;
    vector<vector<int>> pred(n), succ(n);
    for (auto [u, v] : edges) { succ[u].push_back(v); pred[v].push_back(u); }

    // 1) DFS 編號（迭代）
    vector<int> dfn(n, 0), vertex(n + 1, 0), parent(n, -1);
    int timer = 0;
    vector<pair<int,int>> st;
    st.push_back({s, 0});
    while (!st.empty()) {
        int v = st.back().first;
        int idx = st.back().second;
        if (idx == 0) { dfn[v] = ++timer; vertex[timer] = v; }
        if (idx < (int)succ[v].size()) {
            int u = succ[v][idx];
            st.back().second = idx + 1;
            if (dfn[u] == 0) { parent[u] = v; st.push_back({u, 0}); }
        } else {
            st.pop_back();
        }
    }

    // 2) 半支配點 + eval/link
    vector<int> semi = dfn;          // 以 dfn 值表示
    vector<int> label(n);
    iota(label.begin(), label.end(), 0);
    vector<int> ancestor(n, -1);
    vector<vector<int>> bucket(n);

    auto compress = [&](auto&& self, int v) -> void {
        if (ancestor[v] != -1 && ancestor[ancestor[v]] != -1) {
            if (semi[label[ancestor[v]]] < semi[label[v]])
                label[v] = label[ancestor[v]];
            ancestor[v] = ancestor[ancestor[v]];
            self(self, v);
        }
    };
    auto evalNode = [&](int v) -> int {
        if (ancestor[v] == -1) return label[v];
        compress(compress, v);
        return label[v];
    };

    for (int i = timer; i >= 2; --i) {
        int w = vertex[i];
        for (int v : pred[w]) {
            if (dfn[v] == 0) continue;
            int u = evalNode(v);
            if (semi[u] < semi[w]) semi[w] = semi[u];
        }
        bucket[vertex[semi[w]]].push_back(w);
        int p = parent[w];
        if (p != -1) ancestor[w] = p;
        if (p != -1) {
            for (int v : bucket[p]) {
                int u = evalNode(v);
                idom[v] = (semi[u] < semi[v]) ? u : p;
            }
            bucket[p].clear();
        }
    }

    // 3) 第二遍修正
    for (int i = 2; i <= timer; ++i) {
        int w = vertex[i];
        if (idom[w] != vertex[semi[w]]) idom[w] = idom[idom[w]];
    }
    idom[s] = -1;
    for (int v = 0; v < n; ++v) if (dfn[v] == 0) idom[v] = -1;
    return idom;
}

// ------------- 參考實作：迭代支配集合（O(n^3)）-------------
vector<int> bruteIdom(int n, const vector<pair<int,int>>& edges, int s) {
    vector<int> idom(n, -1);
    if (n == 0) return idom;
    vector<vector<int>> pred(n), succ(n);
    for (auto [u, v] : edges) { succ[u].push_back(v); pred[v].push_back(u); }
    vector<char> reach(n, 0);
    vector<int> st = {s};
    while (!st.empty()) {
        int v = st.back(); st.pop_back();
        if (reach[v]) continue;
        reach[v] = 1;
        for (int u : succ[v]) st.push_back(u);
    }
    vector<vector<char>> dom(n, vector<char>(n, 0));
    for (int v = 0; v < n; ++v) if (reach[v]) {
        if (v == s) dom[v][s] = 1;
        else for (int u = 0; u < n; ++u) if (reach[u]) dom[v][u] = 1;
    }
    bool changed = true;
    while (changed) {
        changed = false;
        for (int v = 0; v < n; ++v) {
            if (!reach[v] || v == s) continue;
            vector<int> ps;
            for (int p : pred[v]) if (reach[p]) ps.push_back(p);
            vector<char> nd(n, 0);
            if (ps.empty()) nd[v] = 1;
            else {
                for (int u = 0; u < n; ++u) {
                    bool inAll = true;
                    for (int p : ps) if (!dom[p][u]) { inAll = false; break; }
                    if (inAll) nd[u] = 1;
                }
                nd[v] = 1;
            }
            if (nd != dom[v]) { dom[v] = nd; changed = true; }
        }
    }
    for (int v = 0; v < n; ++v) {
        if (!reach[v] || v == s) { idom[v] = -1; continue; }
        int best = -1, bestSize = -1;
        for (int u = 0; u < n; ++u) {
            if (u != v && dom[v][u]) {
                int sz = 0; for (int x = 0; x < n; ++x) sz += dom[u][x];
                if (sz > bestSize || (sz == bestSize && u > best)) { bestSize = sz; best = u; }
            }
        }
        idom[v] = best;
    }
    return idom;
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
    int n = (int)min(max(val, 0LL), 1000LL);
    string tm = nxt(), ts = nxt();
    long long vm, vs;
    if (tm.empty() || ts.empty() || !parseLL(tm, vm) || !parseLL(ts, vs)) return "";
    int m = (int)min(max(vm, 0LL), 20000LL);
    int s = (int)vs;
    if (!(0 <= s && s < n)) s = 0;
    vector<pair<int,int>> edges;
    for (int i = 0; i < m; ++i) {
        string a = nxt(), b = nxt();
        long long u, v;
        if (a.empty() || b.empty()) break;
        if (!parseLL(a, u) || !parseLL(b, v)) break;
        if (0 <= u && u < n && 0 <= v && v < n) edges.push_back({(int)u, (int)v});
    }
    vector<int> res = dominatorTree(n, edges, s);
    string out;
    for (int i = 0; i < n; ++i) { if (i) out += ' '; out += to_string(res[i]); }
    out += '\n';
    return out;
}

uint32_t dstate = 20261119;
uint32_t drnd(uint32_t m) { dstate = dstate * 1103515245u + 12345u; return dstate % m; }

void runTests() {
    {
        auto r = dominatorTree(3, {{0,1},{1,2}}, 0);
        assert((r == vector<int>{-1,0,1}));
        auto r2 = dominatorTree(4, {{0,1},{0,2},{1,3},{2,3}}, 0);
        assert((r2 == vector<int>{-1,0,0,0}));
    }
    for (int tc = 0; tc < 500; ++tc) {
        int n = (int)(drnd(9) + 1);
        int s = (int)drnd(n);
        vector<pair<int,int>> edges;
        for (int u = 0; u < n; ++u)
            for (int v = 0; v < n; ++v)
                if (u != v && drnd(4) == 0) {
                    edges.push_back({u, v});
                    if (drnd(6) == 0) edges.push_back({u, v});  // 摻入重邊
                }
        assert(dominatorTree(n, edges, s) == bruteIdom(n, edges, s));
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
    string s;
    while (iss >> s) tok.push_back(s);
    bool has = false;
    for (char c : data) if (!isspace((unsigned char)c)) { has = true; break; }
    if (has) cout << runIO(tok);
    else { runTests(); cout << "all tests passed\n"; }
    return 0;
}
