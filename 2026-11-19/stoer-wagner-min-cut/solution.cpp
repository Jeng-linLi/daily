#include <algorithm>
#include <cassert>
#include <cctype>
#include <cstdint>
#include <iostream>
#include <iterator>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>
#ifdef _WIN32
#include <fcntl.h>
#include <io.h>
#endif
using namespace std;

const long long INF = (1LL << 60);

// 無向圖全局最小割，O(n^3)。回傳最小割權值。
long long stoerWagner(int n, const vector<tuple<int,int,long long>>& edges) {
    if (n <= 1) return 0;
    vector<vector<long long>> w(n, vector<long long>(n, 0));
    for (auto [u, v, c] : edges) {
        if (u == v) continue;
        w[u][v] += c;
        w[v][u] += c;
    }
    vector<char> active(n, 1);
    long long best = INF;
    int remaining = n;
    while (remaining > 1) {
        vector<char> visited(n, 0);
        vector<long long> weight(n, 0);
        vector<int> order;
        for (int step = 0; step < remaining; ++step) {
            int pick = -1;
            for (int v = 0; v < n; ++v) {
                if (active[v] && !visited[v]) {
                    if (pick == -1 || weight[v] > weight[pick]) pick = v;
                }
            }
            visited[pick] = 1;
            order.push_back(pick);
            for (int u = 0; u < n; ++u) {
                if (active[u] && !visited[u]) weight[u] += w[pick][u];
            }
        }
        int t = order[order.size() - 1];
        int s = order[order.size() - 2];
        long long cut = weight[t];
        if (cut < best) best = cut;
        // 合併 t 進 s
        for (int u = 0; u < n; ++u) {
            if (active[u] && u != s && u != t) {
                w[s][u] += w[t][u];
                w[u][s] += w[u][t];
            }
        }
        active[t] = 0;
        --remaining;
    }
    return best;
}

// 暴力枚舉所有割（固定點 0 在一側），僅供對拍
long long bruteMinCut(int n, const vector<tuple<int,int,long long>>& edges) {
    if (n <= 1) return 0;
    vector<vector<long long>> w(n, vector<long long>(n, 0));
    for (auto [u, v, c] : edges) {
        if (u != v) { w[u][v] += c; w[v][u] += c; }
    }
    long long ans = INF;
    for (int mask = 0; mask < (1 << (n - 1)) - 1; ++mask) {
        vector<char> in(n, 0);
        in[0] = 1;
        for (int i = 1; i < n; ++i) if (mask >> (i - 1) & 1) in[i] = 1;
        long long cut = 0;
        for (int u = 0; u < n; ++u)
            for (int v = u + 1; v < n; ++v)
                if (in[u] != in[v]) cut += w[u][v];
        if (cut < ans) ans = cut;
    }
    return ans;
}

// ------------- 整數解析（嚴格 ^[+-]?[0-9]+$）-------------
bool parseLL(const string& t, long long& out) {
    if (t.empty()) return false;
    size_t i = 0;
    if (t[0] == '+' || t[0] == '-') i = 1;
    if (i >= t.size()) return false;
    for (size_t k = i; k < t.size(); ++k) if (t[k] < '0' || t[k] > '9') return false;
    out = 0;
    try { out = stoll(t); } catch (...) { return false; }
    return true;
}

string runIO(const vector<string>& tok) {
    size_t pos = 0;
    auto nxt = [&]() -> string {
        if (pos >= tok.size()) return string();
        return tok[pos++];
    };
    long long val;
    string tn = nxt();
    if (tn.empty() || !parseLL(tn, val)) return "";
    int n = (int)min(max(val, 0LL), 200LL);
    string tm = nxt();
    if (tm.empty() || !parseLL(tm, val)) return "";
    int m = (int)min(max(val, 0LL), 20000LL);
    vector<tuple<int,int,long long>> edges;
    for (int i = 0; i < m; ++i) {
        string a = nxt(), b = nxt(), c = nxt();
        long long u, v, w;
        if (a.empty() || b.empty() || c.empty()) break;
        if (!parseLL(a, u) || !parseLL(b, v) || !parseLL(c, w)) break;
        if (0 <= u && u < n && 0 <= v && v < n)
            edges.push_back({(int)u, (int)v, min(max(w, 0LL), 1000000000LL)});
    }
    return to_string(stoerWagner(n, edges)) + "\n";
}

uint32_t st = 20261119;
uint32_t rnd(uint32_t m) {
    st = st * 1103515245u + 12345u;
    return st % m;
}

void runTests() {
    {
        vector<tuple<int,int,long long>> e = {{0,1,2},{0,2,3},{1,2,1},{1,3,4},{2,3,5}};
        long long got = stoerWagner(4, e);
        assert(got == 5);
        assert(bruteMinCut(4, e) == 5);
        assert(stoerWagner(1, {}) == 0);
        assert(stoerWagner(0, {}) == 0);
        assert(stoerWagner(2, {{0,1,7}}) == 7);
        assert(stoerWagner(3, {{0,1,1},{1,2,1}}) == 1);
    }
    for (int tc = 0; tc < 300; ++tc) {
        int n = (int)(rnd(9) + 2);
        vector<tuple<int,int,long long>> edges;
        for (int u = 0; u < n; ++u)
            for (int v = u + 1; v < n; ++v)
                if (rnd(3) == 0) edges.push_back({u, v, (long long)(rnd(20) + 1)});
        assert(stoerWagner(n, edges) == bruteMinCut(n, edges));
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
    if (has) {
        cout << runIO(tok);
    } else {
        runTests();
        cout << "all tests passed\n";
    }
    return 0;
}
