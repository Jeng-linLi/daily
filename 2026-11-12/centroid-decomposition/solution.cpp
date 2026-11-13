// 點分治（Centroid Decomposition）：統計樹上「邊數恰為 k」的點對路徑數。
//
// 命令（IO 模式，從 stdin 讀入）：
//   paths n k u1 v1 u2 v2 ...   （n 個點、n-1 條無向邊，1-based）
//       -> 輸出點對 (u < v) 的個數，使兩點間路徑的邊數恰為 k
//
// stdin 無內容時跑內建斷言測試，輸出 all tests passed。
//
// 核心要點：
//   * 每輪找分量重心，統計經過重心的長度恰 k 的路徑（有序對容斥：
//     全分量 Σ cnt[d]*cnt[k-d] − 每棵子樹 Σ cnt_s[d]*cnt_s[k-d]，
//     k >= 1 時必為偶數，//2 即無序對數），刪重心後對子分量遞迴。
//   * 任意點對在「兩端首次被分開」的那層恰被統計一次，不重不漏。
//   * 複雜度 O(n log n)。

#include <cstdio>
#include <cstdlib>
#include <cctype>
#include <algorithm>
#include <string>
#include <vector>

using namespace std;

typedef long long ll;

static int gN, gK;
static vector<vector<int>> gAdj;
static vector<bool> gRemoved;

// 在未刪除節點中從 start 出發 BFS（可避開 avoid），回傳 (order, parent, dist)
static void bfs(int start, int avoid, vector<int>& order,
                vector<int>& parent, vector<int>& dist) {
    order.clear();
    parent.assign(gN, -1);
    dist.assign(gN, -1);
    dist[start] = 0;
    order.push_back(start);
    for (size_t head = 0; head < order.size(); ++head) {
        int u = order[head];
        for (size_t t = 0; t < gAdj[u].size(); ++t) {
            int v = gAdj[u][t];
            if (!gRemoved[v] && v != avoid && dist[v] < 0) {
                dist[v] = dist[u] + 1;
                parent[v] = u;
                order.push_back(v);
            }
        }
    }
}

static int findCentroid(int start) {
    vector<int> order, parent, dist;
    bfs(start, -1, order, parent, dist);
    vector<int> size(gN, 0);
    for (int t = (int)order.size() - 1; t >= 0; --t) {
        int u = order[t];
        size[u] = 1;
        for (size_t s = 0; s < gAdj[u].size(); ++s) {
            int v = gAdj[u][s];
            if (!gRemoved[v] && v != parent[u]) size[u] += size[v];
        }
    }
    int tot = size[start];
    int c = start;
    while (true) {
        int heavy = -1;
        for (size_t t = 0; t < gAdj[c].size(); ++t) {
            int v = gAdj[c][t];
            if (!gRemoved[v] && v != parent[c] && size[v] * 2 > tot) {
                heavy = v;
                break;
            }
        }
        if (heavy < 0) return c;
        c = heavy;
    }
}

static ll gTotal;

static void decompose(int start) {
    int c = findCentroid(start);
    vector<int> order, parent, distc;
    bfs(c, -1, order, parent, distc);
    vector<ll> cntAll(gK + 1, 0);
    for (size_t t = 0; t < order.size(); ++t) {
        int u = order[t];
        if (distc[u] <= gK) cntAll[distc[u]] += 1;
    }
    ll through = 0;
    for (int d = 0; d <= gK; ++d) through += cntAll[d] * cntAll[gK - d];
    for (size_t t = 0; t < gAdj[c].size(); ++t) {
        int v0 = gAdj[c][t];
        if (gRemoved[v0]) continue;
        vector<int> o2, p2, d2;
        bfs(v0, c, o2, p2, d2);  // 子樹內 BFS（避開 c）
        vector<ll> cntS(gK + 1, 0);
        for (size_t s = 0; s < o2.size(); ++s) {
            int u = o2[s];
            int du = d2[u] + 1;  // 距重心 = 距 v0 + 1
            if (du <= gK) cntS[du] += 1;
        }
        for (int d = 0; d <= gK; ++d) through -= cntS[d] * cntS[gK - d];
    }
    gTotal += through / 2;
    gRemoved[c] = true;
    for (size_t t = 0; t < gAdj[c].size(); ++t) {
        int v = gAdj[c][t];
        if (!gRemoved[v]) decompose(v);
    }
}

static ll countPaths(int n, int k, const vector<pair<int, int>>& edges) {
    if (n <= 1 || k < 1) return 0;
    gN = n;
    gK = k;
    gAdj.assign(n, vector<int>());
    for (size_t t = 0; t < edges.size(); ++t) {
        int u = edges[t].first, v = edges[t].second;
        gAdj[u].push_back(v);
        gAdj[v].push_back(u);
    }
    gRemoved.assign(n, false);
    gTotal = 0;
    decompose(0);
    return gTotal;
}

// O(n^2) 暴力（每點 BFS 全樹），僅供內建測試對拍
static ll countPathsBrute(int n, int k, const vector<pair<int, int>>& edges) {
    if (n <= 1 || k < 1) return 0;
    vector<vector<int>> adj(n);
    for (size_t t = 0; t < edges.size(); ++t) {
        adj[edges[t].first].push_back(edges[t].second);
        adj[edges[t].second].push_back(edges[t].first);
    }
    ll pairs = 0;
    vector<int> dist(n), q;
    for (int s = 0; s < n; ++s) {
        fill(dist.begin(), dist.end(), -1);
        dist[s] = 0;
        q.clear();
        q.push_back(s);
        for (size_t head = 0; head < q.size(); ++head) {
            int u = q[head];
            for (size_t t = 0; t < adj[u].size(); ++t) {
                int v = adj[u][t];
                if (dist[v] < 0) {
                    dist[v] = dist[u] + 1;
                    q.push_back(v);
                }
            }
        }
        for (int t2 = s + 1; t2 < n; ++t2) {
            if (dist[t2] == k) pairs += 1;
        }
    }
    return pairs;
}

// 嚴格整數解析：^[+-]?[0-9]+$，非法 token 視為輸入結束
static bool tryParse(const string& s, ll& v) {
    if (s.empty()) return false;
    size_t p = 0;
    if (s[0] == '+' || s[0] == '-') p = 1;
    if (p >= s.size()) return false;
    for (size_t t = p; t < s.size(); ++t) {
        if (s[t] < '0' || s[t] > '9') return false;
    }
    v = strtoll(s.c_str(), NULL, 10);
    return true;
}

#ifdef _WIN32
#include <io.h>
#include <fcntl.h>
#endif

int main() {
#ifdef _WIN32
    _setmode(_fileno(stdout), _O_BINARY);
#endif
    string data;
    char buf[65536];
    size_t got;
    while ((got = fread(buf, 1, sizeof(buf), stdin)) > 0) data.append(buf, got);

    bool hasInput = false;
    for (char c : data) {
        if (!isspace((unsigned char)c)) {
            hasInput = true;
            break;
        }
    }

    if (hasInput) {
        vector<string> toks;
        {
            string cur;
            for (char c : data) {
                if (isspace((unsigned char)c)) {
                    if (!cur.empty()) {
                        toks.push_back(cur);
                        cur.clear();
                    }
                } else {
                    cur += c;
                }
            }
            if (!cur.empty()) toks.push_back(cur);
        }
        size_t pos = 0;
        auto nxt = [&](string& outTok) -> bool {
            if (pos >= toks.size()) return false;
            outTok = toks[pos++];
            return true;
        };
        string cmd;
        while (nxt(cmd)) {
            if (cmd != "paths") break;
            string t;
            ll n, k;
            if (!nxt(t) || !tryParse(t, n)) break;
            if (!nxt(t) || !tryParse(t, k)) break;
            int nn = (int)n, kk = (int)k;
            vector<pair<int, int>> edges;
            bool ok = true;
            for (int i = 0; i + 1 < nn && ok; ++i) {
                ll u, v;
                if (!nxt(t) || !tryParse(t, u)) ok = false;
                else if (!nxt(t) || !tryParse(t, v)) ok = false;
                else edges.push_back(make_pair((int)u - 1, (int)v - 1));
            }
            if (!ok) break;
            printf("%lld\n", countPaths(nn, kk, edges));
        }
        fflush(stdout);
        return 0;
    }

    // ---- 內建斷言測試 ----
    {
        vector<pair<int, int>> chain = {{0, 1}, {1, 2}, {2, 3}, {3, 4}};
        vector<pair<int, int>> star = {{0, 1}, {0, 2}, {0, 3}, {0, 4}};
        if (countPaths(5, 1, chain) != 4) return 1;
        if (countPaths(5, 2, chain) != 3) return 1;
        if (countPaths(5, 3, chain) != 2) return 1;
        if (countPaths(5, 4, chain) != 1) return 1;
        if (countPaths(5, 5, chain) != 0) return 1;
        if (countPaths(5, 1, star) != 4) return 1;
        if (countPaths(5, 2, star) != 6) return 1;
        if (countPaths(5, 3, star) != 0) return 1;
        if (countPaths(1, 1, {}) != 0) return 1;
        if (countPaths(2, 1, {{0, 1}}) != 1) return 1;
    }
    {
        unsigned int st = 20261112u;
        auto rnd = [&](unsigned int mm) -> unsigned int {
            st = st * 1103515245u + 12345u;
            return (st >> 1) % mm;
        };
        for (int tc = 0; tc < 200; ++tc) {
            int n = (int)rnd(80) + 1;
            vector<pair<int, int>> edges;
            for (int i = 1; i < n; ++i) {
                edges.push_back(make_pair((int)rnd((unsigned int)i), i));
            }
            int k = (int)rnd(12) + 1;
            if (countPaths(n, k, edges) != countPathsBrute(n, k, edges)) return 1;
        }
        int n = 400;
        vector<pair<int, int>> edges;
        for (int i = 0; i + 1 < n; ++i) edges.push_back(make_pair(i, i + 1));
        if (countPaths(n, 7, edges) != n - 7) return 1;
        if (countPaths(n, 399, edges) != 1) return 1;
    }

    printf("all tests passed\n");
    return 0;
}
