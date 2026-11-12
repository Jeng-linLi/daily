// 樹雜湊與樹同構判定（AHU 演算法 / 重心為根的典範雜湊）。
//
// 命令（IO 模式，從 stdin 讀入）：
//   iso n1 m1 u1 v1 ... n2 m2 u2 v2 ...
//       -> 兩棵無標號樹（各 n 個點、n-1 條邊，1-based）
//       -> 同構輸出 1，否則 0
//   hash n m u1 v1 ...
//       -> 輸出該樹的典範雜湊（兩個模數下的值，與 Python 版逐位元一致）
//
// stdin 無內容時跑內建斷言測試，輸出 all tests passed。
//
// 核心要點：
//   * 子樹簽名：hash(u) 由「所有孩子簽名排序後」多項式摺疊而得，
//     與孩子順序無關——這正是無標號樹同構的本質。
//   * 以重心為根：以固定 1 號點為根時，鏈狀樹兩端簽名不同，會把
//     同一棵樹判成異構。重心至多兩個，各算一次取最小即可。
//   * 雙雜湊（兩個質數模）壓低碰撞機率；平手取字典序最小，保證確定性。
//   * 應用：無標號樹計數（A000055）、子樹重複偵測、分子結構比對。

#include <cstdio>
#include <cstdlib>
#include <cctype>
#include <algorithm>
#include <string>
#include <vector>
#include <queue>
#include <set>

using namespace std;

typedef long long ll;
typedef pair<ll, ll> pll;

static const ll M1 = 2147483647LL;  // 2^31 - 1，質數
static const ll M2 = 2147483629LL;  // 質數
static const ll B1 = 1000003;
static const ll B2 = 1000033;

static vector<vector<int>> gAdj;

// 回傳樹的全部重心（1 或 2 個，0-based）
static vector<int> centroids(int n) {
    if (n == 1) return vector<int>(1, 0);
    const vector<vector<int>>& adj = gAdj;
    vector<int> size(n, 0), maxsub(n, 0), parent(n, -1), order;
    vector<bool> visited(n, false);
    vector<int> stack;
    stack.push_back(0);
    visited[0] = true;
    while (!stack.empty()) {
        int u = stack.back();
        stack.pop_back();
        order.push_back(u);
        for (size_t t = 0; t < adj[u].size(); ++t) {
            int v = adj[u][t];
            if (!visited[v]) {
                visited[v] = true;
                parent[v] = u;
                stack.push_back(v);
            }
        }
    }
    for (int t = (int)order.size() - 1; t >= 0; --t) {
        int u = order[t];
        size[u] = 1;
        for (size_t s = 0; s < adj[u].size(); ++s) {
            int v = adj[u][s];
            if (v != parent[u]) size[u] += size[v];
        }
        if (parent[u] >= 0) {
            maxsub[parent[u]] = max(maxsub[parent[u]], size[u]);
        }
    }
    int best = -1;
    vector<int> cents;
    for (int u = 0; u < n; ++u) {
        int up = n - size[u];
        int m = max(maxsub[u], up);
        if (best < 0 || m < best) {
            best = m;
            cents.clear();
            cents.push_back(u);
        } else if (m == best) {
            cents.push_back(u);
        }
    }
    return cents;
}

// 以 u 為根的子樹簽名（孩子簽名排序後摺疊）
static pll hashSub(int u, int p) {
    vector<pll> forms;
    for (size_t t = 0; t < gAdj[u].size(); ++t) {
        int v = gAdj[u][t];
        if (v != p) forms.push_back(hashSub(v, u));
    }
    sort(forms.begin(), forms.end());
    ll h1 = 1, h2 = 1;
    for (size_t t = 0; t < forms.size(); ++t) {
        ll a = forms[t].first, b = forms[t].second;
        h1 = (h1 * B1 + a * 2147483648LL + b) % M1;
        h2 = (h2 * B2 + a * 2147483648LL + b) % M2;
    }
    return make_pair(h1, h2);
}

// 典範雜湊：對每個重心求子樹簽名，取最小
static pll treeHash(int n) {
    vector<int> cents = centroids(n);
    pll best;
    for (size_t t = 0; t < cents.size(); ++t) {
        pll h = hashSub(cents[t], -1);
        if (t == 0 || h < best) best = h;
    }
    return best;
}

// 邊清單（1-based）-> 全域鄰接表
static void buildTree(int n, const vector<pair<int, int>>& edges) {
    gAdj.assign(n, vector<int>());
    for (size_t t = 0; t < edges.size(); ++t) {
        int u = edges[t].first - 1, v = edges[t].second - 1;
        gAdj[u].push_back(v);
        gAdj[v].push_back(u);
    }
}

static bool isomorphic(int n1, const vector<pair<int, int>>& e1,
                       int n2, const vector<pair<int, int>>& e2) {
    if (n1 != n2) return false;
    if (n1 == 0) return true;
    buildTree(n1, e1);
    pll h1 = treeHash(n1);
    buildTree(n2, e2);
    pll h2 = treeHash(n2);
    return h1 == h2;
}

// Prüfer 序列（1-based，長度 n-2）-> 邊清單，僅供內建測試
static vector<pair<int, int>> pruferToEdges(const vector<int>& seq, int n) {
    vector<int> deg(n + 1, 1);
    for (size_t t = 0; t < seq.size(); ++t) deg[seq[t]] += 1;
    priority_queue<int, vector<int>, greater<int>> leaves;
    for (int i = 1; i <= n; ++i) {
        if (deg[i] == 1) leaves.push(i);
    }
    vector<pair<int, int>> edges;
    for (size_t t = 0; t < seq.size(); ++t) {
        int leaf = leaves.top();
        leaves.pop();
        edges.push_back(make_pair(leaf, seq[t]));
        deg[leaf] -= 1;
        deg[seq[t]] -= 1;
        if (deg[seq[t]] == 1) leaves.push(seq[t]);
    }
    vector<int> last;
    for (int i = 1; i <= n; ++i) {
        if (deg[i] == 1) last.push_back(i);
    }
    edges.push_back(make_pair(last[0], last[1]));
    return edges;
}

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
        auto readTree = [&](int& n, vector<pair<int, int>>& edges) -> bool {
            string t;
            ll nn, m;
            if (!nxt(t) || !tryParse(t, nn)) return false;
            if (!nxt(t) || !tryParse(t, m)) return false;
            n = (int)nn;
            edges.clear();
            for (ll i = 0; i < m; ++i) {
                ll u, v;
                if (!nxt(t) || !tryParse(t, u)) return false;
                if (!nxt(t) || !tryParse(t, v)) return false;
                edges.push_back(make_pair((int)u, (int)v));
            }
            return true;
        };
        string cmd;
        while (nxt(cmd)) {
            if (cmd == "iso") {
                int n1, n2;
                vector<pair<int, int>> e1, e2;
                if (!readTree(n1, e1)) break;
                if (!readTree(n2, e2)) break;
                printf("%d\n", isomorphic(n1, e1, n2, e2) ? 1 : 0);
            } else if (cmd == "hash") {
                int n;
                vector<pair<int, int>> e;
                if (!readTree(n, e)) break;
                buildTree(n, e);
                pll h = treeHash(n);
                printf("%lld %lld\n", h.first, h.second);
            } else {
                break;
            }
        }
        fflush(stdout);
        return 0;
    }

    // ---- 內建斷言測試 ----
    {
        // 邊界與手算例
        if (!isomorphic(1, {}, 1, {})) return 1;
        if (!isomorphic(2, {{1, 2}}, 2, {{1, 2}})) return 1;
        if (isomorphic(2, {{1, 2}}, 3, {{1, 2}, {2, 3}})) return 1;
        if (!isomorphic(3, {{1, 2}, {2, 3}}, 3, {{2, 1}, {2, 3}})) return 1;
        if (isomorphic(4, {{1, 2}, {2, 3}, {3, 4}},
                       4, {{1, 2}, {1, 3}, {1, 4}})) return 1;
        // 重心為根的必要性：同一條鏈反轉編號必須同構
        if (!isomorphic(5, {{1, 2}, {2, 3}, {3, 4}, {4, 5}},
                        5, {{2, 1}, {3, 2}, {4, 3}, {5, 4}})) return 1;
    }
    {
        unsigned int st = 20261111u;
        auto rnd = [&](unsigned int mm) -> unsigned int {
            st = st * 1103515245u + 12345u;
            return (st >> 1) % mm;
        };
        auto randomTree = [&](int n) -> vector<pair<int, int>> {
            vector<pair<int, int>> edges;
            for (int v = 2; v <= n; ++v) {
                int u = (int)rnd((unsigned int)(v - 1)) + 1;
                edges.push_back(make_pair(u, v));
            }
            return edges;
        };
        auto relabel = [&](int n, const vector<pair<int, int>>& edges)
            -> vector<pair<int, int>> {
            vector<int> perm(n + 1);
            for (int i = 0; i <= n; ++i) perm[i] = i;
            for (int i = n; i >= 2; --i) {
                int j = (int)rnd((unsigned int)i) + 1;
                swap(perm[i], perm[j]);
            }
            vector<pair<int, int>> out;
            for (size_t t = 0; t < edges.size(); ++t) {
                out.push_back(make_pair(perm[edges[t].first],
                                        perm[edges[t].second]));
            }
            return out;
        };
        // 重標號後雜湊必須不變
        for (int tc = 0; tc < 200; ++tc) {
            int n = (int)rnd(40) + 1;
            vector<pair<int, int>> e = randomTree(n);
            buildTree(n, e);
            pll h1 = treeHash(n);
            vector<pair<int, int>> e2 = relabel(n, e);
            buildTree(n, e2);
            pll h2 = treeHash(n);
            if (h1 != h2) return 1;
        }
        // 無標號樹計數：n=1..7 的雜湊去重數應為 A000055：1,1,1,2,3,6,11
        int unlabeled[7] = {1, 1, 1, 2, 3, 6, 11};
        for (int n = 1; n <= 7; ++n) {
            set<pll> hashes;
            if (n == 1) {
                buildTree(1, {});
                hashes.insert(treeHash(1));
            } else {
                int total = 1;
                for (int t = 0; t < n - 2; ++t) total *= n;
                // Prüfer 序列窮舉（n-2 位、每位 1..n）
                vector<int> seq(n - 2, 1);
                for (int idx = 0; idx < total; ++idx) {
                    int tmp = idx;
                    for (int t = 0; t < n - 2; ++t) {
                        seq[t] = tmp % n + 1;
                        tmp /= n;
                    }
                    vector<pair<int, int>> e = pruferToEdges(seq, n);
                    buildTree(n, e);
                    hashes.insert(treeHash(n));
                }
            }
            if ((int)hashes.size() != unlabeled[n - 1]) return 1;
        }
    }

    printf("all tests passed\n");
    return 0;
}
