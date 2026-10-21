// 弦圖與完美消除序列：MCS / Lex-BFS 辨識，並在弦圖上線性時間求 ω、χ、α
//
// 題意：
//     給一張 n 點 m 邊的**無向簡單圖**，回答若干詢問：
//
//         mode 1   這張圖是不是**弦圖**（chordal / triangulated）？輸出 1 / 0
//         mode 2   輸出 MCS（Maximum Cardinality Search）的訪問順序
//         mode 3   輸出**完美消除序列**（PEO，= MCS 順序的反序）
//         mode 4   ω：最大團的大小
//         mode 5   χ：最小染色數，並給一組具體染色（每個點的顏色，0 起算）
//         mode 6   α：最大獨立集的大小，並給一組具體的最大獨立集（點號升序）
//         mode 7   輸出 PEO 中每個位置的「後續鄰居數」d_i（ω = 1 + max d_i）
//
// 思路：
//     ### 1. 什麼是弦圖
//     「長度 ≥ 4 的迴圈都**有弦**（有一條邊直接連起迴圈上兩個不相鄰的點）」，
//     等價說法是「不存在 induced cycle of length ≥ 4」。樹、區間圖、分裂圖都是弦圖。
//
//     弦圖最強的性質是它**有完美消除序列（PEO）**：
//     一個排列 v1..vn，使得每個 vi 在「後面」（v_{i+1}..v_n）的鄰居**構成一個團**。
//
//     ### 2. 兩個產生 PEO 的演算法（互相對拍）
//     - **MCS**：每次選「已訪問鄰居數最多」的未訪問點（平手取編號最小）。
//     - **Lex-BFS**：維護一個「集合序列」，每次取第一個集合的第一個點，
//       再把每個集合依「是否與該點相鄰」切成兩半（**相鄰的在前**，標籤較大）。
//
//     兩者的**反序**都是 PEO —— 但前提是圖真的是弦圖。
//     所以「把候選順序拿去驗 PEO」本身就是一個**弦圖辨識演算法**：
//     G 是弦圖 ⟺ reverse(MCS order) 是 PEO。
//
//     ### 3. 為什麼要用兩套
//     MCS 與 Lex-BFS 的思路完全不同（一個數鄰居、一個維護字典序分區），
//     卻必須對所有圖給出**相同**的「是否弦圖」答案。這是本題最主要的交叉驗證。
//     另外在 n ≤ 11 時再用「枚舉所有點子集、找 induced cycle」的 2^n 暴力法做第三方確認。
//
//     ### 4. 弦圖上三個 NP-hard 問題全部變線性
//     一般圖上求 ω（最大團）、χ（最小染色）、α（最大獨立集）都是 NP-hard，
//     但拿到 PEO 之後：
//
//         ω = 1 + max_i |{v_i 在 PEO 中後面的鄰居}|         （那個集合本身就是團）
//         χ：依 PEO 的**反序**做貪心染色（每次選最小的可用顏色）→ 恰好 ω 種顏色
//         α：依 PEO 的**正序**做貪心（與已選集合不相鄰就選）→ 最大獨立集
//
//     因為弦圖是**完美圖**（perfect graph），χ = ω；而 α = θ（最小團覆蓋數）。
//     本題會把 χ 與 ω 的相等性當成一條不變式來檢查。
//
//     ### 5. 兩語言一致的細節
//     - 平手規則全部寫死：MCS 取編號最小；Lex-BFS 的每個集合內部維持編號升序、取第一個。
//     - 非弦圖時 ω 的 PEO 公式**不成立**（後續鄰居集合不是團），所以不與暴力值比較；
//       α 仍必定是合法獨立集，因此任何圖都有 α ≤ 暴力值；χ 只檢查「染色合法」與「χ ≥ ω_brute」。
//     - 輸入的邊先正規化：自環、越界、重邊一律丟棄，兩語言才會看到同一張圖。
//
// 應用場景：
//     - 稀疏矩陣 Cholesky 分解：弦圖 = 完美消除順序 = 分解時不產生新的非零元（fill-in）
//     - 貝氏網路 / 概率圖模型：junction tree 演算法必須先把圖「三角化」成弦圖
//     - 資料庫：acyclic hypergraph（α-acyclic）與 join tree，等價於弦圖性質
//     - 約束滿足與排程：區間圖（弦圖的子類）上的教室排課、頻道分配是多項式時間
//     - 電腦圖學：三角化多邊形（polygon triangulation）本質上就是弦圖
//
// 複雜度（n 點 m 邊）：
//     MCS / Lex-BFS（本題的樸素實作）     O(n²)
//     PEO 驗證                            O(n + m)
//     ω / χ / α                           O(n + m)
//     暴力對拍（n ≤ 14 / n ≤ 11）          O(2^n · n)
//     空間                                O(n + m)

#include <algorithm>
#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <deque>
#include <iostream>
#include <random>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#ifdef _WIN32
#include <fcntl.h>
#include <io.h>
#endif

using namespace std;
using ull = unsigned long long;

static const int NMAX = 60;
static const int MMAX = 1200;
static const int QMAX = 20;
static const int BRUTE_N = 14;
static const int BRUTE_CYCLE_N = 11;

// ------------------------------------------------------------------ 圖

class Graph {
public:
    int n;
    vector<set<int>> adj;
    vector<ull> amask;

    explicit Graph(int nn) : n(nn), adj(nn), amask(nn, 0ULL) {}

    void addEdge(int u, int v) {
        if (u == v) return;
        if (u < 0 || u >= n || v < 0 || v >= n) return;
        adj[u].insert(v);
        adj[v].insert(u);
    }

    void finish() {
        for (int i = 0; i < n; ++i) {
            ull m = 0;
            for (int j : adj[i]) m |= (1ULL << j);
            amask[i] = m;
        }
    }
};

// ------------------------------------------------------------ 順序產生

static vector<int> mcsOrder(const Graph& g) {
    int n = g.n;
    vector<char> visited(n, 0);
    vector<int> cnt(n, 0), order;
    for (int k = 0; k < n; ++k) {
        int best = -1;
        for (int v = 0; v < n; ++v) {
            if (!visited[v] && (best < 0 || cnt[v] > cnt[best])) best = v;  // 平手取編號最小
        }
        visited[best] = 1;
        order.push_back(best);
        for (int w : g.adj[best]) cnt[w] += 1;
    }
    return order;
}

static vector<int> lexbfsOrder(const Graph& g) {
    vector<vector<int>> labels;
    if (g.n > 0) {
        vector<int> all(g.n);
        for (int i = 0; i < g.n; ++i) all[i] = i;
        labels.push_back(all);           // 已升序
    }
    vector<int> order;
    vector<char> picked(g.n, 0);
    while (!labels.empty()) {
        vector<int> s = labels[0];
        labels.erase(labels.begin());
        int v = s[0];
        s.erase(s.begin());
        if (!s.empty()) labels.insert(labels.begin(), s);
        picked[v] = 1;
        order.push_back(v);
        vector<vector<int>> nl;
        const set<int>& nb = g.adj[v];
        for (vector<int>& t : labels) {
            vector<int> tt;
            for (int x : t) if (!picked[x]) tt.push_back(x);
            if (tt.empty()) continue;
            vector<int> a, b;            // a = 與 v 相鄰（標籤較大，排前面）
            for (int x : tt) {
                if (nb.count(x)) a.push_back(x);
                else b.push_back(x);
            }
            if (!a.empty()) nl.push_back(a);
            if (!b.empty()) nl.push_back(b);
        }
        labels = nl;
    }
    return order;
}

// ------------------------------------------------------------------ PEO

static vector<int> peoLaterCounts(const Graph& g, const vector<int>& peo) {
    int n = g.n;
    vector<int> pos(n, 0), out(n, 0);
    for (int i = 0; i < n; ++i) pos[peo[i]] = i;
    for (int i = 0; i < n; ++i) {
        int v = peo[i], c = 0;
        for (int w : g.adj[v]) if (pos[w] > i) ++c;
        out[i] = c;
    }
    return out;
}

static bool isPeo(const Graph& g, const vector<int>& peo) {
    int n = g.n;
    vector<int> pos(n, 0);
    for (int i = 0; i < n; ++i) pos[peo[i]] = i;
    for (int i = 0; i < n; ++i) {
        int v = peo[i];
        vector<int> later;
        for (int w : g.adj[v]) if (pos[w] > i) later.push_back(w);
        sort(later.begin(), later.end(), [&](int a, int b) { return pos[a] < pos[b]; });
        for (int a = 0; a < (int)later.size(); ++a)
            for (int b = a + 1; b < (int)later.size(); ++b)
                if (!g.adj[later[a]].count(later[b])) return false;
    }
    return true;
}

static vector<int> colorByPeo(const Graph& g, const vector<int>& peo) {
    int n = g.n;
    vector<int> color(n, -1);
    for (int i = n - 1; i >= 0; --i) {
        int v = peo[i];
        vector<char> used(n + 2, 0);
        for (int w : g.adj[v]) if (color[w] >= 0) used[color[w]] = 1;
        int c = 0;
        while (used[c]) ++c;
        color[v] = c;
    }
    return color;
}

static vector<int> misByPeo(const Graph& g, const vector<int>& peo) {
    vector<int> chosen;
    vector<char> taken(g.n, 0);
    for (int v : peo) {
        bool bad = false;
        for (int w : g.adj[v]) if (taken[w]) { bad = true; break; }
        if (!bad) {
            chosen.push_back(v);
            taken[v] = 1;
        }
    }
    return chosen;
}

// ------------------------------------------------------------- 暴力對照

static int popcnt(ull x) { return __builtin_popcountll(x); }

static int bruteMaxClique(const Graph& g) {
    int n = g.n;
    int best = 0;
    for (ull mask = 1; mask < (1ULL << n); ++mask) {
        bool ok = true;
        ull m = mask;
        while (m) {
            ull b = m & (~m + 1ULL);
            int v = __builtin_ctzll(b);
            if (((g.amask[v] | b) & mask) != mask) { ok = false; break; }
            m -= b;
        }
        if (ok) best = max(best, popcnt(mask));
    }
    return best;
}

static int bruteMIS(const Graph& g) {
    int n = g.n;
    int best = 0;
    for (ull mask = 1; mask < (1ULL << n); ++mask) {
        bool ok = true;
        ull m = mask;
        while (m) {
            ull b = m & (~m + 1ULL);
            int v = __builtin_ctzll(b);
            if ((g.amask[v] & mask) != 0) { ok = false; break; }
            m -= b;
        }
        if (ok) best = max(best, popcnt(mask));
    }
    return best;
}

static bool bruteIsChordal(const Graph& g) {
    int n = g.n;
    for (ull mask = 1; mask < (1ULL << n); ++mask) {
        if (popcnt(mask) < 4) continue;
        bool ok = true;
        ull m = mask;
        while (m) {
            ull b = m & (~m + 1ULL);
            int v = __builtin_ctzll(b);
            if (popcnt(g.amask[v] & mask) != 2) { ok = false; break; }
            m -= b;
        }
        if (!ok) continue;
        // 所有點度數都是 2，再確認連通 ⟹ 是一個 induced cycle
        int start = __builtin_ctzll(mask & (~mask + 1ULL));
        ull seen = 1ULL << start;
        deque<int> q;
        q.push_back(start);
        while (!q.empty()) {
            int x = q.front();
            q.pop_front();
            for (int w : g.adj[x]) {
                if ((mask >> w) & 1ULL && !((seen >> w) & 1ULL)) {
                    seen |= 1ULL << w;
                    q.push_back(w);
                }
            }
        }
        if (seen == mask) return false;
    }
    return true;
}

static bool coloringProper(const Graph& g, const vector<int>& color) {
    for (int v = 0; v < g.n; ++v)
        for (int w : g.adj[v])
            if (color[v] == color[w]) return false;
    return true;
}

static bool independentValid(const Graph& g, const vector<int>& s) {
    set<int> st(s.begin(), s.end());
    for (int v : s)
        for (int w : g.adj[v])
            if (st.count(w)) return false;
    return true;
}

// ------------------------------------------------------------------ 驅動

static string run(int n, const vector<pair<int, int>>& edges, const vector<int>& modes) {
    Graph g(n);
    for (auto& e : edges) g.addEdge(e.first, e.second);
    g.finish();

    vector<int> mcs = mcsOrder(g);
    vector<int> lex = lexbfsOrder(g);
    vector<int> peo(mcs.rbegin(), mcs.rend());
    vector<int> peo2(lex.rbegin(), lex.rend());

    bool chordal = isPeo(g, peo);
    bool chordal2 = isPeo(g, peo2);
    if (chordal != chordal2) {
        fprintf(stderr, "chordal mismatch n=%d\n", n);
        abort();
    }
    if (n <= BRUTE_CYCLE_N) {
        bool bc = bruteIsChordal(g);
        if (chordal != bc) {
            fprintf(stderr, "brute chordal mismatch n=%d\n", n);
            abort();
        }
    }

    vector<int> d = peoLaterCounts(g, peo);
    int omega = d.empty() ? 0 : 1 + *max_element(d.begin(), d.end());
    vector<int> color = colorByPeo(g, peo);
    int chi = 0;
    for (int c : color) chi = max(chi, c + 1);
    vector<int> mis = misByPeo(g, peo);
    int alpha = (int)mis.size();

    if (!coloringProper(g, color)) { fprintf(stderr, "coloring invalid\n"); abort(); }
    if (!independentValid(g, mis)) { fprintf(stderr, "MIS invalid\n"); abort(); }

    if (n <= BRUTE_N) {
        int bo = bruteMaxClique(g);
        int ba = bruteMIS(g);
        // α 一定是合法獨立集，任何圖都 ≤ 暴力值；ω 的 PEO 公式只在弦圖才成立
        if (alpha > ba) { fprintf(stderr, "alpha > brute\n"); abort(); }
        if (chi < bo) { fprintf(stderr, "chi < clique\n"); abort(); }
        if (chordal) {
            if (omega != bo) { fprintf(stderr, "chordal omega\n"); abort(); }
            if (alpha != ba) { fprintf(stderr, "chordal alpha\n"); abort(); }
            if (chi != bo) { fprintf(stderr, "chordal chi\n"); abort(); }
        }
    }

    ostringstream oss;
    for (int m : modes) {
        if (m == 1) {
            oss << (chordal ? "1" : "0") << "\n";
        } else if (m == 2) {
            for (int i = 0; i < n; ++i) oss << (i ? " " : "") << (mcs[i] + 1);
            oss << "\n";
        } else if (m == 3) {
            for (int i = 0; i < n; ++i) oss << (i ? " " : "") << (peo[i] + 1);
            oss << "\n";
        } else if (m == 4) {
            oss << omega << "\n";
        } else if (m == 5) {
            oss << chi;
            for (int v = 0; v < n; ++v) oss << " " << color[v];
            oss << "\n";
        } else if (m == 6) {
            vector<int> s = mis;
            sort(s.begin(), s.end());
            oss << alpha;
            for (int v : s) oss << " " << (v + 1);
            oss << "\n";
        } else if (m == 7) {
            for (int i = 0; i < n; ++i) oss << (i ? " " : "") << d[i];
            oss << "\n";
        } else {
            oss << "error\n";
        }
    }
    return oss.str();
}

static void check(int n, const vector<pair<int, int>>& edges, const vector<int>& modes) {
    run(n, edges, modes);
}

static void runTests() {
    static const int ALL[] = {1, 2, 3, 4, 5, 6, 7};
    vector<int> allModes(ALL, ALL + 7);

    check(4, {{0, 1}, {1, 2}, {2, 3}}, allModes);                              // 路徑
    check(5, {{0, 1}, {0, 2}, {0, 3}, {0, 4}}, allModes);                      // 星形
    check(4, {{0, 1}, {1, 2}, {2, 3}, {3, 0}}, allModes);                      // 4-迴圈（非弦圖）
    check(4, {{0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}}, allModes);      // K4
    check(5, {{0, 1}, {1, 2}, {0, 2}, {2, 3}, {3, 4}, {2, 4}}, allModes);      // 兩個三角形
    check(0, {}, allModes);
    check(1, {}, allModes);
    check(3, {{0, 0}, {0, 1}, {0, 1}, {0, 5}}, allModes);                      // 自環 / 重邊 / 越界
    check(5, {{0, 1}, {1, 2}, {2, 3}, {3, 4}, {4, 0}, {0, 2}}, allModes);      // 5-迴圈 + 一條弦

    mt19937 rng(20261021);
    auto ri = [&](int a, int b) -> int { return (int)(rng() % (unsigned)(b - a + 1)) + a; };

    // 隨機樹（一定弦圖）
    for (int t = 0; t < 60; ++t) {
        int n = ri(1, 12);
        vector<pair<int, int>> edges;
        for (int v = 1; v < n; ++v) edges.push_back({v, ri(0, v - 1)});
        check(n, edges, allModes);
    }

    // 隨機一般圖
    for (int t = 0; t < 120; ++t) {
        int n = ri(0, 12);
        double p = 0.1 * ri(1, 6);
        vector<pair<int, int>> edges;
        for (int u = 0; u < n; ++u)
            for (int v = u + 1; v < n; ++v)
                if (rng() % 1000 < (unsigned)(p * 1000)) edges.push_back({u, v});
        check(n, edges, allModes);
    }

    // 隨機分裂圖：團 + 獨立集（一定弦圖）
    for (int t = 0; t < 40; ++t) {
        int n = ri(2, 12);
        int k = ri(1, n);
        vector<pair<int, int>> edges;
        for (int u = 0; u < k; ++u)
            for (int v = u + 1; v < k; ++v) edges.push_back({u, v});
        for (int v = k; v < n; ++v)
            for (int u = 0; u < k; ++u)
                if (ri(0, 1)) edges.push_back({u, v});
        check(n, edges, allModes);
    }

    // 隨機含雜訊邊（自環 / 越界 / 重邊）
    for (int t = 0; t < 40; ++t) {
        int n = ri(2, 10);
        vector<pair<int, int>> edges;
        int m = ri(1, 12);
        for (int i = 0; i < m; ++i) edges.push_back({ri(0, n + 1), ri(0, n + 1)});
        check(n, edges, allModes);
    }
}

// -- IO 解析 --

static bool isIntTokStrict(const string& s) {
    if (s.empty()) return false;
    size_t i = 0;
    if (s[0] == '+' || s[0] == '-') i = 1;
    if (i >= s.size()) return false;
    for (size_t k = i; k < s.size(); ++k)
        if (s[k] < '0' || s[k] > '9') return false;
    return true;
}

int main() {
#ifdef _WIN32
    _setmode(_fileno(stdout), _O_BINARY);
#endif
    string raw, line;
    ostringstream buf;
    while (getline(cin, line)) buf << line << "\n";
    raw = buf.str();

    bool hasInput = false;
    for (char c : raw) if (!isspace((unsigned char)c)) { hasInput = true; break; }
    if (!hasInput) {
        runTests();
        fputs("all tests passed\n", stdout);
        return 0;
    }

    vector<string> toks;
    {
        string cur;
        for (char c : raw) {
            if (isspace((unsigned char)c)) { if (!cur.empty()) { toks.push_back(cur); cur.clear(); } }
            else cur.push_back(c);
        }
        if (!cur.empty()) toks.push_back(cur);
    }
    size_t pos = 0;
    auto nxt = [&](long long defVal) -> long long {
        if (pos >= toks.size()) return defVal;
        string t = toks[pos++];
        if (!isIntTokStrict(t)) return defVal;
        long long v = strtoll(t.c_str(), nullptr, 10);
        if (v > 1000000000000000LL) v = 1000000000000000LL;
        if (v < -1000000000000000LL) v = -1000000000000000LL;
        return v;
    };

    long long nl = nxt(0);
    if (nl < 0) nl = 0;
    if (nl > NMAX) nl = NMAX;
    long long ml = nxt(0);
    if (ml < 0) ml = 0;
    if (ml > MMAX) ml = MMAX;
    vector<pair<int, int>> edges;
    for (long long i = 0; i < ml; ++i) {
        long long a = nxt(1) - 1;
        long long b = nxt(1) - 1;
        edges.push_back({(int)a, (int)b});   // 越界 / 自環由 addEdge 過濾
    }
    long long ql = nxt(0);
    if (ql < 0) ql = 0;
    if (ql > QMAX) ql = QMAX;
    vector<int> modes;
    for (long long i = 0; i < ql; ++i) modes.push_back((int)nxt(0));

    string out = run((int)nl, edges, modes);
    fwrite(out.data(), 1, out.size(), stdout);
    return 0;
}
