// 矩陣樹定理（Matrix-Tree Theorem / Kirchhoff 定理）：帶權生成樹計數
//
// 題意：
//     給定 n 個點、m 條邊的無向帶權圖（允許重邊、不允許自環；權值可為 0，
//     權 0 的邊等價於不存在）。生成樹的權值 = 樹上所有邊權的乘積，
//     求所有生成樹權值總和，對 998244353 取模。不連通時答案為 0。
//     輸入：第一行 n m，接下來 m 行每行 u v w（1 起算，u != v）。
//     輸出：一個整數。
//
// 思路：
//     拉普拉斯矩陣 L = D - A（D 為帶權度數矩陣、A 為帶權鄰接矩陣，重邊累加）。
//     去掉第 0 行第 0 列的 n-1 階子式 M 滿足 det(M) = 生成樹邊權乘積總和。
//     在模質數 p = 998244353 下高斯消元：找非零主元（找不到 => det = 0，
//     對應不連通）、換行變號、乘主元、用逆元消下方各行。
//     全程整數運算無捨入誤差（實數消元會失真）。
//
// 複雜度：
//     建矩陣 O(n + m)，行列式 O(n^3)，空間 O(n^2)。
//
// 應用場景：
//     網路可靠性（生成樹數 = 冗餘度）、電路等效電阻（Kirchhoff 電阻網）、
//     驗證 Cayley 公式（K_n 有 n^(n-2) 棵生成樹）、Markov chain tree theorem。

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>
#include <algorithm>

#ifdef _WIN32
#include <io.h>
#include <fcntl.h>
#endif

typedef long long ll;

const ll MOD = 998244353;
const int NMAX = 200;
const int MMAX = 2000;
const ll WMAX = 1000000000LL;

static ll power_mod(ll a, ll e, ll mod) {
    ll r = 1 % mod;
    a %= mod;
    while (e > 0) {
        if (e & 1) r = r * a % mod;
        a = a * a % mod;
        e >>= 1;
    }
    return r;
}

// 模 p 高斯消元求行列式（mat 為 k x k，會被就地修改）
static ll det_mod(std::vector<std::vector<ll> > &mat) {
    int k = (int)mat.size();
    ll det = 1;
    for (int col = 0; col < k; ++col) {
        int piv = -1;
        for (int r = col; r < k; ++r) {
            if (mat[r][col] != 0) { piv = r; break; }
        }
        if (piv == -1) return 0;
        if (piv != col) {
            std::swap(mat[col], mat[piv]);
            det = MOD - det;  // 換行變號
        }
        det = det * mat[col][col] % MOD;
        ll inv = power_mod(mat[col][col], MOD - 2, MOD);
        for (int r = col + 1; r < k; ++r) {
            ll f = mat[r][col] * inv % MOD;
            if (f) {
                for (int c = col; c < k; ++c) {
                    mat[r][c] = (mat[r][c] - f * mat[col][c] % MOD + MOD) % MOD;
                }
            }
        }
    }
    return det % MOD;
}

// 矩陣樹定理：回傳帶權生成樹總權值和（mod MOD）
// edges 為 (u, v, w)，u / v 為 0 起算且 u != v
static ll spanning_tree_count(int n,
                              const std::vector<std::vector<ll> > &edges) {
    if (n == 1) return 1;  // 空積：單點圖恰有一棵「空生成樹」
    std::vector<std::vector<ll> > L(n, std::vector<ll>(n, 0));
    for (size_t i = 0; i < edges.size(); ++i) {
        ll u = edges[i][0], v = edges[i][1], w = edges[i][2];
        if (w == 0) continue;
        L[u][u] = (L[u][u] + w) % MOD;
        L[v][v] = (L[v][v] + w) % MOD;
        L[u][v] = (L[u][v] - w % MOD + MOD) % MOD;
        L[v][u] = (L[v][u] - w % MOD + MOD) % MOD;
    }
    // 去掉第 0 行第 0 列
    std::vector<std::vector<ll> > minor(n - 1, std::vector<ll>(n - 1, 0));
    for (int r = 1; r < n; ++r)
        for (int c = 1; c < n; ++c)
            minor[r - 1][c - 1] = L[r][c];
    return det_mod(minor);
}

// ---- 並查集（僅供暴力對拍）----
struct DSU {
    std::vector<int> fa;
    explicit DSU(int n) : fa(n) {
        for (int i = 0; i < n; ++i) fa[i] = i;
    }
    int find(int x) {
        while (fa[x] != x) { fa[x] = fa[fa[x]]; x = fa[x]; }
        return x;
    }
    bool union_xy(int a, int b) {
        int ra = find(a), rb = find(b);
        if (ra == rb) return false;
        fa[ra] = rb;
        return true;
    }
};

// 暴力：枚舉所有 2^m 個子集，取大小恰為 n-1 且構成生成樹者
static ll brute_force(int n, const std::vector<std::vector<ll> > &edges) {
    int m = (int)edges.size();
    ll total = 0;
    for (int mask = 0; mask < (1 << m); ++mask) {
        int bits = 0;
        for (int i = 0; i < m; ++i) if (mask >> i & 1) ++bits;
        if (bits != n - 1) continue;
        DSU dsu(n);
        ll prod = 1;
        bool ok = true;
        for (int i = 0; i < m && ok; ++i) {
            if (mask >> i & 1) {
                prod = prod * edges[i][2] % MOD;
                if (!dsu.union_xy((int)edges[i][0], (int)edges[i][1])) ok = false;
            }
        }
        if (ok && prod) total = (total + prod) % MOD;
    }
    return total;
}

// 自製 LCG（與 Python 版同一序列）
static unsigned int lcg_state = 20261103u;
static unsigned int lcg_next() {
    lcg_state = (unsigned int)((ll)lcg_state * 1103515245LL + 12345LL) & 0x7FFFFFFFu;
    return lcg_state;
}

static void emit(const std::string &s) {
    fwrite(s.data(), 1, s.size(), stdout);
}

static void run_tests() {
    // ---- 固定用例 ----
    std::vector<std::vector<ll> > e;
    // 1) 兩點一邊：帶權計數 = 邊權 = 5
    e.clear(); e.push_back({0, 1, 5});
    if (spanning_tree_count(2, e) != 5) { emit("FAIL 1\n"); return; }
    // 2) 三角形帶權：2*3 + 2*5 + 3*5 = 31
    e.clear(); e.push_back({0, 1, 2}); e.push_back({1, 2, 3}); e.push_back({0, 2, 5});
    if (spanning_tree_count(3, e) != 31) { emit("FAIL 2\n"); return; }
    // 3) 不連通
    e.clear(); e.push_back({0, 1, 1});
    if (spanning_tree_count(3, e) != 0) { emit("FAIL 3\n"); return; }
    // 4) 重邊：權 1 與 3 的平行邊 → 1 + 3 = 4
    e.clear(); e.push_back({0, 1, 1}); e.push_back({0, 1, 3});
    if (spanning_tree_count(2, e) != 4) { emit("FAIL 4\n"); return; }
    // 5) 權 0 的邊等價於不存在
    e.clear(); e.push_back({0, 1, 0});
    if (spanning_tree_count(2, e) != 0) { emit("FAIL 5\n"); return; }
    // 6) 單點圖
    if (spanning_tree_count(1, e) != 1) { emit("FAIL 6\n"); return; }
    // 7) K4：Cayley 公式 4^2 = 16
    e.clear();
    for (int u = 0; u < 4; ++u)
        for (int v = u + 1; v < 4; ++v)
            e.push_back({(ll)u, (ll)v, 1});
    if (spanning_tree_count(4, e) != 16) { emit("FAIL 7\n"); return; }

    // ---- 與暴力對拍 ----
    for (int it = 0; it < 300; ++it) {
        int n = 2 + (int)(lcg_next() % 4);           // 2..5
        int m = (int)(lcg_next() % (unsigned)(n * 2));
        std::vector<std::vector<ll> > edges;
        for (int i = 0; i < m; ++i) {
            int u = (int)(lcg_next() % (unsigned)n);
            int v = (int)(lcg_next() % (unsigned)n);
            if (u == v) v = (v + 1) % n;
            ll w = (ll)(lcg_next() % 4u);
            std::vector<ll> t(3);
            t[0] = u; t[1] = v; t[2] = w;
            edges.push_back(t);
        }
        ll got = spanning_tree_count(n, edges);
        ll want = brute_force(n, edges);
        if (got != want) { emit("FAIL random\n"); return; }
    }
}

// ---- IO 模式 ----
// 嚴格整數解析：^[+-]?[0-9]+$，超出 ±10^18 回 false
static bool parse_ll(const char *s, ll &out) {
    if (!s || !*s) return false;
    int i = 0;
    if (s[0] == '+' || s[0] == '-') i = 1;
    if (!s[i]) return false;
    for (int j = i; s[j]; ++j)
        if (s[j] < '0' || s[j] > '9') return false;
    // 位數上限（<= 19 位數字）
    if (strlen(s) - i > 19) return false;
    long long v = 0;
    for (int j = i; s[j]; ++j) v = v * 10 + (s[j] - '0');
    if (v > 1000000000000000000LL) return false;
    out = (s[0] == '-') ? -v : v;
    return true;
}

int main() {
#ifdef _WIN32
    _setmode(_fileno(stdout), _O_BINARY);
#endif
    // 讀入整個 stdin
    std::string raw;
    {
        char buf[1 << 16];
        size_t n;
        while ((n = fread(buf, 1, sizeof(buf), stdin)) > 0) raw.append(buf, n);
    }
    // 切 token
    std::vector<std::string> toks;
    {
        std::string cur;
        for (size_t i = 0; i <= raw.size(); ++i) {
            char ch = (i < raw.size()) ? raw[i] : ' ';
            if (ch == ' ' || ch == '\n' || ch == '\r' || ch == '\t') {
                if (!cur.empty()) { toks.push_back(cur); cur.clear(); }
            } else {
                cur += ch;
            }
        }
    }
    if (toks.empty()) {
        run_tests();
        emit("all tests passed\n");
        return 0;
    }

    size_t pos = 0;
    auto bad = [&]() { emit("error\n"); };
    ll n_ll, m_ll;
    if (pos >= toks.size() || !parse_ll(toks[pos++].c_str(), n_ll) ||
        n_ll < 1 || n_ll > NMAX) { bad(); return 0; }
    if (pos >= toks.size() || !parse_ll(toks[pos++].c_str(), m_ll) ||
        m_ll < 0 || m_ll > MMAX) { bad(); return 0; }
    int n = (int)n_ll, m = (int)m_ll;

    std::vector<std::vector<ll> > edges;
    for (int i = 0; i < m; ++i) {
        ll u, v, w;
        if (pos >= toks.size() || !parse_ll(toks[pos++].c_str(), u) ||
            pos >= toks.size() || !parse_ll(toks[pos++].c_str(), v) ||
            pos >= toks.size() || !parse_ll(toks[pos++].c_str(), w)) { bad(); return 0; }
        if (!(1 <= u && u <= n && 1 <= v && v <= n) || u == v ||
            !(0 <= w && w <= WMAX)) { bad(); return 0; }
        std::vector<ll> t(3);
        t[0] = u - 1; t[1] = v - 1; t[2] = w;
        edges.push_back(t);
    }

    char buf[32];
    ll ans = spanning_tree_count(n, edges);
    sprintf(buf, "%lld\n", ans);
    emit(buf);
    return 0;
}
