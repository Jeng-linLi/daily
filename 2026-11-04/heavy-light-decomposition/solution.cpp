// 重鏈剖分 Heavy-Light Decomposition（HLD）：樹上路徑修改與查詢
//
// 題意：
//     n 個點的樹（根為 1 號點），每個點有點權。q 筆操作：
//         1 u v x —— u 到 v 路徑上所有點權 +x
//         2 u v   —— 查詢路徑點權最大值
//         3 u v   —— 查詢路徑點權和
//     輸入：n q、n 個點權、n-1 條邊、q 筆操作（1 起算）。
//     輸出：對每筆 2 / 3 操作各輸出一行。
//
// 思路：
//     重輕劃分：子樹最大的兒子為重兒子，重邊連成重鏈。
//     從任意點到根經過的輕邊不超過 O(log n) 條（每過一條輕邊
//     子樹大小至少翻倍）。DFS 時重兒子優先，重鏈內 pos 連續，
//     路徑拆成 O(log n) 段連續區間 → 線段樹（區間加 + 區間和 + 區間 max）。
//     單次操作 O(log^2 n)。
//
// 複雜度：
//     預處理 O(n)，單次操作 O(log^2 n)，空間 O(n)。
//
// 應用場景：
//     組織架構匯報鏈批量調整 / 網路路由路徑加權與瓶頸查詢 /
//     版本依賴樹統計 / 樹上問題轉序列問題的通用框架。

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

const int NMAX = 20000;
const int QMAX = 20000;
const ll NEG_INF = -(1LL << 62);

// ---- 線段樹：區間加、區間和、區間最大值 ----
struct SegTree {
    int n;
    std::vector<ll> sm, mx, lz;

    void init(int n_, const std::vector<ll> &vals) {
        n = n_;
        sm.assign(4 * n, 0);
        mx.assign(4 * n, NEG_INF);
        lz.assign(4 * n, 0);
        if (n > 0) build(1, 0, n - 1, vals);
    }
    void build(int o, int l, int r, const std::vector<ll> &vals) {
        if (l == r) { sm[o] = mx[o] = vals[l]; return; }
        int m = (l + r) >> 1;
        build(o * 2, l, m, vals);
        build(o * 2 + 1, m + 1, r, vals);
        pull(o);
    }
    void pull(int o) {
        sm[o] = sm[o * 2] + sm[o * 2 + 1];
        mx[o] = std::max(mx[o * 2], mx[o * 2 + 1]);
    }
    void apply(int o, int l, int r, ll x) {
        sm[o] += x * (ll)(r - l + 1);
        mx[o] += x;
        lz[o] += x;
    }
    void push(int o, int l, int r) {
        if (lz[o]) {
            int m = (l + r) >> 1;
            apply(o * 2, l, m, lz[o]);
            apply(o * 2 + 1, m + 1, r, lz[o]);
            lz[o] = 0;
        }
    }
    void add(int o, int l, int r, int ql, int qr, ll x) {
        if (qr < l || r < ql) return;
        if (ql <= l && r <= qr) { apply(o, l, r, x); return; }
        push(o, l, r);
        int m = (l + r) >> 1;
        add(o * 2, l, m, ql, qr, x);
        add(o * 2 + 1, m + 1, r, ql, qr, x);
        pull(o);
    }
    ll query_sum(int o, int l, int r, int ql, int qr) {
        if (qr < l || r < ql) return 0;
        if (ql <= l && r <= qr) return sm[o];
        push(o, l, r);
        int m = (l + r) >> 1;
        return query_sum(o * 2, l, m, ql, qr) +
               query_sum(o * 2 + 1, m + 1, r, ql, qr);
    }
    ll query_max(int o, int l, int r, int ql, int qr) {
        if (qr < l || r < ql) return NEG_INF;
        if (ql <= l && r <= qr) return mx[o];
        push(o, l, r);
        int m = (l + r) >> 1;
        return std::max(query_max(o * 2, l, m, ql, qr),
                        query_max(o * 2 + 1, m + 1, r, ql, qr));
    }
};

// ---- HLD（0 起算，根為 0 號點）----
struct HLD {
    int n;
    std::vector<int> par, dep, sz, heavy, head, pos;
    SegTree seg;

    void build(int n_, const std::vector<int> &parents,
               const std::vector<ll> &vals) {
        n = n_;
        par = parents;
        dep.assign(n, 0);
        sz.assign(n, 1);
        heavy.assign(n, -1);
        head.assign(n, 0);
        pos.assign(n, 0);
        std::vector<std::vector<int> > ch(n);
        for (int v = 1; v < n; ++v) ch[par[v]].push_back(v);
        // BFS 定深
        std::vector<int> order;
        order.reserve(n);
        order.push_back(0);
        for (size_t i = 0; i < order.size(); ++i) {
            int u = order[i];
            for (size_t j = 0; j < ch[u].size(); ++j) {
                int w = ch[u][j];
                dep[w] = dep[u] + 1;
                order.push_back(w);
            }
        }
        // 逆序求子樹大小與重兒子
        for (int i = n - 1; i >= 0; --i) {
            int u = order[i], best = -1, bw = -1;
            for (size_t j = 0; j < ch[u].size(); ++j) {
                int w = ch[u][j];
                sz[u] += sz[w];
                if (sz[w] > best) { best = sz[w]; bw = w; }
            }
            heavy[u] = bw;
        }
        // 重兒子優先的迭代 DFS，分配 pos 與 head
        int p = 0;
        std::vector<int> stk;
        stk.reserve(n);
        stk.push_back(0);
        while (!stk.empty()) {
            int u = stk.back();
            stk.pop_back();
            pos[u] = p++;
            head[u] = (par[u] == -1 || heavy[par[u]] != u)
                          ? u : head[par[u]];
            for (size_t j = 0; j < ch[u].size(); ++j)
                if (ch[u][j] != heavy[u]) stk.push_back(ch[u][j]);
            if (heavy[u] != -1) stk.push_back(heavy[u]);
        }
        // 按 pos 重排點權
        std::vector<ll> ordered(n, 0);
        for (int u = 0; u < n; ++u) ordered[pos[u]] = vals[u];
        seg.init(n, ordered);
    }

    // 把 u-v 路徑拆成若干 [l, r] 區間（l <= r，pos 座標）
    std::vector<std::pair<int, int> > jump(int u, int v) const {
        std::vector<std::pair<int, int> > segs;
        while (head[u] != head[v]) {
            if (dep[head[u]] < dep[head[v]]) std::swap(u, v);
            segs.push_back(std::make_pair(pos[head[u]], pos[u]));
            u = par[head[u]];
        }
        if (dep[u] > dep[v]) std::swap(u, v);
        segs.push_back(std::make_pair(pos[u], pos[v]));
        return segs;
    }
    void path_add(int u, int v, ll x) {
        std::vector<std::pair<int, int> > segs = jump(u, v);
        for (size_t i = 0; i < segs.size(); ++i)
            seg.add(1, 0, n - 1, segs[i].first, segs[i].second, x);
    }
    ll path_sum(int u, int v) {
        std::vector<std::pair<int, int> > segs = jump(u, v);
        ll r = 0;
        for (size_t i = 0; i < segs.size(); ++i)
            r += seg.query_sum(1, 0, n - 1, segs[i].first, segs[i].second);
        return r;
    }
    ll path_max(int u, int v) {
        std::vector<std::pair<int, int> > segs = jump(u, v);
        ll r = NEG_INF;
        for (size_t i = 0; i < segs.size(); ++i)
            r = std::max(r, seg.query_max(1, 0, n - 1,
                                          segs[i].first, segs[i].second));
        return r;
    }
};

// 暴力：收集 u-v 路徑上所有點
static void path_nodes(const std::vector<int> &par,
                       const std::vector<int> &dep,
                       int u, int v, std::vector<int> &out) {
    out.clear();
    while (dep[u] > dep[v]) { out.push_back(u); u = par[u]; }
    while (dep[v] > dep[u]) { out.push_back(v); v = par[v]; }
    while (u != v) {
        out.push_back(u);
        out.push_back(v);
        u = par[u];
        v = par[v];
    }
    out.push_back(u);
}

static void emit(const std::string &s) { fwrite(s.data(), 1, s.size(), stdout); }

// 自製 LCG（與 Python 版同一序列）
static unsigned int lcg_state = 20261104u;
static int nxt_mod(int mod) {
    lcg_state = (unsigned int)((ll)lcg_state * 1103515245LL + 12345LL) & 0x7FFFFFFFu;
    return (int)(lcg_state % (unsigned)mod);
}

static void run_tests() {
    // ---- 固定用例：鏈 0-1-2-3-4 ----
    {
        HLD h;
        std::vector<int> par = {-1, 0, 1, 2, 3};
        std::vector<ll> vals = {1, 2, 3, 4, 5};
        h.build(5, par, vals);
        if (h.path_sum(0, 4) != 15) { emit("FAIL 1\n"); return; }
        if (h.path_max(0, 4) != 5) { emit("FAIL 2\n"); return; }
        h.path_add(1, 3, 10);
        if (h.path_sum(0, 4) != 45) { emit("FAIL 3\n"); return; }
        if (h.path_max(2, 4) != 14) { emit("FAIL 4\n"); return; }
        if (h.path_max(0, 1) != 12) { emit("FAIL 5\n"); return; }
    }
    // ---- 星形：0 連 1..4 ----
    {
        HLD s;
        std::vector<int> par = {-1, 0, 0, 0, 0};
        std::vector<ll> vals = {5, 1, 4, 2, 3};
        s.build(5, par, vals);
        if (s.path_sum(1, 2) != 10) { emit("FAIL 6\n"); return; }
        if (s.path_max(1, 3) != 5) { emit("FAIL 7\n"); return; }
        if (s.path_sum(0, 3) != 7) { emit("FAIL 8\n"); return; }
        s.path_add(0, 4, 3);
        if (s.path_sum(0, 4) != 14) { emit("FAIL 9\n"); return; }
        if (s.path_max(0, 4) != 8) { emit("FAIL 10\n"); return; }
    }
    // ---- 與暴力對拍 ----
    for (int it = 0; it < 120; ++it) {
        int n = 2 + nxt_mod(60);              // 2..61
        std::vector<int> par(n), dep(n, 0);
        std::vector<ll> vals(n);
        par[0] = -1;
        for (int i = 0; i < n; ++i) vals[i] = (ll)nxt_mod(201) - 100;
        for (int v = 1; v < n; ++v) par[v] = nxt_mod(v);
        HLD h;
        h.build(n, par, vals);
        dep = h.dep;
        for (int k = 0; k < 30; ++k) {
            int u = nxt_mod(n), v = nxt_mod(n), op = nxt_mod(3);
            if (op == 0) {
                ll x = nxt_mod(21) - 10;
                h.path_add(u, v, x);
                std::vector<int> nodes;
                path_nodes(par, dep, u, v, nodes);
                for (size_t i = 0; i < nodes.size(); ++i)
                    vals[nodes[i]] += x;
            } else {
                ll got = (op == 1) ? h.path_max(u, v) : h.path_sum(u, v);
                std::vector<int> nodes;
                path_nodes(par, dep, u, v, nodes);
                ll want = (op == 1) ? NEG_INF : 0;
                for (size_t i = 0; i < nodes.size(); ++i) {
                    if (op == 1) want = std::max(want, vals[nodes[i]]);
                    else want += vals[nodes[i]];
                }
                if (got != want) { emit("FAIL random\n"); return; }
            }
        }
    }
}

// ---- IO 模式 ----
static bool parse_ll(const char *s, ll &out) {
    if (!s || !*s) return false;
    int i = 0;
    if (s[0] == '+' || s[0] == '-') i = 1;
    if (!s[i]) return false;
    for (int j = i; s[j]; ++j)
        if (s[j] < '0' || s[j] > '9') return false;
    if (strlen(s) - i > 18) return false;
    long long v = 0;
    for (int j = i; s[j]; ++j) v = v * 10 + (s[j] - '0');
    out = (s[0] == '-') ? -v : v;
    return true;
}

int main() {
#ifdef _WIN32
    _setmode(_fileno(stdout), _O_BINARY);
#endif
    std::string raw;
    {
        char buf[1 << 16];
        size_t n;
        while ((n = fread(buf, 1, sizeof(buf), stdin)) > 0) raw.append(buf, n);
    }
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
    auto next_ll = [&](ll &out) -> bool {
        if (pos >= toks.size()) return false;
        return parse_ll(toks[pos++].c_str(), out);
    };

    ll n_ll, q_ll;
    if (!next_ll(n_ll) || n_ll < 1 || n_ll > NMAX) { bad(); return 0; }
    if (!next_ll(q_ll) || q_ll < 0 || q_ll > QMAX) { bad(); return 0; }
    int n = (int)n_ll, q = (int)q_ll;

    std::vector<ll> vals(n);
    for (int i = 0; i < n; ++i) {
        if (!next_ll(vals[i]) || vals[i] > 1000000000LL ||
            vals[i] < -1000000000LL) { bad(); return 0; }
    }
    std::vector<std::vector<int> > adj(n);
    for (int i = 0; i < n - 1; ++i) {
        ll a, b;
        if (!next_ll(a) || !next_ll(b) || a < 1 || a > n || b < 1 || b > n ||
            a == b) { bad(); return 0; }
        adj[a - 1].push_back((int)b - 1);
        adj[b - 1].push_back((int)a - 1);
    }
    // 以 0 為根的有根樹
    std::vector<int> parents(n, -1);
    {
        std::vector<bool> seen(n, false);
        std::vector<int> bfs;
        bfs.push_back(0);
        seen[0] = true;
        for (size_t i = 0; i < bfs.size(); ++i) {
            int u = bfs[i];
            for (size_t j = 0; j < adj[u].size(); ++j) {
                int w = adj[u][j];
                if (!seen[w]) {
                    seen[w] = true;
                    parents[w] = u;
                    bfs.push_back(w);
                }
            }
        }
        for (int i = 0; i < n; ++i)
            if (!seen[i]) { bad(); return 0; }
    }

    HLD h;
    h.build(n, parents, vals);
    for (int k = 0; k < q; ++k) {
        ll op, u, v, x;
        if (!next_ll(op)) { bad(); return 0; }
        if (op == 1) {
            if (!next_ll(u) || !next_ll(v) || !next_ll(x) ||
                u < 1 || u > n || v < 1 || v > n || x > 1000000000LL ||
                x < -1000000000LL) { bad(); return 0; }
            h.path_add((int)u - 1, (int)v - 1, x);
        } else if (op == 2 || op == 3) {
            if (!next_ll(u) || !next_ll(v) ||
                u < 1 || u > n || v < 1 || v > n) { bad(); return 0; }
            ll r = (op == 2) ? h.path_max((int)u - 1, (int)v - 1)
                             : h.path_sum((int)u - 1, (int)v - 1);
            char buf[32];
            sprintf(buf, "%lld\n", r);
            emit(buf);
        } else {
            bad();
            return 0;
        }
    }
    return 0;
}
