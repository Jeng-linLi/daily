// 長鏈剖分（Long-Path Decomposition）：O(1) 求 k 級祖先
//
// 每個節點往「子樹深度最深的兒子」（長兒子）延伸一條長鏈。
// 對每條鏈的鏈頂 top 預處理：
//   down[top][i] : 沿鏈向下走 i 步的節點（0..L，L 為鏈長）
//   up[top][j]   : top 的 j 級祖先（j = 0..min(L, dep[top])）
// 查詢 u 的 k 級祖先（k <= len[u]，len[u] = u 往下最長距離 = L - d）可 O(1)：
//   d = dep[u] - dep[top[u]]
//   k <= d          -> down[top[u]][d - k]
//   d < k <= len[u] -> up[top[u]][k - d]   （k - d <= dep[top]，不會越界）
// k > len[u] 時退化到倍增表（O(log n)）。
//
// IO 格式與 solution.py 完全一致（見同目錄 README.md）。
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>
#ifdef _WIN32
#include <io.h>
#include <fcntl.h>
#endif

static int n_;

// ---------- 內建測試用的同一套邏輯（也供 IO 模式呼叫） ----------
struct LongPath {
    std::vector<int> par, dep, maxlen, top, hson;
    std::vector<std::vector<int>> down, up;   // 以鏈頂編號索引
    std::vector<std::vector<int>> upk;        // 倍增表
    int LOG;

    void build(int n, const std::vector<int>& par_) {
        n_ = n;
        par = par_;
        dep.assign(n, 0);
        maxlen.assign(n, 0);
        for (int i = 1; i < n; ++i) dep[i] = dep[par[i]] + 1;
        for (int i = n - 1; i >= 1; --i) {
            int p = par[i];
            if (maxlen[i] + 1 > maxlen[p]) maxlen[p] = maxlen[i] + 1;
        }
        hson.assign(n, -1);
        for (int i = 1; i < n; ++i) {          // 平手取編號最小的兒子
            int p = par[i];
            if (hson[p] < 0 && maxlen[i] + 1 == maxlen[p]) hson[p] = i;
        }
        top.assign(n, 0);
        for (int i = 1; i < n; ++i)
            top[i] = (i == hson[par[i]]) ? top[par[i]] : i;
        down.assign(n, {});
        up.assign(n, {});
        for (int t = 0; t < n; ++t) {
            if (top[t] != t) continue;
            int L = maxlen[t];
            std::vector<int> chain;
            chain.push_back(t);
            int cur = t;
            while (hson[cur] >= 0) {
                cur = hson[cur];
                chain.push_back(cur);
            }
            down[t] = chain;                    // 長度 L + 1
            int cnt = L < dep[t] ? L : dep[t];
            std::vector<int> anc;
            anc.push_back(t);
            cur = t;
            for (int j = 0; j < cnt; ++j) {
                cur = par[cur];
                anc.push_back(cur);
            }
            up[t] = anc;                        // up[t][j] = t 的 j 級祖先
        }
        LOG = 1;
        while ((1 << LOG) < n) ++LOG;
        if (LOG < 1) LOG = 1;
        upk.assign(LOG, std::vector<int>(n));
        upk[0] = par;
        for (int j = 1; j < LOG; ++j)
            for (int i = 0; i < n; ++i)
                upk[j][i] = upk[j - 1][upk[j - 1][i]];
    }

    int kth(int u, long long k) const {
        if (k == 0) return u;
        if (k > dep[u]) return -1;
        int t = top[u];
        int d = dep[u] - dep[t];
        if (k <= d) return down[t][d - k];
        if (k <= maxlen[u]) return up[t][k - d];
        int v = u, j = 0;
        long long kk = k;
        while (kk) {
            if (kk & 1) v = upk[j][v];
            kk >>= 1;
            ++j;
        }
        return v;
    }
};

// ---------- 內建測試 ----------
struct RNG {
    unsigned int s;
    explicit RNG(unsigned int seed) : s(seed) {}
    int next(int m) {
        s = (s * 1103515245u + 12345u) & 0x7FFFFFFFu;
        return (int)(s % (unsigned int)m);
    }
};

static int maxlenBrute(int u, const std::vector<std::vector<int>>& ch) {
    int best = 0;
    for (size_t i = 0; i < ch[u].size(); ++i) {
        int cand = maxlenBrute(ch[u][i], ch) + 1;
        if (cand > best) best = cand;
    }
    return best;
}

static int topBrute(int u, const std::vector<int>& par,
                    const std::vector<int>& hson) {
    int v = u;
    while (par[v] != v && hson[par[v]] == v) v = par[v];
    return v;
}

static void runTests() {
    RNG rng(20261106u);
    std::vector<std::vector<int>> cases;
    int fixed[][13] = {
        {6, 0, 0, 1, 2, 3, 4, -1, -1, -1, -1, -1, -1},              // 鏈
        {5, 0, 0, 0, 0, 0, -1, -1, -1, -1, -1, -1, -1},             // 星形
        {9, 0, 0, 1, 1, 2, 2, 3, 3, 4, -1, -1, -1},                 // 二叉樹
        {9, 0, 0, 1, 1, 2, 3, 4, 5, 6, -1, -1, -1},                 // 毛毛蟲
        {12, 0, 0, 0, 1, 3, 4, 5, 2, 7, 8, 9, 10},                  // 雙長鏈
    };
    for (int c = 0; c < 5; ++c) {
        std::vector<int> par;
        for (int j = 1; j <= fixed[c][0]; ++j) par.push_back(fixed[c][j]);
        cases.push_back(par);
    }
    for (int t = 0; t < 200; ++t) {
        int n = rng.next(200) + 1;
        std::vector<int> par(n, 0);
        for (int i = 1; i < n; ++i) par[i] = rng.next(i);
        cases.push_back(par);
    }

    for (size_t c = 0; c < cases.size(); ++c) {
        const std::vector<int>& par = cases[c];
        int n = (int)par.size();
        LongPath lp;
        lp.build(n, par);
        std::vector<std::vector<int>> ch(n);
        for (int i = 1; i < n; ++i) ch[par[i]].push_back(i);
        for (int u = 0; u < n; ++u) {
            int v = u, steps = 0;                    // dep 暴力
            while (v != 0) { v = par[v]; ++steps; }
            if (lp.dep[u] != steps) { printf("FAIL dep\n"); exit(1); }
            if (lp.maxlen[u] != maxlenBrute(u, ch)) { printf("FAIL maxlen\n"); exit(1); }
            if (lp.top[u] != topBrute(u, par, lp.hson)) { printf("FAIL top\n"); exit(1); }
            for (int k = 0; k <= lp.dep[u] + 2; ++k) {
                int want = -1;                       // kth 暴力
                if (k <= lp.dep[u]) {
                    want = u;
                    for (int s = 0; s < k; ++s) want = par[want];
                }
                int got = lp.kth(u, k);
                if (got != want) { printf("FAIL kth u=%d k=%d\n", u, k); exit(1); }
            }
        }
    }
    // 鏈長合計 = n
    for (size_t c = 0; c < 20 && c < cases.size(); ++c) {
        const std::vector<int>& par = cases[c];
        int n = (int)par.size();
        LongPath lp;
        lp.build(n, par);
        int total = 0;
        for (int t = 0; t < n; ++t)
            if (lp.top[t] == t) total += (int)lp.down[t].size();
        if (total != n) { printf("FAIL chain total\n"); exit(1); }
    }
    printf("all tests passed\n");
}

int main() {
#ifdef _WIN32
    _setmode(_fileno(stdout), _O_BINARY);
#endif
    // 嘗試讀取 stdin：先讀首個 token 判斷模式
    long long first;
    if (scanf("%lld", &first) != 1) {
        runTests();
        return 0;
    }
    int n = (int)first;
    long long q;
    if (scanf("%lld", &q) != 1) { runTests(); return 0; }
    std::vector<int> par(n, 0);
    for (int i = 1; i < n; ++i)
        if (scanf("%d", &par[i]) != 1) { runTests(); return 0; }
    LongPath lp;
    lp.build(n, par);
    std::string out;
    out.reserve((size_t)q * 8);
    char buf[32];
    for (long long t = 0; t < q; ++t) {
        int op, u;
        if (scanf("%d %d", &op, &u) != 2) break;
        if (op == 1) {
            long long k;
            scanf("%lld", &k);
            sprintf(buf, "%d\n", lp.kth(u, k));
        } else if (op == 2) {
            sprintf(buf, "%d\n", lp.dep[u]);
        } else if (op == 3) {
            sprintf(buf, "%d\n", lp.maxlen[u]);
        } else {
            sprintf(buf, "%d\n", lp.top[u]);
        }
        out += buf;
    }
    fputs(out.c_str(), stdout);
    return 0;
}
