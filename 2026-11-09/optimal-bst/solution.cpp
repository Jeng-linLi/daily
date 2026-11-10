// -*- coding: utf-8 -*-
// 最優二叉搜尋樹 Optimal BST 與 Knuth 四邊形不等式優化
//
// 題意：n 個已排序的鍵 1..n，第 i 個鍵被搜尋的權重為 p[i]。
//      總成本 = Σ p[i] × depth(key i)，求最小總成本與整樹根。
//
// 樸素區間 DP（O(n³)）：
//   dp[i][j] = min_{r=i..j} { dp[i][r-1] + dp[r+1][j] } + w[i][j]
// Knuth 優化（O(n²)）：前綴和權重滿足四邊形不等式 ⟹
//   s[i][j-1] ≤ s[i][j] ≤ s[i+1][j]，根只在窗口內列舉。
//
// IO 模式（stdin 有內容時）：
//   n
//   n 個整數 p[1..n]（1 ≤ n ≤ 3000，0 ≤ p[i] ≤ 10^9）
//   輸出兩行：最小成本、整樹根 s[1][n]。n = 0 輸出 "0 0"。
// 無 stdin 時跑內建斷言（O(n³) 樸素 DP 對拍 + s 單調性檢查），
// 輸出 all tests passed。
#include <cstdio>
#include <vector>
#include <algorithm>

using namespace std;
typedef long long ll;
typedef unsigned long long ull;

static ull g_seed;

static inline unsigned rnd(unsigned n) {
    // 與 Python 版相同的 LCG（mod 2^31），保證兩語言測試序列一致
    g_seed = (g_seed * 1103515245ULL + 12345ULL) & 0x7FFFFFFFULL;
    return (unsigned)(g_seed % n);
}

// ---- Knuth O(n²) 版 ----
// p 為 1-indexed，長度 n+1。回傳 pair(最小成本, s[1][n])。
static pair<ll, int> optbstKnuth(const vector<ll>& p) {
    int n = (int)p.size() - 1;
    if (n == 0) return make_pair(0LL, 0);
    vector<vector<ll>> w(n + 1, vector<ll>(n + 1, 0));
    vector<vector<ll>> dp(n + 2, vector<ll>(n + 2, 0));
    vector<vector<int>> s(n + 1, vector<int>(n + 1, 0));
    for (int i = 1; i <= n; ++i) {
        w[i][i] = p[i];
        dp[i][i] = p[i];
        s[i][i] = i;
    }
    for (int len = 2; len <= n; ++len) {
        for (int i = 1; i + len - 1 <= n; ++i) {
            int j = i + len - 1;
            w[i][j] = w[i][j - 1] + p[j];
            int lo = s[i][j - 1], hi = s[i + 1][j];
            ll best = -1;
            int bestr = lo;
            for (int r = lo; r <= hi; ++r) {
                ll cost = dp[i][r - 1] + dp[r + 1][j];
                if (best < 0 || cost < best) {
                    best = cost;
                    bestr = r;
                }
            }
            dp[i][j] = best + w[i][j];
            s[i][j] = bestr;
        }
    }
    return make_pair(dp[1][n], s[1][n]);
}

// ---- 樸素 O(n³) 版（對拍用，另回傳 dp 表供 _rootCost 使用）----
static pair<ll, vector<vector<ll>>> optbstBruteDP(const vector<ll>& p) {
    int n = (int)p.size() - 1;
    vector<vector<ll>> w(n + 1, vector<ll>(n + 1, 0));
    vector<vector<ll>> dp(n + 2, vector<ll>(n + 2, 0));
    for (int i = 1; i <= n; ++i) {
        w[i][i] = p[i];
        dp[i][i] = p[i];
    }
    for (int len = 2; len <= n; ++len) {
        for (int i = 1; i + len - 1 <= n; ++i) {
            int j = i + len - 1;
            w[i][j] = w[i][j - 1] + p[j];
            ll best = -1;
            for (int r = i; r <= j; ++r) {
                ll cost = dp[i][r - 1] + dp[r + 1][j];
                if (best < 0 || cost < best) best = cost;
            }
            dp[i][j] = best + w[i][j];
        }
    }
    return make_pair(n == 0 ? 0LL : dp[1][n], dp);
}

// 以指定根 r 組整棵樹（左右子樹各取最優）的總成本
static ll rootCost(const vector<ll>& p, int root, const vector<vector<ll>>& dp) {
    int n = (int)p.size() - 1;
    vector<vector<ll>> w(n + 1, vector<ll>(n + 1, 0));
    for (int i = 1; i <= n; ++i) w[i][i] = p[i];
    for (int len = 2; len <= n; ++len)
        for (int i = 1; i + len - 1 <= n; ++i) {
            int j = i + len - 1;
            w[i][j] = w[i][j - 1] + p[j];
        }
    return w[1][n] + dp[1][root - 1] + dp[root + 1][n];
}

// 重建 Knuth 的 s 表並檢查單調性 s[i][j-1] ≤ s[i][j] ≤ s[i+1][j]
static void checkSMonotone(const vector<ll>& p) {
    int n = (int)p.size() - 1;
    vector<vector<ll>> w(n + 1, vector<ll>(n + 1, 0));
    vector<vector<ll>> dp(n + 2, vector<ll>(n + 2, 0));
    vector<vector<int>> s(n + 1, vector<int>(n + 1, 0));
    for (int i = 1; i <= n; ++i) {
        w[i][i] = p[i];
        dp[i][i] = p[i];
        s[i][i] = i;
    }
    for (int len = 2; len <= n; ++len) {
        for (int i = 1; i + len - 1 <= n; ++i) {
            int j = i + len - 1;
            w[i][j] = w[i][j - 1] + p[j];
            int lo = s[i][j - 1], hi = s[i + 1][j];
            ll best = -1;
            int bestr = lo;
            for (int r = lo; r <= hi; ++r) {
                ll cost = dp[i][r - 1] + dp[r + 1][j];
                if (best < 0 || cost < best) {
                    best = cost;
                    bestr = r;
                }
            }
            dp[i][j] = best + w[i][j];
            s[i][j] = bestr;
            if (!(s[i][j - 1] <= s[i][j] && s[i][j] <= s[i + 1][j])) {
                fprintf(stderr, "s monotonicity broken\n");
                exit(1);
            }
        }
    }
}

static void runTests() {
    g_seed = 20261109ULL;

    // 1) 固定用例（期望值由實跑取得）
    {
        vector<ll> p = {0, 2, 3, 1, 4};
        pair<ll, vector<vector<ll>>> br = optbstBruteDP(p);
        pair<ll, int> kn = optbstKnuth(p);
        if (!(br.first == 18 && kn.first == 18 && kn.second == 2)) {
            fprintf(stderr, "fixed case 1 failed: %lld %lld %d\n", br.first, kn.first, kn.second);
            exit(1);
        }
        if (!(rootCost(p, 2, br.second) == 18 && rootCost(p, 3, br.second) == 21)) {
            fprintf(stderr, "fixed rootCost failed\n");
            exit(1);
        }
    }
    {
        vector<ll> p = {0, 5};
        pair<ll, int> kn = optbstKnuth(p);
        if (!(kn.first == 5 && kn.second == 1)) {
            fprintf(stderr, "fixed case 2 failed\n");
            exit(1);
        }
    }
    {
        vector<ll> p = {0, 1, 1, 1};
        pair<ll, vector<vector<ll>>> br = optbstBruteDP(p);
        pair<ll, int> kn = optbstKnuth(p);
        if (!(br.first == 5 && kn.first == 5 && kn.second == 2)) {
            fprintf(stderr, "fixed case 3 failed\n");
            exit(1);
        }
    }

    // 2) 隨機對拍：Knuth 成本 = 樸素 O(n³) 成本；Knuth 根為最優根
    for (int tc = 0; tc < 300; ++tc) {
        int n = 1 + rnd(25);
        vector<ll> p(n + 1, 0);
        for (int i = 1; i <= n; ++i) p[i] = rnd(21);
        pair<ll, vector<vector<ll>>> br = optbstBruteDP(p);
        pair<ll, int> kn = optbstKnuth(p);
        if (br.first != kn.first) {
            fprintf(stderr, "case %d mismatch: %lld vs %lld\n", tc, br.first, kn.first);
            exit(1);
        }
        if (rootCost(p, kn.second, br.second) != br.first) {
            fprintf(stderr, "case %d root not optimal\n", tc);
            exit(1);
        }
    }

    // 3) s 單調性
    for (int tc = 0; tc < 50; ++tc) {
        int n = 2 + rnd(29);
        vector<ll> p(n + 1, 0);
        for (int i = 1; i <= n; ++i) p[i] = rnd(21);
        checkSMonotone(p);
    }
}

int main() {
    // 判斷是否有 stdin 內容：先嘗試讀第一個整數
    int first;
    if (scanf("%d", &first) != 1) {
        runTests();
        printf("all tests passed\n");
        return 0;
    }
    int n = first;
    if (n == 0) {
        printf("0 0\n");
        return 0;
    }
    vector<ll> p(n + 1, 0);
    for (int i = 1; i <= n; ++i) scanf("%lld", &p[i]);
    pair<ll, int> kn = optbstKnuth(p);
    printf("%lld\n", kn.first);
    printf("%d\n", kn.second);
    return 0;
}
