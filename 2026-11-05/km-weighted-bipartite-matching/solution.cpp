// KM 演算法：二分圖最大權完美匹配（Kuhn–Munkres / 匈牙利算法帶權版）
//
// 題意：
//   給定 n × n 的權值矩陣 w[i][j]（可為負），求完美匹配的最大權值和。
//   輸入：第一行 n；接下來 n 行、每行 n 個整數。
//   輸出：一個整數。
//   無 stdin 時執行內建斷言（枚舉排列暴力對拍），輸出 all tests passed。
//
// 思路：
//   頂標（potentials）u[i] + v[j] >= w[i][j]，只走相等子圖的緊繃邊找
//   增廣路；找不到時取最小鬆弛量 delta 調整頂標，相等子圖擴大至少一條邊。
//   每個左部點增廣一次、每次最多 O(n) 輪調整、每輪掃描 O(n)，總 O(n^3)。
//   逐步增廣與頂標調整的細節與 Python 版完全一致（e-maxx 潛能版），
//   保證兩語言對同一輸入輸出一致。

#include <cstdio>
#include <cstdlib>
#include <cctype>
#include <cstdint>
#include <vector>
#include <string>

#ifdef _WIN32
#include <io.h>
#include <fcntl.h>
#endif

static const long long INF = (long long)1 << 62;
static const long long NMAX = 300;
static const long long W_LIM = 1000000000LL;
static const long long VAL_LIM = 1000000000000000000LL;

// KM 主體：w 為 (n+1) x (n+1)，1-indexed 使用，回傳最大權完美匹配總權值
// 可行性不變式：u[i] + v[j] >= w[i][j]；初始 u[i] = max_j w[i][j]、v = 0
static long long kmMaxWeightMatching(const std::vector<std::vector<long long>>& w, int n) {
    std::vector<long long> u(n + 1, 0), v(n + 1, 0);
    std::vector<int> p(n + 1, 0), way(n + 1, 0);   // p[j]: 與右部 j 匹配的左部點
    for (int i = 1; i <= n; ++i) {
        long long mx = -(INF / 4);
        for (int j = 1; j <= n; ++j) mx = mx > w[i][j] ? mx : w[i][j];
        u[i] = mx;
    }
    for (int i = 1; i <= n; ++i) {
        p[0] = i;                                   // 虛擬右部 0 掛上當前左部點
        int j0 = 0;
        std::vector<long long> minv(n + 1, INF);
        std::vector<char> used(n + 1, 0);
        do {
            used[j0] = 1;
            int i0 = p[j0];
            long long delta = INF;
            int j1 = 0;
            for (int j = 1; j <= n; ++j) {
                if (!used[j]) {
                    long long cur = u[i0] + v[j] - w[i0][j];   // 鬆弛量 >= 0
                    if (cur < minv[j]) { minv[j] = cur; way[j] = j0; }
                    if (minv[j] < delta) { delta = minv[j]; j1 = j; }
                }
            }
            // 頂標調整：訪問側 u -= delta、v += delta；未訪問側鬆弛量同步縮小
            for (int j = 0; j <= n; ++j) {
                if (used[j]) {
                    u[p[j]] -= delta;
                    v[j] += delta;
                } else {
                    minv[j] -= delta;
                }
            }
            j0 = j1;
        } while (p[j0] != 0);                       // j0 尚未匹配 → 增廣路打通
        do {                                        // 沿 way 回溯翻轉匹配邊
            int j1 = way[j0];
            p[j0] = p[j1];
            j0 = j1;
        } while (j0);
    }
    long long ans = 0;
    for (int j = 1; j <= n; ++j) ans += w[p[j]][j];
    return ans;
}

// --------------------------------------------------------------- 暴力對拍

// 枚舉右部點的全排列，O(n · n!)
static long long bruteMax(const std::vector<std::vector<long long>>& w, int n,
                          std::vector<int>& perm, int depth, long long cur, long long& best) {
    if (n == 0) return 0;
    if (depth > n) {
        if (cur > best) best = cur;
        return best;
    }
    for (int j = 1; j <= n; ++j) {
        if (!perm[j]) {
            perm[j] = 1;
            bruteMax(w, n, perm, depth + 1, cur + w[depth][j], best);
            perm[j] = 0;
        }
    }
    return best;
}

static unsigned int g_rng = 12345;
static unsigned int rnd() {
    g_rng = g_rng * 1103515245u + 12345u;
    return g_rng;
}

static void runTests() {
    // 固定用例 1：經典指派
    {
        std::vector<std::vector<long long>> w = {
            {0, 0, 0, 0, 0},
            {0, 4, 1, 3, 2},
            {0, 2, 0, 5, 3},
            {0, 3, 2, 2, 3},
            {0, 4, 0, 2, 1}};
        std::vector<int> perm(5, 0);
        long long best = -(VAL_LIM * 2);
        long long want = bruteMax(w, 4, perm, 1, 0, best);
        long long got = kmMaxWeightMatching(w, 4);
        if (got != want) {
            fprintf(stderr, "case1 got=%lld want=%lld\n", got, want);
            exit(1);
        }
    }
    // 固定用例 2：含負權
    {
        std::vector<std::vector<long long>> w = {
            {0, 0, 0, 0},
            {0, -5, -2, 3},
            {0, 1, -7, 2},
            {0, 0, 4, -1}};
        std::vector<int> perm(4, 0);
        long long best = -(VAL_LIM * 2);
        long long want = bruteMax(w, 3, perm, 1, 0, best);
        long long got = kmMaxWeightMatching(w, 3);
        if (got != want) {
            fprintf(stderr, "case2 got=%lld want=%lld\n", got, want);
            exit(1);
        }
    }
    // 固定用例 3：全零矩陣
    {
        std::vector<std::vector<long long>> w = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
        if (kmMaxWeightMatching(w, 2) != 0) { fprintf(stderr, "case3 fail\n"); exit(1); }
    }
    // 固定用例 4：單點
    {
        std::vector<std::vector<long long>> w = {{0, 0}, {0, 7}};
        if (kmMaxWeightMatching(w, 1) != 7) { fprintf(stderr, "case4 fail\n"); exit(1); }
    }
    // 隨機對拍（n ≤ 6）
    for (int it = 0; it < 200; ++it) {
        int n = (int)(rnd() % 6) + 1;
        std::vector<std::vector<long long>> w(n + 1, std::vector<long long>(n + 1, 0));
        for (int i = 1; i <= n; ++i)
            for (int j = 1; j <= n; ++j)
                w[i][j] = (long long)(rnd() % 21) - 10;   // 權值範圍 [-10, 10]
        std::vector<int> perm(n + 1, 0);
        long long best = -(VAL_LIM * 2);
        long long want = bruteMax(w, n, perm, 1, 0, best);
        long long got = kmMaxWeightMatching(w, n);
        if (got != want) {
            fprintf(stderr, "random#%d n=%d got=%lld want=%lld\n", it, n, got, want);
            exit(1);
        }
    }
}

// --------------------------------------------------------------- IO 模式

static std::vector<std::string> g_toks;
static size_t g_pos = 0;

static const std::string* nextTok() {
    if (g_pos >= g_toks.size()) return nullptr;
    return &g_toks[g_pos++];
}

// 嚴格整數解析：^[+-]?[0-9]+$，超界回 nullptr（與 Python 版規則一致）
static bool parseLL(const std::string& s, long long* out) {
    if (s.empty()) return false;
    size_t i = 0;
    if (s[0] == '+' || s[0] == '-') i = 1;
    if (i == s.size()) return false;
    long long val = 0;
    bool neg = false;
    if (s[0] == '-') neg = true;
    for (; i < s.size(); ++i) {
        char c = s[i];
        if (c < '0' || c > '9') return false;
        if (val > VAL_LIM / 2) return false;    // 防溢位，先擋超大值
        val = val * 10 + (c - '0');
    }
    if (val > VAL_LIM) return false;
    *out = neg ? -val : val;
    return true;
}

int main() {
#ifdef _WIN32
    _setmode(_fileno(stdout), _O_BINARY);
#endif
    // 讀全部 stdin
    std::string all;
    {
        char buf[1 << 16];
        size_t r;
        while ((r = fread(buf, 1, sizeof(buf), stdin)) > 0) all.append(buf, r);
    }
    // 分詞
    {
        std::string cur;
        for (size_t i = 0; i <= all.size(); ++i) {
            char c = i < all.size() ? all[i] : ' ';
            if (c == ' ' || c == '\t' || c == '\n' || c == '\r' || c == '\v' || c == '\f') {
                if (!cur.empty()) { g_toks.push_back(cur); cur.clear(); }
            } else {
                cur.push_back(c);
            }
        }
    }
    if (g_toks.empty()) {         // 無輸入 → 內建測試
        runTests();
        printf("all tests passed\n");
        return 0;
    }
    long long n;
    const std::string* t = nextTok();
    if (t == nullptr || !parseLL(*t, &n) || n < 0 || n > NMAX) {
        printf("error\n");
        return 0;
    }
    int nn = (int)n;
    std::vector<std::vector<long long>> w(nn + 1, std::vector<long long>(nn + 1, 0));
    for (int i = 1; i <= nn; ++i) {
        for (int j = 1; j <= nn; ++j) {
            long long x;
            const std::string* tj = nextTok();
            if (tj == nullptr || !parseLL(*tj, &x) || x < -W_LIM || x > W_LIM) {
                printf("error\n");
                return 0;
            }
            w[i][j] = x;
        }
    }
    printf("%lld\n", kmMaxWeightMatching(w, nn));
    return 0;
}
