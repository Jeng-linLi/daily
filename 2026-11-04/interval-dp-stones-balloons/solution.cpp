// 區間 DP：合併石子與戳氣球
//
// 題意：
//     輸入第一行為任務編號 T（1 或 2）：
//     T=1（合併石子）：n 堆石子排成一排，每次只能合併相鄰兩堆，
//        代價 = 兩堆石子數之和，求合成一堆的最小總代價。
//     T=2（戳氣球）：n 個氣球排成一排，戳破第 i 個獲得
//        nums[l-1]*nums[i]*nums[r+1]（l/r 為相鄰未戳破者，越界為 1），
//        求最大收益。
//     輸入：T、n、n 個整數。輸出：一個整數。
//
// 思路：
//     合併石子：dp[i][j] = min_k{ dp[i][k] + dp[k+1][j] } + sum(i..j)，
//     最後一步代價固定為整段區間和（前綴和 O(1)）。
//     戳氣球：把「先戳」改「最後戳」——dp[l][r] 為開區間 (l,r) 內
//     全戳破的最大收益，最後戳 k：dp[l][r] = max_k{ dp[l][k] +
//     dp[k][r] + a[l]*a[k]*a[r] }，兩端補虛擬氣球 1。
//
// 複雜度：兩個任務都是 O(n^3) 時間、O(n^2) 空間。
//
// 應用場景：
//     文件/區塊合併的最小 IO 代價、資料庫歸併成本、
//     「刪除元素收益與左右鄰有關」類問題（廣告位拆除等）、
//     矩陣鏈乘法與多邊形三角剖分（同構）。

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

const int NMAX = 300;
const ll LL_MAX = (ll)4e18;

// 合併石子最小總代價（O(n^3)）
static ll stone_merge_min(const std::vector<ll> &piles) {
    int n = (int)piles.size();
    if (n <= 1) return 0;
    std::vector<ll> pre(n + 1, 0);
    for (int i = 0; i < n; ++i) pre[i + 1] = pre[i] + piles[i];
    std::vector<std::vector<ll> > dp(n, std::vector<ll>(n, 0));
    for (int len = 2; len <= n; ++len) {
        for (int i = 0; i + len - 1 < n; ++i) {
            int j = i + len - 1;
            ll total = pre[j + 1] - pre[i];
            ll best = LL_MAX;
            for (int k = i; k < j; ++k) {
                ll cur = dp[i][k] + dp[k + 1][j];
                if (cur < best) best = cur;
            }
            dp[i][j] = best + total;
        }
    }
    return dp[0][n - 1];
}

// 戳氣球最大收益（開區間 DP，兩端補虛擬氣球 1）
static ll burst_balloons_max(const std::vector<ll> &nums) {
    int n = (int)nums.size();
    std::vector<ll> a(n + 2, 1);
    for (int i = 0; i < n; ++i) a[i + 1] = nums[i];
    int m = n + 2;
    std::vector<std::vector<ll> > dp(m, std::vector<ll>(m, 0));
    for (int len = 2; len < m; ++len) {
        for (int l = 0; l + len < m; ++l) {
            int r = l + len;
            ll best = 0;
            for (int k = l + 1; k < r; ++k) {
                ll cur = dp[l][k] + dp[k][r] + a[l] * a[k] * a[r];
                if (cur > best) best = cur;
            }
            dp[l][r] = best;
        }
    }
    return dp[0][m - 1];
}

// ---- 暴力（對拍用） ----
static ll stone_brute(std::vector<ll> piles) {
    int n = (int)piles.size();
    if (n <= 1) return 0;
    ll best = LL_MAX;
    for (int i = 0; i + 1 < n; ++i) {
        ll merged = piles[i] + piles[i + 1];
        std::vector<ll> rest;
        for (int j = 0; j < n; ++j) {
            if (j == i) rest.push_back(merged);
            else if (j != i + 1) rest.push_back(piles[j]);
        }
        ll cost = merged + stone_brute(rest);
        if (cost < best) best = cost;
    }
    return best;
}

static ll balloons_brute(std::vector<ll> nums) {
    int n = (int)nums.size();
    if (n == 0) return 0;
    ll best = 0;
    for (int i = 0; i < n; ++i) {
        ll left = (i > 0) ? nums[i - 1] : 1;
        ll right = (i + 1 < n) ? nums[i + 1] : 1;
        std::vector<ll> rest;
        for (int j = 0; j < n; ++j)
            if (j != i) rest.push_back(nums[j]);
        ll gain = left * nums[i] * right + balloons_brute(rest);
        if (gain > best) best = gain;
    }
    return best;
}

static void emit(const std::string &s) { fwrite(s.data(), 1, s.size(), stdout); }

// 自製 LCG（與 Python 版同一序列）
static unsigned int lcg_state = 20261104u;
static int nxt_mod(int mod) {
    lcg_state = (unsigned int)((ll)lcg_state * 1103515245LL + 12345LL) & 0x7FFFFFFFu;
    return (int)(lcg_state % (unsigned)mod);
}

static void run_tests() {
    // ---- 合併石子固定用例 ----
    if (stone_merge_min({}) != 0) { emit("FAIL 1\n"); return; }
    if (stone_merge_min({7}) != 0) { emit("FAIL 2\n"); return; }
    if (stone_merge_min({4, 5}) != 9) { emit("FAIL 3\n"); return; }
    if (stone_merge_min({1, 3, 5, 2}) != 22) { emit("FAIL 4\n"); return; }
    if (stone_merge_min({3, 1, 2, 4}) != 19) { emit("FAIL 5\n"); return; }
    // ---- 戳氣球固定用例 ----
    if (burst_balloons_max({}) != 0) { emit("FAIL 6\n"); return; }
    if (burst_balloons_max({5}) != 5) { emit("FAIL 7\n"); return; }
    if (burst_balloons_max({3, 1, 5, 8}) != 167) { emit("FAIL 8\n"); return; }
    if (burst_balloons_max({1, 5}) != 10) { emit("FAIL 9\n"); return; }
    // ---- 對拍 ----
    for (int it = 0; it < 150; ++it) {
        int n = 1 + nxt_mod(10);             // 1..10
        std::vector<ll> piles(n);
        for (int i = 0; i < n; ++i) piles[i] = 1 + nxt_mod(20);
        if (stone_merge_min(piles) != stone_brute(piles)) {
            emit("FAIL stone random\n");
            return;
        }
    }
    for (int it = 0; it < 120; ++it) {
        int n = 1 + nxt_mod(8);              // 1..8
        std::vector<ll> nums(n);
        for (int i = 0; i < n; ++i) nums[i] = nxt_mod(10);
        if (burst_balloons_max(nums) != balloons_brute(nums)) {
            emit("FAIL balloons random\n");
            return;
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
    ll t;
    if (pos >= toks.size() || !parse_ll(toks[pos++].c_str(), t) ||
        (t != 1 && t != 2)) { bad(); return 0; }
    ll n_ll;
    if (pos >= toks.size() || !parse_ll(toks[pos++].c_str(), n_ll) ||
        n_ll < 0 || n_ll > NMAX) { bad(); return 0; }
    int n = (int)n_ll;
    std::vector<ll> arr(n);
    for (int i = 0; i < n; ++i) {
        if (pos >= toks.size() || !parse_ll(toks[pos++].c_str(), arr[i]) ||
            arr[i] < 0 || arr[i] > 1000000000LL) { bad(); return 0; }
    }
    char buf[32];
    ll ans = (t == 1) ? stone_merge_min(arr) : burst_balloons_max(arr);
    sprintf(buf, "%lld\n", ans);
    emit(buf);
    return 0;
}
