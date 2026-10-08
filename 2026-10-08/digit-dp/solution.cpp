// 數位 DP（Digit DP）：數位和計數 / 禁用數字 / 吉利數字 / 數位和總和 / Windy Number / 不要 62
//
// 與 solution.py 完全同構：同樣的演算法、同樣的確定性規則、同樣的輸入輸出格式。
// 編譯：g++ -std=c++17 -O2 -Wall solution.cpp -o solution
// 無 stdin 輸入時運行內置斷言測試並輸出 `all tests passed`。

#include <algorithm>
#include <cassert>
#include <cctype>
#include <cstdlib>
#include <iostream>
#include <random>
#include <set>
#include <sstream>
#include <string>
#include <vector>

using namespace std;

using ll = long long;

static const ll NMAX = (ll)1e15;    // n 的上限（與 Python 一致，避免溢位）
static const int NO_DIGIT = 10;     // 「尚未開始（全是前導零）」的哨兵值

// ---------------------------------------------------------------- 工具

// 嚴格整數規則 ^[+-]?[0-9]+$，與 Python 的 parse_int 完全一致。
static bool tryLL(const string &s, ll &out) {
    if (s.empty()) return false;
    size_t i = 0;
    if (s[0] == '+' || s[0] == '-') i = 1;
    if (i >= s.size()) return false;
    ll val = 0;
    for (; i < s.size(); i++) {
        if (!isdigit((unsigned char)s[i])) return false;
        val = val * 10 + (s[i] - '0');
        if (val > NMAX) val = NMAX;          // 截斷，避免溢位
    }
    out = (s[0] == '-') ? -val : val;
    return true;
}

// 把非負整數拆成十進位數字串（高位在前）。
static vector<int> digitsOf(ll n) {
    string s = to_string(n);
    vector<int> ds;
    for (char c : s) ds.push_back(c - '0');
    return ds;
}

// 十進位數位和。
static int digitSum(ll x) {
    string s = to_string(x);
    int t = 0;
    for (char c : s) t += c - '0';
    return t;
}

// ---------------------------------------------------------------- 數位 DP

// 統計 [0, n] 中數位和等於 k 的數的個數。狀態 (sum, tight)。
static ll countDigitSum(ll n, int k) {
    if (n < 0 || k < 0) return 0;
    vector<int> ds = digitsOf(n);
    int L = (int)ds.size();
    if (k > 9 * L) return 0;
    // dp[sum][tight]，sum ≤ k
    vector<vector<ll>> dp(k + 1, vector<ll>(2, 0));
    dp[0][1] = 1;
    for (int pos = 0; pos < L; pos++) {
        vector<vector<ll>> nd(k + 1, vector<ll>(2, 0));
        for (int s = 0; s <= k; s++) {
            for (int tight = 0; tight < 2; tight++) {
                ll cnt = dp[s][tight];
                if (cnt == 0) continue;
                int lim = tight ? ds[pos] : 9;
                for (int d = 0; d <= lim; d++) {
                    int ns = s + d;
                    if (ns > k) continue;            // 剪枝：數位和已超過目標
                    int ntight = (tight && d == lim) ? 1 : 0;
                    nd[ns][ntight] += cnt;
                }
            }
        }
        dp = nd;
    }
    return dp[k][0] + dp[k][1];
}

// 統計 [1, n] 中每一位都屬於 allowed 的正整數個數（不含 0）。狀態 (started, tight)。
static ll countDigitSet(ll n, const set<int> &allowed) {
    if (n <= 0) return 0;
    vector<int> ds = digitsOf(n);
    ll dp[2][2] = {{0, 0}, {0, 0}};                  // [started][tight]
    dp[0][1] = 1;
    for (int pos = 0; pos < (int)ds.size(); pos++) {
        ll nd[2][2] = {{0, 0}, {0, 0}};
        for (int started = 0; started < 2; started++)
            for (int tight = 0; tight < 2; tight++) {
                ll cnt = dp[started][tight];
                if (cnt == 0) continue;
                int lim = tight ? ds[pos] : 9;
                for (int d = 0; d <= lim; d++) {
                    int ntight = (tight && d == lim) ? 1 : 0;
                    if (!started && d == 0) {
                        nd[0][ntight] += cnt;        // 仍是前導零，數字尚未開始
                    } else if (allowed.count(d)) {
                        nd[1][ntight] += cnt;
                    }
                }
            }
        for (int a = 0; a < 2; a++)
            for (int b = 0; b < 2; b++) dp[a][b] = nd[a][b];
    }
    return dp[1][0] + dp[1][1];
}

// [1, n] 中的 Lucky Number（每一位都是 4 或 7）個數。
static ll countLucky(ll n) { return countDigitSet(n, set<int>{4, 7}); }

// [0, n] 中不含數字 forbid 的數的個數。
// 數字 0 的表示是 "0"：只有當禁用的是 0 時才把 0 排除，其餘情況 0 都要計入。
static ll countAvoidDigit(ll n, int forbid) {
    if (n < 0) return 0;
    set<int> allowed;
    for (int x = 0; x <= 9; x++)
        if (x != forbid) allowed.insert(x);
    return countDigitSet(n, allowed) + (forbid == 0 ? 0 : 1);
}

// Σ_{i=0}^{n} digitsum(i)。狀態只需 tight，每格存 (個數, 數位和總和)。
static ll digitSumTotal(ll n) {
    if (n < 0) return 0;
    vector<int> ds = digitsOf(n);
    ll cnt[2] = {0, 1};                              // [tight] = 個數
    ll tot[2] = {0, 0};                              // [tight] = 數位和總和
    for (int pos = 0; pos < (int)ds.size(); pos++) {
        ll ncnt[2] = {0, 0}, ntot[2] = {0, 0};
        for (int tight = 0; tight < 2; tight++) {
            if (cnt[tight] == 0) continue;
            int lim = tight ? ds[pos] : 9;
            for (int d = 0; d <= lim; d++) {
                int nt = (tight && d == lim) ? 1 : 0;
                ncnt[nt] += cnt[tight];
                ntot[nt] += tot[tight] + (ll)d * cnt[tight];
            }
        }
        for (int t = 0; t < 2; t++) {
            cnt[t] = ncnt[t];
            tot[t] = ntot[t];
        }
    }
    return tot[0] + tot[1];
}

// [1, n] 中的 Windy Number（相鄰位數字差 ≥ 2）個數。狀態 (上一位, tight)。
static ll countWindy(ll n) {
    if (n < 0) return 0;
    vector<int> ds = digitsOf(n);
    ll dp[11][2] = {{0, 0}};                         // [上一位或 NO_DIGIT][tight]
    dp[NO_DIGIT][1] = 1;
    for (int pos = 0; pos < (int)ds.size(); pos++) {
        ll nd[11][2] = {{0, 0}};
        for (int prev = 0; prev <= NO_DIGIT; prev++)
            for (int tight = 0; tight < 2; tight++) {
                ll c = dp[prev][tight];
                if (c == 0) continue;
                int lim = tight ? ds[pos] : 9;
                for (int d = 0; d <= lim; d++) {
                    int ntight = (tight && d == lim) ? 1 : 0;
                    if (prev == NO_DIGIT && d == 0) {
                        nd[NO_DIGIT][ntight] += c;   // 尚未開始，仍是前導零
                    } else {
                        if (prev != NO_DIGIT && abs(d - prev) < 2) continue;
                        nd[d][ntight] += c;
                    }
                }
            }
        for (int a = 0; a <= NO_DIGIT; a++)
            for (int b = 0; b < 2; b++) dp[a][b] = nd[a][b];
    }
    ll ans = 0;
    for (int prev = 0; prev <= 9; prev++) ans += dp[prev][0] + dp[prev][1];
    return ans;
}

// [0, n] 中不含數字 4、且不含連續 "62" 的數的個數（經典題「不要 62」）。
static ll countAuspicious(ll n) {
    if (n < 0) return 0;
    vector<int> ds = digitsOf(n);
    ll dp[2][2] = {{0, 0}, {0, 0}};                  // [上一位是否為 6][tight]
    dp[0][1] = 1;
    for (int pos = 0; pos < (int)ds.size(); pos++) {
        ll nd[2][2] = {{0, 0}, {0, 0}};
        for (int prev6 = 0; prev6 < 2; prev6++)
            for (int tight = 0; tight < 2; tight++) {
                ll c = dp[prev6][tight];
                if (c == 0) continue;
                int lim = tight ? ds[pos] : 9;
                for (int d = 0; d <= lim; d++) {
                    if (d == 4) continue;            // 禁用數字 4
                    if (prev6 && d == 2) continue;   // 禁用連續子串 "62"
                    int ntight = (tight && d == lim) ? 1 : 0;
                    nd[d == 6 ? 1 : 0][ntight] += c;
                }
            }
        for (int a = 0; a < 2; a++)
            for (int b = 0; b < 2; b++) dp[a][b] = nd[a][b];
    }
    return dp[0][0] + dp[0][1] + dp[1][0] + dp[1][1];
}

// ---------------------------------------------------------------- 暴力對拍

static ll bfDigitSum(ll n, int k) {
    ll c = 0;
    for (ll i = 0; i <= n; i++)
        if (digitSum(i) == k) c++;
    return c;
}

static ll bfAvoid(ll n, int forbid) {
    if (forbid < 0 || forbid > 9) return n + 1;      // 禁用值不是一個數字 → 不禁用任何數字
    string f = to_string(forbid);
    ll c = 0;
    for (ll i = 0; i <= n; i++)
        if (to_string(i).find(f) == string::npos) c++;
    return c;
}

static ll bfLucky(ll n) {
    ll c = 0;
    for (ll i = 1; i <= n; i++) {
        string s = to_string(i);
        bool ok = true;
        for (char ch : s)
            if (ch != '4' && ch != '7') ok = false;
        if (ok) c++;
    }
    return c;
}

static ll bfSumDigits(ll n) {
    ll t = 0;
    for (ll i = 0; i <= n; i++) t += digitSum(i);
    return t;
}

static ll bfWindy(ll n) {
    ll c = 0;
    for (ll i = 1; i <= n; i++) {
        string s = to_string(i);
        bool ok = true;
        for (int j = 0; j + 1 < (int)s.size(); j++)
            if (abs(s[j] - s[j + 1]) < 2) ok = false;
        if (ok) c++;
    }
    return c;
}

static ll bfAuspicious(ll n) {
    ll c = 0;
    for (ll i = 0; i <= n; i++) {
        string s = to_string(i);
        if (s.find('4') == string::npos && s.find("62") == string::npos) c++;
    }
    return c;
}

// ---------------------------------------------------------------- IO 與測試

static void runIO(const string &raw) {
    istringstream iss(raw);
    vector<string> toks;
    string tk;
    while (iss >> tk) toks.push_back(tk);
    size_t pos = 0;
    auto nxt = [&]() -> ll {
        ll v = 0;
        if (pos < toks.size()) {
            ll parsed = 0;
            if (tryLL(toks[pos], parsed)) v = parsed;
        }
        pos++;
        return v;
    };

    ll n = nxt(), k = nxt(), d = nxt();
    if (n < 0) n = 0;
    if (k < 0) k = 0;
    if (k > 200) k = 200;                            // 與 Python 的「k > 9L 直接為 0」等價保護

    ostringstream oss;
    oss << countDigitSum(n, (int)k) << "\n"
        << countAvoidDigit(n, (int)d) << "\n"
        << countLucky(n) << "\n"
        << digitSumTotal(n) << "\n"
        << countWindy(n) << "\n"
        << countAuspicious(n) << "\n";
    cout << oss.str();
}

static void runTests() {
    // ---- 退化情形 ----
    assert(countDigitSum(0, 0) == 1);                // 只有數字 0
    assert(countDigitSum(0, 1) == 0);
    assert(countDigitSum(-1, 0) == 0);
    assert(countDigitSum(9, 5) == 1);
    assert(countAvoidDigit(0, 4) == 1);              // [0,0] 不含 4 → 1 個
    assert(countAvoidDigit(0, 0) == 0);              // 0 含數字 0 → 0 個
    assert(countAvoidDigit(-1, 4) == 0);
    assert(countLucky(0) == 0);
    assert(countLucky(3) == 0);
    assert(countLucky(4) == 1);
    assert(countLucky(7) == 2);
    assert(countLucky(47) == 4);                     // 4, 7, 44, 47
    assert(digitSumTotal(0) == 0);
    assert(digitSumTotal(9) == 45);
    assert(countWindy(0) == 0);
    assert(countWindy(9) == 9);                      // 一位數都算 windy
    assert(countAuspicious(0) == 1);                 // 0 不含 4 也不含 62

    // ---- 固定用例 ----
    assert(countDigitSum(100, 1) == 3);              // 1, 10, 100
    assert(countDigitSum(100, 5) == 6);
    assert(countAvoidDigit(100, 4) == 82);           // 0..100 去掉含 4 的
    assert(digitSumTotal(100) == 901);
    assert(countWindy(100) == 73);
    assert(countAuspicious(100) == 81);              // 去掉含 4 與含 62 的

    // ---- 與暴力對拍（n ≤ 3000）----
    vector<ll> ns;
    for (ll i = 0; i <= 119; i++) ns.push_back(i);
    ns.push_back(199);
    ns.push_back(200);
    ns.push_back(999);
    ns.push_back(1000);
    ns.push_back(1999);
    ns.push_back(2000);
    ns.push_back(2999);
    ns.push_back(3000);
    for (ll n : ns) {
        for (int k = 0; k < 29; k++) assert(countDigitSum(n, k) == bfDigitSum(n, k));
        for (int d = 0; d <= 10; d++) assert(countAvoidDigit(n, d) == bfAvoid(n, d));
        assert(countLucky(n) == bfLucky(n));
        assert(digitSumTotal(n) == bfSumDigits(n));
        assert(countWindy(n) == bfWindy(n));
        assert(countAuspicious(n) == bfAuspicious(n));
    }

    // ---- 隨機對拍 ----
    mt19937 rng(20261008);
    for (int iter = 0; iter < 60; iter++) {
        ll n = rng() % 3001;
        int k = (int)(rng() % 31);
        int d = (int)(rng() % 11);
        assert(countDigitSum(n, k) == bfDigitSum(n, k));
        assert(countAvoidDigit(n, d) == bfAvoid(n, d));
        assert(countLucky(n) == bfLucky(n));
        assert(digitSumTotal(n) == bfSumDigits(n));
        assert(countWindy(n) == bfWindy(n));
        assert(countAuspicious(n) == bfAuspicious(n));
    }

    // ---- 大數自洽性（無法暴力，檢查封閉式結果與區間可減性）----
    ll big = (ll)1e15;
    assert(countDigitSum(big, 0) == 1);              // 只有 0
    assert(countDigitSum(big, 9 * 15) == 1);         // 只有 999999999999999
    assert(countLucky(big) == (1LL << 16) - 2);      // 長度 1..15 的 4/7 串
    assert(digitSumTotal(big - 1) == 15LL * 45 * (ll)1e14);
    assert(digitSumTotal(big) == 15LL * 45 * (ll)1e14 + 1);
    {
        set<int> allowed;
        for (int x = 0; x <= 9; x++)
            if (x != 4) allowed.insert(x);
        assert(countAvoidDigit(big, 4) == countDigitSet(big, allowed) + 1);
    }
    // 區間可減性：[L, R] = f(R) - f(L-1)
    {
        ll L0 = 12345, R0 = 98765;
        ll a = countWindy(R0) - countWindy(L0 - 1);
        ll b = 0;
        for (ll i = L0; i <= R0; i++) {
            string s = to_string(i);
            bool ok = true;
            for (int j = 0; j + 1 < (int)s.size(); j++)
                if (abs(s[j] - s[j + 1]) < 2) ok = false;
            if (ok) b++;
        }
        assert(a == b);
        ll c = countAuspicious(R0) - countAuspicious(L0 - 1);
        ll e = 0;
        for (ll i = L0; i <= R0; i++) {
            string s = to_string(i);
            if (s.find('4') == string::npos && s.find("62") == string::npos) e++;
        }
        assert(c == e);
    }
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    string raw, line;
    bool any = false;
    while (getline(cin, line)) {
        raw += line;
        raw += "\n";
        for (char ch : line)
            if (!isspace((unsigned char)ch)) any = true;
    }
    if (any) {
        runIO(raw);
    } else {
        runTests();
        cout << "all tests passed" << endl;
    }
    return 0;
}
