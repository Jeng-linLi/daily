// NTT 快速數論變換：模 998244353 的多項式乘法
//
// 題意：給定兩個係數在模 998244353 意義下的多項式 f、g，
//       求卷積 h(x) = f(x) * g(x) mod 998244353。
//
// 思路（NTT）：
//   998244353 = 119 * 2^23 + 1 是質數，原根 g = 3 的階恰為 p-1，
//   含足夠長的 2 的冪次單位根，可做快速數論變換（精確無浮點誤差）。
//   1. 補零到 len >= n + m - 1 的 2 的冪；
//   2. 位反轉置換 + 蝶形迭代做正變換；
//   3. 點值逐點相乘；
//   4. 用 n^-1 (mod p) 縮放做逆變換。
//
// IO 模式（stdin 有內容時）：
//   n m            （兩多項式的項數，即次數 + 1）
//   n 個整數：f 的係數（低次到高次）
//   m 個整數：g 的係數
//   輸出：n + m - 1 個整數，h 的係數（低次到高次）
//   無 stdin 時跑內建測試（與 O(n*m) 樸素卷積對拍、正逆變換往返校驗）
//   並輸出 all tests passed。
//
// 註：不用 bits/stdc++.h——本機 MinGW 8.1 的 <filesystem> 頭有已知 bug。

#include <cstdio>
#include <cassert>
#include <vector>
#include <string>
#include <sstream>
#include <iostream>
using namespace std;

typedef long long ll;
const ll MOD = 998244353;
const ll G = 3;  // 998244353 的原根

ll power(ll a, ll e, ll mod) {
    ll r = 1;
    a %= mod;
    while (e > 0) {
        if (e & 1) r = r * a % mod;
        a = a * a % mod;
        e >>= 1;
    }
    return r;
}

// 就地迭代 NTT；invert=true 時做逆變換（含除以 n）
void ntt(vector<ll>& a, bool invert = false) {
    int n = (int)a.size();
    // 位反轉置換
    int j = 0;
    for (int i = 1; i < n; ++i) {
        int bit = n >> 1;
        for (; j & bit; bit >>= 1) j ^= bit;
        j |= bit;
        if (i < j) swap(a[i], a[j]);
    }
    // 蝶形：塊長 1 -> n
    for (int length = 2; length <= n; length <<= 1) {
        // 塊內主單位根 w_len = g^((p-1)/length) 或其逆
        ll wLen = power(G, (MOD - 1) / length, MOD);
        if (invert) wLen = power(wLen, MOD - 2, MOD);
        int half = length / 2;
        for (int start = 0; start < n; start += length) {
            ll w = 1;
            for (int k = start; k < start + half; ++k) {
                ll u = a[k];
                ll v = a[k + half] * w % MOD;
                a[k] = (u + v) % MOD;
                a[k + half] = (u - v % MOD + MOD) % MOD;
                w = w * wLen % MOD;
            }
        }
    }
    if (invert) {
        ll invN = power(n, MOD - 2, MOD);
        for (int i = 0; i < n; ++i) a[i] = a[i] * invN % MOD;
    }
}

vector<ll> polyMul(const vector<ll>& f, const vector<ll>& g) {
    int need = (int)f.size() + (int)g.size() - 1;
    int size = 1;
    while (size < need) size <<= 1;
    vector<ll> a(size, 0), b(size, 0);
    for (size_t i = 0; i < f.size(); ++i) a[i] = f[i];
    for (size_t i = 0; i < g.size(); ++i) b[i] = g[i];
    ntt(a);
    ntt(b);
    for (int i = 0; i < size; ++i) a[i] = a[i] * b[i] % MOD;
    ntt(a, true);
    a.resize(need);
    return a;
}

// ---------- 測試 ----------

vector<ll> bruteMul(const vector<ll>& f, const vector<ll>& g) {
    vector<ll> res(f.size() + g.size() - 1, 0);
    for (size_t i = 0; i < f.size(); ++i)
        for (size_t j = 0; j < g.size(); ++j)
            res[i + j] = (res[i + j] + f[i] * g[j]) % MOD;
    return res;
}

static unsigned long long g_seed = 20261107ULL;
static int rnd(int mv) {
    g_seed = g_seed * 1103515245ULL + 12345ULL;
    return (int)((g_seed >> 1) % (unsigned long long)mv);
}

void runAllTests() {
    vector<vector<ll>> fixedF = {
        {1, 1},
        {1, 2, 3, 4},
        {0},
        {5},
        {1, 1, 1, 1, 1, 1, 1, 1},  // 長度恰為 2 的冪
        {2, 0, 3},
    };
    vector<vector<ll>> fixedG = {
        {1, 1},
        {5, 6, 7},
        {0},
        {7},
        {1, 1},
        {0, 4, 0, 5},
    };
    for (size_t t = 0; t < fixedF.size(); ++t)
        assert(polyMul(fixedF[t], fixedG[t]) == bruteMul(fixedF[t], fixedG[t]));

    // 隨機對拍 + 正逆變換往返校驗
    for (int t = 0; t < 200; ++t) {
        int n = rnd(50) + 1;
        int m = rnd(50) + 1;
        vector<ll> f(n), g(m);
        for (int i = 0; i < n; ++i) f[i] = rnd(MOD);
        for (int i = 0; i < m; ++i) g[i] = rnd(MOD);
        assert(polyMul(f, g) == bruteMul(f, g));
        // 往返：NTT -> INTT 應還原原係數
        int size = 1;
        while (size < n) size <<= 1;
        vector<ll> a(size, 0);
        for (int i = 0; i < n; ++i) a[i] = f[i];
        ntt(a);
        ntt(a, true);
        for (int i = 0; i < n; ++i) assert(a[i] == f[i]);
    }

    // 大規模隨機冒煙（長度 2^13 變換）：抽查 30 個係數與樸素卷積一致
    for (int t = 0; t < 5; ++t) {
        int n = 3000;
        vector<ll> f(n), g(n);
        for (int i = 0; i < n; ++i) f[i] = rnd(MOD);
        for (int i = 0; i < n; ++i) g[i] = rnd(MOD);
        vector<ll> h = polyMul(f, g);
        for (int c = 0; c < 30; ++c) {
            int k = rnd(2 * n - 1);
            ll s = 0;
            int lo = k - n + 1 > 0 ? k - n + 1 : 0;
            for (int i = lo; i <= k && i < n; ++i) s += f[i] * g[k - i] % MOD;
            assert(h[k] == s % MOD);
        }
    }

    printf("all tests passed\n");
}

int main() {
    // 一次性讀入全部 stdin；無內容時跑內建測試
    std::stringstream buf;
    buf << cin.rdbuf();
    string data = buf.str();
    if (data.find_first_not_of(" \t\r\n") == string::npos) {
        runAllTests();
        return 0;
    }
    std::istringstream in(data);
    int n, m;
    in >> n >> m;
    vector<ll> f(n), g(m);
    for (int i = 0; i < n; ++i) { ll x; in >> x; f[i] = ((x % MOD) + MOD) % MOD; }
    for (int i = 0; i < m; ++i) { ll x; in >> x; g[i] = ((x % MOD) + MOD) % MOD; }
    vector<ll> h = polyMul(f, g);
    string out;
    char tmp[24];
    for (size_t i = 0; i < h.size(); ++i) {
        if (i) out += ' ';
        snprintf(tmp, sizeof(tmp), "%lld", h[i]);
        out += tmp;
    }
    out += '\n';
    fputs(out.c_str(), stdout);
    return 0;
}
