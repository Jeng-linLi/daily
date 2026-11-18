// Berlekamp-Massey 最短線性遞迴 + Kitamasa 快速求第 n 項（mod 質數 p）。
//
// IO 模式（stdin）：每行 p k n，接著 k 個 s0..s(k-1)，可讀多筆直到 EOF。
//   輸出 2 行：L c1..cL / 第 n 項 (0-based) mod p。
//   p 須為質數。無 stdin 時跑內建對拍，輸出 all tests passed。
// 複雜度：BM O(k²)，Kitamasa O(L² log n)。

#include <cassert>
#include <cstdio>
#include <string>
#include <vector>
#include <sstream>
#include <iostream>
#include <cstdint>

#ifdef _WIN32
#include <fcntl.h>
#include <io.h>
#endif

using namespace std;

typedef long long ll;
typedef unsigned long long ull;

static ll modmul(ll a, ll b, ll p) {
    return (ll)((__int128)a * b % p);
}

static ll modpow(ll a, ll e, ll p) {
    ll r = 1 % p;
    while (e > 0) {
        if (e & 1) r = modmul(r, a, p);
        a = modmul(a, a, p);
        e >>= 1;
    }
    return r;
}

// ----------------------------- Berlekamp-Massey -----------------------------

// 回傳最短遞迴係數 c[1..L]，使 s[n] = Σ c[i]·s[n-i] (mod p)。
static vector<ll> berlekamp_massey(const vector<ll> &s, ll p) {
    int n = (int)s.size();
    vector<ll> C(1, 1), B(1, 1);
    int L = 0, m = 1;
    ll b = 1;
    for (int N = 0; N < n; N++) {
        ll d = s[N] % p;
        for (size_t i = 1; i < C.size(); i++)
            d = (d + modmul(C[i], s[N - (int)i], p)) % p;
        if (d == 0) { m++; continue; }
        vector<ll> T = C;
        ll coef = modmul(d, modpow(b, p - 2, p), p);
        if ((int)C.size() < (int)B.size() + m)
            C.resize(B.size() + m, 0);
        for (size_t i = 0; i < B.size(); i++)
            C[i + m] = (C[i + m] - modmul(coef, B[i], p)) % p;
        if (2 * L <= N) {
            L = N + 1 - L;
            B = T;
            b = d;
            m = 1;
        } else m++;
    }
    vector<ll> c;
    for (int i = 1; i <= L; i++) c.push_back((p - C[i]) % p);
    return c;
}

// ------------------------------- Kitamasa -------------------------------

// (a·b) mod (x^L - c1 x^(L-1) - ... - cL)。c 為 1-index 係數（c[0]=c1）。
static vector<ll> polyMulMod(const vector<ll> &a, const vector<ll> &b,
                             const vector<ll> &c, ll p) {
    int L = (int)c.size();
    vector<ll> tmp(2 * L - 1, 0);
    for (size_t i = 0; i < a.size(); i++) {
        if (a[i] == 0) continue;
        for (size_t j = 0; j < b.size(); j++)
            tmp[i + j] = (tmp[i + j] + modmul(a[i], b[j], p)) % p;
    }
    for (int d = (int)tmp.size() - 1; d >= L; d--) {
        ll val = tmp[d];
        if (val == 0) continue;
        for (int i = 1; i <= L; i++)
            tmp[d - i] = (tmp[d - i] + modmul(val, c[i - 1], p)) % p;
    }
    return vector<ll>(tmp.begin(), tmp.begin() + L);
}

// 求序列第 n 項（0-based）。s 為前 L 項，c 為遞迴係數。
static ll kitamasa(const vector<ll> &s, const vector<ll> &c, ull n, ll p) {
    int L = (int)c.size();
    if (L == 0) return 0;
    if (n < (ull)L) return s[n] % p;
    vector<ll> result(L, 0), base(L, 0);
    result[0] = 1;
    if (L == 1) base[0] = c[0] % p;
    else base[1] = 1;
    ull e = n;
    while (e > 0) {
        if (e & 1) result = polyMulMod(result, base, c, p);
        base = polyMulMod(base, base, c, p);
        e >>= 1;
    }
    ll ans = 0;
    for (int i = 0; i < L; i++) ans = (ans + modmul(result[i], s[i], p)) % p;
    return ans;
}

// ----------------------------- 暴力參考 -----------------------------

// 暴力找最短遞迴（僅小 p、小 L）。回傳最短長度。
static int bruteFirstFit(const vector<ll> &s, ll p) {
    int k = (int)s.size();
    for (int L = 0; L <= k / 2; L++) {
        if (L == 0) {
            bool allz = true;
            for (ll x : s) if (x % p != 0) { allz = false; break; }
            if (allz) return 0;
            continue;
        }
        // 枚舉所有係數向量（p^L）
        vector<int> idx(L, 0);
        bool found = false;
        while (true) {
            bool ok = true;
            for (int t = L; t < k && ok; t++) {
                ll v = 0;
                for (int i = 0; i < L; i++) v = (v + modmul(idx[i], s[t - 1 - i], p)) % p;
                if (v != s[t] % p) ok = false;
            }
            if (ok) { found = true; break; }
            int pos = L - 1;
            while (pos >= 0 && idx[pos] == p - 1) { idx[pos] = 0; pos--; }
            if (pos < 0) break;
            idx[pos]++;
        }
        if (found) return L;
    }
    return k / 2;
}

static ll dpNth(const vector<ll> &s, const vector<ll> &c, ull n, ll p) {
    int L = (int)c.size();
    if (L == 0) return 0;
    if (n < s.size()) return s[n] % p;
    vector<ll> seq = s;
    while (seq.size() <= n) {
        ll nxt = 0;
        for (int i = 0; i < L; i++)
            nxt = (nxt + modmul(c[i], seq[seq.size() - 1 - i], p)) % p;
        seq.push_back(nxt);
    }
    return seq[n] % p;
}

// ------------------------------- IO 模式 -------------------------------

static string run_io(const vector<string> &toks) {
    size_t pos = 0;
    string out;
    auto nxt = [&]() -> const string* {
        if (pos >= toks.size()) return nullptr;
        return &toks[pos++];
    };
    while (true) {
        const string *a = nxt(), *b = nxt(), *d = nxt();
        if (!a || !b || !d) break;
        ll p, k;
        ull n;
        try {
            p = stoll(*a); k = stoll(*b); n = stoull(*d);
        } catch (...) { break; }
        if (p < 2 || p > 2000000000LL || k < 0 || k > 5000) break;
        vector<ll> s;
        bool ok = true;
        for (ll i = 0; i < k; i++) {
            const string *v = nxt();
            if (!v) { ok = false; break; }
            try { s.push_back((stoll(*v) % p + p) % p); }
            catch (...) { ok = false; break; }
        }
        if (!ok) break;
        vector<ll> c = berlekamp_massey(s, p);
        ll ans = kitamasa(s, c, n, p);
        ostringstream line;
        line << c.size();
        for (ll x : c) line << ' ' << x;
        line << '\n' << ans << '\n';
        out += line.str();
    }
    return out;
}

// ------------------------------ 內建測試 ------------------------------

static ull g_state = 20261118ULL;
static int rnd(int m) {
    g_state = (g_state * 1103515245ULL + 12345ULL) & 0x7FFFFFFFULL;
    return (int)(g_state % (unsigned)m);
}

static void run_tests() {
    const ll P = 1000000007LL;
    // Fibonacci：c = [1,1]
    {
        vector<ll> fib(12);
        fib[0] = fib[1] = 1;
        for (int i = 2; i < 12; i++) fib[i] = (fib[i-1] + fib[i-2]) % P;
        vector<ll> c = berlekamp_massey(fib, P);
        assert(c.size() == 2 && c[0] == 1 && c[1] == 1);
        assert(kitamasa(fib, c, 10, P) == 89);
        assert(kitamasa(fib, c, 100, P) == dpNth(fib, c, 100, P));
    }
    // s[n] = 3 s[n-1] - 2 s[n-2]
    {
        vector<ll> seq(14);
        seq[0] = 2; seq[1] = 3;
        for (int i = 2; i < 14; i++) seq[i] = ((3 * seq[i-1] - 2 * seq[i-2]) % P + P) % P;
        vector<ll> c = berlekamp_massey(seq, P);
        assert(c.size() == 2 && c[0] == 3 && c[1] == P - 2);
        assert(kitamasa(seq, c, 500, P) == dpNth(seq, c, 500, P));
    }
    // 隨機對拍
    const ll smallp = 7;
    for (int tc = 0; tc < 300; tc++) {
        int L = rnd(4) + 1;
        vector<ll> ctrue(L);
        for (int i = 0; i < L; i++) ctrue[i] = rnd((int)smallp);
        if (ctrue[L-1] == 0) ctrue[L-1] = 1;
        vector<ll> s(L);
        for (int i = 0; i < L; i++) s[i] = rnd((int)smallp);
        for (int i = 0; i < 2 * L + 8; i++) {
            ll nxt = 0;
            for (int j = 0; j < L; j++)
                nxt = (nxt + modmul(ctrue[j], s[s.size()-1-j], smallp)) % smallp;
            s.push_back(nxt);
        }
        vector<ll> got = berlekamp_massey(s, smallp);
        assert((int)got.size() <= L + 1);
        for (int t = (int)got.size(); t < (int)s.size(); t++) {
            ll v = 0;
            for (size_t i = 0; i < got.size(); i++)
                v = (v + modmul(got[i], s[t-1-(int)i], smallp)) % smallp;
            assert(v == s[t] % smallp);
        }
        assert((int)got.size() == bruteFirstFit(s, smallp));
        ull ns[] = {0, 1, (ull)L, (ull)(2*L+3)};
        for (ull n : ns)
            assert(kitamasa(s, got, n, smallp) == dpNth(s, got, n, smallp));
    }
    // 全零
    assert(berlekamp_massey({0,0,0,0}, P).empty());
    assert(berlekamp_massey({}, P).empty());
}

int main() {
#ifdef _WIN32
    _setmode(_fileno(stdout), _O_BINARY);
#endif
    ostringstream buf;
    buf << cin.rdbuf();
    string input = buf.str();
    if (!input.empty() && input.find_first_not_of(" \t\r\n") != string::npos) {
        istringstream iss(input);
        vector<string> toks;
        string x;
        while (iss >> x) toks.push_back(x);
        fputs(run_io(toks).c_str(), stdout);
    } else {
        run_tests();
        fputs("all tests passed\n", stdout);
    }
    return 0;
}
