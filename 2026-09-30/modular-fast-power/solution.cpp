// 快速冪與模運算：二進制冪 / 矩陣快速冪（Fibonacci）/ 乘法逆元
// 編譯：g++ -std=c++17 -O2 -Wall solution.cpp -o solution && ./solution
//
// 思路：
//   **二進制冪** 把指數看成二進位，底數反覆平方（a → a^2 → a^4 → …），
//   當前位爲 1 時把底數乘進答案，只需 O(log b) 次乘法。
//   **矩陣快速冪** 把 F(n) 的線性遞推寫成 [F(n),F(n-1)]^T = M·[F(n-1),F(n-2)]^T，
//   於是用 M^n 的 [0][1] 位置拿到 F(n)，O(log n) 解決線性遞推。
//   **乘法逆元**用擴展歐幾裡得解 x·a + y·m = gcd(a,m)，gcd = 1 時 x 即逆元；
//   模數爲質數時費馬小定理 a^(p-2) 同樣可用，測試裏兩者互相對拍。
//
//   所有乘法都走 mulMod 用 __int128 中轉，避免模數接近 10^18 時 long long 溢位。
//   行爲與 Python 版逐字節一致。
//
// 輸入（空白分隔）：n q k / q 行「a b m」/ p / x1..xk
// 輸出：F(n) mod 1e9+7 / q 行 a^b mod m / 一行 k 個逆元（無逆元輸出 -1）

#include <array>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <random>
#include <sstream>
#include <string>
#include <vector>

using namespace std;

using i64 = long long;
using i128 = __int128;

static const i64 MOD = 1000000007LL;

// ---------- 基本運算 ----------

// 安全的模乘：先把乘積放進 128 位再取模，模數接近 10^18 也不會溢位
static i64 mulMod(i64 a, i64 b, i64 mod) {
    return (i64)((i128)a * (i128)b % mod);
}

// 二進制冪：a^b mod mod，O(log b)
static i64 modPow(i64 a, i64 b, i64 mod) {
    if (mod == 1) return 0;
    a %= mod;
    if (a < 0) a += mod;
    i64 res = 1 % mod;
    while (b > 0) {
        if (b & 1) res = mulMod(res, a, mod);
        a = mulMod(a, a, mod);
        b >>= 1;
    }
    return res;
}

// 樸素冪：O(b)，只在測試裏當基準
static i64 slowPow(i64 a, i64 b, i64 mod) {
    a %= mod;
    if (a < 0) a += mod;
    i64 res = 1 % mod;
    for (i64 i = 0; i < b; ++i) res = mulMod(res, a, mod);
    return res;
}

using Mat = array<array<i64, 2>, 2>;

static Mat matMul(const Mat& A, const Mat& B, i64 mod) {
    Mat C{};
    for (int i = 0; i < 2; ++i)
        for (int j = 0; j < 2; ++j) {
            i128 s = 0;
            for (int k = 0; k < 2; ++k) s += (i128)A[i][k] * B[k][j];
            C[i][j] = (i64)(s % mod);
        }
    return C;
}

static Mat matPow(Mat base, long long e, i64 mod) {
    Mat res{};
    res[0][0] = 1 % mod; res[0][1] = 0;
    res[1][0] = 0;       res[1][1] = 1 % mod;
    for (int i = 0; i < 2; ++i)
        for (int j = 0; j < 2; ++j) base[i][j] %= mod;
    while (e > 0) {
        if (e & 1) res = matMul(res, base, mod);
        base = matMul(base, base, mod);
        e >>= 1;
    }
    return res;
}

// F(0)=0, F(1)=1，矩陣快速冪 O(log n)
static i64 fib(long long n, i64 mod = MOD) {
    if (n <= 0) return 0;
    Mat M{};
    M[0][0] = 1; M[0][1] = 1;
    M[1][0] = 1; M[1][1] = 0;
    return matPow(M, n, mod)[0][1];
}

// 擴展歐幾裡得：(g, x, y) 滿足 x*a + y*b = g = gcd(a, b)
static i64 extGcd(i64 a, i64 b, i64& x, i64& y) {
    i64 oldR = a, r = b, oldS = 1, s = 0, oldT = 0, t = 1;
    while (r != 0) {
        i64 q = oldR / r;
        i64 nr = oldR - q * r; oldR = r; r = nr;
        i64 ns = oldS - q * s; oldS = s; s = ns;
        i64 nt = oldT - q * t; oldT = t; t = nt;
    }
    x = oldS; y = oldT;
    return oldR;
}

static i64 modInverse(i64 a, i64 mod) {
    if (mod <= 1) return -1;
    a %= mod;
    if (a < 0) a += mod;
    if (a == 0) return -1;
    i64 x, y;
    i64 g = extGcd(a, mod, x, y);
    if (g != 1) return -1;
    x %= mod;
    if (x < 0) x += mod;
    return x;
}

static i64 modInverseFermat(i64 a, i64 p) {
    if (p <= 1 || a % p == 0) return -1;
    return modPow(a, p - 2, p);
}

static bool isPrime(i64 x) {
    if (x < 2) return false;
    for (i64 d = 2; d * d <= x; ++d)
        if (x % d == 0) return false;
    return true;
}

// ---------- 輸入輸出 ----------

static void runIo(const string& data) {
    istringstream iss(data);
    vector<i64> tk;
    i64 v;
    while (iss >> v) tk.push_back(v);
    size_t pos = 0;
    auto take = [&](i64 def) -> i64 {
        if (pos < tk.size()) return tk[pos++];
        return def;
    };

    long long n = take(0);
    long long q = take(0);
    long long k = take(0);

    ostringstream out;
    out << fib(n, MOD) << '\n';
    for (long long i = 0; i < q; ++i) {
        i64 a = take(0), b = take(0), m = take(1);
        out << modPow(a, b, m) << '\n';
    }
    i64 p = take(1);
    for (long long i = 0; i < k; ++i) {
        if (i) out << ' ';
        out << modInverse(take(0), p);
    }
    out << '\n';
    cout << out.str();
}

// ---------- 內置測試 ----------

static void runTests() {
    // ---- 快速冪 ----
    assert(modPow(2, 10, 1000) == 24);
    assert(modPow(2, 0, 7) == 1);
    assert(modPow(0, 0, 7) == 1);
    assert(modPow(0, 5, 7) == 0);
    assert(modPow(12345, 1, 97) == 12345 % 97);
    assert(modPow(5, 3, 1) == 0);
    assert(modPow(-2, 3, 7) == 6);          // (-2)^3 = -8 ≡ 6 (mod 7)
    assert(modPow(3, 100, MOD) == 886041711);   // 與 Python 版 pow(3, 100, 1e9+7) 一致

    mt19937_64 rng(20260930);
    for (int t = 0; t < 3000; ++t) {
        i64 a = (i64)(rng() % 1001) - 500;
        i64 b = (i64)(rng() % 201);
        i64 m = (i64)(rng() % 500) + 1;
        if (b <= 60) assert(modPow(a, b, m) == slowPow(a, b, m));
        // 與 Python 版一致的性質：a^(b+c) = a^b * a^c
        i64 c = (i64)(rng() % 50);
        assert(modPow(a, b + c, m) == mulMod(modPow(a, b, m), modPow(a, c, m), m));
    }
    // 超大指數
    assert(modPow(2, 1000000000LL, MOD) == 140625001);  // 2^(10^9) mod 1e9+7

    // ---- Fibonacci ----
    // 迭代版（取模）當基準，對拍矩陣快速冪
    for (i64 mod : {(i64)MOD, (i64)1000, (i64)1000003}) {
        vector<i64> seq(91);
        seq[0] = 0; seq[1] = 1 % mod;
        for (int i = 2; i <= 90; ++i) seq[i] = (seq[i - 1] + seq[i - 2]) % mod;
        for (int i = 0; i <= 90; ++i) assert(fib(i, mod) == seq[i]);
    }
    assert(fib(0) == 0 && fib(1) == 1 && fib(2) == 1 && fib(10) == 55);
    assert(fib(50, MOD) == 586268941);
    assert(fib(90, MOD) == 2880067194370816120LL % MOD);
    assert(fib(200, 1000) == (fib(199, 1000) + fib(198, 1000)) % 1000);
    assert(fib(40, 1000) == 102334155 % 1000);

    // ---- 逆元 ----
    assert(modInverse(3, 7) == 5);
    assert(modInverse(1, 100) == 1);
    assert(modInverse(2, 4) == -1);
    assert(modInverse(0, 7) == -1);
    assert(modInverse(5, 1) == -1);
    assert(modInverseFermat(3, 7) == 5);

    for (int t = 0; t < 2000; ++t) {
        i64 p = (i64)(rng() % 90) + 2;
        while (!isPrime(p)) ++p;
        i64 x = (i64)(rng() % p);
        i64 inv = modInverse(x, p);
        if (x % p == 0) {
            assert(inv == -1);
            assert(modInverseFermat(x, p) == -1);
        } else {
            assert(inv != -1);
            assert(mulMod(x, inv, p) == 1);
            assert(modInverseFermat(x, p) == inv);
        }
    }
    for (int t = 0; t < 2000; ++t) {
        i64 m = (i64)(rng() % 59) + 2;
        i64 x = (i64)(rng() % m);
        i64 inv = modInverse(x, m);
        i64 g, xx, yy;
        g = extGcd(x, m, xx, yy);
        if (g != 1) assert(inv == -1);
        else assert(mulMod(x, inv, m) == 1);
    }

    cout << "all tests passed" << '\n';
}

int main() {
    string data, line;
    bool hasInput = false;
    while (getline(cin, line)) {
        data += line;
        data += '\n';
        if (!line.empty()) hasInput = true;
    }
    if (hasInput) runIo(data);
    else runTests();
    return 0;
}
