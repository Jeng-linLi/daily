// FWT 快速沃爾什變換：XOR / OR / AND 三種位運算卷積
//
// 題意：
//     給定整數 k（0 <= k <= 16）與兩個長度為 2^k 的序列 A、B，
//     以及一個運算 op ∈ {XOR, OR, AND}，求它們的**位運算卷積**：
//
//         XOR 卷積   C[t] = Σ_{i^j = t} A[i]·B[j]
//         OR  卷積   C[t] = Σ_{i|j = t} A[i]·B[j]
//         AND 卷積   C[t] = Σ_{i&j = t} A[i]·B[j]
//
//     所有計算對 MOD = 998244353 取模；輸入的負數先化歸到 [0, MOD)。
//     輸出長度 2^k 的序列 C。
//
// 思路：
//     ### 1. 位運算卷積 = 布爾格上的多項式乘法
//     把 A、B 看成「以布爾格 {0,1}^k 為下標、以逐位元運算為乘法」的
//     多項式。三種運算各對應一個交換半群，卷積就是這個半群上的
//     多項式乘法 —— 直接算是 O(4^k)。
//
//     ### 2. FWT：把卷積變成逐點乘
//     每種位運算都有一個對應的線性變換（zeta 變換），把序列送到
//     「下標和」的域，卷積就退化為逐點相乘，再逆變換回來：
//
//         XOR：蝶形 (a, b) → (a+b, a−b)            逆變換同形，最後乘 2^{-k}
//         OR ：zeta  (a, b) → (a, a+b)             逆變換 (a, b) → (a, b−a)
//         AND：zeta  (a, b) → (a+b, b)             逆變換 (a, b) → (a−b, b)
//
//     三種變換都是「按位長度倍增」的原地蝶形，總複雜度 O(k·2^k)。
//     推導方式完全相同：找出 (a±b) 這樣的線性組合，使得
//     transform(A)[s] = Σ_{i ⊆ s 或 i^j=s} A[i]（依運算而定）。
//
//     ### 3. 為什麼 OR 與 AND 的zeta方向看起來是「反的」
//     OR 的 zeta 是「把低位併進高位」(a, b) → (a, a+b)：s = i|j 涵蓋
//     所有 i ⊆ s。AND 則是「把高位併進低位」(a, b) → (a+b, b)：
//     s = i&j 涵蓋所有 s ⊆ i。兩者互為對偶（把 0/1 互換、下標取反即可
//     互相轉化），寫反是最常見的 bug，所以測試裡對兩種方向都做
//     暴力對拍。
//
//     ### 4. 模運算
//     MOD = 998244353 是質數，2^{-k} 用費馬小定理求逆元。
//     負數輸入先化歸到 [0, MOD)，C++ 與 Python 用同一條規則
//     ((v % MOD) + MOD) % MOD，保證輸出逐位元組一致。
//
// 複雜度（n = 2^k）：
//     三次變換 + 逐點乘      時間 O(k·n)
//     空間                  O(n)
//
// 應用場景：
//     - 多項式乘法：OR / AND 卷積是「按位受限」的多項式乘法，
//       集合冪級數（set power series）的基礎運算
//     - 子集計數：統計「兩組條件的 OR/AND 恰好等於某掩碼」的方案數，
//       如特徵組合、權限位彙總
//     - 期望與機率：擲骰子/抽卡模型中用 XOR 卷積快速求分佈的卷積
//       （k 個 bit 的狀態轉移）
//     - 圖論：傳遞閉包的布爾版（OR 卷積 + 重複平方法）可加速
//       某些位元矩陣運算
//     - 競賽常青題：BZOJ 4263 / CF 1119F 一類的 FWT 入門題

#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <optional>
#include <random>
#include <sstream>
#include <string>
#include <vector>

#ifdef _WIN32
#include <fcntl.h>
#include <io.h>
#endif

namespace {

constexpr int KMAX = 16;
constexpr long long MOD = 998244353LL;

void fail(const std::string& msg) {
    std::cerr << "test failed: " << msg << std::endl;
    std::exit(1);
}

long long norm(long long v) {
    return ((v % MOD) + MOD) % MOD;
}

// a, b ∈ [0, MOD) 時的減法，結果保持在 [0, MOD)
long long submod(long long a, long long b) {
    long long r = (a - b) % MOD;
    return r < 0 ? r + MOD : r;
}

long long powerMod(long long a, long long e) {
    long long r = 1;
    a %= MOD;
    while (e > 0) {
        if (e & 1) {
            r = r * a % MOD;
        }
        a = a * a % MOD;
        e >>= 1;
    }
    return r;
}

// ---------------------------------------------------------------- FWT

void fwtXor(std::vector<long long>& a, bool inverse) {
    int n = static_cast<int>(a.size());
    for (int step = 1; step < n; step *= 2) {
        for (int i = 0; i < n; i += step * 2) {
            for (int j = i; j < i + step; ++j) {
                long long u = a[j];
                long long v = a[j + step];
                a[j] = (u + v) % MOD;
                a[j + step] = submod(u, v);
            }
        }
    }
    if (inverse) {
        long long invN = powerMod(n, MOD - 2);
        for (int i = 0; i < n; ++i) {
            a[i] = a[i] * invN % MOD;
        }
    }
}

void fwtOr(std::vector<long long>& a, bool inverse) {
    int n = static_cast<int>(a.size());
    for (int step = 1; step < n; step *= 2) {
        for (int i = 0; i < n; i += step * 2) {
            for (int j = i; j < i + step; ++j) {
                long long u = a[j];
                long long v = a[j + step];
                if (inverse) {
                    a[j] = u % MOD;
                    a[j + step] = submod(v, u);
                } else {
                    a[j] = u % MOD;
                    a[j + step] = (u + v) % MOD;
                }
            }
        }
    }
}

void fwtAnd(std::vector<long long>& a, bool inverse) {
    int n = static_cast<int>(a.size());
    for (int step = 1; step < n; step *= 2) {
        for (int i = 0; i < n; i += step * 2) {
            for (int j = i; j < i + step; ++j) {
                long long u = a[j];
                long long v = a[j + step];
                if (inverse) {
                    a[j] = submod(u, v);
                    a[j + step] = v % MOD;
                } else {
                    a[j] = (u + v) % MOD;
                    a[j + step] = v % MOD;
                }
            }
        }
    }
}

void fwt(std::vector<long long>& a, const std::string& op, bool inverse) {
    if (op == "XOR") {
        fwtXor(a, inverse);
    } else if (op == "OR") {
        fwtOr(a, inverse);
    } else {
        fwtAnd(a, inverse);
    }
}

std::vector<long long> convolve(std::vector<long long> a,
                                std::vector<long long> b,
                                const std::string& op) {
    fwt(a, op, false);
    fwt(b, op, false);
    for (size_t i = 0; i < a.size(); ++i) {
        a[i] = a[i] * b[i] % MOD;
    }
    fwt(a, op, true);
    return a;
}

std::vector<long long> bruteConvolve(const std::vector<long long>& a,
                                     const std::vector<long long>& b,
                                     const std::string& op) {
    int n = static_cast<int>(a.size());
    std::vector<long long> c(n, 0);
    for (int i = 0; i < n; ++i) {
        if (a[i] == 0) {
            continue;
        }
        for (int j = 0; j < n; ++j) {
            if (b[j] == 0) {
                continue;
            }
            int t;
            if (op == "XOR") {
                t = i ^ j;
            } else if (op == "OR") {
                t = i | j;
            } else {
                t = i & j;
            }
            c[t] = (c[t] + a[i] * b[j]) % MOD;
        }
    }
    return c;
}

std::vector<long long> solve(int k, const std::vector<long long>& a,
                             const std::vector<long long>& b,
                             const std::string& op) {
    std::vector<long long> fa(a.size()), fb(b.size());
    for (size_t i = 0; i < a.size(); ++i) {
        fa[i] = norm(a[i]);
        fb[i] = norm(b[i]);
    }
    return convolve(fa, fb, op);
}

// ---------------------------------------------------------------- 測試

void expectVec(const std::vector<long long>& got,
               const std::vector<long long>& want, const std::string& tag) {
    if (got != want) {
        std::string msg = tag + "\n got:";
        for (long long v : got) msg += " " + std::to_string(v);
        msg += "\nwant:";
        for (long long v : want) msg += " " + std::to_string(v);
        fail(msg);
    }
}

void runTests() {
    // ---- 固定用例（手工驗算：A=[1,2], B=[3,4]）----
    expectVec(solve(1, {1, 2}, {3, 4}, "XOR"), {11, 10}, "fixed xor");
    expectVec(solve(1, {1, 2}, {3, 4}, "OR"), {3, 18}, "fixed or");
    expectVec(solve(1, {1, 2}, {3, 4}, "AND"), {13, 8}, "fixed and");
    // k = 0：長度 1，三種運算都歸結為 A[0]·B[0]
    expectVec(solve(0, {7}, {9}, "XOR"), {63}, "fixed k0 xor");
    expectVec(solve(0, {7}, {9}, "OR"), {63}, "fixed k0 or");
    expectVec(solve(0, {7}, {9}, "AND"), {63}, "fixed k0 and");
    // 取模邊界：MOD-1 ≡ -1，(-1)·(-1) = 1
    long long m1 = MOD - 1;
    expectVec(solve(1, {m1, m1}, {m1, m1}, "XOR"), {2, 2}, "fixed mod xor");
    expectVec(solve(1, {m1, m1}, {m1, m1}, "OR"), {1, 3}, "fixed mod or");
    expectVec(solve(1, {m1, m1}, {m1, m1}, "AND"), {3, 1}, "fixed mod and");
    // 負數輸入：先化歸到 [0, MOD)
    expectVec(solve(1, {-1, 0}, {-1, 0}, "XOR"), {1, 0}, "fixed neg xor");
    {
        std::vector<long long> a = {-1, -2, -3, -4};
        std::vector<long long> b = {-5, -6, -7, -8};
        std::vector<long long> na(4), nb(4);
        for (int i = 0; i < 4; ++i) {
            na[i] = norm(a[i]);
            nb[i] = norm(b[i]);
        }
        for (const std::string& op : {"XOR", "OR", "AND"}) {
            expectVec(solve(2, a, b, op), bruteConvolve(na, nb, op),
                      "fixed neg " + op);
        }
    }
    // 全零 / 全相同
    expectVec(solve(2, {0, 0, 0, 0}, {1, 2, 3, 4}, "XOR"), {0, 0, 0, 0},
              "fixed zero");
    {
        std::vector<long long> a(4, 5), b(4, 5);
        expectVec(solve(2, a, b, "AND"), bruteConvolve(a, b, "AND"),
                  "fixed same");
    }

    // ---- 變換往返：inverse(transform(x)) == x ----
    std::mt19937 rng(20261026);
    for (int tc = 0; tc < 30; ++tc) {
        int k = static_cast<int>(rng() % 9);
        int n = 1 << k;
        std::vector<long long> x(n);
        for (int i = 0; i < n; ++i) {
            x[i] = static_cast<long long>(rng()) % MOD;
        }
        for (const std::string& op : {"XOR", "OR", "AND"}) {
            std::vector<long long> y = x;
            fwt(y, op, false);
            fwt(y, op, true);
            if (y != x) {
                fail("roundtrip " + op);
            }
        }
    }

    // ---- 隨機對拍：FWT vs 暴力雙重循環 ----
    for (int tc = 0; tc < 120; ++tc) {
        int k = static_cast<int>(rng() % 8);
        int n = 1 << k;
        std::vector<long long> a(n), b(n);
        for (int i = 0; i < n; ++i) {
            a[i] = static_cast<long long>(rng()) % MOD;
            b[i] = static_cast<long long>(rng()) % MOD;
        }
        std::vector<long long> na(n), nb(n);
        for (int i = 0; i < n; ++i) {
            na[i] = norm(a[i]);
            nb[i] = norm(b[i]);
        }
        for (const std::string& op : {"XOR", "OR", "AND"}) {
            expectVec(solve(k, a, b, op), bruteConvolve(na, nb, op),
                      "random " + op);
        }
    }

    // ---- 稀疏性壓力：大量 0 與少量大值 ----
    for (int tc = 0; tc < 20; ++tc) {
        int k = 8;
        int n = 1 << k;
        std::vector<long long> a(n, 0), b(n, 0);
        for (int t = 0; t < 5; ++t) {
            a[rng() % n] = static_cast<long long>(rng()) % MOD;
            b[rng() % n] = static_cast<long long>(rng()) % MOD;
        }
        for (const std::string& op : {"XOR", "OR", "AND"}) {
            expectVec(solve(k, a, b, op), bruteConvolve(a, b, op),
                      "sparse " + op);
        }
    }
}

// ---------------------------------------------------------------- IO 模式

std::optional<long long> parseInt(const std::string& s) {
    // 與 Python 版一致：可選正負號 + 純數字，空串/非數字/越界皆 nullopt
    if (s.empty()) {
        return std::nullopt;
    }
    size_t i = 0;
    if (s[0] == '+' || s[0] == '-') {
        i = 1;
    }
    if (i == s.size()) {
        return std::nullopt;
    }
    for (size_t j = i; j < s.size(); ++j) {
        if (s[j] < '0' || s[j] > '9') {
            return std::nullopt;
        }
    }
    long long val = 0;
    for (size_t j = i; j < s.size(); ++j) {
        val = val * 10 + (s[j] - '0');
        if (val > 1000000000000000000LL) {
            return std::nullopt;
        }
    }
    if (s[0] == '-') {
        val = -val;
    }
    return val;
}

void runIo(const std::string& raw) {
    std::istringstream in(raw);
    std::vector<std::string> toks;
    std::string t;
    while (in >> t) {
        toks.push_back(t);
    }
    size_t pos = 0;
    auto nxt = [&]() -> const std::string* {
        if (pos >= toks.size()) {
            return nullptr;
        }
        return &toks[pos++];
    };

    auto bad = []() {
        std::cout << "error\n";
    };

    const std::string* p = nxt();
    std::optional<long long> k = p ? parseInt(*p) : std::nullopt;
    if (!k || !(*k >= 0 && *k <= KMAX)) {
        bad();
        return;
    }
    int n = 1 << static_cast<int>(*k);
    std::vector<long long> a, b;
    for (int idx = 0; idx < 2; ++idx) {
        std::vector<long long>& seq = (idx == 0) ? a : b;
        for (int i = 0; i < n; ++i) {
            p = nxt();
            std::optional<long long> v = p ? parseInt(*p) : std::nullopt;
            if (!v) {
                bad();
                return;
            }
            seq.push_back(*v);
        }
    }
    p = nxt();
    if (!p || (*p != "XOR" && *p != "OR" && *p != "AND")) {
        bad();
        return;
    }
    std::string op = *p;
    std::vector<long long> res = solve(static_cast<int>(*k), a, b, op);
    for (size_t i = 0; i < res.size(); ++i) {
        if (i) {
            std::cout << ' ';
        }
        std::cout << res[i];
    }
    std::cout << "\n";
}

}  // namespace

int main() {
#ifdef _WIN32
    _setmode(_fileno(stdout), _O_BINARY);
#endif
    std::string raw((std::istreambuf_iterator<char>(std::cin)),
                    std::istreambuf_iterator<char>());
    bool blank = true;
    for (size_t i = 0; i < raw.size(); ++i) {
        unsigned char u = static_cast<unsigned char>(raw[i]);
        if (u > ' ') {
            blank = false;
            break;
        }
    }
    if (!blank) {
        runIo(raw);
    } else {
        runTests();
        std::cout << "all tests passed" << std::endl;
    }
    return 0;
}
