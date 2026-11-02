// 莫比烏斯反演（Möbius Inversion）：求 gcd 恰好為 g 的數對數量
//
// 題意：
//     給定 n 個正整數 a[1..n] 與 q 個查詢 g_1..g_q。
//     對每個查詢 g，求有多少個**無序數對** (i, j)（1 <= i < j <= n）
//     滿足 gcd(a[i], a[j]) == g。
//     限制：1 <= n, q <= 100000，1 <= a[i] <= 100000，1 <= g <= 100000。
//
// 思路：
//     ### 1. 從容斥講起
//     設 cnt[d] = 滿足 d | a[i] 的元素個數。
//     則「兩個數都是 d 的倍數」的數對數是 C(cnt[d], 2) = cnt[d](cnt[d]-1)/2。
//
//     我們要的是 gcd **恰好等於** g，而不是「g 是它們的公因數」。
//     容斥的做法是：
//         答案(g) = Σ_{k>=1} μ(k) * C(cnt[g*k], 2)
//     其中 μ 是莫比烏斯函數：
//         μ(1) = 1
//         μ(n) = 0              若 n 含平方因子
//         μ(n) = (-1)^k         若 n 是 k 個相異質數的乘積
//     這個式子的直觀意義是「對『g 是公因數』做帶權容斥」，
//     權重就是 μ。
//
//     ### 2. 莫比烏斯函數的定義式
//     莫比烏斯反演定理：若 f(n) = Σ_{d|n} g(d)，則 g(n) = Σ_{d|n} μ(d) * f(n/d)。
//     反過來，若 g(n) = Σ_{d|n} f(d)，則 f(n) = Σ_{d|n} μ(d) * g(n/d)。
//     套到本題：
//         f(g) = #{數對 : gcd 恰好為 g}
//         g(d) = C(cnt[d], 2) = #{數對 : d 同時整除兩數}
//     因為「gcd 是 d 的倍數」等價於「d 同時整除兩數」，所以 g(d) = C(cnt[d],2)，
//     而 f(g) = Σ_{k>=1} μ(k) * g(g*k)，正是上面的容斥式。
//
//     ### 3. 怎麼算 cnt[d]
//     直接對每個 d 去數有幾個 a[i] 是 d 的倍數，O(n * maxA) 太慢。
//     用「倒數篩」：對每個 d from 1 to M，cnt[d] = freq[d] + freq[2d] + freq[3d] + ...
//     總共 Σ_{d=1}^{M} M/d = O(M log M)。
//
//     ### 4. 預先算出所有 g 的答案
//     對固定 g，Σ_{k>=1} μ(k) * C(cnt[g*k], 2) 的 k 範圍是 1..M/g。
//     我們對**每個 g**（1..最大查詢值）都預先算好：
//         總運算量 = Σ_{g=1}^{M} (M/g) = O(M log M)
//     有了 ans[1..M] 之後，每個查詢 O(1) 回應。
//
//     ### 5. 邊界
//     cnt[d] 可能為 0 或 1，此時 C(cnt,2) = 0（沒有數對），
//     計算時必須先判 cnt >= 2，避免負數乘積。
//
// 複雜度：
//     線性篩求 μ            O(M)
//     倍數列舉求 cnt[]     O(M log M)
//     對所有 g 求 ans[]     O(M log M)
//     每個查詢              O(1)
//     總計                  時間 O(M log M + n + q)，空間 O(M + n)
//
// 應用場景（Applications）：
//     - **數論題的批量求解**：一次反演算出所有 g 的答案，之後查詢 O(1)
//     - **資訊安全**：RSA 的 φ(n) 計算、互質數對統計（金鑰生成前的檢查）
//     - **生物統計**：族群中「親緣關係恰為某值」的配對數量
//     - **資料清洗**：找出「共同特徵恰好為某類」的記錄對數
//     - **競賽題型**：數位反演、歐拉函數求和、Divisor Sum 的各種變體

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <numeric>
#include <optional>
#include <string>
#include <vector>

#ifdef _WIN32
#include <fcntl.h>
#include <io.h>
#endif

namespace {

constexpr int NMAX = 100000;
constexpr int QMAX = 100000;
constexpr int VAL_MAX = 100000;
constexpr long long PARSE_LIMIT = 1000000000000000000LL;

// 線性篩求莫比烏斯函數 μ[1..m]
std::vector<int> mobiusSieve(int m) {
    std::vector<int> mu(static_cast<std::size_t>(m) + 1, 0);
    if (m >= 1) mu[1] = 1;
    std::vector<char> isComposite(static_cast<std::size_t>(m) + 1, 0);
    std::vector<int> primes;
    for (int i = 2; i <= m; ++i) {
        if (!isComposite[static_cast<std::size_t>(i)]) {
            primes.push_back(i);
            mu[static_cast<std::size_t>(i)] = -1;
        }
        for (int p : primes) {
            long long x = 1LL * i * p;
            if (x > m) break;
            isComposite[static_cast<std::size_t>(x)] = 1;
            if (i % p == 0) {
                mu[static_cast<std::size_t>(x)] = 0;
                break;
            }
            mu[static_cast<std::size_t>(x)] = -mu[static_cast<std::size_t>(i)];
        }
    }
    return mu;
}

// 計算 ans[g] = gcd 恰好為 g 的無序數對數量（g 從 0 到 maxQuery）
std::vector<long long> countPairsByGcd(const std::vector<int>& values, int maxQuery) {
    std::vector<long long> ans(static_cast<std::size_t>(maxQuery) + 1, 0);
    int n = static_cast<int>(values.size());
    if (n < 2 || maxQuery <= 0) return ans;

    int m = maxQuery;
    for (int v : values) m = std::max(m, v);

    // ---- 1. 線性篩求 μ ----
    std::vector<int> mu = mobiusSieve(m);
    // ---- 2. 頻率表 ----
    std::vector<int> freq(static_cast<std::size_t>(m) + 1, 0);
    for (int v : values) freq[static_cast<std::size_t>(v)] += 1;
    // ---- 3. 倍數列舉求 cnt[d] = 有多少個元素是 d 的倍數 ----
    std::vector<int> cnt(static_cast<std::size_t>(m) + 1, 0);
    for (int d = 1; d <= m; ++d) {
        int s = 0;
        for (int v = d; v <= m; v += d) s += freq[static_cast<std::size_t>(v)];
        cnt[static_cast<std::size_t>(d)] = s;
    }
    // ---- 4. C(cnt[d], 2) ----
    std::vector<long long> pairs(static_cast<std::size_t>(m) + 1, 0);
    for (int d = 1; d <= m; ++d) {
        long long c = cnt[static_cast<std::size_t>(d)];
        if (c >= 2) pairs[static_cast<std::size_t>(d)] = c * (c - 1) / 2;
    }
    // ---- 5. 反演：ans[g] = Σ_k μ(k) * pairs[g*k] ----
    for (int g = 1; g <= maxQuery; ++g) {
        long long total = 0;
        int limit = m / g;
        for (int k = 1; k <= limit; ++k) {
            total += 1LL * mu[static_cast<std::size_t>(k)] *
                     pairs[static_cast<std::size_t>(1LL * g * k)];
        }
        ans[static_cast<std::size_t>(g)] = total;
    }
    return ans;
}

long long gcdLL(long long a, long long b) {
    while (b != 0) {
        long long t = a % b;
        a = b;
        b = t;
    }
    return a;
}

// 參考實作：直接枚舉所有數對算 gcd（O(n^2 log V)，僅供對拍）
long long bruteForce(const std::vector<int>& values, int g) {
    int n = static_cast<int>(values.size());
    long long cnt = 0;
    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            if (gcdLL(values[static_cast<std::size_t>(i)],
                      values[static_cast<std::size_t>(j)]) == g) {
                ++cnt;
            }
        }
    }
    return cnt;
}

// 簡易 LCG，兩種語言共用同一組常數以確保亂數序列一致
struct Lcg {
    unsigned long long s;
    explicit Lcg(unsigned long long seed) : s(seed) {}
    unsigned long long next() {
        s = (s * 6364136223846793005ULL + 1442695040888963407ULL) &
            0xFFFFFFFFFFFFFFFFULL;
        return s;
    }
    int range(int lo, int hi) {  // [lo, hi]
        return lo + static_cast<int>(
                          next() % static_cast<unsigned long long>(hi - lo + 1));
    }
};

void fail(const char* msg) {
    std::cerr << msg << "\n";
    std::exit(1);
}

void checkEq(long long got, long long want, const char* msg) {
    if (got != want) fail(msg);
}

void runTests() {
    // ---- 固定用例 ----
    {
        std::vector<int> v{6, 10};
        std::vector<long long> a = countPairsByGcd(v, 10);
        checkEq(a[1], 0, "case 1a failed");
        checkEq(a[2], 1, "case 1b failed");
    }
    {
        std::vector<int> v{9, 28};
        std::vector<long long> a = countPairsByGcd(v, 10);
        checkEq(a[1], 1, "case 2 failed");
    }
    {
        std::vector<int> v{6, 10, 15};
        std::vector<long long> a = countPairsByGcd(v, 15);
        checkEq(a[1], 0, "case 3a failed");
        checkEq(a[2], 1, "case 3b failed");
        checkEq(a[3], 1, "case 3c failed");
        checkEq(a[5], 1, "case 3d failed");
    }
    {
        std::vector<int> v{7, 7, 7};
        std::vector<long long> a = countPairsByGcd(v, 7);
        checkEq(a[7], 3, "case 4 failed");
    }
    {
        std::vector<int> v{5};
        std::vector<long long> a = countPairsByGcd(v, 5);
        checkEq(a[1], 0, "case 5 failed");
    }
    {
        std::vector<int> v;
        std::vector<long long> a = countPairsByGcd(v, 5);
        checkEq(a[1], 0, "case 6 failed");
    }
    {
        std::vector<int> v{2, 4, 8};
        std::vector<long long> a = countPairsByGcd(v, 8);
        checkEq(a[2], 2, "case 7a failed");
        checkEq(a[4], 1, "case 7b failed");
        checkEq(a[1], 0, "case 7c failed");
    }
    {
        std::vector<int> v{1, 2, 3};
        std::vector<long long> a = countPairsByGcd(v, 3);
        checkEq(a[1], 3, "case 8a failed");
        checkEq(a[2], 0, "case 8b failed");
        checkEq(a[3], 0, "case 8c failed");
    }
    {
        std::vector<int> v{12, 18, 24};
        std::vector<long long> a = countPairsByGcd(v, 24);
        checkEq(a[6], 2, "case 9a failed");
        checkEq(a[12], 1, "case 9b failed");
        checkEq(a[1], 0, "case 9c failed");
    }
    {
        std::vector<int> v{4, 4, 4, 4, 4};
        std::vector<long long> a = countPairsByGcd(v, 4);
        checkEq(a[4], 10, "case 10a failed");
        checkEq(a[1], 0, "case 10b failed");
        checkEq(a[2], 0, "case 10c failed");
    }
    {
        std::vector<int> v{100000, 99999};
        std::vector<long long> a = countPairsByGcd(v, 100000);
        checkEq(a[1], 1, "case 11a failed");
        checkEq(a[100000], 0, "case 11b failed");
    }
    {
        std::vector<int> v{8, 9, 10, 11};
        std::vector<long long> a = countPairsByGcd(v, 11);
        checkEq(a[1], 5, "case 12a failed");
        checkEq(a[2], 1, "case 12b failed");
    }

    // ---- 與參考實作對拍 ----
    {
        Lcg rng(20261102ULL);
        for (int iter = 0; iter < 200; ++iter) {
            int n = rng.range(0, 12);
            std::vector<int> values;
            for (int i = 0; i < n; ++i) values.push_back(rng.range(1, 24));
            int maxq = rng.range(1, 24);
            std::vector<long long> a = countPairsByGcd(values, maxq);
            for (int g = 1; g <= maxq; ++g) {
                checkEq(a[static_cast<std::size_t>(g)], bruteForce(values, g),
                        "random mismatch");
            }
        }
    }
    {
        Lcg rng(99887766ULL);
        for (int iter = 0; iter < 80; ++iter) {
            int n = rng.range(0, 10);
            std::vector<int> values;
            for (int i = 0; i < n; ++i) values.push_back(rng.range(1, 200));
            int maxq = rng.range(1, 200);
            std::vector<long long> a = countPairsByGcd(values, maxq);
            for (int g = 1; g <= maxq; ++g) {
                checkEq(a[static_cast<std::size_t>(g)], bruteForce(values, g),
                        "mid random mismatch");
            }
        }
    }
    {
        Lcg rng(13572468ULL);
        for (int iter = 0; iter < 60; ++iter) {
            int v = rng.range(1, 30);
            int n = rng.range(2, 15);
            std::vector<int> values(static_cast<std::size_t>(n), v);
            int maxq = std::max(v, rng.range(1, 30));
            std::vector<long long> a = countPairsByGcd(values, maxq);
            long long wantTotal = 1LL * n * (n - 1) / 2;
            checkEq(a[static_cast<std::size_t>(v)], wantTotal, "same value failed");
            for (int g = 1; g <= maxq; ++g) {
                if (g != v) checkEq(a[static_cast<std::size_t>(g)], 0, "same value others");
            }
        }
    }
    {
        Lcg rng(24681357ULL);
        for (int iter = 0; iter < 40; ++iter) {
            int base = rng.range(2, 12);
            int n = rng.range(2, 12);
            std::vector<int> values;
            int mx = 0;
            for (int i = 0; i < n; ++i) {
                int v = base * rng.range(1, 8);
                values.push_back(v);
                mx = std::max(mx, v);
            }
            int maxq = mx + 5;
            std::vector<long long> a = countPairsByGcd(values, maxq);
            long long total = 0;
            for (int g = 1; g <= maxq; ++g) {
                long long want = bruteForce(values, g);
                checkEq(a[static_cast<std::size_t>(g)], want, "multiple base failed");
                total += want;
            }
            checkEq(total, 1LL * n * (n - 1) / 2, "total pairs failed");
        }
    }
    // 檢查莫比烏斯函數本身
    {
        std::vector<int> mu = mobiusSieve(100);
        checkEq(mu[1], 1, "mu 1");
        checkEq(mu[2], -1, "mu 2");
        checkEq(mu[3], -1, "mu 3");
        checkEq(mu[4], 0, "mu 4");
        checkEq(mu[5], -1, "mu 5");
        checkEq(mu[6], 1, "mu 6");
        checkEq(mu[7], -1, "mu 7");
        checkEq(mu[8], 0, "mu 8");
        checkEq(mu[9], 0, "mu 9");
        checkEq(mu[10], 1, "mu 10");
        checkEq(mu[12], 0, "mu 12");
        checkEq(mu[30], -1, "mu 30");
        // Σ_{d|n} μ(d) = [n == 1]
        for (int n = 1; n <= 100; ++n) {
            int s = 0;
            for (int d = 1; d * d <= n; ++d) {
                if (n % d == 0) {
                    s += mu[static_cast<std::size_t>(d)];
                    if (d * d != n) s += mu[static_cast<std::size_t>(n / d)];
                }
            }
            checkEq(s, n == 1 ? 1 : 0, "mu sum identity");
        }
    }
}

// ---------------------------------------------------------------- IO 模式

std::optional<long long> parseInt(const std::string& tok) {
    if (tok.empty()) return std::nullopt;
    std::size_t i = 0;
    if (tok[0] == '+' || tok[0] == '-') i = 1;
    if (i == tok.size()) return std::nullopt;
    for (std::size_t j = i; j < tok.size(); ++j) {
        if (tok[j] < '0' || tok[j] > '9') return std::nullopt;
    }
    __int128 val = 0;
    bool neg = (tok[0] == '-');
    for (std::size_t j = i; j < tok.size(); ++j) {
        val = val * 10 + (tok[j] - '0');
        if (val > (static_cast<__int128>(1) << 100)) return std::nullopt;
    }
    if (neg) val = -val;
    if (val > static_cast<__int128>(PARSE_LIMIT)) return std::nullopt;
    if (val < -static_cast<__int128>(PARSE_LIMIT)) return std::nullopt;
    return static_cast<long long>(val);
}

void runIo(const std::string& raw) {
    std::vector<std::string> toks;
    {
        std::size_t i = 0;
        while (i < raw.size()) {
            while (i < raw.size() && static_cast<unsigned char>(raw[i]) <= ' ') ++i;
            if (i >= raw.size()) break;
            std::size_t j = i;
            while (j < raw.size() && static_cast<unsigned char>(raw[j]) > ' ') ++j;
            toks.push_back(raw.substr(i, j - i));
            i = j;
        }
    }
    std::size_t pos = 0;
    auto nxt = [&toks, &pos]() -> const std::string* {
        if (pos >= toks.size()) return nullptr;
        return &toks[pos++];
    };
    auto bad = []() { std::cout << "error\n"; };

    const std::string* pt = nxt();
    std::optional<long long> nOpt = std::nullopt;
    if (pt) nOpt = parseInt(*pt);
    if (!nOpt || *nOpt < 1 || *nOpt > NMAX) {
        bad();
        return;
    }
    const std::string* pt2 = nxt();
    std::optional<long long> qOpt = std::nullopt;
    if (pt2) qOpt = parseInt(*pt2);
    if (!qOpt || *qOpt < 0 || *qOpt > QMAX) {
        bad();
        return;
    }
    int n = static_cast<int>(*nOpt);
    int q = static_cast<int>(*qOpt);

    std::vector<int> values;
    values.reserve(n);
    for (int i = 0; i < n; ++i) {
        const std::string* t = nxt();
        std::optional<long long> v = std::nullopt;
        if (t) v = parseInt(*t);
        if (!v || *v < 1 || *v > VAL_MAX) {
            bad();
            return;
        }
        values.push_back(static_cast<int>(*v));
    }

    std::vector<int> queries;
    queries.reserve(q);
    int maxq = 1;
    for (int i = 0; i < q; ++i) {
        const std::string* t = nxt();
        std::optional<long long> g = std::nullopt;
        if (t) g = parseInt(*t);
        if (!g || *g < 1 || *g > VAL_MAX) {
            bad();
            return;
        }
        queries.push_back(static_cast<int>(*g));
        maxq = std::max(maxq, static_cast<int>(*g));
    }

    std::vector<long long> ans = countPairsByGcd(values, maxq);
    for (int g : queries) {
        std::cout << ans[static_cast<std::size_t>(g)] << "\n";
    }
}

}  // namespace

int main() {
#ifdef _WIN32
    _setmode(_fileno(stdout), _O_BINARY);
#endif
    std::string raw((std::istreambuf_iterator<char>(std::cin)),
                    std::istreambuf_iterator<char>());
    bool blank = true;
    for (std::size_t i = 0; i < raw.size(); ++i) {
        if (static_cast<unsigned char>(raw[i]) > ' ') {
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
