// Miller-Rabin 質數測試與 Pollard-Rho 因數分解
//
// 題意：
//     給定若干個非負整數（0 ≤ n ≤ 10^12），對每一個輸出：
//       prime   n 是否為質數（1 / 0）
//       lpf     最大的質因數（n < 2 時為 0）
//       omega   相異質因數個數 ω(n)
//       Omega   質因數總數（含重複）Ω(n)
//       phi     尤拉函數 φ(n)
//       sigma   因數和 σ(n)
//       tau     因數個數 τ(n)
//       fac     質因數分解，寫成 p^e 用 * 串起來（n ≤ 1 時為 -）
//     以上每一項都與**獨立的對照算法**互相驗證（試除法分解、枚舉因數、
//     暴力數互質個數、∑_{d|n} φ(d) = n 這個恆等式）。
//
// 思路：
//     ### 費馬小定理為什麼不夠
//     若 p 是質數，則 a^(p-1) ≡ 1 (mod p)。反過來「a^(n-1) ≡ 1 就說 n 是質數」
//     會被 Carmichael 數（例如 561 = 3·11·17）騙到：它們對所有與 n 互質的 a 都滿足費馬式。
//
//     ### Miller-Rabin：把「平方剩餘」的條件也加進來
//     把 n-1 寫成 d·2^s（d 為奇數）。若 n 是質數，則對任意 a 都有
//         a^d ≡ 1 (mod n)  或  存在某個 r ∈ [0, s) 使 a^(d·2^r) ≡ -1 (mod n)
//     （因為在質數的乘法群裡 x² ≡ 1 只有 x ≡ ±1 兩個解，所以「一路平方到 1」的過程
//       中最後一個不是 1 的值必定是 -1。）
//     任一底數不通過 → n 一定是合數（單邊錯誤）。本題用前 12 個質數當底數
//     {2,3,5,7,11,13,17,19,23,29,31,37}：這組對所有 n < 3.3×10^24 都是**確定性**的，
//     所以對 n ≤ 10^12 而言本程式的質數測試沒有任何機率成分。
//
//     ### Pollard-Rho：在模 n 的世界裡走一趟函數圖，撞到「同餘」就拿到因數
//     這是本倉庫 2026-10-17 `functional-graph` 那一題的直接應用：
//     迭代 x ← (x² + c) mod n，值域只有 n 個值，路徑一定形成 ρ 形狀（尾巴 + 環）。
//     若 p | n，則「模 p 的序列」也會撞環；用 Floyd 龜兔法維護 d = gcd(|x - y|, n)，
//     一旦 x ≡ y (mod p) 但 x ≢ y (mod n)，就得到 1 < d < n，即 n 的非平凡因數。
//     期望步數 O(√p)（p 為 n 最小的質因數），比試除法的 O(√n) 快非常多。
//     d == n 表示這條參數 c 失敗，換一個 c 重來即可。
//
//     ### 因數分解 → φ / σ / τ
//     分解出 n = Π p_i^e_i 後用積性公式
//         φ(n) = Π p_i^(e_i-1)·(p_i - 1)
//         σ(n) = Π (p_i^(e_i+1) - 1)/(p_i - 1)
//         τ(n) = Π (e_i + 1)
//     而程式另外「枚舉所有因數」把 τ 與 σ 直接算出來對拍，並用 ∑_{d|n} φ(d) = n 驗證自洽。
//
// 應用場景：
//     - 密碼學：RSA 需要兩個大質數（用 Miller-Rabin 產生），而 RSA 的安全性正是建立
//       在大整數「難以分解」上，Pollard-Rho 是暴力分解的主力之一
//     - 競賽題：判斷大數是否為質數、求最大質因數、計算 φ / σ / τ、
//       處理「n 的所有因數」（枚舉因數比枚舉 1..n 快得多）
//     - 工程：雜湊表大小取質數、亂數週期長度設計、區塊鏈與零知識證明裡的模運算
//     - ρ 形思路也出現在：生日攻擊、雜湊碰撞搜尋、離散對數的 BSGS
//
// 複雜度：
//     Miller-Rabin   O(k log³ n)（k = 12 個底數）
//     Pollard-Rho    期望 O(√p) 次迭代找出一個因數（p 為最小質因數）
//     完整分解       期望 O(n^1/4 log n) 等級，對 10^12 是瞬間完成
//     空間           O(log n)

#include <algorithm>
#include <cassert>
#include <cerrno>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <iterator>
#include <numeric>
#include <random>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#ifdef _WIN32
#include <fcntl.h>
#include <io.h>
#endif

namespace {

const long long NMAX = 1000000000000LL;   // IO 模式下單一數值的上限（10^12）
const int KMAX = 20;                       // 最多 20 個數值
const long long BRUTE_TRIAL = 1000000;     // ≤ 此值就跑試除法分解 + 枚舉因數對拍
const long long BRUTE_PHI = 20000;         // ≤ 此值才用「數互質個數」對拍 φ

const int MR_BASES[12] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37};

const uint64_t SEED_MAIN = 20261017ULL;
const uint64_t SEED_CHECK = 987654321ULL;   // 第二次獨立分解用的種子

using int64 = long long;

// ------------------------------------------------------- 基本工具

int64 mulMod(int64 a, int64 b, int64 mod) {
    return static_cast<int64>(static_cast<unsigned __int128>(static_cast<uint64_t>(a)) *
                              static_cast<uint64_t>(b) % static_cast<uint64_t>(mod));
}

int64 powMod(int64 a, int64 e, int64 mod) {
    int64 res = 1 % mod;
    a %= mod;
    while (e > 0) {
        if (e & 1) res = mulMod(res, a, mod);
        a = mulMod(a, a, mod);
        e >>= 1;
    }
    return res;
}

int64 gcd64(int64 a, int64 b) {
    while (b) {
        int64 t = a % b;
        a = b;
        b = t;
    }
    return a < 0 ? -a : a;
}

// ------------------------------------------------------- 確定性亂數（LCG）

struct Rng {
    uint64_t s;
    explicit Rng(uint64_t seed) : s(seed) {}
    uint64_t next() {
        s = s * 6364136223846793005ULL + 1442695040888963407ULL;
        return s;
    }
};

// ------------------------------------------------------- Miller-Rabin

bool isPrime(int64 n) {
    if (n < 2) return false;
    for (int p : MR_BASES) {
        if (n % p == 0) return n == p;
    }
    int64 d = n - 1;
    int s = 0;
    while (d % 2 == 0) {
        d /= 2;
        ++s;
    }
    for (int a : MR_BASES) {
        int64 x = powMod(a, d, n);
        if (x == 1 || x == n - 1) continue;
        bool witness = true;               // 沒有出現 -1 → a 是「n 為合數」的見證
        for (int i = 0; i < s - 1; ++i) {
            x = mulMod(x, x, n);
            if (x == n - 1) {
                witness = false;
                break;
            }
        }
        if (witness) return false;
    }
    return true;
}

bool isPrimeBrute(int64 n) {
    if (n < 2) return false;
    for (int64 d = 2; d * d <= n; ++d) {
        if (n % d == 0) return false;
    }
    return true;
}

// ------------------------------------------------------- Pollard-Rho

int64 pollardRho(int64 n, Rng& rng) {
    if (n % 2 == 0) return 2;
    if (n % 3 == 0) return 3;
    while (true) {
        int64 c = static_cast<int64>(rng.next() % static_cast<uint64_t>(n - 1)) + 1;
        int64 x = static_cast<int64>(rng.next() % static_cast<uint64_t>(n - 1)) + 1;
        int64 y = x;
        int64 d = 1;
        while (d == 1) {
            x = (mulMod(x, x, n) + c) % n;
            y = (mulMod(y, y, n) + c) % n;
            y = (mulMod(y, y, n) + c) % n;
            int64 diff = x >= y ? x - y : y - x;
            d = gcd64(diff, n);
        }
        if (d != n) return d;
    }
}

std::vector<int64> factorize(int64 n, uint64_t seed) {
    std::vector<int64> out;
    if (n < 2) return out;
    Rng rng(seed);
    std::vector<int64> stack(1, n);
    while (!stack.empty()) {
        int64 v = stack.back();
        stack.pop_back();
        if (v == 1) continue;
        if (isPrime(v)) {
            out.push_back(v);
            continue;
        }
        int64 d = pollardRho(v, rng);
        stack.push_back(d);
        stack.push_back(v / d);
    }
    std::sort(out.begin(), out.end());
    return out;
}

std::vector<int64> trialFactor(int64 n) {
    std::vector<int64> res;
    if (n < 2) return res;
    for (int64 d = 2; d * d <= n; ++d) {
        while (n % d == 0) {
            res.push_back(d);
            n /= d;
        }
    }
    if (n > 1) res.push_back(n);
    return res;
}

// ------------------------------------------------------- 因數與積性函數

std::vector<std::pair<int64, int>> groupFactors(const std::vector<int64>& fac) {
    std::vector<std::pair<int64, int>> out;
    for (int64 p : fac) {
        if (!out.empty() && out.back().first == p) {
            out.back().second += 1;
        } else {
            out.emplace_back(p, 1);
        }
    }
    return out;
}

std::vector<int64> allDivisors(const std::vector<int64>& fac) {
    std::vector<int64> divs(1, 1);
    for (auto& pe : groupFactors(fac)) {
        int64 p = pe.first;
        int e = pe.second;
        std::vector<int64> nxt;
        int64 mul = 1;
        for (int k = 0; k <= e; ++k) {
            for (int64 d : divs) nxt.push_back(d * mul);
            mul *= p;
        }
        divs.swap(nxt);
    }
    std::sort(divs.begin(), divs.end());
    return divs;
}

// 枚舉所有因數並順便算出每個因數的 φ：[(因數, φ(因數))]
std::vector<std::pair<int64, int64>> divisorsWithPhi(const std::vector<std::pair<int64, int>>& g) {
    std::vector<std::pair<int64, int64>> res(1, {1, 1});
    for (auto& pe : g) {
        int64 p = pe.first;
        int e = pe.second;
        std::vector<std::pair<int64, int64>> nxt;
        for (auto& cur : res) {
            int64 d = cur.first, phid = cur.second;
            int64 pk = 1;
            for (int k = 0; k <= e; ++k) {
                if (k == 0) {
                    nxt.emplace_back(d, phid);
                } else {
                    pk *= p;
                    nxt.emplace_back(d * pk, phid * ((p - 1) * (pk / p)));
                }
            }
        }
        res.swap(nxt);
    }
    return res;
}

int64 phiFrom(const std::vector<int64>& fac) {
    int64 res = 1;
    for (auto& pe : groupFactors(fac)) {
        int64 p = pe.first;
        int e = pe.second;
        int64 pw = 1;
        for (int i = 0; i < e - 1; ++i) pw *= p;
        res *= (p - 1) * pw;
    }
    return res;
}

int64 phiBrute(int64 n) {
    if (n < 1) return 0;
    int64 cnt = 0;
    for (int64 i = 1; i <= n; ++i) {
        if (gcd64(i, n) == 1) ++cnt;
    }
    return cnt;
}

std::pair<int64, int64> bruteDivisors(int64 n) {
    if (n < 1) return {0, 0};
    int64 tau = 0, sigma = 0;
    for (int64 d = 1; d * d <= n; ++d) {
        if (n % d == 0) {
            ++tau;
            sigma += d;
            int64 q = n / d;
            if (q != d) {
                ++tau;
                sigma += q;
            }
        }
    }
    return {tau, sigma};
}

// ------------------------------------------------------- 輸出

std::string fmtFactors(const std::vector<int64>& fac) {
    if (fac.empty()) return "-";
    std::ostringstream oss;
    bool first = true;
    for (auto& pe : groupFactors(fac)) {
        if (!first) oss << "*";
        first = false;
        oss << pe.first << "^" << pe.second;
    }
    return oss.str();
}

std::string run(const std::vector<int64>& nums) {
    std::ostringstream oss;
    int ok = 1;
    for (int64 n : nums) {
        std::vector<int64> fac = factorize(n, SEED_MAIN);
        std::vector<int64> fac2 = factorize(n, SEED_CHECK);      // 換一組亂數種子再分解一次
        if (fac != fac2) ok = 0;
        if (!std::is_sorted(fac.begin(), fac.end())) ok = 0;
        int64 prod = 1;
        for (int64 p : fac) {
            prod *= p;
            if (!isPrime(p)) ok = 0;
        }
        if (n >= 2 && prod != n) ok = 0;
        if (n < 2 && !fac.empty()) ok = 0;

        int prime = isPrime(n) ? 1 : 0;
        if (2 <= n && n <= BRUTE_TRIAL && prime != (isPrimeBrute(n) ? 1 : 0)) ok = 0;

        std::vector<std::pair<int64, int>> g = groupFactors(fac);
        int64 omega = static_cast<int64>(g.size());
        int64 omegaTotal = static_cast<int64>(fac.size());
        int64 lpf = g.empty() ? 0 : g.back().first;
        int64 phi = n >= 1 ? phiFrom(fac) : 0;
        int64 tauFormula = 1;
        for (auto& pe : g) tauFormula *= pe.second + 1;
        if (n < 1) tauFormula = 0;

        int64 tau = 0, sigma = 0;
        if (n >= 1) {
            std::vector<std::pair<int64, int64>> dp = divisorsWithPhi(g);
            tau = static_cast<int64>(dp.size());
            sigma = 0;
            int64 sumPhi = 0;
            for (auto& cur : dp) {
                sigma += cur.first;
                sumPhi += cur.second;
            }
            if (tau != tauFormula) ok = 0;
            if (sumPhi != n) ok = 0;                             // ∑_{d|n} φ(d) = n
            if (n <= BRUTE_TRIAL) {
                std::pair<int64, int64> bd = bruteDivisors(n);
                if (tau != bd.first || sigma != bd.second) ok = 0;
                if (fac != trialFactor(n)) ok = 0;
            }
            if (1 <= n && n <= BRUTE_PHI && phi != phiBrute(n)) ok = 0;
        }

        oss << n << " " << prime << " " << lpf << " " << omega << " " << omegaTotal << " "
            << phi << " " << sigma << " " << tau << " " << fmtFactors(fac) << "\n";
    }
    oss << ok << "\n";
    return oss.str();
}

// ------------------------------------------------------- 內置測試

void runTests() {
    // ---- 固定用例：質數 ----
    const int64 primes[] = {2, 3, 5, 7, 11, 13, 97, 1009, 7919, 104729, 999983, 1000003};
    for (int64 p : primes) {
        assert(isPrime(p));
        assert(factorize(p, SEED_MAIN) == std::vector<int64>(1, p));
    }
    const int64 comps[] = {0, 1, 4, 9, 25, 561, 1105, 1729, 2465, 2821, 6601, 8911};
    for (int64 c : comps) assert(!isPrime(c));

    // Carmichael 數：對所有與 n 互質的 a 都滿足費馬式，但確實是合數
    assert(!isPrime(561));
    {
        std::vector<int64> want = {3, 11, 17};
        assert(factorize(561, SEED_MAIN) == want);
    }

    // ---- 固定用例：分解 ----
    assert(factorize(1, SEED_MAIN).empty());
    assert(factorize(0, SEED_MAIN).empty());
    {
        std::vector<int64> want = {2, 2, 3};
        assert(factorize(12, SEED_MAIN) == want);
        std::vector<int64> w2 = {2, 2, 5, 5};
        assert(factorize(100, SEED_MAIN) == w2);
        std::vector<int64> w3(10, 2);
        assert(factorize(1024, SEED_MAIN) == w3);
        std::vector<int64> w4;
        for (int i = 0; i < 12; ++i) w4.push_back(2);
        for (int i = 0; i < 12; ++i) w4.push_back(5);
        assert(factorize(NMAX, SEED_MAIN) == w4);
        std::vector<int64> w5 = {999983, 1000003};            // 兩個大質數的乘積
        assert(factorize(999985999949LL, SEED_MAIN) == w5);
    }

    // ---- 固定用例：積性函數 ----
    {
        std::vector<int64> f = {2, 2, 3};
        assert(phiFrom(f) == 4);                              // φ(12) = 4
        std::vector<int64> f7 = {7};
        assert(phiFrom(f7) == 6);
        assert(phiFrom(std::vector<int64>()) == 1);           // φ(1) = 1
        assert(bruteDivisors(12) == std::make_pair(6LL, 28LL));
        assert(bruteDivisors(1) == std::make_pair(1LL, 1LL));
        assert(phiBrute(12) == 4);
        assert(phiBrute(1) == 1);
        std::vector<int64> want = {1, 2, 3, 4, 6, 12};
        assert(allDivisors(f) == want);
    }

    // 整行輸出
    {
        std::string a = run(std::vector<int64>(1, 12));
        assert(a.substr(0, a.find('\n')) == "12 0 3 2 3 4 28 6 2^2*3^1");
        std::string b = run(std::vector<int64>(1, 7));
        assert(b.substr(0, b.find('\n')) == "7 1 7 1 1 6 8 2 7^1");
        std::string c = run(std::vector<int64>(1, 1));
        assert(c.substr(0, c.find('\n')) == "1 0 0 0 0 1 1 1 -");
        std::string d = run(std::vector<int64>(1, 0));
        assert(d.substr(0, d.find('\n')) == "0 0 0 0 0 0 0 0 -");
    }

    // ---- 隨機對拍 ----
    {
        std::mt19937 rng(20261017);
        const int64 semis[8] = {2, 3, 5, 7, 997, 10007, 99991, 999983};
        for (int it = 0; it < 400; ++it) {
            int kind = static_cast<int>(rng() % 6);
            int64 n;
            if (kind == 0) {
                n = static_cast<int64>(rng() % 61);
            } else if (kind == 1) {
                n = static_cast<int64>(rng() % 20001);
            } else if (kind == 2) {
                n = static_cast<int64>(rng() % (BRUTE_TRIAL + 1));
            } else if (kind == 3) {
                n = static_cast<int64>((static_cast<uint64_t>(rng()) << 20 | rng()) % static_cast<uint64_t>(NMAX + 1));
            } else if (kind == 4) {
                int64 a = semis[rng() % 8], b = semis[rng() % 8];
                n = std::min(a * b, NMAX);
            } else {
                int64 base = 2 + static_cast<int64>(rng() % 59);
                int64 pw = 2 + static_cast<int64>(rng() % 5);
                n = 1;
                for (int i = 0; i < pw; ++i) n *= base;
                n = std::min(n, NMAX);
            }
            std::string out = run(std::vector<int64>(1, n));
            assert(out.size() >= 2 && out[out.size() - 2] == '1');
        }
    }

    // ---- 一批數字同時處理 ----
    {
        std::vector<int64> batch = {0, 1, 2, 561, 1024, 1000000, 999983, NMAX};
        std::string out = run(batch);
        assert(out.size() >= 2 && out[out.size() - 2] == '1');
    }
}

// ------------------------------------------------------- IO

bool isInt(const std::string& t) {
    if (t.empty()) return false;
    size_t i = 0;
    if (t[0] == '+' || t[0] == '-') i = 1;
    if (i >= t.size()) return false;
    for (; i < t.size(); ++i) {
        if (t[i] < '0' || t[i] > '9') return false;
    }
    return true;
}

long long clampLL(const std::string& t, long long lo, long long hi, long long def) {
    if (!isInt(t)) return def;
    errno = 0;
    char* end = nullptr;
    long long v = std::strtoll(t.c_str(), &end, 10);
    if (errno == ERANGE || end != t.c_str() + t.size()) return def;
    if (v < lo) v = lo;
    if (v > hi) v = hi;
    return v;
}

void runIo(const std::string& raw) {
    std::istringstream iss(raw);
    std::vector<std::string> toks;
    std::string t;
    while (iss >> t) toks.push_back(t);
    auto parseAt = [&](size_t i, long long def) -> long long {
        if (i >= toks.size()) return def;
        return clampLL(toks[i], -1000000000000000000LL, 1000000000000000000LL, def);
    };
    long long k = parseAt(0, 0);
    if (k < 0) k = 0;
    if (k > KMAX) k = KMAX;
    std::vector<int64> nums;
    for (long long i = 1; i <= k; ++i) {
        long long v = parseAt(static_cast<size_t>(i), 0);
        if (v < 0) v = 0;
        if (v > NMAX) v = NMAX;
        nums.push_back(v);
    }
    std::cout << run(nums);
}

}  // namespace

int main() {
#ifdef _WIN32
    _setmode(_fileno(stdout), _O_BINARY);
#endif
    std::string raw((std::istreambuf_iterator<char>(std::cin)), std::istreambuf_iterator<char>());
    bool blank = true;
    for (char c : raw) {
        unsigned char u = static_cast<unsigned char>(c);
        if (u > ' ') { blank = false; break; }
    }
    if (!blank) {
        runIo(raw);
    } else {
        runTests();
        std::cout << "all tests passed" << std::endl;
    }
    return 0;
}
