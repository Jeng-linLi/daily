// 中國剩餘定理（CRT）與擴展中國剩餘定理（exCRT）
//
// 題意：
//     一口氣解四組同餘問題，每組都同時給出「高效演算法」與「獨立實作的暴力版本」互相對拍。
//
//     A. **同餘方程組**：給 n 組 (a_i, m_i)，求所有滿足 x ≡ a_i (mod m_i) 的 x。
//        分別用「互質版 CRT」（乘積 + 逆元）與「一般版 exCRT」（兩兩合併）各解一次。
//     B. **線性同餘方程** a·x ≡ b (mod m)：求 gcd、解的個數、最小非負解。
//     C. **CRT 還原**：給定 x 與兩個模數 p、q，先取 (x mod p, x mod q)，
//        再用 exCRT 還原出 x mod lcm(p, q)，跟直接算的 x mod lcm(p, q) 比對。
//
// 思路：
//     ### 1. 擴展歐幾里得（exgcd）
//     迭代版即可：維護 (r, s, t) 使 r = a·s + b·t，不斷做帶餘除法直到 r = 0。
//     限制輸入為**非負**，就完全避開「Python 的 // 是 floor、C++ 的 / 是 truncate」
//     在負數上會分歧的老坑。
//
//     ### 2. 合併兩個同餘式
//     寫 x = a₁ + m₁·t，則 m₁·t ≡ a₂ - a₁ (mod m₂)。令 g = gcd(m₁, m₂)，
//     (a₂-a₁) 不是 g 的倍數 → **無解**；否則同除以 g：
//         (m₁/g)·t ≡ (a₂-a₁)/g  (mod m₂/g)
//     此時 m₁/g 與 m₂/g 互質，逆元存在，回代得 x₀ = a₁ + m₁·t (mod lcm)。
//     逐條合併下去就是 exCRT，模數不互質也能用。
//
//     ### 3. 互質版 CRT（孫子定理的標準形式）
//     M = ∏ m_i，M_i = M / m_i，inv_i = M_i^{-1} mod m_i，則
//         x = Σ a_i · M_i · inv_i  (mod M)
//     因為 M_i·inv_i ≡ 1 (mod m_i) 且 ≡ 0 (mod m_j, j≠i)，每一項只在自己那一維生效。
//     這個形式**只在模數兩兩互質時成立**，非互質必須回到第 2 步。
//
//     ### 4. 線性同餘 a·x ≡ b (mod m)
//     g = gcd(a, m)；b 不是 g 的倍數就無解；否則
//         x₀ = (b/g) · (a/g)^{-1}  (mod m/g)
//     在 [0, m) 裡恰好有 g 個解：x₀, x₀ + m/g, …
//
//     ### 5. 溢位處理（跨語言一致）
//     Python 無限精度、C++ 的 int64 會繞回，兩邊必須用同一套截斷規則：
//     合併模數時先做**除法檢查**「M > LIMIT / (m/g) 就判定溢位」（回傳 -2 -2），
//     乘法則在 C++ 用 __int128 承接。LIMIT 定為 10^18。
//
//     ### 6. 暴力對拍
//     lcm 不超過 200000 時從 0 枚舉到 lcm-1 找最小解；
//     線性同餘則在 m ≤ 200000 時枚舉 [0, m)。這能抓出「公式對但正規化錯」的 bug。
//
// 應用場景：
//     - **RSA 的 CRT 加速**：分別在模 p 與模 q 下解密，最後 CRT 合併，速度約 4 倍
//     - **大數運算**：把巨大整數拆成幾個小模數下的運算，最後 CRT 還原（NTT 多模合併）
//     - **秘密共享（Asmuth-Bloom 門限方案）**：每人拿到不同模數下的餘數，
//       湊齊足夠多人才能還原祕密，少一個都不行
//     - **排程 / 週期對齊**：「每 3 天、每 5 天、每 7 天一次」何時同時發生
//     - 競賽題裡的「韓信點兵」、模數非互質的計數題
//
// 複雜度：
//     exgcd / 逆元            O(log min(a,b)) 時間 / O(1) 空間
//     exCRT 合併 n 條          O(n log M) 時間 / O(n) 空間
//     互質版 CRT               O(n log M) 時間 / O(n) 空間
//     線性同餘                 O(log m) 時間 / O(1) 空間
//     暴力對拍                 O(lcm) 或 O(m) 時間

#include <algorithm>
#include <cerrno>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#ifdef _WIN32
#include <fcntl.h>
#include <io.h>
#endif

namespace {

using i64 = long long;
using i128 = __int128;

const i64 LIMIT = 1000000000000000000LL;   // 合併模數上限（10^18）
const i64 BRUTE_LIMIT = 200000;            // lcm / m 不超過這個值才跑暴力
const int NMAX = 12;                       // 同餘方程組條數上限
const i64 MMAX = 100000;                   // 單一模數上限
const i64 AMAX = 1000000000;               // 餘數 / 係數上限
const i64 XMAX = LIMIT;                    // 還原用的大數上限
const i64 PQMAX = 1000000;                 // 還原用的兩個模數上限

// ---------------------------------------------------------------- 數論基礎

// 迭代版擴展歐幾里得：回傳 {g, s, t} 使 a·s + b·t = g。要求 a, b ≥ 0。
struct Exg {
    i64 g, s, t;
};

Exg exgcd(i64 a, i64 b) {
    i64 old_r = a, r = b;
    i64 old_s = 1, s = 0;
    i64 old_t = 0, t = 1;
    while (r != 0) {
        i64 q = old_r / r;
        i64 nr = old_r - q * r;
        i64 ns = old_s - q * s;
        i64 nt = old_t - q * t;
        old_r = r;
        old_s = s;
        old_t = t;
        r = nr;
        s = ns;
        t = nt;
    }
    if (old_r < 0) {
        old_r = -old_r;
        old_s = -old_s;
        old_t = -old_t;
    }
    return Exg{old_r, old_s, old_t};
}

i64 gcd(i64 a, i64 b) {
    while (b) {
        i64 t = a % b;
        a = b;
        b = t;
    }
    return a;
}

// a 在模 m 下的乘法逆元。m == 1 回傳 0；不存在回傳 -1。
i64 modInverse(i64 a, i64 m) {
    if (m == 1) return 0;
    Exg e = exgcd(a % m, m);
    if (e.g != 1) return -1;
    i64 r = e.s % m;
    if (r < 0) r += m;
    return r;
}

// 合併 x ≡ a1 (mod m1) 與 x ≡ a2 (mod m2)。回傳 {狀態, r, lcm}。
// 狀態：0 成功 / 1 無解 / 2 溢位。
struct MergeRes {
    int st;
    i64 r;
    i64 M;
};

MergeRes merge2(i64 a1, i64 m1, i64 a2, i64 m2) {
    if (m1 <= 0 || m2 <= 0) return MergeRes{1, -1, -1};
    a1 %= m1;
    a2 %= m2;
    i64 g = gcd(m1, m2);
    i64 diff = a2 - a1;
    if (diff % g != 0) return MergeRes{1, -1, -1};
    i64 step = m2 / g;
    if (m1 > LIMIT / step) return MergeRes{2, -2, -2};
    i64 M = m1 * step;
    i64 m2g = m2 / g;
    i64 inv = modInverse((m1 / g) % m2g, m2g);
    if (inv < 0) return MergeRes{1, -1, -1};
    i64 t = static_cast<i64>((static_cast<i128>((diff / g) % m2g) * inv) % m2g);
    i64 r = static_cast<i64>((static_cast<i128>(a1) + static_cast<i128>(m1 % M) * t) % M);
    if (r < 0) r += M;
    return MergeRes{0, r, M};
}

// ---------------------------------------------------------------- A. 同餘方程組

struct SysRes {
    int st;
    i64 r;
    i64 M;
};

SysRes excrt(const std::vector<std::pair<i64, i64>>& eqs) {
    i64 r = 0, M = 1;
    for (const auto& e : eqs) {
        MergeRes res = merge2(r, M, e.first % e.second, e.second);
        if (res.st != 0) return SysRes{res.st, res.r, res.M};
        r = res.r;
        M = res.M;
    }
    return SysRes{0, r % M, M};
}

SysRes crtCoprime(const std::vector<std::pair<i64, i64>>& eqs) {
    int n = static_cast<int>(eqs.size());
    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            if (gcd(eqs[i].second, eqs[j].second) != 1) return SysRes{1, -1, -1};
        }
    }
    i64 M = 1;
    for (const auto& e : eqs) {
        if (M > LIMIT / e.second) return SysRes{2, -2, -2};
        M *= e.second;
    }
    i64 r = 0;
    for (const auto& e : eqs) {
        i64 m = e.second;
        if (m == 0) return SysRes{1, -1, -1};
        i64 Mi = M / m;
        i64 inv = modInverse(Mi % m, m);
        if (inv < 0) return SysRes{1, -1, -1};
        i64 term = static_cast<i64>((static_cast<i128>(e.first % m) * (Mi % M) % M) * inv % M);
        r = (r + term) % M;
    }
    return SysRes{0, r % M, M};
}

// 所有模數的 lcm；超過 LIMIT 回傳 -1。
i64 lcmAll(const std::vector<i64>& mods) {
    i64 L = 1;
    for (i64 m : mods) {
        i64 g = gcd(L, m);
        i64 step = m / g;
        if (L > LIMIT / step) return -1;
        L *= step;
    }
    return L;
}

// 暴力枚舉 [0, lcm)。回傳 {狀態, 最小解}。狀態 0 有解 / 1 無解 / 2 太大不跑。
std::pair<int, i64> bruteSystem(const std::vector<std::pair<i64, i64>>& eqs) {
    std::vector<i64> mods;
    for (const auto& e : eqs) mods.push_back(e.second);
    i64 L = lcmAll(mods);
    if (L < 0 || L > BRUTE_LIMIT) return std::make_pair(2, -1LL);
    for (i64 x = 0; x < L; ++x) {
        bool good = true;
        for (const auto& e : eqs) {
            if (e.second > 0 && x % e.second != e.first % e.second) {
                good = false;
                break;
            }
        }
        if (good) return std::make_pair(0, x);
    }
    return std::make_pair(1, -1LL);
}

// ---------------------------------------------------------------- B. 線性同餘

struct LinRes {
    i64 g, cnt, xmin;
};

LinRes solveLinear(i64 a, i64 b, i64 m) {
    if (m == 0) return LinRes{0, 0, -1};
    a %= m;
    b %= m;
    i64 g = gcd(a, m);
    if (b % g != 0) return LinRes{g, 0, -1};
    i64 m2 = m / g;
    if (m2 == 1) return LinRes{g, g, 0};
    i64 inv = modInverse((a / g) % m2, m2);
    if (inv < 0) return LinRes{g, 0, -1};
    i64 x0 = static_cast<i64>(static_cast<i128>((b / g) % m2) * inv % m2);
    return LinRes{g, g, x0};
}

// 暴力枚舉 [0, m)。m 太大回傳 {-1,-1,-1}。
LinRes bruteLinear(i64 a, i64 b, i64 m) {
    if (m == 0 || m > BRUTE_LIMIT) return LinRes{-1, -1, -1};
    a %= m;
    b %= m;
    i64 g = gcd(a, m);
    i64 cnt = 0, first = -1;
    for (i64 x = 0; x < m; ++x) {
        if (static_cast<i64>(static_cast<i128>(a) * x % m) == b) {
            ++cnt;
            if (first < 0) first = x;
        }
    }
    return LinRes{g, cnt, first};
}

// ---------------------------------------------------------------- C. CRT 還原

std::pair<i64, i64> crtRestore(i64 x, i64 p, i64 q) {
    if (p <= 0 || q <= 0) return std::make_pair(-1LL, -1LL);
    MergeRes res = merge2(x % p, p, x % q, q);
    if (res.st != 0) return std::make_pair(-2LL, -2LL);
    return std::make_pair(res.r, res.M);
}

// ---------------------------------------------------------------- 主流程

std::string run(const std::vector<std::pair<i64, i64>>& eqs, i64 a, i64 b, i64 m,
                i64 x, i64 p, i64 q) {
    bool ok = true;

    SysRes s1 = excrt(eqs);
    SysRes s2 = crtCoprime(eqs);

    if (s1.st == 0) {
        for (const auto& e : eqs) {
            if (e.second > 0 && s1.r % e.second != e.first % e.second) {
                ok = false;
                break;
            }
        }
    }
    std::pair<int, i64> bs = bruteSystem(eqs);
    if (bs.first == 0) {
        if (s1.st != 0 || s1.r != bs.second) ok = false;
    } else if (bs.first == 1) {
        if (s1.st != 1) ok = false;
    }
    if (s2.st == 0) {
        if (s1.st != 0 || s1.r != s2.r || s1.M != s2.M) ok = false;
    }

    LinRes lr = solveLinear(a, b, m);
    if (lr.cnt > 0 && m > 0) {
        if (static_cast<i64>(static_cast<i128>(a % m) * lr.xmin % m) != b % m) ok = false;
    }
    LinRes bl = bruteLinear(a, b, m);
    if (bl.g != -1) {
        if (lr.g != bl.g || lr.cnt != bl.cnt || lr.xmin != bl.xmin) ok = false;
    }

    std::pair<i64, i64> rs = crtRestore(x, p, q);
    if (rs.second > 0) {
        if (rs.first != x % rs.second) ok = false;
    } else {
        ok = false;
    }

    std::ostringstream oss;
    oss << s1.r << ' ' << s1.M << '\n';
    oss << s2.r << ' ' << s2.M << '\n';
    oss << lr.g << ' ' << lr.cnt << ' ' << lr.xmin << '\n';
    oss << rs.first << ' ' << rs.second << '\n';
    oss << (ok ? '1' : '0') << '\n';
    return oss.str();
}

// ---------------------------------------------------------------- 內置測試

struct Lcg {
    unsigned long long s;
    explicit Lcg(unsigned long long seed) : s(seed) {}
    unsigned long long next() {
        s = s * 6364136223846793005ULL + 1442695040888963407ULL;
        return (s >> 32) & 0xFFFFFFFFULL;
    }
    long long below(long long n) { return static_cast<long long>(next() % static_cast<unsigned long long>(n)); }
};

void expect(const std::vector<std::pair<i64, i64>>& eqs, i64 a, i64 b, i64 m, i64 x, i64 p,
            i64 q, int line, const std::string& want) {
    std::string got = run(eqs, a, b, m, x, p, q);
    if (got != want) {
        std::cout << "FAIL line " << line << "\n--- got ---\n" << got
                  << "--- want ---\n" << want << std::endl;
        std::exit(1);
    }
}

void runTests() {
    expect({{2, 3}, {3, 5}, {2, 7}}, 1, 1, 1, 0, 1, 1, 1, "23 105\n23 105\n1 1 0\n0 1\n1\n");
    expect({{2, 4}, {4, 6}}, 1, 1, 1, 0, 1, 1, 2, "10 12\n-1 -1\n1 1 0\n0 1\n1\n");
    expect({{1, 4}, {2, 6}}, 1, 1, 1, 0, 1, 1, 3, "-1 -1\n-1 -1\n1 1 0\n0 1\n1\n");
    expect({{0, 1}}, 6, 4, 10, 0, 1, 1, 4, "0 1\n0 1\n2 2 4\n0 1\n1\n");
    expect({{0, 1}}, 6, 5, 10, 0, 1, 1, 5, "0 1\n0 1\n2 0 -1\n0 1\n1\n");
    expect({{0, 1}}, 1, 1, 1, 123456789, 97, 1000, 6, "0 1\n0 1\n1 1 0\n72789 97000\n1\n");
    expect({}, 1, 1, 1, 0, 1, 1, 7, "0 1\n0 1\n1 1 0\n0 1\n1\n");

    // ---- 數論基礎 ----
    {
        Exg e = exgcd(0, 0);
        if (e.g != 0 || e.s != 1 || e.t != 0) {
            std::cout << "FAIL exgcd(0,0)" << std::endl;
            std::exit(1);
        }
        e = exgcd(35, 15);
        if (e.g != 5 || e.s != 1 || e.t != -2) {
            std::cout << "FAIL exgcd(35,15)" << std::endl;
            std::exit(1);
        }
        if (modInverse(3, 11) != 4 || modInverse(3, 12) != -1 || modInverse(0, 1) != 0) {
            std::cout << "FAIL modInverse" << std::endl;
            std::exit(1);
        }
        if (gcd(0, 0) != 0) {
            std::cout << "FAIL gcd(0,0)" << std::endl;
            std::exit(1);
        }
    }

    // ---- 隨機對拍（只驗證 ok 旗標） ----
    {
        Lcg rng(20261019ULL);
        for (int t = 0; t < 600; ++t) {
            int n = static_cast<int>(rng.below(9));
            std::vector<std::pair<i64, i64>> eqs;
            for (int i = 0; i < n; ++i) {
                eqs.push_back({rng.below(51), rng.below(30) + 1});
            }
            i64 a = rng.below(201);
            i64 b = rng.below(201);
            i64 m = rng.below(300) + 1;
            i64 x = rng.below(1000001);
            i64 p = rng.below(200) + 1;
            i64 q = rng.below(200) + 1;
            std::string out = run(eqs, a, b, m, x, p, q);
            size_t pos = out.rfind('\n', out.size() - 2);
            std::string last = (pos == std::string::npos) ? out
                                                          : out.substr(pos + 1, out.size() - 2 - pos);
            if (last != "1") {
                std::cout << "FAIL random t=" << t << "\n" << out << std::endl;
                std::exit(1);
            }
        }
    }

    // ---- 故意構造「一定有解」的系統 ----
    {
        Lcg rng(31337ULL);
        for (int t = 0; t < 300; ++t) {
            int n = static_cast<int>(rng.below(8)) + 1;
            std::vector<i64> mods;
            for (int i = 0; i < n; ++i) mods.push_back(rng.below(40) + 1);
            i64 L = lcmAll(mods);
            if (L < 0 || L > LIMIT) continue;
            i64 trueX = rng.below(L);
            std::vector<std::pair<i64, i64>> eqs;
            for (i64 mm : mods) eqs.push_back({trueX % mm, mm});
            std::string out = run(eqs, 1, 1, 1, trueX, mods[0], mods.back());
            std::ostringstream want;
            want << (trueX % L) << ' ' << L << '\n';
            if (out.compare(0, want.str().size(), want.str()) != 0) {
                std::cout << "FAIL constructed t=" << t << "\n" << out << " want " << want.str()
                          << std::endl;
                std::exit(1);
            }
        }
    }
}

// ---------------------------------------------------------------- IO

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

i64 clampLL(const std::string& t, i64 lo, i64 hi, i64 def) {
    if (!isInt(t)) return def;
    errno = 0;
    char* end = nullptr;
    i64 v = std::strtoll(t.c_str(), &end, 10);
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
    size_t idx = 0;
    auto nextInt = [&](i64 def) -> i64 {
        if (idx >= toks.size()) return def;
        std::string s = toks[idx++];
        return clampLL(s, -1000000000000000000LL, 1000000000000000000LL, def);
    };

    i64 n = nextInt(0);
    if (n < 0) n = 0;
    if (n > NMAX) n = NMAX;
    std::vector<std::pair<i64, i64>> eqs;
    for (i64 i = 0; i < n; ++i) {
        i64 a = nextInt(0);
        if (a < 0) a = 0;
        if (a > AMAX) a = AMAX;
        i64 m = nextInt(1);
        if (m < 1) m = 1;
        if (m > MMAX) m = MMAX;
        eqs.push_back({a, m});
    }

    i64 a = nextInt(1);
    if (a < 0) a = 0;
    if (a > AMAX) a = AMAX;
    i64 b = nextInt(0);
    if (b < 0) b = 0;
    if (b > AMAX) b = AMAX;
    i64 m = nextInt(1);
    if (m < 1) m = 1;
    if (m > MMAX) m = MMAX;

    i64 x = nextInt(0);
    if (x < 0) x = 0;
    if (x > XMAX) x = XMAX;
    i64 p = nextInt(1);
    if (p < 1) p = 1;
    if (p > PQMAX) p = PQMAX;
    i64 q = nextInt(1);
    if (q < 1) q = 1;
    if (q > PQMAX) q = PQMAX;

    std::cout << run(eqs, a, b, m, x, p, q);
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
