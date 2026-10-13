// 高斯消元：GF(p) 上的線性方程組 + GF(2) 上的 XOR 線性基
//
// 題意：
//     給 n×m 矩陣 A、右端向量 b、質數模數 p、一組非負整數 X 與整數 k，輸出
//     A 在模 p 下的秩與行列式、方程組 A·x ≡ b (mod p) 的解分類與一組解、
//     A 是否可逆、X 的 XOR 線性基大小與最大異或值、以及第 k 小可組出的異或值。
//
// 思路：
//     ### 高斯消元在做什麼
//     三種不改動解集合的列操作：① 換列（行列式變號）② 一列乘非零常數 ③ 一列加另一列的倍數。
//     用它們把矩陣化成簡化列階梯型（RREF），每個主元所在行只有該主元非零，解就能直接讀出。
//
//     ### 為什麼用模質數域 GF(p)
//     實數上做消元有浮點誤差，兩個語言不可能輸出一致。改到模質數 p 的有限域後
//     所有運算都是精確整數：除法變成乘逆元，p 為質數時 a⁻¹ ≡ a^(p-2) (mod p)（費馬小定理），
//     每個非零元素都可逆，主元搜尋只要找「非零」即可。
//
//     ### 行列式
//     記列交換次數 s 與每個主元 pᵢ（正規化之前的值）。三種操作對行列式的影響是
//     ×(-1)、×(1/pᵢ)、×1，最後得到 RREF（可逆時為 I，行列式 1），故
//         det(A) = (-1)^s · Π pᵢ  (mod p)
//     若秩 < n（方陣不滿秩）則行列式為 0。
//
//     ### 解的分類
//     增廣矩陣一起消元，比較 rank(A) 與 rank([A|b])：
//       rank(A) < rank([A|b]) → 出現「0 = 非零」的列 → 無解
//       rank(A) = m          → 唯一解
//       否則                 → 無窮多解，自由變數 m - rank 個，GF(p) 上共 p^(m-rank) 組
//     本題自由變數一律補 0，保證兩個語言輸出一致。
//
//     ### XOR 線性基 = GF(2) 上的高斯消元
//     每個數字看成 63 維 GF(2) 向量，異或即加法。basis[i] 是最高位為 i 的基底向量，
//     插入時從高位往下消，與「拿主元消掉其他列」完全同構。
//       最大異或值：從高位貪心，res = max(res, res ^ basis[i])
//       第 k 小    ：先把基底重整（basis[i] 的低位不含其它基底的最高位），
//                   再把非零基底由小到大排成 d，把 k-1 的位元當選擇向量即可。
//
// 複雜度：RREF O(n·m·min(n,m))（n,m ≤ 40）；線性基插入 O(63)；空間 O(n·m)。
//
// 輸入（stdin，空白分隔）：
//     n m p
//     a_11 ... a_1m        （n 列）
//     b_1 ... b_n
//     c                    （線性基數字個數）
//     x_1 ... x_c
//     k
//     非數字 token 視為輸入結束；讀不滿的那一列/那一組直接丟棄；n,m,c ≤ 40；
//     p 提升到「不小於輸入值的最小質數」（p < 2 時取 2）。
// 輸出（stdout，6 行）：
//     1: rank det          2: solkind free_cnt   3: 一組解（無解為空行）
//     4: invertible        5: basis_size max_xor 6: kth_xor（越界為 -1）
// 無 stdin 時跑內置斷言測試並輸出 all tests passed。

#include <cassert>
#include <cerrno>
#include <cstdlib>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <map>
#include <random>
#include <set>
#include <string>
#include <vector>

#ifdef _WIN32
#include <fcntl.h>
#include <io.h>
#endif

namespace {

constexpr int MAXNM = 40;
constexpr long long P_MAX = 2147483647LL;      // 2^31 - 1，本身是質數
constexpr int MAXB = 62;                       // 線性基最高位元
constexpr unsigned long long MASK = (1ULL << 63) - 1;

using Matrix = std::vector<std::vector<long long>>;

// ---------------------------------------------------------------- 質數與模運算
bool isPrime(long long v) {
    if (v < 2) return false;
    if (v % 2 == 0) return v == 2;
    for (long long d = 3; d * d <= v; d += 2)
        if (v % d == 0) return false;
    return true;
}

long long toPrime(long long v) {
    if (v < 2) return 2;
    while (!isPrime(v)) ++v;
    return v;
}

long long norm(long long v, long long p) {
    long long r = v % p;
    return r < 0 ? r + p : r;
}

long long modPow(long long a, long long e, long long p) {
    long long res = 1 % p;
    a = norm(a, p);
    while (e > 0) {
        if (e & 1) res = (res * a) % p;
        a = (a * a) % p;
        e >>= 1;
    }
    return res;
}

long long modInv(long long a, long long p) { return modPow(a, p - 2, p); }

// ---------------------------------------------------------------- 高斯消元
struct RrefResult {
    int rank;
    long long det;
    Matrix aug;                                  // n × (m+1)
};

RrefResult rref(const Matrix& a, const std::vector<long long>& b, int n, int m, long long p) {
    Matrix aug(n, std::vector<long long>(m + 1, 0));
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < m; ++j) aug[i][j] = norm(a[i][j], p);
        aug[i][m] = norm(b[i], p);
    }
    int rank = 0;
    long long det = 1 % p;
    long long sign = 1;
    int col = 0;
    while (col < m && rank < n) {
        int piv = -1;
        for (int i = rank; i < n; ++i)
            if (aug[i][col] % p != 0) { piv = i; break; }
        if (piv < 0) { ++col; continue; }
        if (piv != rank) {
            std::swap(aug[piv], aug[rank]);
            sign = -sign;
        }
        long long pivot = norm(aug[rank][col], p);
        det = (det * pivot) % p;
        long long inv = modInv(pivot, p);
        for (int j = 0; j <= m; ++j) aug[rank][j] = (aug[rank][j] * inv) % p;
        for (int i = 0; i < n; ++i) {
            if (i == rank || aug[i][col] % p == 0) continue;
            long long f = norm(aug[i][col], p);
            for (int j = 0; j <= m; ++j)
                aug[i][j] = norm(aug[i][j] - f * aug[rank][j] % p, p);
        }
        ++rank;
        ++col;
    }
    long long d = 0;
    if (n == m) d = (rank == n) ? norm(sign * det, p) : 0;
    return {rank, d, aug};
}

struct SolveResult {
    int rank;
    long long det;
    int kind;                                    // 0 無解 / 1 唯一 / 2 無窮多
    int freeCnt;
    std::vector<long long> x;
};

SolveResult solveLinear(const Matrix& a, const std::vector<long long>& b, int n, int m,
                        long long p) {
    RrefResult rr = rref(a, b, n, m, p);
    for (int i = rr.rank; i < n; ++i)
        if (rr.aug[i][m] % p != 0) return {rr.rank, rr.det, 0, m - rr.rank, {}};

    std::vector<long long> x(m, 0);
    int col = 0, r = 0;
    while (col < m && r < rr.rank) {
        while (col < m && rr.aug[r][col] % p == 0) ++col;
        if (col >= m) break;
        x[col] = norm(rr.aug[r][m], p);
        ++r;
        ++col;
    }
    if (rr.rank == m) return {rr.rank, rr.det, 1, 0, x};
    return {rr.rank, rr.det, 2, m - rr.rank, x};
}

bool matrixInverse(const Matrix& a, int n, long long p, Matrix& out) {
    Matrix aug(n, std::vector<long long>(2 * n, 0));
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) aug[i][j] = norm(a[i][j], p);
        aug[i][n + i] = 1 % p;
    }
    for (int col = 0; col < n; ++col) {
        int piv = -1;
        for (int i = col; i < n; ++i)
            if (aug[i][col] % p != 0) { piv = i; break; }
        if (piv < 0) return false;
        std::swap(aug[col], aug[piv]);
        long long inv = modInv(norm(aug[col][col], p), p);
        for (int j = 0; j < 2 * n; ++j) aug[col][j] = (aug[col][j] * inv) % p;
        for (int i = 0; i < n; ++i) {
            if (i == col || aug[i][col] % p == 0) continue;
            long long f = norm(aug[i][col], p);
            for (int j = 0; j < 2 * n; ++j) aug[i][j] = norm(aug[i][j] - f * aug[col][j] % p, p);
        }
    }
    out.assign(n, std::vector<long long>(n, 0));
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j) out[i][j] = norm(aug[i][n + j], p);
    return true;
}

// ---------------------------------------------------------------- XOR 線性基
class LinearBasis {
public:
    LinearBasis() : basis_(MAXB + 1, 0) {}

    bool insert(unsigned long long x) {
        x &= MASK;
        for (int i = MAXB; i >= 0; --i) {
            if (((x >> i) & 1ULL) == 0) continue;
            if (basis_[i] == 0) {
                basis_[i] = x;
                return true;
            }
            x ^= basis_[i];
        }
        return false;
    }

    int size() const {
        int s = 0;
        for (unsigned long long v : basis_) if (v != 0) ++s;
        return s;
    }

    unsigned long long maxXor() const {
        unsigned long long res = 0;
        for (int i = MAXB; i >= 0; --i)
            if (basis_[i] != 0 && (res ^ basis_[i]) > res) res ^= basis_[i];
        return res;
    }

    std::vector<unsigned long long> reduceBasis() const {
        std::vector<unsigned long long> b = basis_;
        for (int i = 0; i <= MAXB; ++i) {
            if (b[i] == 0) continue;
            for (int j = i - 1; j >= 0; --j)
                if (b[j] != 0 && ((b[i] >> j) & 1ULL)) b[i] ^= b[j];
        }
        std::vector<unsigned long long> d;
        for (unsigned long long v : b) if (v != 0) d.push_back(v);
        return d;
    }

    long long kthSmallest(long long k) const {
        std::vector<unsigned long long> d = reduceBasis();
        if (d.size() >= 63) return -1;                 // 2^63 組以上一定超出 long long 範圍
        unsigned long long total = 1ULL << d.size();
        if (k < 1 || static_cast<unsigned long long>(k) > total) return -1;
        unsigned long long kk = static_cast<unsigned long long>(k) - 1;
        unsigned long long ans = 0;
        for (std::size_t i = 0; i < d.size(); ++i)
            if ((kk >> i) & 1ULL) ans ^= d[i];
        return static_cast<long long>(ans);
    }

private:
    std::vector<unsigned long long> basis_;
};

// ---------------------------------------------------------------- 暴力對拍
int bfRank(const Matrix& a, int n, int m, long long p) {
    Matrix mat(n, std::vector<long long>(m, 0));
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < m; ++j) mat[i][j] = norm(a[i][j], p);
    int row = 0;
    for (int col = m - 1; col >= 0 && row < n; --col) {
        int piv = -1;
        for (int i = row; i < n; ++i)
            if (mat[i][col] % p != 0) { piv = i; break; }
        if (piv < 0) continue;
        std::swap(mat[row], mat[piv]);
        long long inv = modInv(norm(mat[row][col], p), p);
        for (int j = 0; j < m; ++j) mat[row][j] = (mat[row][j] * inv) % p;
        for (int i = 0; i < n; ++i) {
            if (i == row || mat[i][col] % p == 0) continue;
            long long f = norm(mat[i][col], p);
            for (int j = 0; j < m; ++j) mat[i][j] = norm(mat[i][j] - f * mat[row][j] % p, p);
        }
        ++row;
    }
    return row;
}

long long bfDet(const Matrix& a, int n, long long p) {
    std::vector<int> perm(n);
    for (int i = 0; i < n; ++i) perm[i] = i;
    long long total = 0;
    do {
        int sign = 1;
        std::vector<char> seen(n, 0);
        for (int i = 0; i < n; ++i)
            if (!seen[i]) {
                int j = i, len = 0;
                while (!seen[j]) { seen[j] = 1; j = perm[j]; ++len; }
                if (len % 2 == 0) sign = -sign;
            }
        long long prod = 1 % p;
        for (int i = 0; i < n; ++i) prod = (prod * norm(a[i][perm[i]], p)) % p;
        total = norm(total + sign * prod, p);
    } while (std::next_permutation(perm.begin(), perm.end()));
    return total;
}

void bfSolutionsEnum(int idx, int m, long long p, const Matrix& a,
                     const std::vector<long long>& b, int n,
                     std::vector<long long>& cur, long long& cnt,
                     std::vector<long long>& first) {
    if (idx == m) {
        for (int i = 0; i < n; ++i) {
            long long s = 0;
            for (int j = 0; j < m; ++j) s = (s + norm(a[i][j], p) * cur[j]) % p;
            if (norm(s, p) != norm(b[i], p)) return;
        }
        ++cnt;
        if (first.empty()) first = cur;
        return;
    }
    for (long long v = 0; v < p; ++v) {
        cur[idx] = v;
        bfSolutionsEnum(idx + 1, m, p, a, b, n, cur, cnt, first);
    }
    cur[idx] = 0;
}

bool smallEnough(int m, long long p, long long limit) {
    long long acc = 1;
    for (int i = 0; i < m; ++i) {
        if (acc > limit / p) return false;
        acc *= p;
    }
    return true;
}

std::vector<unsigned long long> bfXorValues(const std::vector<unsigned long long>& xs) {
    std::set<unsigned long long> vals;
    vals.insert(0);
    for (unsigned long long x : xs) {
        std::vector<unsigned long long> add;
        for (unsigned long long v : vals) add.push_back(v ^ x);
        for (unsigned long long v : add) vals.insert(v);
    }
    return std::vector<unsigned long long>(vals.begin(), vals.end());
}

// ---------------------------------------------------------------- 輸出
std::string run(const Matrix& a, const std::vector<long long>& b, int n, int m, long long p,
                const std::vector<unsigned long long>& xs, long long k) {
    SolveResult sr = solveLinear(a, b, n, m, p);
    int invertible = (n == m && sr.rank == n) ? 1 : 0;

    LinearBasis lb;
    for (unsigned long long x : xs) lb.insert(x);

    std::string out;
    out += std::to_string(sr.rank) + " " + std::to_string(sr.det) + "\n";
    out += std::to_string(sr.kind) + " " + std::to_string(sr.freeCnt) + "\n";
    if (sr.kind != 0) {
        for (int j = 0; j < m; ++j) {
            if (j) out += " ";
            out += std::to_string(sr.x[j]);
        }
    }
    out += "\n";
    out += std::to_string(invertible) + "\n";
    out += std::to_string(lb.size()) + " " + std::to_string(lb.maxXor()) + "\n";
    out += std::to_string(lb.kthSmallest(k));
    return out;
}

// ---------------------------------------------------------------- 測試
void checkSystem(const Matrix& a, const std::vector<long long>& b, int n, int m, long long p) {
    SolveResult sr = solveLinear(a, b, n, m, p);
    assert(sr.rank == bfRank(a, n, m, p));
    if (n == m && n <= 6) assert(sr.det == bfDet(a, n, p));

    if (sr.kind != 0) {
        assert(static_cast<int>(sr.x.size()) == m);
        for (int i = 0; i < n; ++i) {
            long long s = 0;
            for (int j = 0; j < m; ++j) s = (s + norm(a[i][j], p) * sr.x[j]) % p;
            assert(norm(s, p) == norm(b[i], p));
        }
    }
    if (sr.kind == 1) assert(sr.rank == m && sr.freeCnt == 0);
    if (sr.kind == 2) assert(sr.freeCnt == m - sr.rank && sr.freeCnt > 0);
    if (sr.kind == 0) assert(sr.x.empty());

    if (n == m) {
        Matrix inv;
        if (sr.rank == n) {
            assert(matrixInverse(a, n, p, inv));
            for (int i = 0; i < n; ++i)
                for (int j = 0; j < n; ++j) {
                    long long s = 0;
                    for (int t = 0; t < n; ++t) s = (s + norm(a[i][t], p) * inv[t][j]) % p;
                    assert(norm(s, p) == (i == j ? 1 % p : 0));
                }
        } else {
            assert(n == 0 || !matrixInverse(a, n, p, inv));
        }
    }

    if (smallEnough(m, p, 3000)) {
        long long cnt = 0;
        std::vector<long long> first;
        std::vector<long long> cur(m, 0);
        bfSolutionsEnum(0, m, p, a, b, n, cur, cnt, first);
        if (sr.kind == 0) {
            assert(cnt == 0);
        } else if (sr.kind == 1) {
            assert(cnt == 1);
            for (int j = 0; j < m; ++j) assert(sr.x[j] == first[j]);
        } else {
            long long expect = 1;
            for (int i = 0; i < m - sr.rank; ++i) expect *= p;
            assert(cnt == expect);
            assert(!first.empty());
        }
    }
}

void checkBasis(const std::vector<unsigned long long>& xs) {
    LinearBasis lb;
    for (unsigned long long x : xs) lb.insert(x);
    std::vector<unsigned long long> vals = bfXorValues(xs);
    assert(static_cast<long long>(vals.size()) == (1LL << lb.size()));
    assert(lb.maxXor() == (vals.empty() ? 0ULL : vals.back()));
    for (std::size_t i = 0; i < vals.size(); ++i) {
        long long k = static_cast<long long>(i) + 1;
        assert(lb.kthSmallest(k) == static_cast<long long>(vals[i]));
    }
    assert(lb.kthSmallest(0) == -1);
    assert(lb.kthSmallest(static_cast<long long>(vals.size()) + 1) == -1);
}

void runTests() {
    std::mt19937 rng(20261013);

    // ---- 固定用例（與 Python 版逐字節一致）----
    assert(toPrime(0) == 2 && toPrime(1) == 2 && toPrime(2) == 2);
    assert(toPrime(9) == 11 && toPrime(100) == 101);
    assert(toPrime(1000000007) == 1000000007);

    assert(run({{1, 2}, {3, 4}}, {5, 11}, 2, 2, 7, {3, 5}, 3) == "2 5\n1 0\n1 2\n1\n2 6\n5");
    assert(run({{1, 1}, {2, 2}}, {1, 2}, 2, 2, 7, {7, 3}, 2) == "1 0\n2 1\n1 0\n0\n2 7\n3");
    assert(run({{1, 1}, {2, 2}}, {1, 3}, 2, 2, 7, {1}, 1) == "1 0\n0 1\n\n0\n1 1\n0");

    checkSystem({{1, 2}, {3, 4}}, {5, 11}, 2, 2, 7);
    checkSystem({{1, 1}, {2, 2}}, {1, 2}, 2, 2, 7);
    checkSystem({{1, 1}, {2, 2}}, {1, 3}, 2, 2, 7);

    checkBasis({3, 5});
    checkBasis({7, 3});
    checkBasis({1, 2, 3});
    checkBasis({});

    // ---- 隨機對拍：小矩陣 + 全枚舉 ----
    const long long smallP[5] = {2, 3, 5, 7, 11};
    for (int iter = 0; iter < 300; ++iter) {
        long long p = smallP[rng() % 5];
        int n = static_cast<int>(rng() % 5);
        int m = static_cast<int>(rng() % 4);
        Matrix a(static_cast<std::size_t>(n), std::vector<long long>(m));
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < m; ++j)
                a[i][j] = static_cast<long long>(rng() % 13) - 6;
        std::vector<long long> b(n);
        for (int i = 0; i < n; ++i) b[i] = static_cast<long long>(rng() % 13) - 6;
        checkSystem(a, b, n, m, p);
    }

    // ---- 隨機對拍：方陣行列式與萊布尼茲公式 ----
    const long long bigP[4] = {7, 13, 101, 1000000007};
    for (int iter = 0; iter < 150; ++iter) {
        long long p = bigP[rng() % 4];
        int n = 1 + static_cast<int>(rng() % 5);
        Matrix a(static_cast<std::size_t>(n), std::vector<long long>(n));
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j)
                a[i][j] = static_cast<long long>(rng() % 19) - 9;
        std::vector<long long> b(n);
        for (int i = 0; i < n; ++i) b[i] = static_cast<long long>(rng() % 19) - 9;
        SolveResult sr = solveLinear(a, b, n, n, p);
        assert(sr.det == bfDet(a, n, p));
        checkSystem(a, b, n, n, p);
    }

    // ---- 隨機對拍：中矩陣 ----
    const long long midP[3] = {101, 1009, 1000000007};
    for (int iter = 0; iter < 120; ++iter) {
        long long p = midP[rng() % 3];
        int n = 1 + static_cast<int>(rng() % 12);
        int m = 1 + static_cast<int>(rng() % 12);
        Matrix a(static_cast<std::size_t>(n), std::vector<long long>(m));
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < m; ++j)
                a[i][j] = static_cast<long long>(rng() % 41) - 20;
        std::vector<long long> b(n);
        for (int i = 0; i < n; ++i) b[i] = static_cast<long long>(rng() % 41) - 20;
        SolveResult sr = solveLinear(a, b, n, m, p);
        assert(sr.rank == bfRank(a, n, m, p));
        if (sr.kind != 0)
            for (int i = 0; i < n; ++i) {
                long long s = 0;
                for (int j = 0; j < m; ++j) s = (s + norm(a[i][j], p) * sr.x[j]) % p;
                assert(norm(s, p) == norm(b[i], p));
            }
        if (n == m && sr.rank == n) {
            Matrix inv;
            assert(matrixInverse(a, n, p, inv));
            for (int i = 0; i < n; ++i)
                for (int j = 0; j < n; ++j) {
                    long long s = 0;
                    for (int t = 0; t < n; ++t) s = (s + norm(a[i][t], p) * inv[t][j]) % p;
                    assert(norm(s, p) == (i == j ? 1 % p : 0));
                }
        }
    }

    // ---- 隨機對拍：線性基 ----
    for (int iter = 0; iter < 200; ++iter) {
        int c = static_cast<int>(rng() % 13);
        std::vector<unsigned long long> xs;
        for (int i = 0; i < c; ++i) {
            int bits = 1 + static_cast<int>(rng() % 20);
            unsigned long long v = 0;
            for (int j = 0; j < bits; ++j) v = (v << 1) | (rng() & 1ULL);
            xs.push_back(v);
        }
        checkBasis(xs);
    }

    // ---- 邊界 ----
    assert(run({}, {}, 0, 0, 7, {}, 1) == "0 1\n1 0\n\n1\n0 0\n0");
    assert(run({}, {}, 0, 3, 7, {}, 1) == "0 0\n2 3\n0 0 0\n0\n0 0\n0");
    assert(run({{0, 0}}, {5}, 1, 2, 7, {}, 1) == "0 0\n0 2\n\n0\n0 0\n0");

    LinearBasis lb;
    assert(lb.size() == 0 && lb.maxXor() == 0 && lb.kthSmallest(1) == 0);
    assert(lb.kthSmallest(2) == -1);
    LinearBasis lb2;
    lb2.insert(MASK);
    assert(lb2.maxXor() == MASK);
    assert(lb2.insert(MASK) == false);
    LinearBasis lb3;
    for (int i = 0; i < 63; ++i) lb3.insert(1ULL << i);
    assert(lb3.size() == 63);
    assert(lb3.maxXor() == MASK);
}

// ---------------------------------------------------------------- 解析
bool isSpaceCh(unsigned char c) {
    return c == ' ' || c == '\t' || c == '\n' || c == '\r' || c == '\v' || c == '\f';
}

std::vector<std::string> splitWs(const std::string& raw) {
    std::vector<std::string> out;
    std::string cur;
    for (unsigned char c : raw) {
        if (isSpaceCh(c)) {
            if (!cur.empty()) { out.push_back(cur); cur.clear(); }
        } else {
            cur.push_back(static_cast<char>(c));
        }
    }
    if (!cur.empty()) out.push_back(cur);
    return out;
}

bool isNumToken(const std::string& tok) {
    if (tok.empty()) return false;
    std::size_t i = (tok[0] == '+' || tok[0] == '-') ? 1 : 0;
    if (i >= tok.size()) return false;
    for (; i < tok.size(); ++i)
        if (tok[i] < '0' || tok[i] > '9') return false;
    return true;
}

long long parseLLClamp(const std::string& tok) {
    errno = 0;
    char* end = nullptr;
    long long v = std::strtoll(tok.c_str(), &end, 10);
    if (end == tok.c_str()) return 0;
    return v;                       // 溢位時 strtoll 已飽和到 LLONG_MAX / LLONG_MIN
}

class Reader {
public:
    explicit Reader(const std::vector<std::string>& tok) : tok_(tok) {}

    bool readInt(long long& out) {
        if (pos_ >= tok_.size() || !isNumToken(tok_[pos_])) return false;
        out = parseLLClamp(tok_[pos_]);
        ++pos_;
        return true;
    }

private:
    const std::vector<std::string>& tok_;
    std::size_t pos_ = 0;
};

void runIo(const std::string& raw) {
    std::vector<std::string> tok = splitWs(raw);
    Reader rd(tok);
    auto clamp = [](long long v, long long lo, long long hi) {
        return v < lo ? lo : (v > hi ? hi : v);
    };

    long long rawN, rawM, rawP;
    int n = rd.readInt(rawN) ? static_cast<int>(clamp(rawN, 0, MAXNM)) : 0;
    int m = rd.readInt(rawM) ? static_cast<int>(clamp(rawM, 0, MAXNM)) : 0;
    long long p = toPrime(rd.readInt(rawP) ? clamp(rawP, 0, P_MAX) : 2);

    Matrix a;
    for (int i = 0; i < n; ++i) {
        std::vector<long long> row;
        bool full = true;
        for (int j = 0; j < m; ++j) {
            long long v;
            if (!rd.readInt(v)) { full = false; break; }
            row.push_back(norm(v, p));
        }
        if (!full) break;                    // 讀不滿 → 整列丟棄
        a.push_back(row);
    }
    n = static_cast<int>(a.size());

    std::vector<long long> b;
    for (int i = 0; i < n; ++i) {
        long long v;
        if (!rd.readInt(v)) break;
        b.push_back(norm(v, p));
    }
    if (static_cast<int>(b.size()) < n) {
        n = static_cast<int>(b.size());
        a.resize(static_cast<std::size_t>(n));
    }

    long long rawC;
    int c = rd.readInt(rawC) ? static_cast<int>(clamp(rawC, 0, MAXNM)) : 0;
    std::vector<unsigned long long> xs;
    for (int i = 0; i < c; ++i) {
        long long v;
        if (!rd.readInt(v)) break;
        xs.push_back(static_cast<unsigned long long>(v) & MASK);
    }

    long long k = 0;
    long long got;
    if (rd.readInt(got)) k = got;

    std::string out = run(a, b, n, m, p, xs, k);
    out += "\n";
    std::cout << out;
}

}  // namespace

int main() {
#ifdef _WIN32
    // Windows 文字模式會把 \n 轉成 \r\n，與 Python 的二進位輸出對不起來
    _setmode(_fileno(stdout), _O_BINARY);
#endif
    std::string raw((std::istreambuf_iterator<char>(std::cin)),
                    std::istreambuf_iterator<char>());
    bool blank = true;
    for (unsigned char c : raw)
        if (!isSpaceCh(c)) { blank = false; break; }
    if (!blank) {
        runIo(raw);
    } else {
        runTests();
        std::cout << "all tests passed" << std::endl;
    }
    return 0;
}
