"""快速冪與模運算：二進制冪 / 矩陣快速冪（Fibonacci）/ 乘法逆元

題意：
    給定若干組詢問，要求高效計算：
      1. a^b mod m —— 樸素做法要做 b 次乘法，b = 10^9 時直接超時；
      2. F(n) mod MOD —— 第 n 個 Fibonacci 數（F(0)=0, F(1)=1），
         用矩陣快速冪把線性遞推從 O(n) 降到 O(log n)；
      3. x 在模 p 下的乘法逆元，即滿足 x * inv ≡ 1 (mod p) 的數。

思路：
    **二進制冪（快速冪）** 的核心是把指數 b 看成二進位：
        a^b = a^(2^0 * b0) * a^(2^1 * b1) * ... ，b_k ∈ {0, 1}
    於是一邊讓底數反覆平方（a → a^2 → a^4 → ...），
    一邊在 b 的當前位爲 1 時把底數乘進答案，共 O(log b) 次乘法。
    這是所有「冪 / 遞推加速」題的地基。

    **矩陣快速冪**：F(n) 的一階線性遞推可以寫成向量乘矩陣
        [F(n), F(n-1)]^T = [[1,1],[1,0]] · [F(n-1), F(n-2)]^T
    所以 [F(n), F(n-1)]^T = M^(n-1) · [F(1), F(0)]^T，
    取結果矩陣的 [0][1] 位置就是 F(n)（本實作直接用 M^n，M^n = [[F(n+1),F(n)],[F(n),F(n-1)]]）。
    任何 k 階線性遞推（Tribonacci、a_n = 3a_{n-1} - 2a_{n-3} ...）都能照搬這個寫法。

    **乘法逆元**：擴展歐幾裡得演算法解 x·a + y·m = gcd(a, m)，
    當 gcd = 1 時 x 就是 a 的逆元（適用於任意互質的 a, m）；
    若 m 是質數，還能用費馬小定理 a^(m-2) ≡ a^(-1) (mod m) 一行搞定。
    本實作兩個版本都寫了，測試裏互相對拍（費馬版僅在模數爲質數時可用）。

輸入格式（stdin，數字按空白分隔即可）：
    n q k
    a1 b1 m1
    a2 b2 m2
    …… （共 q 行，每行一次快速冪詢問）
    p
    x1 x2 ... xk         （要求逆元的 k 個數，模數統一爲 p）
輸出格式（stdout）：
    第 1 行：F(n) mod 1000000007
    第 2 … q+1 行：每組詢問的 a^b mod m（m = 1 時結果必爲 0）
    最後 1 行：k 個逆元，空格分隔；逆元不存在時輸出 -1（k = 0 時輸出空行）
無 stdin 輸入時運行內置斷言測試並輸出 `all tests passed`。
"""

import random
import sys
from typing import List, Tuple

MOD = 1_000_000_007

Matrix = List[List[int]]


def mat_mul(A: Matrix, B: Matrix, mod: int = MOD) -> Matrix:
    """2×2 矩陣乘法（結果對 mod 取模）。"""
    a00 = (A[0][0] * B[0][0] + A[0][1] * B[1][0]) % mod
    a01 = (A[0][0] * B[0][1] + A[0][1] * B[1][1]) % mod
    a10 = (A[1][0] * B[0][0] + A[1][1] * B[1][0]) % mod
    a11 = (A[1][0] * B[0][1] + A[1][1] * B[1][1]) % mod
    return [[a00, a01], [a10, a11]]


def mat_pow(A: Matrix, e: int, mod: int = MOD) -> Matrix:
    """矩陣快速冪：A^e，O(log e) 次矩陣乘法。"""
    # 單位矩陣當初始答案
    res = [[1 % mod, 0], [0, 1 % mod]]
    base = [[A[0][0] % mod, A[0][1] % mod], [A[1][0] % mod, A[1][1] % mod]]
    while e > 0:
        if e & 1:
            res = mat_mul(res, base, mod)
        base = mat_mul(base, base, mod)
        e >>= 1
    return res


def fib(n: int, mod: int = MOD) -> int:
    """第 n 個 Fibonacci 數（F(0)=0, F(1)=1），矩陣快速冪 O(log n)。"""
    if n < 0:
        raise ValueError("n must be non-negative")
    if n == 0:
        return 0
    # M^n = [[F(n+1), F(n)], [F(n), F(n-1)]]，取 [0][1]
    return mat_pow([[1, 1], [1, 0]], n, mod)[0][1]


def mod_pow(a: int, b: int, mod: int) -> int:
    """二進制冪：a^b mod mod，O(log b)。要求 b >= 0、mod >= 1。"""
    if mod <= 0:
        raise ValueError("mod must be positive")
    if b < 0:
        raise ValueError("exponent must be non-negative")
    if mod == 1:
        return 0
    res = 1 % mod
    a %= mod
    while b > 0:
        if b & 1:
            res = res * a % mod
        a = a * a % mod
        b >>= 1
    return res


def slow_pow(a: int, b: int, mod: int) -> int:
    """樸素冪：O(b)，只在測試裏當基準對拍用。"""
    res = 1 % mod
    a %= mod
    for _ in range(b):
        res = res * a % mod
    return res


def ext_gcd(a: int, b: int) -> Tuple[int, int, int]:
    """擴展歐幾裡得：回傳 (g, x, y) 使 x*a + y*b = g = gcd(a, b)。a, b 需非負。"""
    old_r, r = a, b
    old_s, s = 1, 0
    old_t, t = 0, 1
    while r != 0:
        q = old_r // r
        old_r, r = r, old_r - q * r
        old_s, s = s, old_s - q * s
        old_t, t = t, old_t - q * t
    return old_r, old_s, old_t


def mod_inverse(a: int, mod: int) -> int:
    """擴展歐幾裡得求逆元：a * inv ≡ 1 (mod mod)。不存在（gcd != 1）時回傳 -1。"""
    if mod <= 1:
        return -1
    a %= mod
    if a == 0:
        return -1
    g, x, _ = ext_gcd(a, mod)
    if g != 1:
        return -1
    return x % mod


def mod_inverse_fermat(a: int, p: int) -> int:
    """費馬小定理求逆元：a^(p-2) mod p，僅當 p 爲質數且 a 不被 p 整除時正確。"""
    if p <= 1 or a % p == 0:
        return -1
    return mod_pow(a, p - 2, p)


def is_prime(x: int) -> bool:
    """試除法判質數，只用於測試裏挑質數模數。"""
    if x < 2:
        return False
    d = 2
    while d * d <= x:
        if x % d == 0:
            return False
        d += 1
    return True


# ---------------------------------------------------------------- 輸入輸出


def run_io(data: str) -> None:
    """按題目格式解析 stdin 並輸出結果。"""
    it = iter(data.split())
    nxt = lambda: int(next(it))  # noqa: E731

    def take(default: int = 0) -> int:
        try:
            return nxt()
        except StopIteration:
            return default

    n = take()
    q = take()
    k = take()

    out: List[str] = [str(fib(n))]
    for _ in range(q):
        a = take()
        b = take()
        m = take(1)
        out.append(str(mod_pow(a, b, m)))

    p = take(1)
    invs = [str(mod_inverse(take(), p)) for _ in range(k)]
    out.append(" ".join(invs))

    sys.stdout.write("\n".join(out) + "\n")


# ---------------------------------------------------------------- 內置測試


def run_tests() -> None:
    # ---- 快速冪 ----
    assert mod_pow(2, 10, 1000) == 24          # 2^10 = 1024
    assert mod_pow(2, 0, 7) == 1               # 任何數的 0 次方都是 1
    assert mod_pow(0, 0, 7) == 1
    assert mod_pow(0, 5, 7) == 0
    assert mod_pow(12345, 1, 97) == 12345 % 97
    assert mod_pow(5, 3, 1) == 0               # 模 1 恆爲 0
    assert mod_pow(-2, 3, 7) == ((-2) ** 3) % 7  # 負底數先取模
    assert mod_pow(3, 100, 1_000_000_007) == pow(3, 100, 1_000_000_007)

    random.seed(20260930)
    for _ in range(3000):
        a = random.randint(-500, 500)
        b = random.randint(0, 200)
        m = random.randint(1, 500)
        assert mod_pow(a, b, m) == pow(a, b, m)      # 與內建 pow 對拍
        if b <= 60:
            assert mod_pow(a, b, m) == slow_pow(a, b, m)

    # 超大指數（樸素做法不可能算完，只有快速冪能秒出）
    big = mod_pow(2, 10**9, 1_000_000_007)
    assert big == pow(2, 10**9, 1_000_000_007)

    # ---- Fibonacci ----
    seq = [0, 1]
    while len(seq) <= 90:
        seq.append(seq[-1] + seq[-2])
    for i, want in enumerate(seq):
        assert fib(i, 10**30 + 1) == want      # 模數夠大時與迭代版完全一致
    assert fib(0) == 0 and fib(1) == 1 and fib(2) == 1 and fib(10) == 55
    assert fib(50, 1_000_000_007) == 586268941   # F(50) = 12586269025，取模後 586268941
    assert fib(90, 1_000_000_007) == seq[90] % 1_000_000_007
    assert fib(200, 1000) == (fib(199, 1000) + fib(198, 1000)) % 1000  # 遞推自洽
    # 大 n 也要和「矩陣連乘 n-1 次」一致（用小模數避免溢出）
    assert fib(40, 1000) == seq[40] % 1000

    # ---- 逆元 ----
    assert mod_inverse(3, 7) == 5               # 3 * 5 = 15 ≡ 1 (mod 7)
    assert mod_inverse(1, 100) == 1
    assert mod_inverse(2, 4) == -1              # gcd(2,4) = 2，無逆元
    assert mod_inverse(0, 7) == -1
    assert mod_inverse(5, 1) == -1
    assert mod_inverse_fermat(3, 7) == 5

    for _ in range(2000):
        p = random.choice([7, 11, 13, 97, 101, 1009, 10007])
        assert is_prime(p)
        x = random.randint(0, p - 1)
        inv = mod_inverse(x, p)
        if x % p == 0:
            assert inv == -1
            assert mod_inverse_fermat(x, p) == -1
        else:
            assert inv != -1
            assert (x * inv) % p == 1                    # 定義驗證
            assert mod_inverse_fermat(x, p) == inv       # 費馬版與擴歐版一致
    # 合數模數下只有互質的數才有逆元
    for _ in range(2000):
        m = random.randint(2, 60)
        x = random.randint(0, m - 1)
        inv = mod_inverse(x, m)
        g, _, _ = ext_gcd(x, m)
        if g != 1:
            assert inv == -1
        else:
            assert (x * inv) % m == 1

    print("all tests passed")


if __name__ == "__main__":
    raw = sys.stdin.read() if not sys.stdin.isatty() else ""
    if raw.strip():
        run_io(raw)
    else:
        run_tests()
