"""Miller-Rabin 質數測試與 Pollard-Rho 因數分解

題意：
    給定若干個非負整數（0 ≤ n ≤ 10^12），對每一個輸出：

    | 欄位 | 說明 |
    |---|---|
    | `prime` | n 是否為質數（1 / 0） |
    | `lpf`   | 最大的質因數（n < 2 時為 0） |
    | `omega` | 相異質因數個數 ω(n) |
    | `Omega` | 質因數總數（含重複）Ω(n) |
    | `phi`   | 尤拉函數 φ(n)：1..n 中與 n 互質的數的個數 |
    | `sigma` | 因數和 σ(n) |
    | `tau`   | 因數個數 τ(n) |
    | `fac`   | 質因數分解，寫成 `p^e` 用 `*` 串起來（n ≤ 1 時為 `-`） |

    以上每一項都與**獨立的對照算法**互相驗證（試除法分解、枚舉因數、
    暴力數互質個數、∑_{d|n} φ(d) = n 這個恆等式）。

思路：
    ### 費馬小定理為什麼不夠
    若 p 是質數，則 a^(p-1) ≡ 1 (mod p)。反過來「a^(n-1) ≡ 1 就說 n 是質數」
    會被 Carmichael 數（例如 561 = 3·11·17）騙到：它們對所有與 n 互質的 a 都滿足費馬式。

    ### Miller-Rabin：把「平方剩餘」的條件也加進來
    把 n-1 寫成 d·2^s（d 為奇數）。若 n 是質數，則對任意 a 都有
        a^d ≡ 1 (mod n)  或  存在某個 r ∈ [0, s) 使 a^(d·2^r) ≡ -1 (mod n)
    （因為在質數的乘法群裡，x² ≡ 1 只有 x ≡ ±1 兩個解，
      所以「一路平方到 1」的過程中最後一個不是 1 的值必定是 -1。）
    對一組固定的底數 a 測試上述條件，全部通過就**很可能**是質數；
    任一個不通過則 n 一定是合數（單邊錯誤）。
    本題用前 12 個質數當底數 {2,3,5,7,11,13,17,19,23,29,31,37}：
    這組底數對所有 n < 3.3×10^24 都是**確定性**的（已由數學家驗證），
    所以對 n ≤ 10^12 而言，本程式的質數測試沒有任何機率成分。

    ### Pollard-Rho：隨機走一趟函數圖，撞到「同餘」就拿到因數
    這是本倉庫 2026-10-17 `functional-graph` 那一題的直接應用：
    在模 n 的世界裡迭代 x ← (x² + c) mod n，因為值域只有 n 個，
    路徑一定會形成 ρ 形狀（一條尾巴接一個環）。
    若 p | n，則「模 p 的序列」也會撞環；用 Floyd 龜兔法維護
        d = gcd(|x - y|, n)
    一旦 x ≡ y (mod p) 但 x ≢ y (mod n)，就得到 1 < d < n，即 n 的一個非平凡因數。
    期望步數約 O(√p)（p 是 n 最小的質因數），比試除法的 O(√n) 快非常多。
    d == n 表示這條參數 c 失敗，換一個 c 重來即可。

    ### 因數分解 → φ / σ / τ
    分解出 n = Π p_i^e_i 後直接用積性公式：
        φ(n) = Π p_i^(e_i-1) · (p_i - 1)
        σ(n) = Π (p_i^(e_i+1) - 1) / (p_i - 1)
        τ(n) = Π (e_i + 1)
    這些公式本身也會被「枚舉所有因數」的獨立做法交叉驗證。

應用場景：
    - 密碼學：RSA 需要兩個大質數（用 Miller-Rabin 產生）、
      而 RSA 的安全性正是建立在大整數「難以分解」上，Pollard-Rho 是暴力分解的主力之一
    - 競賽題：判斷大數是否為質數、求最大質因數、計算 φ / σ / τ、
      處理「n 的所有因數」（枚舉因數比枚舉 1..n 快得多）
    - 工程：雜湊表大小取質數、亂數週期長度設計、區塊鏈與零知識證明裡的模運算
    - Pollard-Rho 的 ρ 形思路也出現在：生日攻擊、雜湊碰撞搜尋、離散對數的 BSGS

複雜度：
    Miller-Rabin      O(k log³ n)（k 為底數個數，此處 k = 12）
    Pollard-Rho       期望 O(√p) 次迭代找出一個因數（p 為最小質因數）
    完整分解         期望 O(n^1/4 log n) 等級，實務上對 10^12 是瞬間完成
    空間              O(log n)
"""

import random
import re
import sys
from typing import Dict, List, Tuple

INT_RE = re.compile(r"^[+-]?[0-9]+$")

NMAX = 10 ** 12      # IO 模式下單一數值的上限
KMAX = 20            # 最多 20 個數值
BRUTE_TRIAL = 10 ** 6    # ≤ 這個值就跑試除法分解 + 枚舉因數對拍
BRUTE_PHI = 20000        # ≤ 這個值才用「數互質個數」對拍 φ

MR_BASES = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37)
MASK64 = (1 << 64) - 1

SEED_MAIN = 20261017
SEED_CHECK = 987654321          # 第二次獨立分解用的種子

# ---------------------------------------------------------------- 確定性亂數（LCG）


class Rng:
    """64 位元 LCG；Python 與 C++ 行為完全一致（都靠遮罩 / 自然溢位截斷）。"""

    __slots__ = ("s",)

    def __init__(self, seed: int) -> None:
        self.s = seed & MASK64

    def next(self) -> int:
        self.s = (self.s * 6364136223846793005 + 1442695040888963407) & MASK64
        return self.s


# ---------------------------------------------------------------- Miller-Rabin


def is_prime(n: int) -> bool:
    """確定性 Miller-Rabin（對 n < 3.3×10^24 成立）。"""
    if n < 2:
        return False
    for p in MR_BASES:
        if n % p == 0:
            return n == p
    d = n - 1
    s = 0
    while d % 2 == 0:
        d //= 2
        s += 1
    for a in MR_BASES:
        x = pow(a, d, n)
        if x == 1 or x == n - 1:
            continue
        witness = True                       # 尚未出現 -1 → a 是「n 為合數」的見證
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1:
                witness = False
                break
        if witness:
            return False
    return True


def is_prime_brute(n: int) -> bool:
    """試除法：只用於小數對拍。"""
    if n < 2:
        return False
    d = 2
    while d * d <= n:
        if n % d == 0:
            return False
        d += 1
    return True


# ---------------------------------------------------------------- Pollard-Rho


def pollard_rho(n: int, rng: Rng) -> int:
    """回傳 n 的一個非平凡因數（1 < d < n）。"""
    if n % 2 == 0:
        return 2
    if n % 3 == 0:
        return 3
    while True:
        c = rng.next() % (n - 1) + 1
        x = rng.next() % (n - 1) + 1
        y = x
        d = 1
        while d == 1:
            x = (x * x + c) % n
            y = (y * y + c) % n
            y = (y * y + c) % n
            diff = x - y if x >= y else y - x
            d = _gcd(diff, n)
        if d != n:
            return d


def _gcd(a: int, b: int) -> int:
    while b:
        a, b = b, a % b
    return a


def factorize(n: int, seed: int) -> List[int]:
    """回傳 n 的質因數（含重複，遞增排序）；n < 2 回傳空串列。"""
    if n < 2:
        return []
    rng = Rng(seed)
    out: List[int] = []
    stack = [n]
    while stack:
        v = stack.pop()
        if v == 1:
            continue
        if is_prime(v):
            out.append(v)
            continue
        d = pollard_rho(v, rng)
        stack.append(d)
        stack.append(v // d)
    out.sort()
    return out


def trial_factor(n: int) -> List[int]:
    """試除法分解（對照組）。"""
    if n < 2:
        return []
    res = []
    d = 2
    while d * d <= n:
        while n % d == 0:
            res.append(d)
            n //= d
        d += 1
    if n > 1:
        res.append(n)
    return res


# ---------------------------------------------------------------- 因數與積性函數


def group_factors(fac: List[int]) -> List[Tuple[int, int]]:
    out: List[Tuple[int, int]] = []
    for p in fac:
        if out and out[-1][0] == p:
            out[-1] = (p, out[-1][1] + 1)
        else:
            out.append((p, 1))
    return out


def all_divisors(fac: List[int]) -> List[int]:
    """枚舉所有因數（遞增排序）。"""
    divs = [1]
    for p, e in group_factors(fac):
        nxt = []
        mul = 1
        for _ in range(e + 1):
            for d in divs:
                nxt.append(d * mul)
            mul *= p
        divs = nxt
    divs.sort()
    return divs


def divisors_with_phi(g: List[Tuple[int, int]]) -> List[Tuple[int, int]]:
    """枚舉所有因數並順便算出每個因數的 φ，回傳 [(因數, φ(因數))]。

    φ 也用指數向量算：d = Π p_i^k_i 時 φ(d) = Π_{k_i ≥ 1} p_i^(k_i-1)·(p_i-1)。
    這樣「枚舉因數」就不必對每個因數重新做一次分解。
    """
    res: List[Tuple[int, int]] = [(1, 1)]
    for p, e in g:
        nxt: List[Tuple[int, int]] = []
        for d, phid in res:
            pk = 1
            for k in range(e + 1):
                if k == 0:
                    nxt.append((d, phid))
                else:
                    pk *= p
                    nxt.append((d * pk, phid * (p - 1) * pk // p))
        res = nxt
    return res


def phi_from(fac: List[int]) -> int:
    """φ(n) = Π p^(e-1)·(p-1)。"""
    res = 1
    for p, e in group_factors(fac):
        res *= (p - 1) * p ** (e - 1)
    return res


def phi_brute(n: int) -> int:
    """暴力數「1..n 中與 n 互質的個數」（只用於小數對拍）。"""
    if n < 1:
        return 0
    cnt = 0
    for i in range(1, n + 1):
        if _gcd(i, n) == 1:
            cnt += 1
    return cnt


def brute_divisors(n: int) -> Tuple[int, int]:
    """回傳 (τ(n), σ(n))，用試除枚舉（對照組）。"""
    if n < 1:
        return 0, 0
    tau = 0
    sigma = 0
    d = 1
    while d * d <= n:
        if n % d == 0:
            tau += 1
            sigma += d
            q = n // d
            if q != d:
                tau += 1
                sigma += q
        d += 1
    return tau, sigma


# ---------------------------------------------------------------- 輸出


def fmt_factors(fac: List[int]) -> str:
    if not fac:
        return "-"
    return "*".join("%d^%d" % (p, e) for p, e in group_factors(fac))


def run(nums: List[int]) -> str:
    lines = []
    ok = 1
    for n in nums:
        fac = factorize(n, SEED_MAIN)
        fac2 = factorize(n, SEED_CHECK)          # 換一組亂數種子再分解一次
        if fac != fac2:
            ok = 0
        if fac != sorted(fac):
            ok = 0
        prod = 1
        for p in fac:
            prod *= p
            if not is_prime(p):
                ok = 0
        if n >= 2 and prod != n:
            ok = 0
        if n < 2 and fac:
            ok = 0

        prime = 1 if is_prime(n) else 0
        if 2 <= n <= BRUTE_TRIAL and prime != (1 if is_prime_brute(n) else 0):
            ok = 0

        g = group_factors(fac)
        omega = len(g)
        omega_total = len(fac)
        lpf = g[-1][0] if g else 0
        phi = phi_from(fac) if n >= 1 else 0
        tau_formula = 1
        for _p, e in g:
            tau_formula *= e + 1
        tau_formula = tau_formula if n >= 1 else 0

        if n >= 1:
            dp = divisors_with_phi(g)
            tau = len(dp)
            sigma = sum(d for d, _ in dp)
            if tau != tau_formula:
                ok = 0
            # ∑_{d|n} φ(d) = n：尤拉函數的經典恆等式，獨立驗證「分解」與「φ」自洽
            if sum(ph for _, ph in dp) != n:
                ok = 0
            if n <= BRUTE_TRIAL:
                bt, bs = brute_divisors(n)
                if (tau, sigma) != (bt, bs):
                    ok = 0
                if fac != trial_factor(n):
                    ok = 0
            if 1 <= n <= BRUTE_PHI and phi != phi_brute(n):
                ok = 0
        else:
            tau = 0
            sigma = 0

        lines.append("%d %d %d %d %d %d %d %d %s" % (
            n, prime, lpf, omega, omega_total, phi, sigma, tau, fmt_factors(fac)))
    lines.append(str(ok))
    return "\n".join(lines)


# ---------------------------------------------------------------- 內置測試


def run_tests() -> None:
    # ---- 固定用例：質數 ----
    for p in (2, 3, 5, 7, 11, 13, 97, 1009, 7919, 104729, 999983, 1000003):
        assert is_prime(p), p
        assert factorize(p, SEED_MAIN) == [p]
    for c in (0, 1, 4, 9, 25, 561, 1105, 1729, 2465, 2821, 6601, 8911):
        assert not is_prime(c), c

    # Carmichael 數：對所有與 n 互質的 a 都滿足費馬式，但確實是合數
    assert 561 == 3 * 11 * 17
    assert not is_prime(561)
    assert factorize(561, SEED_MAIN) == [3, 11, 17]

    # ---- 固定用例：分解 ----
    assert factorize(1, SEED_MAIN) == []
    assert factorize(0, SEED_MAIN) == []
    assert factorize(12, SEED_MAIN) == [2, 2, 3]
    assert factorize(100, SEED_MAIN) == [2, 2, 5, 5]
    assert factorize(1024, SEED_MAIN) == [2] * 10
    assert factorize(10 ** 12, SEED_MAIN) == [2] * 12 + [5] * 12
    assert factorize(999985999949, SEED_MAIN) == [999983, 1000003]   # 兩個大質數的乘積

    # ---- 固定用例：積性函數 ----
    assert phi_from([2, 2, 3]) == 4            # φ(12) = 4
    assert phi_from([7]) == 6
    assert phi_from([]) == 1                   # φ(1) = 1
    assert brute_divisors(12) == (6, 28)
    assert brute_divisors(1) == (1, 1)
    assert phi_brute(12) == 4
    assert phi_brute(1) == 1
    assert all_divisors([2, 2, 3]) == [1, 2, 3, 4, 6, 12]

    # 因數個數與因數和的經典值
    assert run([12]).split("\n")[0] == "12 0 3 2 3 4 28 6 2^2*3^1"
    assert run([7]).split("\n")[0] == "7 1 7 1 1 6 8 2 7^1"
    assert run([1]).split("\n")[0] == "1 0 0 0 0 1 1 1 -"
    assert run([0]).split("\n")[0] == "0 0 0 0 0 0 0 0 -"

    # ---- 隨機對拍 ----
    rng = random.Random(20261017)
    for _ in range(400):
        kind = rng.randint(0, 5)
        if kind == 0:
            n = rng.randint(0, 60)
        elif kind == 1:
            n = rng.randint(0, 20000)
        elif kind == 2:
            n = rng.randint(0, BRUTE_TRIAL)
        elif kind == 3:
            n = rng.randint(0, NMAX)
        elif kind == 4:                       # 半質數（兩個質數相乘）
            a = rng.choice([2, 3, 5, 7, 997, 10007, 99991, 999983])
            b = rng.choice([2, 3, 5, 7, 997, 10007, 99991, 999983])
            n = min(a * b, NMAX)
        else:                                 # 高次冪
            n = rng.randint(2, 60) ** rng.randint(2, 6)
            n = min(n, NMAX)
        out = run([n])
        assert out.split("\n")[-1] == "1", (n, out)

    # ---- 一批數字同時處理 ----
    batch = [0, 1, 2, 561, 1024, 1000000, 999983, 10 ** 12]
    assert run(batch).split("\n")[-1] == "1"


# ---------------------------------------------------------------- IO


def emit(text: str) -> None:
    sys.stdout.buffer.write(text.encode("latin-1"))
    sys.stdout.buffer.flush()


def run_io(raw: str) -> None:
    toks = raw.split()

    def parse(t: str, default: int) -> int:
        if not INT_RE.match(t):
            return default
        v = int(t)
        return max(-(10 ** 18), min(10 ** 18, v))

    k = parse(toks[0], 0) if len(toks) >= 1 else 0
    k = max(0, min(k, KMAX))
    nums = []
    for i in range(1, 1 + k):
        v = parse(toks[i], 0) if i < len(toks) else 0
        nums.append(max(0, min(v, NMAX)))
    emit(run(nums) + "\n")


def main() -> None:
    sys.setrecursionlimit(1000000)
    raw = sys.stdin.buffer.read().decode("latin-1")
    if raw.strip():
        run_io(raw)
    else:
        run_tests()
        emit("all tests passed\n")


if __name__ == "__main__":
    main()
