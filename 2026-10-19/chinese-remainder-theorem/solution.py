"""中國剩餘定理（CRT）與擴展中國剩餘定理（exCRT）

題意：
    一口氣解四組同餘問題，每組都同時給出「高效演算法」與「獨立實作的暴力版本」互相對拍。

    A. **同餘方程組**：給 n 組 (a_i, m_i)，求所有滿足 x ≡ a_i (mod m_i) 的 x。
       分別用「互質版 CRT」（乘積 + 逆元）與「一般版 exCRT」（兩兩合併）各解一次。
    B. **線性同餘方程** a·x ≡ b (mod m)：求 gcd、解的個數、最小非負解。
    C. **CRT 還原**：給定 x 與兩個模數 p、q，先取 (x mod p, x mod q)，
       再用 exCRT 還原出 x mod lcm(p, q)，跟直接算的 x mod lcm(p, q) 比對。

思路：
    ### 1. 擴展歐幾里得（exgcd）
    迭代版即可：維護 (r, s, t) 使 r = a·s + b·t，不斷做帶餘除法直到 r = 0。
    限制輸入為**非負**，就完全避開「Python 的 // 是 floor、C++ 的 / 是 truncate」
    在負數上會分歧的老坑。

    ### 2. 合併兩個同餘式
    要同時滿足 x ≡ a₁ (mod m₁) 與 x ≡ a₂ (mod m₂)，寫 x = a₁ + m₁·t：
        m₁·t ≡ a₂ - a₁  (mod m₂)
    令 g = gcd(m₁, m₂)。若 (a₂ - a₁) 不能被 g 整除 → **無解**。
    否則整個式子同除以 g：
        (m₁/g)·t ≡ (a₂-a₁)/g  (mod m₂/g)
    此時 m₁/g 與 m₂/g 互質，逆元存在，直接得到 t，再回代
        x₀ = a₁ + m₁·t  (mod lcm(m₁, m₂))
    逐條合併下去就是 exCRT，模數不互質也能用。

    ### 3. 互質版 CRT（孫子定理的標準形式）
    M = ∏ m_i，對每個 i 令 M_i = M / m_i、inv_i = M_i^{-1} mod m_i，則
        x = Σ a_i · M_i · inv_i  (mod M)
    因為 M_i·inv_i ≡ 1 (mod m_i) 且 ≡ 0 (mod m_j, j≠i)，每一項「只在自己那一維生效」。
    這個形式漂亮但**只在模數兩兩互質時成立**，非互質必須回到第 2 步。

    ### 4. 線性同餘 a·x ≡ b (mod m)
    令 g = gcd(a, m)。b 不是 g 的倍數就無解；否則同除以 g 後 a/g 在模 m/g 下可逆，
        x₀ = (b/g) · (a/g)^{-1}  (mod m/g)
    在 [0, m) 裡恰好有 **g 個**解：x₀, x₀ + m/g, x₀ + 2·m/g, …

    ### 5. 溢位處理（跨語言一致）
    Python 是無限精度、C++ 的 int64 會繞回，兩邊必須用同一套截斷規則：
    合併模數時先做**除法檢查** `M > LIMIT / (m/g)` 就判定溢位（回傳 -2 -2），
    乘法則在 C++ 用 `__int128` 承接。LIMIT 定為 10^18。

    ### 6. 暴力對拍
    當所有模數的 lcm 不超過 200000 時，直接從 0 枚舉到 lcm-1 找最小解；
    線性同餘則在 m ≤ 200000 時枚舉 [0, m)。這能抓出「公式對但正規化錯」的 bug。

應用場景：
    - **RSA 的 CRT 加速**：分別在模 p 與模 q 下做解密，最後用 CRT 合併，速度約 4 倍
    - **大數運算**：把一個巨大整數拆成幾個小模數下的運算，最後 CRT 還原（NTT 多模合併）
    - **秘密共享（Asmuth-Bloom 門限方案）**：每人拿到一個不同模數下的餘數，
      湊齊足夠多人才能用 CRT 還原祕密，少一個都不行
    - **排程 / 週期對齊**：「每 3 天一次、每 5 天一次、每 7 天一次」何時同時發生
    - 競賽題裡的「韓信點兵」、模數非互質的計數題

複雜度：
    exgcd / 逆元            O(log min(a,b)) 時間 / O(1) 空間
    exCRT 合併 n 條          O(n log M) 時間 / O(n) 空間
    互質版 CRT               O(n log M) 時間 / O(n) 空間
    線性同餘                 O(log m) 時間 / O(1) 空間
    暴力對拍                 O(lcm) 或 O(m) 時間
"""

import random
import re
import sys
from typing import List, Optional, Tuple

INT_RE = re.compile(r"^[+-]?[0-9]+$")

LIMIT = 10 ** 18            # 合併模數的上限，超過就標記為溢位
BRUTE_LIMIT = 200000        # lcm 或 m 不超過這個值才跑暴力對拍
NMAX = 12                   # 同餘方程組的條數上限
MMAX = 10 ** 5              # 單一模數上限
AMAX = 10 ** 9              # 餘數 / 係數上限
XMAX = 10 ** 18             # 還原用的大數上限
PQMAX = 10 ** 6             # 還原用的兩個模數上限


# ---------------------------------------------------------------- 數論基礎


def exgcd(a: int, b: int) -> Tuple[int, int, int]:
    """回傳 (g, s, t) 使 a·s + b·t = g = gcd(a, b)。要求 a, b ≥ 0。"""
    old_r, r = a, b
    old_s, s = 1, 0
    old_t, t = 0, 1
    while r != 0:
        q = old_r // r
        old_r, r = r, old_r - q * r
        old_s, s = s, old_s - q * s
        old_t, t = t, old_t - q * t
    if old_r < 0:
        old_r, old_s, old_t = -old_r, -old_s, -old_t
    return old_r, old_s, old_t


def gcd(a: int, b: int) -> int:
    while b:
        a, b = b, a % b
    return a


def mod_inverse(a: int, m: int) -> Optional[int]:
    """a 在模 m 下的乘法逆元；不存在（或 m == 1）回傳 None / 0。"""
    if m == 1:
        return 0
    g, s, _t = exgcd(a % m, m)
    if g != 1:
        return None
    return s % m


def merge(a1: int, m1: int, a2: int, m2: int) -> Tuple[int, int, int]:
    """合併 x ≡ a1 (mod m1) 與 x ≡ a2 (mod m2)。回傳 (狀態, r, lcm)。
    狀態：0 成功 / 1 無解 / 2 溢位（lcm 超過 LIMIT）。"""
    if m1 <= 0 or m2 <= 0:
        return (1, -1, -1)
    a1 %= m1
    a2 %= m2
    g = gcd(m1, m2)
    diff = a2 - a1
    if diff % g != 0:
        return (1, -1, -1)
    step = m2 // g
    # 溢位檢查：用除法，避免真的乘出超過 LIMIT 的數
    if m1 > LIMIT // step:
        return (2, -2, -2)
    M = m1 * step
    m2g = m2 // g
    inv = mod_inverse((m1 // g) % m2g, m2g)
    if inv is None:
        return (1, -1, -1)
    t = ((diff // g) % m2g) * inv % m2g
    r = (a1 + (m1 % M) * t) % M
    return (0, r, M)


# ---------------------------------------------------------------- A. 同餘方程組


def excrt(eqs: List[Tuple[int, int]]) -> Tuple[int, int, int]:
    """回傳 (狀態, r, M)。狀態：0 有解 / 1 無解 / 2 溢位。"""
    r, M = 0, 1
    for a, m in eqs:
        st, nr, nM = merge(r, M, a % m, m)
        if st != 0:
            return (st, nr, nM)
        r, M = nr, nM
    return (0, r % M, M)


def crt_coprime(eqs: List[Tuple[int, int]]) -> Tuple[int, int, int]:
    """互質版 CRT。回傳 (狀態, r, M)。狀態：0 有解 / 1 模數不兩兩互質 / 2 溢位。"""
    n = len(eqs)
    for i in range(n):
        for j in range(i + 1, n):
            if gcd(eqs[i][1], eqs[j][1]) != 1:
                return (1, -1, -1)
    M = 1
    for _a, m in eqs:
        if M > LIMIT // m:
            return (2, -2, -2)
        M *= m
    r = 0
    for a, m in eqs:
        if m == 0:
            return (1, -1, -1)
        Mi = M // m
        inv = mod_inverse(Mi % m, m)
        if inv is None:
            return (1, -1, -1)
        r = (r + (a % m) * Mi % M * inv) % M
    return (0, r % M, M)


def lcm_all(mods: List[int]) -> Optional[int]:
    """所有模數的 lcm；超過 LIMIT 回傳 None。"""
    L = 1
    for m in mods:
        g = gcd(L, m)
        step = m // g
        if L > LIMIT // step:
            return None
        L *= step
    return L


def brute_system(eqs: List[Tuple[int, int]]) -> Tuple[int, int]:
    """暴力：枚舉 [0, lcm)。回傳 (狀態, 最小解)。狀態 0 有解 / 1 無解 / 2 太大不跑。"""
    mods = [m for _a, m in eqs]
    L = lcm_all(mods) if mods else 1
    if L is None or L > BRUTE_LIMIT:
        return (2, -1)
    for x in range(L):
        good = True
        for a, m in eqs:
            if m == 0:
                good = False
                break
            if x % m != a % m:
                good = False
                break
        if good:
            return (0, x)
    return (1, -1)


# ---------------------------------------------------------------- B. 線性同餘


def solve_linear(a: int, b: int, m: int) -> Tuple[int, int, int]:
    """解 a·x ≡ b (mod m)。回傳 (g, 解的個數, 最小非負解)。無解時個數 0、解 -1。"""
    if m == 0:
        return (0, 0, -1)
    g = gcd(a % m, m)
    a %= m
    b %= m
    if b % g != 0:
        return (g, 0, -1)
    m2 = m // g
    if m2 == 1:
        return (g, g, 0)
    inv = mod_inverse((a // g) % m2, m2)
    if inv is None:
        return (g, 0, -1)
    x0 = ((b // g) % m2) * inv % m2
    return (g, g, x0)


def brute_linear(a: int, b: int, m: int) -> Tuple[int, int, int]:
    """暴力枚舉 [0, m)。回傳 (g, 個數, 最小解)；m 太大回傳 (-1, -1, -1)。"""
    if m == 0 or m > BRUTE_LIMIT:
        return (-1, -1, -1)
    g = gcd(a % m, m)
    cnt = 0
    first = -1
    for x in range(m):
        if (a % m) * x % m == b % m:
            cnt += 1
            if first < 0:
                first = x
    return (g, cnt, first)


# ---------------------------------------------------------------- C. CRT 還原


def crt_restore(x: int, p: int, q: int) -> Tuple[int, int]:
    """從 (x mod p, x mod q) 還原。回傳 (還原值, 模數)；失敗回傳 (-1, -1)。"""
    if p <= 0 or q <= 0:
        return (-1, -1)
    st, r, M = merge(x % p, p, x % q, q)
    if st != 0:
        return (-2, -2)          # 同一個 x 產生的方程組一定有解，只可能是溢位
    return (r, M)


# ---------------------------------------------------------------- 主流程


def run(eqs: List[Tuple[int, int]], lin: Tuple[int, int, int],
        rest: Tuple[int, int, int]) -> str:
    ok = True

    # ---- A. 同餘方程組 ----
    st1, r1, M1 = excrt(eqs)
    st2, r2, M2 = crt_coprime(eqs)

    if st1 == 0:
        for a, m in eqs:
            if m > 0 and r1 % m != a % m:
                ok = False
                break
    bst, bx = brute_system(eqs)
    if bst == 0:
        if st1 != 0 or r1 != bx:
            ok = False
    elif bst == 1:
        if st1 != 1:            # 暴力說無解，exCRT 也必須說無解
            ok = False
    if st2 == 0:
        if st1 != 0 or r1 != r2 or M1 != M2:
            ok = False

    # ---- B. 線性同餘 ----
    a, b, m = lin
    g, cnt, xmin = solve_linear(a, b, m)
    if cnt > 0:
        if m > 0 and (a % m) * xmin % m != b % m:
            ok = False
    bg, bcnt, bxmin = brute_linear(a, b, m)
    if bg != -1:
        if (g, cnt, xmin) != (bg, bcnt, bxmin):
            ok = False

    # ---- C. CRT 還原 ----
    x, p, q = rest
    rv, rm = crt_restore(x, p, q)
    if rm > 0:
        direct = x % rm
        if rv != direct:
            ok = False
    else:
        ok = False

    lines = []
    lines.append(f"{r1} {M1}")
    lines.append(f"{r2} {M2}")
    lines.append(f"{g} {cnt} {xmin}")
    lines.append(f"{rv} {rm}")
    lines.append("1" if ok else "0")
    return "\n".join(lines) + "\n"


# ---------------------------------------------------------------- 內置測試


def run_tests() -> None:
    # ---- 固定用例 1：孫子定理原題 x ≡ 2 (mod 3), 3 (mod 5), 2 (mod 7) ----
    out = run([(2, 3), (3, 5), (2, 7)], (1, 1, 1), (0, 1, 1))
    assert out.split("\n")[0] == "23 105", out
    assert out.split("\n")[1] == "23 105", out

    # ---- 固定用例 2：非互質但有解 x ≡ 2 (mod 4), x ≡ 4 (mod 6) ----
    out = run([(2, 4), (4, 6)], (1, 1, 1), (0, 1, 1))
    assert out.split("\n")[0] == "10 12", out
    assert out.split("\n")[1] == "-1 -1", out      # 非互質 → 互質版不適用

    # ---- 固定用例 3：無解 x ≡ 1 (mod 4), x ≡ 2 (mod 6) ----
    out = run([(1, 4), (2, 6)], (1, 1, 1), (0, 1, 1))
    assert out.split("\n")[0] == "-1 -1", out
    assert out.split("\n")[-2] == "1", out

    # ---- 固定用例 4：線性同餘 6x ≡ 4 (mod 10) ----
    out = run([(0, 1)], (6, 4, 10), (0, 1, 1))
    assert out.split("\n")[2] == "2 2 4", out

    # ---- 固定用例 5：線性同餘無解 6x ≡ 5 (mod 10) ----
    out = run([(0, 1)], (6, 5, 10), (0, 1, 1))
    assert out.split("\n")[2] == "2 0 -1", out

    # ---- 固定用例 6：CRT 還原 ----
    out = run([(0, 1)], (1, 1, 1), (123456789, 97, 1000))
    assert out.split("\n")[3] == "72789 97000", out

    # ---- 空方程組 ----
    out = run([], (1, 1, 1), (0, 1, 1))
    assert out.split("\n")[0] == "0 1", out
    assert out.split("\n")[-2] == "1", out

    # ---- 數論基礎 ----
    assert exgcd(0, 0) == (0, 1, 0)
    assert exgcd(35, 15) == (5, 1, -2), exgcd(35, 15)
    assert mod_inverse(3, 11) == 4
    assert mod_inverse(3, 12) is None
    assert mod_inverse(0, 1) == 0
    assert gcd(0, 0) == 0

    # ---- 隨機對拍 ----
    rng = random.Random(20261019)
    for _ in range(600):
        n = rng.randint(0, 8)
        eqs = [(rng.randint(0, 50), rng.randint(1, 30)) for _ in range(n)]
        a = rng.randint(0, 200)
        b = rng.randint(0, 200)
        m = rng.randint(1, 300)
        x = rng.randint(0, 10 ** 6)
        p = rng.randint(1, 200)
        q = rng.randint(1, 200)
        out = run(eqs, (a, b, m), (x, p, q))
        assert out.split("\n")[-2] == "1", (eqs, a, b, m, x, p, q, out)

    # ---- 故意構造「一定有解」的系統 ----
    for _ in range(300):
        n = rng.randint(1, 8)
        mods = [rng.randint(1, 40) for _ in range(n)]
        L = lcm_all(mods)
        if L is None or L > LIMIT:
            continue
        true_x = rng.randint(0, max(0, L - 1))
        eqs = [(true_x % m, m) for m in mods]
        out = run(eqs, (1, 1, 1), (true_x, mods[0], mods[-1]))
        assert out.split("\n")[0] == f"{true_x % L} {L}", (eqs, L, out)


# ---------------------------------------------------------------- IO


def emit(text: str) -> None:
    sys.stdout.buffer.write(text.encode("latin-1"))
    sys.stdout.buffer.flush()


def run_io(raw: str) -> None:
    toks = raw.split()
    idx = 0

    def next_int(default: int) -> int:
        nonlocal idx
        if idx >= len(toks):
            return default
        t = toks[idx]
        idx += 1
        if not INT_RE.match(t):
            return default
        try:
            v = int(t)
        except ValueError:
            return default
        return max(-(10 ** 18), min(10 ** 18, v))

    n = next_int(0)
    if n < 0:
        n = 0
    if n > NMAX:
        n = NMAX
    eqs: List[Tuple[int, int]] = []
    for _ in range(n):
        a = max(0, min(next_int(0), AMAX))
        m = next_int(1)
        if m < 1:
            m = 1
        if m > MMAX:
            m = MMAX
        eqs.append((a, m))

    a = max(0, min(next_int(1), AMAX))
    b = max(0, min(next_int(0), AMAX))
    m = next_int(1)
    if m < 1:
        m = 1
    if m > MMAX:
        m = MMAX
    lin = (a, b, m)

    x = max(0, min(next_int(0), XMAX))
    p = next_int(1)
    if p < 1:
        p = 1
    if p > PQMAX:
        p = PQMAX
    q = next_int(1)
    if q < 1:
        q = 1
    if q > PQMAX:
        q = PQMAX

    emit(run(eqs, lin, (x, p, q)))


def main() -> None:
    raw = sys.stdin.buffer.read().decode("latin-1")
    if raw.strip():
        run_io(raw)
    else:
        run_tests()
        emit("all tests passed\n")


if __name__ == "__main__":
    main()
