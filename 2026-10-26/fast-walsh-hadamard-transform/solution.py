"""FWT 快速沃爾什變換：XOR / OR / AND 三種位運算卷積

題意：
    給定整數 k（0 <= k <= 16）與兩個長度為 2^k 的序列 A、B，
    以及一個運算 op ∈ {XOR, OR, AND}，求它們的**位運算卷積**：

        XOR 卷積   C[t] = Σ_{i^j = t} A[i]·B[j]
        OR  卷積   C[t] = Σ_{i|j = t} A[i]·B[j]
        AND 卷積   C[t] = Σ_{i&j = t} A[i]·B[j]

    所有計算對 MOD = 998244353 取模；輸入的負數先化歸到 [0, MOD)。
    輸出長度 2^k 的序列 C。

思路：
    ### 1. 位運算卷積 = 布爾格上的多項式乘法
    把 A、B 看成「以布爾格 {0,1}^k 為下標、以逐位元運算為乘法」的
    多項式。三種運算各對應一個交換半群，卷積就是這個半群上的
    多項式乘法 —— 直接算是 O(4^k)。

    ### 2. FWT：把卷積變成逐點乘
    每種位運算都有一個對應的線性變換（zeta 變換），把序列送到
    「下標和」的域，卷積就退化為逐點相乘，再逆變換回來：

        XOR：蝶形 (a, b) → (a+b, a−b)            逆變換同形，最後乘 2^{-k}
        OR ：zeta  (a, b) → (a, a+b)             逆變換 (a, b) → (a, b−a)
        AND：zeta  (a, b) → (a+b, b)             逆變換 (a, b) → (a−b, b)

    三種變換都是「按位長度倍增」的原地蝶形，總複雜度 O(k·2^k)。
    推導方式完全相同：找出 (a±b) 這樣的線性組合，使得
    transform(A)[s] = Σ_{i ⊆ s 或 i^j=s} A[i]（依運算而定）。

    ### 3. 為什麼 OR 與 AND 的zeta方向看起來是「反的」
    OR 的 zeta 是「把低位併進高位」(a, b) → (a, a+b)：s = i|j 涵蓋
    所有 i ⊆ s。AND 則是「把高位併進低位」(a, b) → (a+b, b)：
    s = i&j 涵蓋所有 s ⊆ i。兩者互為對偶（把 0/1 互換、下標取反即可
    互相轉化），寫反是最常見的 bug，所以測試裡對兩種方向都做
    暴力對拍。

    ### 4. 模運算
    MOD = 998244353 是質數，2^{-k} 用費馬小定理求逆元。
    負數輸入先化歸到 [0, MOD)，C++ 與 Python 用同一條規則
    ((v % MOD) + MOD) % MOD，保證輸出逐位元組一致。

複雜度（n = 2^k）：
    三次變換 + 逐點乘      時間 O(k·n)
    空間                  O(n)

應用場景：
    - 多項式乘法：OR / AND 卷積是「按位受限」的多項式乘法，
      集合冪級數（set power series）的基礎運算
    - 子集計數：統計「兩組條件的 OR/AND 恰好等於某掩碼」的方案數，
      如特徵組合、權限位彙總
    - 期望與機率：擲骰子/抽卡模型中用 XOR 卷積快速求分佈的卷積
      （k 個 bit 的狀態轉移）
    - 圖論：傳遞閉包的布爾版（OR 卷積 + 重複平方法）可加速
      某些位元矩陣運算
    - 競賽常青題：BZOJ 4263 / CF 1119F 一類的 FWT 入門題
"""

import random
import sys
from typing import List, Tuple

KMAX = 16
MOD = 998244353


def _norm(v: int) -> int:
    return ((v % MOD) + MOD) % MOD


def fwt_xor(a: List[int], inverse: bool) -> None:
    n = len(a)
    step = 1
    while step < n:
        for i in range(0, n, step * 2):
            for j in range(i, i + step):
                u = a[j]
                v = a[j + step]
                a[j] = (u + v) % MOD
                a[j + step] = (u - v) % MOD
        step *= 2
    if inverse:
        inv_n = pow(n, MOD - 2, MOD)
        for i in range(n):
            a[i] = a[i] * inv_n % MOD


def fwt_or(a: List[int], inverse: bool) -> None:
    n = len(a)
    step = 1
    while step < n:
        for i in range(0, n, step * 2):
            for j in range(i, i + step):
                u = a[j]
                v = a[j + step]
                if inverse:
                    a[j] = u % MOD
                    a[j + step] = (v - u) % MOD
                else:
                    a[j] = u % MOD
                    a[j + step] = (u + v) % MOD
        step *= 2


def fwt_and(a: List[int], inverse: bool) -> None:
    n = len(a)
    step = 1
    while step < n:
        for i in range(0, n, step * 2):
            for j in range(i, i + step):
                u = a[j]
                v = a[j + step]
                if inverse:
                    a[j] = (u - v) % MOD
                    a[j + step] = v % MOD
                else:
                    a[j] = (u + v) % MOD
                    a[j + step] = v % MOD
        step *= 2


FWT = {"XOR": fwt_xor, "OR": fwt_or, "AND": fwt_and}


def convolve(a: List[int], b: List[int], op: str) -> List[int]:
    fwt = FWT[op]
    fa = a[:]
    fb = b[:]
    fwt(fa, False)
    fwt(fb, False)
    for i in range(len(fa)):
        fa[i] = fa[i] * fb[i] % MOD
    fwt(fa, True)
    return fa


def brute_convolve(a: List[int], b: List[int], op: str) -> List[int]:
    n = len(a)
    c = [0] * n
    for i in range(n):
        if a[i] == 0:
            continue
        for j in range(n):
            if b[j] == 0:
                continue
            if op == "XOR":
                t = i ^ j
            elif op == "OR":
                t = i | j
            else:
                t = i & j
            c[t] = (c[t] + a[i] * b[j]) % MOD
    return c


def solve(k: int, a: List[int], b: List[int], op: str) -> List[int]:
    fa = [_norm(v) for v in a]
    fb = [_norm(v) for v in b]
    return convolve(fa, fb, op)


def run_tests() -> None:
    # ---- 固定用例（手工驗算：A=[1,2], B=[3,4]）----
    assert solve(1, [1, 2], [3, 4], "XOR") == [11, 10]
    assert solve(1, [1, 2], [3, 4], "OR") == [3, 18]
    assert solve(1, [1, 2], [3, 4], "AND") == [13, 8]
    # k = 0：長度 1，三種運算都歸結為 A[0]·B[0]
    assert solve(0, [7], [9], "XOR") == [63]
    assert solve(0, [7], [9], "OR") == [63]
    assert solve(0, [7], [9], "AND") == [63]
    # 取模邊界：MOD-1 ≡ -1，(-1)·(-1) = 1
    m1 = MOD - 1
    assert solve(1, [m1, m1], [m1, m1], "XOR") == [2, 2]
    assert solve(1, [m1, m1], [m1, m1], "OR") == [1, 3]
    assert solve(1, [m1, m1], [m1, m1], "AND") == [3, 1]
    # 負數輸入：先化歸到 [0, MOD)
    assert solve(1, [-1, 0], [-1, 0], "XOR") == [1, 0]
    assert solve(2, [-1, -2, -3, -4], [-5, -6, -7, -8], "OR") == \
        brute_convolve([_norm(-1), _norm(-2), _norm(-3), _norm(-4)],
                       [_norm(-5), _norm(-6), _norm(-7), _norm(-8)], "OR")
    # 全零 / 全相同
    assert solve(2, [0, 0, 0, 0], [1, 2, 3, 4], "XOR") == [0, 0, 0, 0]
    assert solve(2, [5, 5, 5, 5], [5, 5, 5, 5], "AND") == \
        brute_convolve([5] * 4, [5] * 4, "AND")

    # ---- 變換往返：inverse(transform(x)) == x ----
    rng = random.Random(20261026)
    for _ in range(30):
        k = rng.randint(0, 8)
        n = 1 << k
        x = [rng.randrange(MOD) for _ in range(n)]
        for op in ("XOR", "OR", "AND"):
            y = x[:]
            FWT[op](y, False)
            FWT[op](y, True)
            assert y == x, (op, k)

    # ---- 隨機對拍：FWT vs 暴力雙重循環 ----
    for tc in range(120):
        k = rng.randint(0, 7)
        n = 1 << k
        a = [rng.randrange(MOD) for _ in range(n)]
        b = [rng.randrange(MOD) for _ in range(n)]
        for op in ("XOR", "OR", "AND"):
            got = solve(k, a, b, op)
            want = brute_convolve([_norm(v) for v in a],
                                  [_norm(v) for v in b], op)
            assert got == want, (op, k, a, b)

    # ---- 稀疏性壓力：大量 0 與少量大值 ----
    for tc in range(20):
        k = 8
        n = 1 << k
        a = [0] * n
        b = [0] * n
        for _ in range(5):
            a[rng.randrange(n)] = rng.randrange(MOD)
            b[rng.randrange(n)] = rng.randrange(MOD)
        for op in ("XOR", "OR", "AND"):
            assert solve(k, a, b, op) == brute_convolve(a, b, op), (op, tc)


# ---------------------------------------------------------------- IO 模式

OPS = {"XOR", "OR", "AND"}


def _parse_int(tok: str):
    """整數字串解析；越界（超出 ±10^18）或非純整數回 None。"""
    if not tok:
        return None
    i = 0
    if tok[0] in "+-":
        i = 1
    if i == len(tok):
        return None
    for ch in tok[i:]:
        if ch < "0" or ch > "9":
            return None
    val = int(tok)
    if val > 10 ** 18:
        return None
    return val


def run_io(raw: str) -> None:
    toks = raw.split()
    pos = 0

    def nxt():
        nonlocal pos
        if pos >= len(toks):
            return None
        t = toks[pos]
        pos += 1
        return t

    def bad() -> None:
        emit("error\n")

    k = _parse_int(nxt() or "")
    if k is None or not (0 <= k <= KMAX):
        bad()
        return
    n = 1 << k
    a: List[int] = []
    for _ in range(n):
        v = _parse_int(nxt() or "")
        if v is None:
            bad()
            return
        a.append(v)
    b: List[int] = []
    for _ in range(n):
        v = _parse_int(nxt() or "")
        if v is None:
            bad()
            return
        b.append(v)
    op = nxt()
    if op not in OPS:
        bad()
        return
    res = solve(k, a, b, op)
    emit(" ".join(map(str, res)) + "\n")


def emit(s: str) -> None:
    sys.stdout.buffer.write(s.encode("latin-1"))


def main() -> None:
    raw = sys.stdin.buffer.read().decode("latin-1")
    if raw.strip():
        run_io(raw)
    else:
        run_tests()
        emit("all tests passed\n")


if __name__ == "__main__":
    main()
