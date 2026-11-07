# -*- coding: utf-8 -*-
"""NTT 快速數論變換：模 998244353 的多項式乘法

題意：給定兩個係數在模 998244353 意義下的多項式 f、g，
      求卷積 h(x) = f(x) * g(x) mod 998244353。

思路（NTT）：
  998244353 = 119 * 2^23 + 1 是質數，且原根 g = 3 的階恰為 p-1，
  含足夠長的 2 的冪次單位根，可做快速數論變換（精確無浮點誤差）。
  1. 補零到 len >= n + m - 1 的 2 的冪；
  2. 位反轉置換 + 蝶形迭代做正變換；
  3. 點值逐點相乘；
  4. 用 n^-1 (mod p) 縮放做逆變換。
  蝶形中旋轉因子用預處理的冪表（w^k = 3^(k * (p-1)/len)），
  自底向上合併，長度 1 -> 2 -> 4 -> ... -> len。

IO 模式（stdin 有內容時）：
  n m            （兩多項式的項數，即次數 + 1）
  n 個整數：f 的係數（低次到高次，0 <= 係數 < 998244353）
  m 個整數：g 的係數
  輸出：n + m - 1 個整數，h 的係數（低次到高次）
  無 stdin 時跑內建測試（與 O(n*m) 樸素卷積對拍、正逆變換往返校驗）
  並輸出 all tests passed。
"""
import sys

MOD = 998244353
G = 3  # 998244353 的原根


def ntt(a, invert=False):
    """就地迭代 NTT。invert=True 時做逆變換（不除 n，由呼叫端縮放）。"""
    n = len(a)
    # 位反轉置換
    j = 0
    for i in range(1, n):
        bit = n >> 1
        while j & bit:
            j ^= bit
            bit >>= 1
        j |= bit
        if i < j:
            a[i], a[j] = a[j], a[i]
    # 蝶形：塊長 1 -> n
    length = 2
    while length <= n:
        # 塊內主單位根 w_len = g^((p-1)/length) 或其逆
        w_len = pow(G, (MOD - 1) // length, MOD)
        if invert:
            w_len = pow(w_len, MOD - 2, MOD)
        for start in range(0, n, length):
            w = 1
            half = length // 2
            for k in range(start, start + half):
                u = a[k]
                v = a[k + half] * w % MOD
                a[k] = (u + v) % MOD
                a[k + half] = (u - v) % MOD
                w = w * w_len % MOD
        length <<= 1
    if invert:
        inv_n = pow(n, MOD - 2, MOD)
        for i in range(n):
            a[i] = a[i] * inv_n % MOD


def poly_mul(f, g):
    """NTT 卷積。"""
    need = len(f) + len(g) - 1
    size = 1
    while size < need:
        size <<= 1
    a = list(f) + [0] * (size - len(f))
    b = list(g) + [0] * (size - len(g))
    ntt(a)
    ntt(b)
    for i in range(size):
        a[i] = a[i] * b[i] % MOD
    ntt(a, invert=True)
    return a[:need]


# ---------- 測試 ----------

def brute_mul(f, g):
    """O(n*m) 樸素卷積，作為對拍基準。"""
    res = [0] * (len(f) + len(g) - 1)
    for i, x in enumerate(f):
        for j, y in enumerate(g):
            res[i + j] = (res[i + j] + x * y) % MOD
    return res


def run_tests():
    seed = 20261107

    def rnd(mv):
        nonlocal seed
        seed = (seed * 1103515245 + 12345) & 0x7FFFFFFF
        return seed % mv

    # 固定用例：(1+x) * (1+x) = 1 + 2x + x^2 等
    fixed = [
        ([1, 1], [1, 1]),
        ([1, 2, 3, 4], [5, 6, 7]),
        ([0], [0]),                          # 全零多項式
        ([5], [7]),                          # 常數 × 常數
        ([1, 1, 1, 1, 1, 1, 1, 1], [1, 1]),  # 長度恰為 2 的冪
        ([2, 0, 3], [0, 4, 0, 5]),
    ]
    for f, g in fixed:
        assert poly_mul(f, g) == brute_mul(f, g), (f, g)

    # 隨機對拍 + 正逆變換往返校驗
    for _ in range(200):
        n = rnd(50) + 1
        m = rnd(50) + 1
        f = [rnd(MOD) for _ in range(n)]
        g = [rnd(MOD) for _ in range(m)]
        assert poly_mul(f, g) == brute_mul(f, g), (n, m)
        # 往返：NTT -> INTT 應還原原係數
        size = 1
        while size < n:
            size <<= 1
        a = list(f) + [0] * (size - n)
        ntt(a)
        ntt(a, invert=True)
        assert a[:n] == f, (n, f)

    # 大規模隨機冒煙（驗證長度 2^13 的變換正確性）
    for _ in range(5):
        n = 3000
        f = [rnd(MOD) for _ in range(n)]
        g = [rnd(MOD) for _ in range(n)]
        h = poly_mul(f, g)
        # 抽查 30 個係數與樸素卷積一致
        for _ in range(30):
            k = rnd(2 * n - 1)
            s = 0
            lo = max(0, k - n + 1)
            for i in range(lo, min(n, k + 1)):
                s += f[i] * g[k - i]
            assert h[k] == s % MOD, (k, h[k])

    print("all tests passed")


def main():
    data = sys.stdin.buffer.read().split()
    if not data:
        run_tests()
        return
    it = iter(data)
    n = int(next(it))
    m = int(next(it))
    f = [int(next(it)) % MOD for _ in range(n)]
    g = [int(next(it)) % MOD for _ in range(m)]
    h = poly_mul(f, g)
    sys.stdout.write(" ".join(map(str, h)) + "\n")


if __name__ == "__main__":
    main()
