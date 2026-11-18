#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Berlekamp-Massey 最短線性遞迴 + Kitamasa 快速求第 n 項（mod 質數 p）。

IO 模式（stdin）：
  第一行 p k n；第二行 k 個整數 s0..s(k-1)。可讀多筆測資直到 EOF。
  每筆輸出 2 行：
    1) L c1 c2 ... cL   —— 最短遞迴長度 L 與係數，使 s[t] = Σ c_i·s[t-i]
    2) 第 n 項（0-based）mod p
  p 須為質數（Kitamasa 需要求 c_L 的逆元）。
無 stdin 時跑內建對拍（vs 直接 DP 展開 + 暴力最短遞迴搜尋），輸出 all tests passed。

複雜度：BM O(k²)，Kitamasa O(L² log n)。
"""

import sys

# ----------------------------- Berlekamp-Massey -----------------------------

def berlekamp_massey(s, p):
    """回傳最短線性遞迴係數 c[1..L]，使 s[n] = Σ_{i=1..L} c[i]·s[n-i] (mod p)。"""
    n = len(s)
    C = [1]           # C[0]=1, C[i] = -c[i]
    B = [1]
    L = 0
    m = 1
    b = 1
    for N in range(n):
        # d = s[N] + Σ_{i=1..L} C[i]·s[N-i]
        d = s[N] % p
        for i in range(1, len(C)):
            d = (d + C[i] * s[N - i]) % p
        if d == 0:
            m += 1
            continue
        T = C[:]
        coef = d * pow(b, p - 2, p) % p
        if len(C) < len(B) + m:
            C.extend([0] * (len(B) + m - len(C)))
        for i in range(len(B)):
            C[i + m] = (C[i + m] - coef * B[i]) % p
        if 2 * L <= N:
            L = N + 1 - L
            B = T
            b = d
            m = 1
        else:
            m += 1
    # c[i] = -C[i] mod p
    return [(-C[i]) % p for i in range(1, L + 1)]


# ------------------------------- Kitamasa -------------------------------

def poly_mul_mod(a, b, c, p):
    """(a·b) mod (x^L - c1 x^(L-1) - ... - cL)。c 為 1-index 係數列表。"""
    L = len(c)
    tmp = [0] * (2 * L - 1)
    for i in range(len(a)):
        if a[i] == 0:
            continue
        for j in range(len(b)):
            tmp[i + j] = (tmp[i + j] + a[i] * b[j]) % p
    # 從高次往下降：x^d (d>=L) = Σ_{i=1..L} c[i]·x^(d-i)
    for d in range(len(tmp) - 1, L - 1, -1):
        val = tmp[d]
        if val == 0:
            continue
        for i in range(1, L + 1):
            tmp[d - i] = (tmp[d - i] + val * c[i - 1]) % p
    return tmp[:L]


def kitamasa(s, c, n, p):
    """求序列第 n 項（0-based）。s 為前 L 項，c 為遞迴係數（長度 L）。"""
    L = len(c)
    if L == 0:
        return 0
    if n < L:
        return s[n] % p
    # 多項式 x^n mod 特徵多項式
    result = [1] + [0] * (L - 1)   # 多項式 1
    base = [0] * L
    if L == 1:
        base[0] = c[0] % p          # 多項式 x ≡ c1 (mod x - c1)
    else:
        base[1] = 1                 # 多項式 x
    e = n
    while e > 0:
        if e & 1:
            result = poly_mul_mod(result, base, c, p)
        base = poly_mul_mod(base, base, c, p)
        e >>= 1
    return sum(result[i] * s[i] for i in range(L)) % p


# ----------------------------- 暴力參考 -----------------------------

def brute_first_fit(s, p):
    """暴力找最短遞迴（僅小 k）：對 L=0..，解線性方程組（枚舉係數）。"""
    k = len(s)
    for L in range(0, k // 2 + 1):
        # 檢查是否存在 c[1..L] 使所有 t>=L 滿足；L=0 表示全為 0
        if L == 0:
            if all(x % p == 0 for x in s):
                return 0
            continue
        # 用前 L 個方程解 c（高斯消元），再驗證其餘
        import itertools
        found = None
        # 直接枚舉 c（p 小、L 小才可行）
        for coeffs in itertools.product(range(p), repeat=L):
            ok = True
            for t in range(L, k):
                v = sum(coeffs[i] * s[t - 1 - i] for i in range(L)) % p
                if v != s[t] % p:
                    ok = False
                    break
            if ok:
                found = list(coeffs)
                break
        if found is not None:
            return L
    return k // 2


def dp_nth(s, c, n, p):
    """直接 DP 展開到第 n 項（n 小時用）。"""
    L = len(c)
    if L == 0:
        return 0
    seq = list(s)
    while len(seq) <= n:
        nxt = sum(c[i] * seq[-1 - i] for i in range(L)) % p
        seq.append(nxt)
    return seq[n] % p


# ------------------------------- IO 模式 -------------------------------

def run_io(data):
    toks = data.split()
    pos = 0
    out = []

    def nxt():
        nonlocal pos
        if pos >= len(toks):
            return None
        t = toks[pos]
        pos += 1
        return t

    while True:
        a, b, d = nxt(), nxt(), nxt()
        if a is None or b is None or d is None:
            break
        try:
            p, k, n = int(a), int(b), int(d)
        except ValueError:
            break
        if not (2 <= p <= 2000000000 and 0 <= k <= 5000 and 0 <= n <= 10**18):
            break
        s = []
        ok = True
        for _ in range(k):
            v = nxt()
            if v is None:
                ok = False
                break
            try:
                s.append(int(v) % p)
            except ValueError:
                ok = False
                break
        if not ok:
            break
        c = berlekamp_massey(s, p)
        ans = kitamasa(s, c, n, p)
        out.append(' '.join(map(str, [len(c)] + c)))
        out.append(str(ans))
    return ''.join(line + '\n' for line in out)


# ------------------------------ 內建測試 ------------------------------

def run_tests():
    P = 1000000007
    # 固定情境：s = 1,1,2,3,5,8,...（Fibonacci），遞迴 c=[1,1]
    fib = [1, 1]
    for _ in range(20):
        fib.append((fib[-1] + fib[-2]) % P)
    c = berlekamp_massey(fib[:12], P)
    assert c == [1, 1], c
    assert kitamasa(fib[:12], c, 10, P) == 89
    assert kitamasa(fib[:12], c, 100, P) == dp_nth(fib[:12], c, 100, P)

    # s = 2,3,5,7,11,13... 用已知二階遞迴 s[n]=3 s[n-1]-2 s[n-2]
    seq = [2, 3]
    for _ in range(30):
        seq.append((3 * seq[-1] - 2 * seq[-2]) % P)
    c2 = berlekamp_massey(seq[:14], P)
    assert c2 == [3, (P - 2) % P], c2
    assert kitamasa(seq[:14], c2, 500, P) == dp_nth(seq[:14], c2, 500, P)

    # 隨機對拍：隨機係數生成序列，BM 還原，Kitamasa vs DP
    state = 20261118

    def rnd(m):
        nonlocal state
        state = (state * 1103515245 + 12345) & 0x7FFFFFFF
        return state % m

    smallp = 7
    for tc in range(300):
        L = rnd(4) + 1
        ctrue = [rnd(smallp) for _ in range(L)]
        if ctrue[-1] == 0:
            ctrue[-1] = 1
        s = [rnd(smallp) for _ in range(L)]
        for _ in range(2 * L + 8):
            s.append(sum(ctrue[i] * s[-1 - i] for i in range(L)) % smallp)
        got = berlekamp_massey(s, smallp)
        # BM 找回的遞迴長度 ≤ L，且能重現所有項
        assert len(got) <= L + 1, (got, L)
        for t in range(len(got), len(s)):
            v = sum(got[i] * s[t - 1 - i] for i in range(len(got))) % smallp
            assert v == s[t] % smallp, (got, s, t)
        # 最短性：暴力確認不存在更短的
        assert len(got) == brute_first_fit(s, smallp), (got, s, brute_first_fit(s, smallp))
        for n in (0, 1, L, 2 * L + 3):
            assert kitamasa(s, got, n, smallp) == dp_nth(s, got, n, smallp)

    # 全零序列 -> L=0
    assert berlekamp_massey([0, 0, 0, 0], P) == []
    assert berlekamp_massey([], P) == []


def main():
    data = sys.stdin.read()
    if data.strip():
        sys.stdout.buffer.write(run_io(data).encode('utf-8'))
    else:
        run_tests()
        sys.stdout.buffer.write(b'all tests passed\n')


if __name__ == '__main__':
    main()
