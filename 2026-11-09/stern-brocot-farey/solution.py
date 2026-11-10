#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Stern-Brocot 樹與 Farey 數列。

三種命令（IO 模式，從 stdin 讀入）：
  farey n        -> 輸出 F_n 的全部既約分數（升冪，空格分隔）
  path p q       -> 輸出 Stern-Brocot 樹中從 1/1 走到 p/q 的壓縮路徑（如 L3 R5）
  cmp a b c d    -> 精確比較 a/b 與 c/d，輸出 -1 / 0 / 1

stdin 無內容時跑內建斷言測試，輸出 all tests passed。

核心要點：
  * Farey 相鄰項遞推：當前項 a/b、前一項 c/d，則後繼 = (k*a - c)/(k*b - d)，
    其中 k = floor((n + d) / b)。O(|F_n|) 而非 O(n^2 log n) 排序。
  * Farey 鄰域性質：相鄰兩項 a/b < c/d 滿足 b*c - a*d = 1。
  * |F_n| = 1 + sum_{m=1..n} phi(m)。
  * Stern-Brocot walk：區間 [a/b, c/d] 從哨兵 0/1 與 1/0 出發，
    mediant = (a+c)/(b+d)；p/q 偏小走 L、偏大走 R。由於每步 mediant
    分母嚴格遞增，步數上界為 p+q（約簡後）。用整數交叉相乘比較，
    全程精確無誤差。
"""

import sys
from math import gcd


def farey(n):
    """Farey 數列 F_n：所有分母 <= n 的既約分數，升冪返回。"""
    if n <= 0:
        return []
    # 初始兩項 0/1, 1/n，c/d 為前一項
    res = [(0, 1), (1, n)]
    c, d = 0, 1
    a, b = 1, n
    while (a, b) != (1, 1):
        k = (n + d) // b
        a, b, c, d = k * a - c, k * b - d, a, b
        res.append((a, b))
    return res


def stern_brocot_path(p, q):
    """從根 mediant 1/1 出發到 p/q 的 L/R 路徑（未壓縮的字元列表）。

    p=0 或 q=0（哨兵 0/1 / 1/0）以及 p/q = 1/1（根本身）回傳空路徑。
    """
    if p <= 0 or q <= 0:
        return []
    g = gcd(p, q)
    p //= g
    q //= g
    if p == q:
        return []
    a, b, c, d = 0, 1, 1, 0  # 區間 [a/b, c/d]，c/d = 1/0 = 無窮大
    steps = []
    while True:
        mn, md = a + c, b + d  # mediant
        lhs, rhs = p * md, q * mn
        if lhs == rhs:
            break
        if lhs < rhs:  # p/q < mediant -> 左子樹
            steps.append('L')
            c, d = mn, md
        else:  # p/q > mediant -> 右子樹
            steps.append('R')
            a, b = mn, md
    return steps


def rle(steps):
    """把 L/R 序列壓縮成連續段，如 L L R R R -> 'L2 R3'。"""
    out = []
    i = 0
    while i < len(steps):
        j = i
        while j < len(steps) and steps[j] == steps[i]:
            j += 1
        out.append('%s%d' % (steps[i], j - i))
        i = j
    return ' '.join(out)


def cmp_frac(a, b, c, d):
    """精確比較 a/b 與 c/d：交叉相乘 a*d vs b*c。"""
    lhs, rhs = a * d, b * c
    if lhs < rhs:
        return -1
    if lhs > rhs:
        return 1
    return 0


def run_io(data):
    """IO 模式：逐命令處理，參數不足或遇到未知命令即停止。"""
    tokens = data.split()
    pos = 0
    out = []

    def nxt():
        nonlocal pos
        if pos >= len(tokens):
            return None
        t = tokens[pos]
        pos += 1
        return t

    def nxt_int():
        t = nxt()
        if t is None:
            return None
        try:
            return int(t)
        except ValueError:
            return None

    while True:
        cmd = nxt()
        if cmd is None:
            break
        if cmd == 'farey':
            n = nxt_int()
            if n is None:
                break
            out.append(' '.join('%d/%d' % x for x in farey(n)))
        elif cmd == 'path':
            p = nxt_int()
            q = nxt_int()
            if p is None or q is None:
                break
            out.append(rle(stern_brocot_path(p, q)))
        elif cmd == 'cmp':
            a = nxt_int()
            b = nxt_int()
            c = nxt_int()
            d = nxt_int()
            if None in (a, b, c, d):
                break
            out.append(str(cmp_frac(a, b, c, d)))
        else:
            break
    return ''.join(line + '\n' for line in out)


def run_tests():
    state = 20261109

    def rnd(m):
        nonlocal state
        state = (state * 1103515245 + 12345) & 0x7FFFFFFF
        return state % m

    def phi(m):
        return sum(1 for x in range(1, m + 1) if gcd(x, m) == 1)

    # 1) Farey 與暴力枚舉排序一致
    for n in range(1, 41):
        brute = sorted({(a, b) for b in range(1, n + 1)
                        for a in range(b + 1) if gcd(a, b) == 1},
                       key=lambda t: t[0] / t[1])
        f = farey(n)
        assert f == brute, (n, f[:8], brute[:8])

    # 2) 鄰域性質：相鄰兩項 b*c - a*d = 1
    for n in range(1, 41):
        f = farey(n)
        for (a, b), (c, d) in zip(f, f[1:]):
            assert b * c - a * d == 1, (n, a, b, c, d)

    # 3) |F_n| = 1 + sum(phi(m))
    for n in range(1, 61):
        assert len(farey(n)) == 1 + sum(phi(m) for m in range(1, n + 1))

    # 4) path 重放：終點 mediant 恰為 p/q；步數 <= p+q；rle 可無損還原
    for p in range(1, 61):
        for q in range(1, 61):
            g = gcd(p, q)
            pp, qq = p // g, q // g
            steps = stern_brocot_path(p, q)
            a, b, c, d = 0, 1, 1, 0
            for s in steps:
                if s == 'L':
                    c, d = a + c, b + d
                else:
                    a, b = a + c, b + d
            assert (a + c) * qq == (b + d) * pp, (p, q)
            assert len(steps) <= pp + qq, (p, q, len(steps))
            dec = []
            for run in (rle(steps).split() if steps else []):
                dec.extend([run[0]] * int(run[1:]))
            assert dec == steps, (p, q, rle(steps))

    # 5) 對稱性：p/q 與 q/p 的路徑 L/R 互換
    for p in range(1, 41):
        for q in range(1, 41):
            sp = stern_brocot_path(p, q)
            sq = stern_brocot_path(q, p)
            assert ['R' if s == 'L' else 'L' for s in sp] == sq, (p, q)

    # 6) cmp 與 Fraction（含無窮哨兵 1/0）一致
    from fractions import Fraction
    for _ in range(400):
        b = rnd(51)
        a = rnd(50) + 1 if b == 0 else rnd(51)
        d = rnd(51)
        c = rnd(50) + 1 if d == 0 else rnd(51)
        got = cmp_frac(a, b, c, d)
        v1 = float('inf') if b == 0 else Fraction(a, b)
        v2 = float('inf') if d == 0 else Fraction(c, d)
        exp = (v1 > v2) - (v1 < v2)
        assert got == exp, (a, b, c, d, got, exp)

    # 7) 固定輸出用例（實跑結果作斷言，防回歸）
    assert ' '.join('%d/%d' % x for x in farey(5)) == '0/1 1/5 1/4 1/3 2/5 1/2 3/5 2/3 3/4 4/5 1/1'
    assert rle(stern_brocot_path(5, 3)) == 'R1 L1 R1'
    assert rle(stern_brocot_path(3, 5)) == 'L1 R1 L1'
    assert rle(stern_brocot_path(1, 1)) == ''
    assert rle(stern_brocot_path(0, 7)) == ''
    assert cmp_frac(2, 4, 1, 2) == 0
    assert cmp_frac(22, 7, 3, 1) == 1
    assert cmp_frac(1, 0, 999999, 1) == 1


def main():
    data = sys.stdin.read()
    if data.strip():
        sys.stdout.write(run_io(data))
    else:
        run_tests()
        sys.stdout.write('all tests passed\n')


if __name__ == '__main__':
    main()
