#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""插頭 DP（輪廓線 DP）：n×m 棋盤 1×2 骨牌完全覆蓋計數。

命令（IO 模式，從 stdin 讀入）：
  count n m   -> 輸出 n×m 棋盤被 1×2 骨牌完全覆蓋的方案數（n, m <= 12）

stdin 無內容時跑內建斷言測試，輸出 all tests passed。

核心要點（輪廓線 / 插頭 DP 的入門形態）：
  * 逐格掃描，狀態是長度 m 的 bitmask：進入格子 (i, j) 時
      - bit k (k < j)  = 第 i 行第 k 列向下伸的插頭（伸給 i+1 行）
      - bit k (k >= j) = 上方伸下來的插頭（(i-1,k) 豎放下伸 = (i,k) 的上插頭；
        (i,k-1) 橫放右伸 = (i,k) 的左插頭，編碼在同一個 bit 上）
    骨牌每格度數為 2，同一輪廓位至多一個插頭，單 bit 即可。
  * 轉移：輪廓位 j 已有插頭 -> 本格已被覆蓋，消費該插頭；
    無插頭 -> 本格必須被覆蓋：豎放（bit j 置 1，伸向下一行）
    或橫放（要求 bit j+1 為 0，置 bit j+1 = 1 給右鄰）。
  * 最後一格無法豎放/橫放時該狀態自然死亡；終態要求 mask == 0
    （殘留插頭意味著骨牌伸出棋盤）。
  * 複雜度 O(n·m·2^m)，12×12 在 long long 內約 5.8e12 種。
"""

import sys


def count_domino(n, m):
    """n×m 棋盤 1×2 骨牌完全覆蓋方案數（輪廓線逐格 DP）。"""
    if n <= 0 or m <= 0:
        return 0
    dp = {0: 1}
    for i in range(n):
        for j in range(m):
            ndp = {}
            for mask, cnt in dp.items():
                if (mask >> j) & 1:
                    # 本格已被左鄰/上鄰的骨牌佔住：消費插頭
                    nmask = mask & ~(1 << j)
                    ndp[nmask] = ndp.get(nmask, 0) + cnt
                else:
                    # 豎放：(i,j)-(i+1,j)，向 bit j 伸出
                    if i + 1 < n:
                        nmask = mask | (1 << j)
                        ndp[nmask] = ndp.get(nmask, 0) + cnt
                    # 橫放：(i,j)-(i,j+1)，bit j+1 必須空閒
                    if j + 1 < m and not (mask >> (j + 1)) & 1:
                        nmask = mask | (1 << (j + 1))
                        ndp[nmask] = ndp.get(nmask, 0) + cnt
            dp = ndp
    return dp.get(0, 0)


def brute_count(n, m):
    """暴力 DFS 填格計數（僅用於小棋盤對拍）。"""
    cells = n * m
    used = [False] * cells

    def dfs(pos):
        while pos < cells and used[pos]:
            pos += 1
        if pos == cells:
            return 1
        total = 0
        r, c = divmod(pos, m)
        used[pos] = True
        # 豎放
        if r + 1 < n and not used[pos + m]:
            used[pos + m] = True
            total += dfs(pos + 1)
            used[pos + m] = False
        # 橫放
        if c + 1 < m and not used[pos + 1]:
            used[pos + 1] = True
            total += dfs(pos + 1)
            used[pos + 1] = False
        used[pos] = False
        return total

    return dfs(0)


def run_io(data):
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
        if cmd == 'count':
            n = nxt_int()
            m = nxt_int()
            if n is None or m is None:
                break
            out.append(str(count_domino(n, m)))
        else:
            break
    return ''.join(line + '\n' for line in out)


def run_tests():
    # 1) 已知數列：1×m、2×n = Fib(n+1)、3×2k 數列
    assert count_domino(1, 1) == 0
    assert count_domino(1, 2) == 1
    assert count_domino(2, 1) == 1
    assert count_domino(1, 7) == 0
    assert count_domino(1, 8) == 1
    fib = [1, 1]
    for k in range(2, 13):
        fib.append(fib[-1] + fib[-2])
    for k in range(1, 11):
        assert count_domino(2, k) == fib[k], (k, count_domino(2, k))
    seq3 = [3, 11, 41, 153, 571]  # 3×2, 3×4, 3×6, 3×8, 3×10
    for idx, want in enumerate(seq3):
        assert count_domino(3, 2 * (idx + 1)) == want, (idx, want)

    # 2) 對稱性 count(n,m) == count(m,n)
    for n in range(1, 8):
        for m in range(1, 8):
            assert count_domino(n, m) == count_domino(m, n), (n, m)

    # 3) 格數為奇數必為 0
    for n in range(1, 7):
        for m in range(1, 7):
            if (n * m) % 2 == 1:
                assert count_domino(n, m) == 0, (n, m)

    # 4) 與暴力 DFS 對拍（小棋盤隨機 + 全覆蓋）
    state = 20261110

    def rnd(mm):
        nonlocal state
        state = (state * 1103515245 + 12345) & 0x7FFFFFFF
        return state % mm

    for _ in range(60):
        n = rnd(4) + 1
        m = rnd(4) + 1
        if n * m > 14:
            continue
        assert count_domino(n, m) == brute_count(n, m), (n, m)

    # 5) 大數值固定回歸（對照 OEIS A004003：n×n 骨牌鋪法數）
    assert count_domino(8, 8) == 12988816
    assert count_domino(10, 10) == 258584046368
    assert count_domino(12, 12) == 53060477521960000
    assert count_domino(7, 5) == 0  # 奇數格數（35 格）


def main():
    data = sys.stdin.read()
    if data.strip():
        sys.stdout.write(run_io(data))
    else:
        run_tests()
        sys.stdout.write('all tests passed\n')


if __name__ == '__main__':
    main()
