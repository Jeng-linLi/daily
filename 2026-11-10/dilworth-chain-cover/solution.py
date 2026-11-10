#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Dilworth 定理：二維偏序集的最小鏈覆蓋 = 最長反鏈。

命令（IO 模式，從 stdin 讀入）：
  chain n x1 y1 x2 y2 ... xn yn
      -> 輸出「最小鏈覆蓋 最長反鏈」兩個數（依 Dilworth 定理恆相等）

偏序定義：(a, b) ≼ (c, d) 當且僅當 a <= c 且 b <= d（支配序）。
  鏈 = 兩兩可比的點集；反鏈 = 兩兩互不支配的點集。
stdin 無內容時跑內建斷言測試，輸出 all tests passed。

核心要點：
  * 去重後按 a 升序排序，則反鏈 = b 的最長嚴格遞減子序列（LDS，O(n^2) DP）。
  * 最小鏈覆蓋轉化為 DAG 最小路徑覆蓋：點 i 拆成出/入兩側，邊 i->j 存在
    當 i < j 且 b_i <= b_j；最小鏈覆蓋 = n − 最大匹配（匈牙利算法求）。
  * 兩個完全獨立的算法（LDS vs 匹配）結果恆相等——這正是 Dilworth 定理
    （最小鏈覆蓋 = 最長反鏈）的構造性驗證。
"""

import sys


def dedup_sort(points):
    """去重 + 按 (a, b) 升序排序；去重後 a 互異。"""
    return sorted(set(points))


def longest_antichain(pts):
    """最長反鏈：按 a 升序後 b 的最長嚴格遞減子序列。"""
    n = len(pts)
    if n == 0:
        return 0
    dp = [1] * n  # dp[i] = 以 pts[i] 結尾的 LDS 長度
    for i in range(n):
        for j in range(i):
            if pts[j][1] > pts[i][1]:  # 嚴格遞減
                dp[i] = max(dp[i], dp[j] + 1)
    return max(dp)


def hungarian(n, adj):
    """二分圖最大匹配（匈牙利算法，左側 n 個點，adj[i] = 右側鄰居列表）。"""
    match_r = [-1] * n  # 右側點 -> 左側點

    def try_k(i, vis):
        for j in adj[i]:
            if vis[j]:
                continue
            vis[j] = True
            if match_r[j] == -1 or try_k(match_r[j], vis):
                match_r[j] = i
                return True
        return False

    cnt = 0
    for i in range(n):
        if try_k(i, [False] * n):
            cnt += 1
    return cnt


def min_chain_cover(pts):
    """最小鏈覆蓋 = n − 最大匹配（DAG 最小路徑覆蓋，點鏈版本）。"""
    n = len(pts)
    if n == 0:
        return 0
    adj = [[] for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            if pts[i][1] <= pts[j][1]:  # a 已升序，i < j 且 b_i <= b_j 即支配
                adj[i].append(j)
    return n - hungarian(n, adj)


def solve(points):
    pts = dedup_sort(points)
    return min_chain_cover(pts), longest_antichain(pts)


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
        if cmd == 'chain':
            n = nxt_int()
            if n is None:
                break
            pts = []
            ok = True
            for _ in range(n):
                x = nxt_int()
                y = nxt_int()
                if x is None or y is None:
                    ok = False
                    break
                pts.append((x, y))
            if not ok:
                break
            cover, anti = solve(pts)
            out.append('%d %d' % (cover, anti))
        else:
            break
    return ''.join(line + '\n' for line in out)


def run_tests():
    # 1) 全鏈（已排序支配鏈）-> 覆蓋 1、反鏈 1
    assert solve([(1, 1), (2, 2), (3, 3)]) == (1, 1)
    # 2) 全反鏈（a 增 b 嚴格減）-> 覆蓋 = 反鏈 = n
    assert solve([(1, 5), (2, 4), (3, 3), (4, 2), (5, 1)]) == (5, 5)
    # 3) 手驗混合例：(1,1)(1,3) 相同 a 成鏈；(2,2) 與 (1,1) 可比、與 (1,3) 不可比
    #    點集 {(1,1),(1,3),(2,2),(3,4)}：反鏈 {(1,3),(2,2)} 大小 2，覆蓋 2
    assert solve([(1, 1), (1, 3), (2, 2), (3, 4)]) == (2, 2)
    # 4) 重複點：去重後一致
    assert solve([(1, 1), (1, 1), (2, 2), (1, 1)]) == solve([(1, 1), (2, 2)]) == (1, 1)
    # 5) 空集
    assert solve([]) == (0, 0)

    state = 20261110

    def rnd(m):
        nonlocal state
        state = (state * 1103515245 + 12345) & 0x7FFFFFFF
        return state % m

    # 6) Dilworth 定理驗證：LDS（反鏈）與匹配（覆蓋）兩獨立算法恆相等
    for _ in range(200):
        n = rnd(20) + 1
        pts = [(rnd(15) + 1, rnd(15) + 1) for _ in range(n)]
        cover, anti = solve(pts)
        assert cover == anti, (pts, cover, anti)

    # 7) LDS 與暴力子集枚舉對拍（兩兩互不支配的最大子集）
    def brute_anti(pts):
        n = len(pts)
        best = 0
        for mask in range(1 << n):
            sel = [pts[i] for i in range(n) if mask >> i & 1]
            ok = all(not (p[0] <= q[0] and p[1] <= q[1]) and
                     not (q[0] <= p[0] and q[1] <= p[1])
                     for i, p in enumerate(sel) for q in sel[i + 1:])
            if ok:
                best = max(best, len(sel))
        return best

    for _ in range(50):
        n = rnd(11) + 2
        pts = dedup_sort([(rnd(6), rnd(6)) for _ in range(n)])
        assert longest_antichain(pts) == brute_anti(pts), pts

    # 8) 匈牙利與暴力最大匹配對拍
    def brute_match(n, adj):
        best = 0
        rmatch = [-1] * n

        def go(i):
            nonlocal best
            if i == n:
                best = max(best, sum(1 for x in rmatch if x != -1))
                return
            go(i + 1)
            for j in adj[i]:
                if rmatch[j] == -1:
                    rmatch[j] = i
                    go(i + 1)
                    rmatch[j] = -1

        go(0)
        return best

    for _ in range(50):
        n = rnd(9) + 1
        pts = dedup_sort([(rnd(8), rnd(8)) for _ in range(n)])
        nn = len(pts)
        adj = [[] for _ in range(nn)]
        for i in range(nn):
            for j in range(i + 1, nn):
                if pts[i][1] <= pts[j][1]:
                    adj[i].append(j)
        assert hungarian(nn, adj) == brute_match(nn, adj), pts


def main():
    data = sys.stdin.read()
    if data.strip():
        sys.stdout.write(run_io(data))
    else:
        run_tests()
        sys.stdout.write('all tests passed\n')


if __name__ == '__main__':
    main()
