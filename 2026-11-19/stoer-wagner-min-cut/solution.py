#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Stoer–Wagner 全局最小割（Global Minimum Cut）。

在一個**無向帶權圖**中，找出權值和最小的割：把頂點分成兩堆，使跨越兩堆的邊權和最小。
Stoer–Wagner 演算法用 n-1 輪「最大權相鄰生長」來合併頂點，時間 O(n^3)。

## 應用
- 網路可靠性：衡量把通訊網路切成兩半所需的最小頻寬（最脆弱環節）。
- 影像分割：把像素視為圖、相似度視為邊權，最小割即最佳前景/背景切分。
- VLSI / 社群劃分：把電路或社交網路切成盡量少連線的兩塊。

## IO（stdin）
    n m
    u1 v1 w1
    ...          # 0-based 頂點，無向邊，權為正整數；重邊累加
    um vm wm
輸出一行：最小割的權值和（圖只有 0/1 個點時輸出 0）。

stdin 無內容時跑內建隨機對拍（與 2^(n-1) 枚舉所有割的暴力法比較），輸出 all tests passed。
"""

import sys
import re

_INT = re.compile(r'^[+-]?[0-9]+$')
INF = 1 << 60


def stoer_wagner(n, edges):
    """回傳無向圖的全局最小割權值。edges: [(u, v, w), ...]（0-based）。"""
    if n <= 1:
        return 0
    # 鄰接矩陣（重邊累加，自環忽略）
    w = [[0] * n for _ in range(n)]
    for u, v, c in edges:
        if u == v:
            continue
        w[u][v] += c
        w[v][u] += c

    active = [True] * n          # 是否仍未被合併
    best = INF
    remaining = n

    while remaining > 1:
        visited = [False] * n
        weight = [0] * n         # weight[x] = 目前已選集合 A 到 x 的邊權和
        order = []
        for _ in range(remaining):
            # 在未被合併、未選入 A 的點中取 weight 最大者（平手取編號最小，確保兩語言一致）
            pick = -1
            for v in range(n):
                if active[v] and not visited[v]:
                    if pick == -1 or weight[v] > weight[pick]:
                        pick = v
            visited[pick] = True
            order.append(pick)
            for u in range(n):
                if active[u] and not visited[u]:
                    weight[u] += w[pick][u]

        t = order[-1]            # 最後加入者
        s = order[-2]            # 倒數第二個加入者
        cut = weight[t]          # A(含 s) 與 {t} 之間的割
        if cut < best:
            best = cut

        # 把 t 合併進 s：s 吸收 t 到所有其它仍活躍點的邊權
        for u in range(n):
            if active[u] and u != s and u != t:
                w[s][u] += w[t][u]
                w[u][s] += w[u][t]
        active[t] = False
        remaining -= 1

    return best


def brute_min_cut(n, edges):
    """暴力：固定頂點 0 必在一側，列舉其餘 2^(n-1) 種劃分。僅供小圖對拍。"""
    if n <= 1:
        return 0
    w = [[0] * n for _ in range(n)]
    for u, v, c in edges:
        if u != v:
            w[u][v] += c
            w[v][u] += c
    ans = INF
    others = list(range(1, n))
    for mask in range((1 << len(others)) - 1):  # 排除全在一側的無效劃分
        in_side = [False] * n
        in_side[0] = True
        for i, v in enumerate(others):
            if mask >> i & 1:
                in_side[v] = True
        cut = 0
        for u in range(n):
            for v in range(u + 1, n):
                if in_side[u] != in_side[v]:
                    cut += w[u][v]
        if cut < ans:
            ans = cut
    return ans


def run_io(data):
    toks = data.split()
    pos = 0

    def nxt():
        nonlocal pos
        if pos >= len(toks):
            return None
        t = toks[pos]
        pos += 1
        return t

    tn = nxt()
    if tn is None or not _INT.match(tn):
        return ''
    n = min(max(int(tn), 0), 200)
    tm = nxt()
    if tm is None or not _INT.match(tm):
        return ''
    m = min(max(int(tm), 0), 20000)
    edges = []
    for _ in range(m):
        a, b, c = nxt(), nxt(), nxt()
        if a is None or b is None or c is None:
            break
        if not (_INT.match(a) and _INT.match(b) and _INT.match(c)):
            break
        u, v, w = int(a), int(b), int(c)
        if not (0 <= u < n and 0 <= v < n):
            continue
        edges.append((u, v, min(max(w, 0), 10**9)))
    return str(stoer_wagner(n, edges)) + '\n'


def run_tests():
    # 固定情境（先實跑確認）
    e = [(0, 1, 2), (0, 2, 3), (1, 2, 1), (1, 3, 4), (2, 3, 5)]
    assert stoer_wagner(4, e) == 5
    assert brute_min_cut(4, e) == 5
    assert stoer_wagner(1, []) == 0
    assert stoer_wagner(0, []) == 0
    assert stoer_wagner(2, [(0, 1, 7)]) == 7
    assert stoer_wagner(3, [(0, 1, 1), (1, 2, 1)]) == 1

    # 隨機對拍
    state = 20261119

    def rnd(m):
        nonlocal state
        state = (state * 1103515245 + 12345) & 0x7FFFFFFF
        return state % m

    for tc in range(300):
        n = rnd(9) + 2                       # 2..10
        edges = []
        for u in range(n):
            for v in range(u + 1, n):
                if rnd(3) == 0:
                    edges.append((u, v, rnd(20) + 1))
        got = stoer_wagner(n, edges)
        want = brute_min_cut(n, edges)
        assert got == want, (tc, n, edges, got, want)


def main():
    data = sys.stdin.read()
    if data.strip():
        sys.stdout.buffer.write(run_io(data).encode('utf-8'))
    else:
        run_tests()
        sys.stdout.buffer.write(b'all tests passed\n')


if __name__ == '__main__':
    main()
