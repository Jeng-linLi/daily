#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Steiner Tree（斯坦纳树）— 最小連通終端集子圖

題意：
    給定一張帶非負邊權的無向圖，以及 k 個「終端節點」(terminals)。
    求一棵連通子圖（樹），使得所有終端都被連通，且總邊權和最小。
    這棵樹可以額外包含非終端節點（所謂 Steiner 點），以換取更小權重。

思路（狀壓 DP + Dijkstra）：
    把 k 個終端編號 0..k-1，定義
        dp[mask][v] = 連通「mask 所指終端集合」、且經過頂點 v 的最小權值。
    1) 初始化：dp[1<<i][term[i]] = 0（第 i 個終端自身即可連通自己）。
    2) 枚舉 mask（由小到大）：
       a. 子樹合併：在頂點 v 把兩個互不相交的子集合併
          dp[mask][v] = min_{sub} dp[sub][v] + dp[mask^sub][v]
          （sub 為 mask 的真子集；sub 與 mask^sub 都已先算完）
       b. 邊鬆弛：以 dp[mask][·] 為初始距離跑 Dijkstra，
          沿邊 (v,u,w) 更新 dp[mask][u] = min(..., dp[mask][v]+w)
    3) 答案 = min_v dp[(1<<k)-1][v]；若全為無窮則圖不連通返回 -1。

時間複雜度：O(3^k · n + 2^k · m · log n)（子集合併 3^k、每個 mask 一次 Dijkstra）。
空間複雜度：O(2^k · n)。

輸入格式（stdin 模式）：
    n m k
    接下來 m 行：u v w   （0-based，無向邊）
    最後 1 行：k 個終端頂點
    輸出：最小總權值，不可連通則 -1。
無 stdin 時跑內建斷言並印出 `all tests passed`。
"""

import sys
import heapq

INF = float("inf")


def steiner_tree(n, edges, terminals):
    """回傳連通所有 terminals 的最小總權值；不可連通回傳 INF。"""
    k = len(terminals)
    adj = [[] for _ in range(n)]
    for u, v, w in edges:
        adj[u].append((v, w))
        adj[v].append((u, w))

    full = 1 << k
    # dp[mask][v]
    dp = [[INF] * n for _ in range(full)]

    # 初始化：每個終端自身可連通自己（代價 0）
    for i, t in enumerate(terminals):
        dp[1 << i][t] = 0

    # 按 mask 由小到大處理（保證子集合併時子狀態已算完）
    for mask in range(1, full):
        # (a) 子樹合併：在頂點 v 把兩個互斥子集拼起來
        sub = (mask - 1) & mask
        while sub:
            other = mask ^ sub
            for v in range(n):
                if dp[sub][v] < INF and dp[other][v] < INF:
                    cand = dp[sub][v] + dp[other][v]
                    if cand < dp[mask][v]:
                        dp[mask][v] = cand
            sub = (sub - 1) & mask

        # (b) 邊鬆弛：以當前 dp[mask] 為初始距離跑 Dijkstra
        pq = [(dp[mask][v], v) for v in range(n) if dp[mask][v] < INF]
        heapq.heapify(pq)
        while pq:
            d, v = heapq.heappop(pq)
            if d > dp[mask][v]:
                continue
            for to, w in adj[v]:
                nd = d + w
                if nd < dp[mask][to]:
                    dp[mask][to] = nd
                    heapq.heappush(pq, (nd, to))

    return min(dp[full - 1])


def solve_io(text=None):
    if text is None:
        text = sys.stdin.read()
    data = text.strip().split()
    if not data:
        return
    it = iter(data)
    n = int(next(it))
    m = int(next(it))
    k = int(next(it))
    edges = []
    for _ in range(m):
        u = int(next(it))
        v = int(next(it))
        w = int(next(it))
        edges.append((u, v, w))
    terminals = [int(next(it)) for _ in range(k)]
    ans = steiner_tree(n, edges, terminals)
    print(-1 if ans >= INF else ans)


def _brute(n, edges, terminals):
    """指數級暴力（窮舉邊子集求最小連通終端權值），僅供測試比對。"""
    m = len(edges)
    best = INF
    # 建無向鄰接表
    adj = [[] for _ in range(n)]
    for u, v, w in edges:
        adj[u].append(v)
        adj[v].append(u)
    from itertools import combinations

    def connected_with(used_edges):
        # 用 used_edges（邊索引集合）建圖，檢查 terminals 是否連通
        parent = list(range(n))

        def find(x):
            while parent[x] != x:
                parent[x] = parent[parent[x]]
                x = parent[x]
            return x

        for ei in used_edges:
            u, v, _ = edges[ei]
            pu, pv = find(u), find(v)
            if pu != pv:
                parent[pu] = pv
        roots = {find(t) for t in terminals}
        return len(roots) == 1

    for r in range(1, m + 1):
        for combo in combinations(range(m), r):
            if connected_with(combo):
                w = sum(edges[ei][2] for ei in combo)
                if w < best:
                    best = w
    return best


def _run_tests():
    # 測試 1：基礎題（預期 5）
    n = 5
    edges = [
        (0, 1, 1),
        (1, 2, 1),
        (0, 3, 4),
        (3, 4, 1),
        (4, 2, 4),
    ]
    terminals = [2, 3, 4]
    ans = steiner_tree(n, edges, terminals)
    assert ans == 5, ans
    assert _brute(n, edges, terminals) == 5

    # 測試 2：不連通（預期 -1 / INF）
    n = 4
    edges = [(0, 1, 1), (2, 3, 2)]
    terminals = [0, 3]
    ans = steiner_tree(n, edges, terminals)
    assert ans >= INF, ans
    assert _brute(n, edges, terminals) >= INF

    # 測試 3：終端恰為全部節點的路徑（預期最小生成樹 = 所有邊和）
    n = 4
    edges = [(0, 1, 2), (1, 2, 3), (2, 3, 5)]
    terminals = [0, 1, 2, 3]
    ans = steiner_tree(n, edges, terminals)
    assert ans == 10, ans
    assert _brute(n, edges, terminals) == 10

    # 測試 4：需要 Steiner 點（預期 4）
    # 0-1(1) 1-2(1) 2-3(1) 3-4(1) 4-5(1)，終端 {0,2,5} -> 取 0-1-2 (2) + 2-3-4-5(3) = 5?
    # 改成：0-1(1),1-2(1),0-3(1),3-4(1),4-2(1),終端{2,4} 直接 4-2=1
    n = 5
    edges = [(0, 1, 1), (1, 2, 1), (0, 3, 1), (3, 4, 1), (4, 2, 1)]
    terminals = [2, 4]
    ans = steiner_tree(n, edges, terminals)
    assert ans == 1, ans  # 邊 4-2 直接連
    assert _brute(n, edges, terminals) == 1

    # 測試 5：用 Steiner 點省錢（預期 4）
    # 終端 {0,3}，若只走 0-1-2-3 = 1+5+1=7；走 Steiner 點繞路：0-1(1)1-4(1)4-2(1)2-3(1)=4
    n = 5
    edges = [(0, 1, 1), (1, 2, 5), (2, 3, 1),
             (1, 4, 1), (4, 2, 1)]
    terminals = [0, 3]
    ans = steiner_tree(n, edges, terminals)
    assert ans == 4, ans
    assert _brute(n, edges, terminals) == 4

    print("all tests passed")


if __name__ == "__main__":
    data = sys.stdin.read()
    if data.strip():
        solve_io(data)        # 有輸入走 IO 模式
    else:
        _run_tests()          # 無 stdin（如 </dev/null）跑內建斷言
