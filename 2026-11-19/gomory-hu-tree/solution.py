#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Gomory–Hu 樹（等價流樹 / Equivalent Flow Tree）。

對一個 n 點無向帶權圖，構造一棵 n-1 邊的樹，滿足：任意兩點 s,t 在**原圖**的最小割，
等於它們在樹上路徑的最小邊權。因此 n-1 次最小割就能回答所有點對的最小割。

## 應用
- 大規模點對最大流查詢：預先 O(n·maxflow) 建樹，之後每次查詢 O(log n)。
- 網路可靠性分析、層次聚類（最小割樹狀圖）。

## IO（stdin）
    n m q
    u1 v1 w1
    ...              # 0-based 無向邊，權為正整數
    um vm wm
    s1 t1            # q 個查詢
    ...
    sq tq
輸出 q 行：每行 s,t 之間的最小割值（等於樹上路徑最小邊；s==t 輸出 0）。

stdin 無內容時跑內建測試：與 Edmonds-Karp 逐查詢暴力最小割對拍，輸出 all tests passed。
"""

import sys
import re
from collections import deque

_INT = re.compile(r'^[+-]?[0-9]+$')
INF = 1 << 60


def maxflow_mincut_side(n, adj, s, t):
    """Edmonds-Karp 求 s->t 最大流，回傳 (flow, s 側的可達點集合)。adj: 剩餘容量鄰接表。"""
    if s == t:
        return 0, None
    cap = [row[:] for row in adj]
    flow = 0
    while True:
        prev = [-1] * n
        prev[s] = -2
        q = deque([s])
        while q and prev[t] == -1:
            v = q.popleft()
            for u, c in enumerate(cap[v]):
                if c > 0 and prev[u] == -1:
                    prev[u] = v
                    q.append(u)
                    if u == t:
                        break
        if prev[t] == -1:
            break
        # 找瓶頸
        aug = INF
        v = t
        while v != s:
            p = prev[v]
            aug = min(aug, cap[p][v])
            v = p
        v = t
        while v != s:
            p = prev[v]
            cap[p][v] -= aug
            cap[v][p] += aug
            v = p
        flow += aug
    # s 側 = 剩餘網路中從 s 可達
    seen = [False] * n
    seen[s] = True
    q = deque([s])
    while q:
        v = q.popleft()
        for u, c in enumerate(cap[v]):
            if c > 0 and not seen[u]:
                seen[u] = True
                q.append(u)
    return flow, seen


def gomory_hu(n, edges):
    """構造 Gomory-Hu 樹，回傳 (parent, weight)。parent[0] = -1。"""
    adj = [[0] * n for _ in range(n)]
    for u, v, w in edges:
        if u == v:
            continue
        adj[u][v] += w
        adj[v][u] += w

    parent = [0] * n
    parent[0] = -1
    for s in range(1, n):
        t = parent[s]
        flow, side = maxflow_mincut_side(n, adj, s, t)
        # 更新父親：目前掛在 t 下、且落在 s 側的點改掛到 s
        for i in range(s + 1, n):
            if parent[i] == t and side is not None and side[i]:
                parent[i] = s
        # 記錄邊權 (s, t)；此處用一維 weight[s] 存 s 與其父的邊權
        # 稍後統一
    return parent, adj


def build_tree_weights(n, edges):
    """完整 Gomory-Hu：回傳 parent[] 與 wpar[]（wpar[v] = v-parent[v] 邊權）。"""
    adj = [[0] * n for _ in range(n)]
    for u, v, wt in edges:
        if u == v:
            continue
        adj[u][v] += wt
        adj[v][u] += wt

    parent = [-1] * n
    wpar = [0] * n
    # 用並查集維護「目前樹中被 s 側含蓋」的合併，簡化父親改掛
    # 採用標準 Gomory-Hu：parent[i] 初始值 0
    parent = [0] * n
    parent[0] = -1
    for s in range(1, n):
        t = parent[s]
        flow, side = maxflow_mincut_side(n, adj, s, t)
        wpar[s] = flow
        for i in range(s + 1, n):
            if parent[i] == t and side is not None and side[i]:
                parent[i] = s
    return parent, wpar


def path_min(n, parent, wpar, s, t):
    """在 Gomory-Hu 樹上求 s-t 路徑最小邊權。"""
    if s == t:
        return 0
    # 提升到同一深度
    def anc_depth(v):
        d = 0
        while parent[v] != -1:
            v = parent[v]
            d += 1
        return d

    ds, dt = anc_depth(s), anc_depth(t)
    a, b = s, t
    best = INF
    while ds > dt:
        best = min(best, wpar[a])
        a = parent[a]
        ds -= 1
    while dt > ds:
        best = min(best, wpar[b])
        b = parent[b]
        dt -= 1
    while a != b:
        best = min(best, wpar[a])
        best = min(best, wpar[b])
        a = parent[a]
        b = parent[b]
    return best


def solve_queries(n, edges, queries):
    parent, wpar = build_tree_weights(n, edges)
    out = []
    for s, t in queries:
        if s == t:
            out.append(0)
        else:
            out.append(path_min(n, parent, wpar, s, t))
    return out


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
    n = min(max(int(tn), 0), 60)
    tm, tq = nxt(), nxt()
    if tm is None or tq is None or not _INT.match(tm) or not _INT.match(tq):
        return ''
    m = min(max(int(tm), 0), 2000)
    q = min(max(int(tq), 0), 2000)
    edges = []
    for _ in range(m):
        a, b, c = nxt(), nxt(), nxt()
        if a is None or b is None or c is None:
            break
        if not (_INT.match(a) and _INT.match(b) and _INT.match(c)):
            break
        u, v, w = int(a), int(b), int(c)
        if 0 <= u < n and 0 <= v < n:
            edges.append((u, v, min(max(w, 0), 10**8)))
    queries = []
    for _ in range(q):
        a, b = nxt(), nxt()
        if a is None or b is None:
            break
        if not (_INT.match(a) and _INT.match(b)):
            break
        s, t = int(a), int(b)
        if not (0 <= s < n and 0 <= t < n):
            continue
        queries.append((s, t))
    if n == 0:
        return ''.join('0\n' for _ in queries)
    res = solve_queries(n, edges, queries)
    return ''.join(str(x) + '\n' for x in res)


def run_tests():
    # 固定情境：三角形 0-1(2),1-2(3),0-2(4)
    # 任意兩點最小割：0-1 = min(2, 3+4)=2? 實際 {0}|{1,2} 割 = 2+4=6；
    #   {0,1}|{2} 割 = 3+4=7；{0,2}|{1} 割 = 2+3=5 → 0-1 割 = min(6,7)=6
    e = [(0, 1, 2), (1, 2, 3), (0, 2, 4)]
    p, wp = build_tree_weights(3, e)
    # 對所有點對與 Edmonds-Karp 對拍
    adj = [[0] * 3 for _ in range(3)]
    for u, v, w in e:
        adj[u][v] += w
        adj[v][u] += w
    for s in range(3):
        for t in range(3):
            want, _ = maxflow_mincut_side(3, adj, s, t)
            assert path_min(3, p, wp, s, t) == want, (s, t, path_min(3, p, wp, s, t), want)

    state = 20261119

    def rnd(k):
        nonlocal state
        state = (state * 1103515245 + 12345) & 0x7FFFFFFF
        return state % k

    for tc in range(120):
        n = rnd(7) + 2                        # 2..8
        edges = []
        for u in range(n):
            for v in range(u + 1, n):
                if rnd(2) == 0:
                    edges.append((u, v, rnd(15) + 1))
        adj = [[0] * n for _ in range(n)]
        for u, v, w in edges:
            adj[u][v] += w
            adj[v][u] += w
        p, wp = build_tree_weights(n, edges)
        for s in range(n):
            for t in range(n):
                want, _ = maxflow_mincut_side(n, adj, s, t)
                got = path_min(n, p, wp, s, t)
                assert got == want, (tc, n, edges, s, t, got, want)


def main():
    data = sys.stdin.read()
    if data.strip():
        sys.stdout.buffer.write(run_io(data).encode('utf-8'))
    else:
        run_tests()
        sys.stdout.buffer.write(b'all tests passed\n')


if __name__ == '__main__':
    main()
