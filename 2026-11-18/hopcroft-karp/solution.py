#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Hopcroft-Karp 二分圖最大匹配（BFS 分層 + DFS 找極大不相交增廣路）。

IO 模式（stdin）：
  第一行 nL nR m，接著 m 行 u v（0-based，u∈[0,nL), v∈[0,nR)）。
  可讀多筆測資直到 EOF。每筆輸出 4 行：
    1) 最大匹配數
    2) 匹配對，格式 u-v，依 u 升冪排序，空白分隔；空則輸出 -
    3) 最小點覆蓋大小（Kőnig 定理）
    4) 最小點覆蓋，格式 Lu / Rv（L 表左部、R 表右部），先 L 後 R 各自升冪；空則 -
無 stdin 時跑內建對拍（vs Kuhn O(VE) + 2^n 暴力最小點覆蓋），輸出 all tests passed。

複雜度 O(E·sqrt(V))，space O(V+E)。
"""

import sys
from collections import deque

# ----------------------------- Hopcroft-Karp -----------------------------

def hopcroft_karp(nL, nR, adj):
    """adj[u] = 右部鄰居列表。回傳 pairU, pairV（比對對象，-1 表示未匹配）。"""
    INF = 10**9
    pairU = [-1] * nL
    pairV = [-1] * nR
    dist = [0] * nL

    def bfs():
        q = deque()
        found = False
        for u in range(nL):
            if pairU[u] == -1:
                dist[u] = 0
                q.append(u)
            else:
                dist[u] = -1
        while q:
            u = q.popleft()
            for v in adj[u]:
                pu = pairV[v]
                if pu == -1:
                    found = True
                elif dist[pu] == -1:
                    dist[pu] = dist[u] + 1
                    q.append(pu)
        return found

    sys.setrecursionlimit(1000000)

    def dfs(u):
        for v in adj[u]:
            pu = pairV[v]
            if pu == -1 or (dist[pu] == dist[u] + 1 and dfs(pu)):
                pairU[u] = v
                pairV[v] = u
                return True
        dist[u] = -1
        return False

    matching = 0
    while bfs():
        for u in range(nL):
            if pairU[u] == -1 and dfs(u):
                matching += 1
    return pairU, pairV


def min_vertex_cover(nL, nR, adj, pairU, pairV):
    """Kőnig：從未匹配左點出發交錯 DFS，覆蓋 = (左部扣掉 Z) 聯集 (右部交 Z)。"""
    seenL = [False] * nL
    seenR = [False] * nR
    q = deque()
    for u in range(nL):
        if pairU[u] == -1:
            seenL[u] = True
            q.append(u)
    while q:
        u = q.popleft()
        for v in adj[u]:
            if not seenR[v]:
                seenR[v] = True
                pu = pairV[v]
                if pu != -1 and not seenL[pu]:
                    seenL[pu] = True
                    q.append(pu)
    left = [u for u in range(nL) if not seenL[u]]
    right = [v for v in range(nR) if seenR[v]]
    return left, right


# ----------------------------- 暴力參考 -----------------------------

def kuhn(nL, nR, adj):
    pairV = [-1] * nR

    def try_k(u, seen):
        for v in adj[u]:
            if seen[v]:
                continue
            seen[v] = True
            if pairV[v] == -1 or try_k(pairV[v], seen):
                pairV[v] = u
                return True
        return False

    res = 0
    for u in range(nL):
        if try_k(u, [False] * nR):
            res += 1
    return res


def brute_min_cover_size(nL, nR, edges):
    """2^(nL+nR) 枚舉最小點覆蓋（僅供小圖對拍）。"""
    n = nL + nR
    best = n
    for mask in range(1 << n):
        k = bin(mask).count('1')
        if k >= best:
            continue
        ok = True
        for u, v in edges:
            if not ((mask >> u) & 1) and not ((mask >> (nL + v)) & 1):
                ok = False
                break
        if ok:
            best = k
    return best


# ------------------------------- 格式化 -------------------------------

def format_result(nL, nR, adj, pairU, pairV, left, right):
    out = []
    out.append(str(sum(1 for u in range(nL) if pairU[u] != -1)))
    pairs = sorted((u, pairU[u]) for u in range(nL) if pairU[u] != -1)
    out.append(' '.join('%d-%d' % (u, v) for u, v in pairs) if pairs else '-')
    out.append(str(len(left) + len(right)))
    nodes = ['L%d' % u for u in left] + ['R%d' % v for v in right]
    out.append(' '.join(nodes) if nodes else '-')
    return '\n'.join(out) + '\n'


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
        a, b, c = nxt(), nxt(), nxt()
        if a is None or b is None or c is None:
            break
        try:
            nL, nR, m = int(a), int(b), int(c)
        except ValueError:
            break
        if not (0 <= nL <= 2000 and 0 <= nR <= 2000 and 0 <= m <= 100000):
            break
        adj = [[] for _ in range(nL)]
        ok = True
        for _ in range(m):
            x, y = nxt(), nxt()
            if x is None or y is None:
                ok = False
                break
            try:
                u, v = int(x), int(y)
            except ValueError:
                ok = False
                break
            if 0 <= u < nL and 0 <= v < nR:
                adj[u].append(v)
        if not ok:
            break
        for u in range(nL):
            adj[u].sort()
            adj[u] = list(dict.fromkeys(adj[u]))
        pu, pv = hopcroft_karp(nL, nR, adj)
        left, right = min_vertex_cover(nL, nR, adj, pu, pv)
        out.append(format_result(nL, nR, adj, pu, pv, left, right).rstrip('\n'))
    return ''.join(line + '\n' for line in out)


# ------------------------------ 內建測試 ------------------------------

def run_tests():
    # 固定情境（輸出已實跑確認）
    adj = [[1, 0], [0], [2]]
    nL, nR = 3, 3
    pu, pv = hopcroft_karp(nL, nR, adj)
    assert sum(1 for u in range(nL) if pu[u] != -1) == 3
    left, right = min_vertex_cover(nL, nR, adj, pu, pv)
    sample = format_result(nL, nR, adj, pu, pv, left, right)
    assert sample == '3\n0-1 1-0 2-2\n3\nL0 L1 L2\n', repr(sample)

    # 隨機對拍：匹配數 vs Kuhn；最小點覆蓋 vs 2^n 暴力
    state = 20261118

    def rnd(m):
        nonlocal state
        state = (state * 1103515245 + 12345) & 0x7FFFFFFF
        return state % m

    for tc in range(500):
        nL, nR = rnd(5) + 1, rnd(5) + 1
        adj = [[] for _ in range(nL)]
        edges = []
        for u in range(nL):
            for v in range(nR):
                if rnd(3) == 0:
                    adj[u].append(v)
                    edges.append((u, v))
        pu, pv = hopcroft_karp(nL, nR, adj)
        match = sum(1 for u in range(nL) if pu[u] != -1)
        assert match == kuhn(nL, nR, adj), (tc, nL, nR, adj, match)
        # 匹配是合法的
        for u in range(nL):
            if pu[u] != -1:
                assert pv[pu[u]] == u and pu[u] in adj[u]
        left, right = min_vertex_cover(nL, nR, adj, pu, pv)
        cover = len(left) + len(right)
        assert cover == match, (tc, cover, match)   # Kőnig
        assert cover == brute_min_cover_size(nL, nR, edges), (tc, edges, cover)
        # 覆蓋確實蓋住每條邊
        for u, v in edges:
            assert (u in left) or (v in right)

    # 邊界：空圖、單邊、完全二分圖
    assert hopcroft_karp(0, 0, []) == ([], [])
    assert sum(1 for u in range(2) if hopcroft_karp(2, 2, [[0], [1]])[0][u] != -1) == 2
    full = [[0, 1, 2], [0, 1, 2], [0, 1, 2]]
    pu, _ = hopcroft_karp(3, 3, full)
    assert all(x != -1 for x in pu)


def main():
    data = sys.stdin.read()
    if data.strip():
        sys.stdout.buffer.write(run_io(data).encode('utf-8'))
    else:
        run_tests()
        sys.stdout.buffer.write(b'all tests passed\n')


if __name__ == '__main__':
    main()
