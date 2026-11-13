#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""點分治（Centroid Decomposition）：統計樹上「邊數恰為 k」的點對路徑數。

命令（IO 模式，從 stdin 讀入）：
  paths n k u1 v1 u2 v2 ...   （n 個點、n-1 條無向邊，1-based）
      -> 輸出點對 (u < v) 的個數，使兩點間路徑的邊數恰為 k

stdin 無內容時跑內建斷言測試，輸出 all tests passed。

核心要點：
  * 每輪在當前連通分量裡找**重心**（刪去後最大子樹 <= 分量一半），
    統計「經過重心的路徑」中長度恰為 k 的點對，然後刪掉重心、
    對每個剩餘子分量遞迴。
  * 正確性：任意點對 (u, v) 的路徑，在「u、v 首次被分開」的那一層
    恰好經過該層的重心，且只被統計一次——分層天然不重不漏。
  * 經過重心的路徑計數（有序對，dist(u,c) + dist(v,c) = k）：
      through = Σ_d cnt_all[d] * cnt_all[k-d] − Σ_子樹 Σ_d cnt_s[d] * cnt_s[k-d]
    減掉「兩端落在同一棵子樹」的部分；k >= 1 時結果必為偶數（每個
    無序對恰被計兩次），除以 2 即答案。
  * 重心查找：BFS 求 subtree size，從起點不斷走向「超過一半」的重孩子，
    收斂點即重心。分量每輪至少減半，遞迴深度 O(log n)。
  * 複雜度 O(n log n)。
"""

import sys

INF = float('inf')


def count_paths(n, k, edges):
    """回傳樹上距離（邊數）恰為 k 的無序點對數。"""
    if n <= 1 or k < 1:
        return 0
    adj = [[] for _ in range(n)]
    for (u, v) in edges:
        adj[u].append(v)
        adj[v].append(u)
    removed = [False] * n
    total = 0

    def bfs(start, avoid):
        """在未刪除的節點中從 start 出發 BFS（可避開 avoid）。
        回傳 (order, parent, dist)。"""
        order = [start]
        parent = [-1] * n
        dist = [-1] * n
        dist[start] = 0
        head = 0
        while head < len(order):
            u = order[head]
            head += 1
            for v in adj[u]:
                if not removed[v] and v != avoid and dist[v] < 0:
                    dist[v] = dist[u] + 1
                    parent[v] = u
                    order.append(v)
        return order, parent, dist

    def find_centroid(start):
        order, parent, _ = bfs(start, -1)
        size = [0] * n
        for u in reversed(order):
            size[u] = 1
            for v in adj[u]:
                if not removed[v] and v != parent[u]:
                    size[u] += size[v]
        tot = size[start]
        c = start
        while True:
            heavy = -1
            for v in adj[c]:
                if not removed[v] and v != parent[c] and size[v] * 2 > tot:
                    heavy = v
                    break
            if heavy < 0:
                return c
            c = heavy

    def decompose(start):
        nonlocal total
        c = find_centroid(start)
        order, _, distc = bfs(c, -1)
        # 經過重心的有序對：Σ cnt_all[d]*cnt_all[k-d] − Σ 子樹 Σ cnt_s[d]*cnt_s[k-d]
        cnt_all = [0] * (k + 1)
        for u in order:
            if distc[u] <= k:
                cnt_all[distc[u]] += 1
        through = 0
        for d in range(k + 1):
            through += cnt_all[d] * cnt_all[k - d]
        for v0 in adj[c]:
            if removed[v0]:
                continue
            o2, _, d2 = bfs(v0, c)  # 子樹內 BFS（避開 c）
            cnt_s = [0] * (k + 1)
            for u in o2:
                du = d2[u] + 1  # 距重心 = 距 v0 + 1
                if du <= k:
                    cnt_s[du] += 1
            for d in range(k + 1):
                through -= cnt_s[d] * cnt_s[k - d]
        total += through // 2
        removed[c] = True
        for v in adj[c]:
            if not removed[v]:
                decompose(v)

    decompose(0)
    return total


def count_paths_brute(n, k, edges):
    """每個點 BFS 全樹的 O(n^2) 暴力，僅供內建測試對拍。"""
    if n <= 1 or k < 1:
        return 0
    adj = [[] for _ in range(n)]
    for (u, v) in edges:
        adj[u].append(v)
        adj[v].append(u)
    pairs = 0
    for s in range(n):
        dist = [-1] * n
        dist[s] = 0
        q = [s]
        head = 0
        while head < len(q):
            u = q[head]
            head += 1
            for v in adj[u]:
                if dist[v] < 0:
                    dist[v] = dist[u] + 1
                    q.append(v)
        for t in range(s + 1, n):
            if dist[t] == k:
                pairs += 1
    return pairs


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
        # 嚴格整數：^[+-]?[0-9]+$
        body = t[1:] if t[0] in '+-' else t
        if not body or not body.isdigit():
            return None
        return int(t)

    while True:
        cmd = nxt()
        if cmd is None:
            break
        if cmd == 'paths':
            n = nxt_int()
            k = nxt_int()
            if n is None or k is None:
                break
            edges = []
            ok = True
            for _ in range(max(0, n - 1)):
                u = nxt_int()
                v = nxt_int()
                if u is None or v is None:
                    ok = False
                    break
                edges.append((u - 1, v - 1))
            if not ok:
                break
            out.append(str(count_paths(n, k, edges)))
        else:
            break
    return ''.join(line + '\n' for line in out)


def run_tests():
    # 1) 固定小例：鏈 0-1-2-3-4、星形（中心 0 連 1..4）
    chain = [(0, 1), (1, 2), (2, 3), (3, 4)]
    star = [(0, 1), (0, 2), (0, 3), (0, 4)]
    assert count_paths(5, 1, chain) == 4
    assert count_paths(5, 2, chain) == 3
    assert count_paths(5, 3, chain) == 2
    assert count_paths(5, 4, chain) == 1
    assert count_paths(5, 5, chain) == 0
    assert count_paths(5, 1, star) == 4
    assert count_paths(5, 2, star) == 6   # C(4,2)，任兩片葉
    assert count_paths(5, 3, star) == 0
    assert count_paths(1, 1, []) == 0
    assert count_paths(2, 1, [(0, 1)]) == 1

    # 2) 與 O(n^2) 暴力對拍：隨機樹 + 隨機 k
    state = 20261112

    def rnd(mm):
        nonlocal state
        state = (state * 1103515245 + 12345) & 0x7FFFFFFF
        return state % mm

    for _ in range(200):
        n = rnd(80) + 1
        edges = [(rnd(i), i) for i in range(1, n)]  # 每點隨機連向前面的點
        k = rnd(12) + 1
        got = count_paths(n, k, edges)
        want = count_paths_brute(n, k, edges)
        assert got == want, (n, k, edges, got, want)

    # 3) 鏈上大 n（考驗迭代 BFS / 重心查找）
    n = 400
    edges = [(i, i + 1) for i in range(n - 1)]
    assert count_paths(n, 7, edges) == n - 7
    assert count_paths(n, 399, edges) == 1


def main():
    data = sys.stdin.read()
    if data.strip():
        # 二進位輸出：避免 Windows 文字模式把 \n 轉成 \r\n（與 C++ 版不一致）
        sys.stdout.buffer.write(run_io(data).encode('utf-8'))
    else:
        run_tests()
        sys.stdout.buffer.write(b'all tests passed\n')


if __name__ == '__main__':
    main()
