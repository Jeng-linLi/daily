#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""樹雜湊與樹同構判定（AHU 演算法 / 重心為根的典範雜湊）。

命令（IO 模式，從 stdin 讀入）：
  iso n1 m1 u1 v1 ... n2 m2 u2 v2 ...
      -> 兩棵無標號樹（各 n 個點、n-1 條邊，1-based）
      -> 同構輸出 1，否則 0
  hash n m u1 v1 ...
      -> 輸出該樹的典範雜湊（兩個模數下的值，與 C++ 版逐位元一致）

stdin 無內容時跑內建斷言測試，輸出 all tests passed。

核心要點：
  * 子樹簽名：hash(u) 由「所有孩子簽名排序後」多項式摺疊而得，
    與孩子順序無關——這正是無標號樹同構的本質。
  * 以重心為根：普通以 1 號點為根時，鏈狀樹的兩端會算出不同簽名，
    把同一棵樹判成異構。重心至多兩個，對每個重心各算一次取最小即可。
    （重心的性質：刪去後最大子樹不超過 n/2。）
  * 雙雜湊（兩個質數模）壓低碰撞機率；平手一律取字典序最小，保證確定性。
  * 應用：無標號樹計數（A000055）、子樹重複偵測、化学分子比对。
"""

import sys

M1 = 2147483647   # 2^31 - 1，質數
M2 = 2147483629   # 質數
B1 = 1000003
B2 = 1000033


def centroids(n, adj):
    """回傳樹的全部重心（1 或 2 個）。adj 為 0-based 鄰接表。"""
    if n == 1:
        return [0]
    size = [0] * n
    maxsub = [0] * n
    order = []
    parent = [-1] * n
    stack = [0]
    visited = [False] * n
    visited[0] = True
    while stack:
        u = stack.pop()
        order.append(u)
        for v in adj[u]:
            if not visited[v]:
                visited[v] = True
                parent[v] = u
                stack.append(v)
    for u in reversed(order):
        size[u] = 1
        for v in adj[u]:
            if v != parent[u]:
                size[u] += size[v]
        if parent[u] >= 0:
            maxsub[parent[u]] = max(maxsub[parent[u]], size[u])
    best = None
    cents = []
    for u in range(n):
        up = n - size[u]
        m = max(maxsub[u], up)
        if best is None or m < best:
            best = m
            cents = [u]
        elif m == best:
            cents.append(u)
    return cents


def tree_hash(n, adj):
    """典範雜湊：對每個重心求子樹簽名，取最小。回傳 (h1, h2)。"""
    best = None
    for c in centroids(n, adj):
        h = hash_sub(c, -1, adj)
        if best is None or h < best:
            best = h
    return best


def hash_sub(u, p, adj):
    """以 u 為根的子樹簽名（雙雜湊 (h1, h2)，孩子簽名排序後摺疊）。"""
    forms = []
    for v in adj[u]:
        if v != p:
            forms.append(hash_sub(v, u, adj))
    forms.sort()
    h1 = 1
    h2 = 1
    for (a, b) in forms:
        h1 = (h1 * B1 + a * 2147483648 + b) % M1
        h2 = (h2 * B2 + a * 2147483648 + b) % M2
    return (h1, h2)


def build_tree(n, edges):
    """邊清單（1-based）-> 0-based 鄰接表。"""
    adj = [[] for _ in range(n)]
    for (u, v) in edges:
        adj[u - 1].append(v - 1)
        adj[v - 1].append(u - 1)
    return adj


def isomorphic(n1, e1, n2, e2):
    if n1 != n2:
        return False
    if n1 == 0:
        return True
    return tree_hash(n1, build_tree(n1, e1)) == tree_hash(n2, build_tree(n2, e2))


def prufer_to_edges(seq, n):
    """Prüfer 序列 -> 邊清單（1-based）。seq 長度 n-2，元素 1..n。"""
    deg = [1] * (n + 1)
    for x in seq:
        deg[x] += 1
    edges = []
    ptr = 1
    leaf_heap = []  # 直接線性掃描即可（測試用，不求效能）
    leaves = [i for i in range(1, n + 1) if deg[i] == 1]
    import heapq
    heapq.heapify(leaves)
    for x in seq:
        leaf = heapq.heappop(leaves)
        edges.append((leaf, x))
        deg[leaf] -= 1
        deg[x] -= 1
        if deg[x] == 1:
            heapq.heappush(leaves, x)
    last = [i for i in range(1, n + 1) if deg[i] == 1]
    edges.append((last[0], last[1]))
    return edges


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

    def read_tree():
        n = nxt_int()
        if n is None:
            return None, None
        m = nxt_int()
        if m is None:
            return None, None
        edges = []
        ok = True
        for _ in range(m):
            u = nxt_int()
            v = nxt_int()
            if u is None or v is None:
                ok = False
                break
            edges.append((u, v))
        if not ok:
            return None, None
        return n, edges

    while True:
        cmd = nxt()
        if cmd is None:
            break
        if cmd == 'iso':
            n1, e1 = read_tree()
            if n1 is None:
                break
            n2, e2 = read_tree()
            if n2 is None:
                break
            out.append('1' if isomorphic(n1, e1, n2, e2) else '0')
        elif cmd == 'hash':
            n, e = read_tree()
            if n is None:
                break
            h1, h2 = tree_hash(n, build_tree(n, e))
            out.append('%d %d' % (h1, h2))
        else:
            break
    return ''.join(line + '\n' for line in out)


def run_tests():
    import itertools

    # 1) 邊界與手算例
    assert isomorphic(1, [], 1, [])
    assert isomorphic(2, [(1, 2)], 2, [(1, 2)])
    assert not isomorphic(2, [(1, 2)], 3, [(1, 2), (2, 3)])
    # 鏈 3 點 vs 星形 3 點（同構）
    assert isomorphic(3, [(1, 2), (2, 3)], 3, [(2, 1), (2, 3)])
    # 鏈 4 點 vs 星形 4 點（不同構）
    assert not isomorphic(4, [(1, 2), (2, 3), (3, 4)],
                          4, [(1, 2), (1, 3), (1, 4)])
    # 重心為根的必要性：把同一條鏈反轉編號，必須判同構
    assert isomorphic(5, [(1, 2), (2, 3), (3, 4), (4, 5)],
                      5, [(2, 1), (3, 2), (4, 3), (5, 4)])

    # 2) 重標號（隨機置換）後雜湊必須不變
    state = 20261111

    def rnd(mm):
        nonlocal state
        state = (state * 1103515245 + 12345) & 0x7FFFFFFF
        return state % mm

    def random_tree(n):
        edges = []
        for v in range(2, n + 1):
            u = rnd(v - 1) + 1
            edges.append((u, v))
        return edges

    def relabel(n, edges):
        perm = list(range(1, n + 1))
        for i in range(n - 1, 0, -1):
            j = rnd(i + 1)
            perm[i], perm[j] = perm[j], perm[i]
        return [(perm[u - 1], perm[v - 1]) for (u, v) in edges]

    for _ in range(200):
        n = rnd(40) + 1
        e = random_tree(n)
        h1 = tree_hash(n, build_tree(n, e))
        e2 = relabel(n, e)
        h2 = tree_hash(n, build_tree(n, e2))
        assert h1 == h2, (n, e, e2)

    # 3) 小 n 與暴力同構（全置換）交叉驗證：hash 相等 <=> 同構
    def brute_iso(n1, e1, n2, e2):
        if n1 != n2:
            return False
        adj2 = set()
        for (u, v) in e2:
            adj2.add((u, v))
            adj2.add((v, u))
        for perm in itertools.permutations(range(1, n1 + 1)):
            ok = True
            for (u, v) in e1:
                if (perm[u - 1], perm[v - 1]) not in adj2:
                    ok = False
                    break
            if ok:
                return True
        return False

    for _ in range(150):
        n = rnd(6) + 1
        e1 = random_tree(n)
        e2 = random_tree(n) if rnd(2) else relabel(n, e1)
        got = isomorphic(n, e1, n, e2)
        want = brute_iso(n, e1, n, e2)
        assert got == want, (n, e1, e2, got, want)

    # 4) 無標號樹計數：n=1..7 全部 Prüfer 序列去重後的雜湊數
    #    應等於 OEIS A000055：1, 1, 1, 2, 3, 6, 11
    unlabeled = [1, 1, 1, 2, 3, 6, 11]
    for n in range(1, 8):
        hashes = set()
        if n == 1:
            hashes.add(tree_hash(1, build_tree(1, [])))
        else:
            for seq in itertools.product(range(1, n + 1), repeat=n - 2):
                e = prufer_to_edges(list(seq), n)
                hashes.add(tree_hash(n, build_tree(n, e)))
        assert len(hashes) == unlabeled[n - 1], (n, len(hashes))


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
