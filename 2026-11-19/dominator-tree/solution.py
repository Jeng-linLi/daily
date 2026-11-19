#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""支配樹（Dominator Tree）與 Lengauer–Tarjan 演算法。

在有向圖中（從入口 s 出發），若「所有從 s 到 v 的路徑都必經 u」，則 u 支配 v。
immediate dominator（idom[v]）是離 v 最近的嚴格支配點，所有 idom 邊構成一棵樹。

## 應用
- 編譯器控制流分析：SSA（靜態單一賦值）形式、支配前向樹、循環自然循環識別。
- 程式最佳化：不變數代碼外提（LICM）、可用表達式分析都建立在支配關係上。

## IO（stdin）
    n m s
    u1 v1
    ...          # 0-based，有向邊 u->v，s 為入口
    um vm
輸出一行：n 個整數 idom[0..n-1]，idom[i] = i 或不可達者輸出 -1，其餘為直接支配點。

stdin 無內容時跑內建測試：與 O(n^3) 迭代支配集合演算法對拍，輸出 all tests passed。
"""

import sys
import re

_INT = re.compile(r'^[+-]?[0-9]+$')


def dominator_tree(n, edges, s):
    """Lengauer–Tarjan，回傳 idom 陣列（idom[s] = -1）。"""
    if n == 0:
        return []
    pred = [[] for _ in range(n)]
    succ = [[] for _ in range(n)]
    for u, v in edges:
        succ[u].append(v)
        pred[v].append(u)

    # --- 1) DFS 編號（迭代版，避免深度遞迴） ---
    dfn = [0] * n            # dfs number，1-based；0 表示未訪問
    vertex = [0] * (n + 1)   # vertex[k] = dfn 為 k 的點
    parent = [-1] * n        # DFS 樹父節點
    timer = 0
    stack = [(s, 0)]
    while stack:
        v, idx = stack[-1]
        if idx == 0:
            timer += 1
            dfn[v] = timer
            vertex[timer] = v
        if idx < len(succ[v]):
            u = succ[v][idx]
            stack[-1] = (v, idx + 1)
            if dfn[u] == 0:
                parent[u] = v
                stack.append((u, 0))
        else:
            stack.pop()

    # --- 2) 半支配點 + eval/link（帶路徑壓縮的並查集） ---
    semi = dfn[:]           # semi[v] 以 dfn 值表示，初值為自身 dfn
    label = list(range(n))
    ancestor = [-1] * n      # 0 表示無；此處 -1 表示無
    idom = [-1] * n
    bucket = [[] for _ in range(n)]

    def compress(v):
        # 路徑壓縮：沿 ancestor 鏈把 label 更新為 semi 最小者（須遞迴到鏈底）
        if ancestor[v] != -1 and ancestor[ancestor[v]] != -1:
            if semi[label[ancestor[v]]] < semi[label[v]]:
                label[v] = label[ancestor[v]]
            ancestor[v] = ancestor[ancestor[v]]
            compress(v)

    def eval_node(v):
        if ancestor[v] == -1:
            return label[v]
        compress(v)
        return label[v]

    for i in range(timer, 1, -1):
        w = vertex[i]
        # 計算 semi[w]
        for v in pred[w]:
            if dfn[v] == 0:
                continue
            u = eval_node(v)
            if semi[u] < semi[w]:
                semi[w] = semi[u]
        bucket[vertex[semi[w]]].append(w)
        # link: w 的父接到 DFS 樹父
        p = parent[w]
        if p != -1:
            ancestor[w] = p
        # 處理 bucket[p]
        if p != -1:
            for v in bucket[p]:
                u = eval_node(v)
                idom[v] = u if semi[u] < semi[v] else p
            bucket[p] = []

    # --- 3) 第二遍：修正 idom ---
    for i in range(2, timer + 1):
        w = vertex[i]
        if idom[w] != vertex[semi[w]]:
            idom[w] = idom[idom[w]]

    idom[s] = -1
    # 不可達點輸出 -1
    for v in range(n):
        if dfn[v] == 0:
            idom[v] = -1
    return idom


def brute_idom(n, edges, s):
    """迭代支配集合演算法（O(n^3)），僅供對拍。"""
    if n == 0:
        return []
    pred = [[] for _ in range(n)]
    for u, v in edges:
        pred[v].append(u)
    # 可達點
    reach = set()
    stack = [s]
    while stack:
        v = stack.pop()
        if v in reach:
            continue
        reach.add(v)
        for u in succ_of(edges, v):
            stack.append(u)
    nodes = sorted(reach)
    full = frozenset(nodes)
    dom = {}
    for v in nodes:
        dom[v] = full if v != s else frozenset([s])
    changed = True
    while changed:
        changed = False
        for v in nodes:
            if v == s:
                continue
            ps = [p for p in pred[v] if p in reach]
            if not ps:
                new = frozenset([v])
            else:
                inter = set(dom[ps[0]])
                for p in ps[1:]:
                    inter &= dom[p]
                inter.add(v)
                new = frozenset(inter)
            if new != dom[v]:
                dom[v] = new
                changed = True
    idom = [-1] * n
    for v in nodes:
        if v == s:
            idom[v] = -1
            continue
        # 嚴格支配點中 |dom| 最大者即 idom
        strict = [u for u in dom[v] if u != v]
        if strict:
            best = max(strict, key=lambda u: (len(dom[u]), -u))
            idom[v] = best
    return idom


def succ_of(edges, v):
    out = []
    for u, w in edges:
        if u == v:
            out.append(w)
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
    n = min(max(int(tn), 0), 1000)
    tm = nxt()
    ts = nxt()
    if tm is None or ts is None or not _INT.match(tm) or not _INT.match(ts):
        return ''
    m = min(max(int(tm), 0), 20000)
    s = int(ts)
    if not (0 <= s < n):
        s = 0
    edges = []
    for _ in range(m):
        a, b = nxt(), nxt()
        if a is None or b is None:
            break
        if not (_INT.match(a) and _INT.match(b)):
            break
        u, v = int(a), int(b)
        if 0 <= u < n and 0 <= v < n:
            edges.append((u, v))
    res = dominator_tree(n, edges, s)
    return ' '.join(str(x) for x in res) + '\n'


def run_tests():
    # 固定情境：簡單鏈 0->1->2，idom = [-1,0,1]
    assert dominator_tree(3, [(0, 1), (1, 2)], 0) == [-1, 0, 1]
    # 菱形 0->1,0->2,1->3,2->3：3 的 idom 是 0
    r = dominator_tree(4, [(0, 1), (0, 2), (1, 3), (2, 3)], 0)
    assert r == [-1, 0, 0, 0], r

    state = 20261119

    def rnd(k):
        nonlocal state
        state = (state * 1103515245 + 12345) & 0x7FFFFFFF
        return state % k

    for tc in range(500):
        n = rnd(9) + 1                       # 1..9
        s = rnd(n)
        edges = []
        for u in range(n):
            for v in range(n):
                if u != v and rnd(4) == 0:
                    edges.append((u, v))
                    if rnd(6) == 0:          # 摻入重邊，考驗 eval/link 的路徑壓縮
                        edges.append((u, v))
        got = dominator_tree(n, edges, s)
        want = brute_idom(n, edges, s)
        assert got == want, (tc, n, s, edges, got, want)


def main():
    data = sys.stdin.read()
    if data.strip():
        sys.stdout.buffer.write(run_io(data).encode('utf-8'))
    else:
        run_tests()
        sys.stdout.buffer.write(b'all tests passed\n')


if __name__ == '__main__':
    main()
