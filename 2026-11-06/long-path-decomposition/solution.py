# -*- coding: utf-8 -*-
"""長鏈剖分（Long-Path Decomposition）：O(1) 求 k 級祖先

核心思想：
  每個節點往「子樹深度最深的兒子」（長兒子）延伸一條長鏈，
  把整棵樹剖成若干條互不相交的鏈。對每條鏈的鏈頂 top 預處理：
    down[top][i] : 沿鏈向下走 i 步的節點（0..L，L 為鏈長）
    up[top][j]   : top 的 j 級祖先（j = 0..min(L, dep[top])）
  查詢 u 的 k 級祖先（k <= len[u]，len[u] 為 u 往下最長距離 = L - d）可 O(1)：
    d = dep[u] - dep[top[u]]
    k <= d        -> down[top[u]][d - k]
    d < k <= len[u] -> up[top[u]][k - d]   （k - d <= dep[top]，不會越界）
  k > len[u] 時退化到倍增表（O(log n)）。
  另附 dep / len（往下最長距離）/ top（鏈頂）查詢。

IO 模式（stdin 有內容時）：
  n q
  p_1 ... p_{n-1}        （節點 1..n-1 的父節點，保證 p_i < i，根為 0）
  q 行查詢：
    1 u k  -> u 的 k 級祖先（不存在輸出 -1）
    2 u    -> dep[u]
    3 u    -> len[u]（u 到子樹最深葉子的邊數）
    4 u    -> u 所在鏈的鏈頂編號
無 stdin 時跑內建斷言並輸出 all tests passed。
"""
import sys

INF = 1 << 60


def build(n, par):
    """建長鏈剖分結構，回傳 (dep, maxlen, top, down, up, upk)。"""
    dep = [0] * n
    for i in range(1, n):
        dep[i] = dep[par[i]] + 1
    # maxlen[u]：u 往下沿最長路徑到葉子的邊數
    maxlen = [0] * n
    for i in range(n - 1, 0, -1):  # 子節點編號大於父節點，倒序保證先算兒子
        p = par[i]
        if maxlen[i] + 1 > maxlen[p]:
            maxlen[p] = maxlen[i] + 1
    # 長兒子：平手取編號最小者（掃描順序固定，兩語言一致）
    hson = [-1] * n
    for i in range(1, n):
        p = par[i]
        if hson[p] < 0 and maxlen[i] + 1 == maxlen[p]:
            hson[p] = i
    # 鏈頂
    top = [0] * n
    for i in range(1, n):
        top[i] = top[par[i]] if i == hson[par[i]] else i
    # 每條鏈的 down / up 陣列
    down = [None] * n
    up = [None] * n
    for t in range(n):
        if top[t] != t:
            continue
        L = maxlen[t]
        chain = [t]
        cur = t
        while hson[cur] >= 0:
            cur = hson[cur]
            chain.append(cur)
        down[t] = chain  # 長度 L + 1
        cnt = min(L, dep[t])
        anc = [t]
        cur = t
        for _ in range(cnt):
            cur = par[cur]
            anc.append(cur)
        up[t] = anc  # up[t][j] = t 的 j 級祖先
    # 倍增表（k > len[u] 時的退化路徑）
    LOG = max(1, n.bit_length())
    upk = [par[:]]
    for j in range(1, LOG):
        prev = upk[-1]
        upk.append([prev[prev[i]] for i in range(n)])
    return dep, maxlen, top, down, up, upk


def kth(u, k, dep, maxlen, top, down, up, upk):
    if k == 0:
        return u
    if k > dep[u]:
        return -1
    t = top[u]
    d = dep[u] - dep[t]
    if k <= d:
        return down[t][d - k]
    if k <= maxlen[u]:
        return up[t][k - d]
    # 倍增退化
    v = u
    j = 0
    while k:
        if k & 1:
            v = upk[j][v]
        k >>= 1
        j += 1
    return v


def kth_brute(u, k, par, dep):
    """暴力：沿父指標爬 k 步。"""
    if k > dep[u]:
        return -1
    v = u
    for _ in range(k):
        v = par[v]
    return v


def maxlen_brute(u, ch):
    """暴力 DFS 求往下最長距離。"""
    best = 0
    for c in ch[u]:
        cand = maxlen_brute(c, ch) + 1
        if cand > best:
            best = cand
    return best


def top_brute(u, par, hson):
    """暴力：沿父指標爬，直到當前點不是父親的長兒子。"""
    v = u
    while par[v] != v and hson[par[v]] == v:
        v = par[v]
    return v


def run_tests():
    seed = 20261106

    def rnd(m):
        nonlocal seed
        seed = (seed * 1103515245 + 12345) & 0x7FFFFFFF
        return seed % m

    cases = []
    # 固定樹：鏈 / 星形 / 毛毛蟲 / 雙鏈蜘蛛
    cases.append([0, 0, 1, 2, 3, 4])            # 鏈 0-1-2-3-4-5
    cases.append([0, 0, 0, 0, 0])               # 星形
    cases.append([0, 0, 1, 1, 2, 2, 3, 3, 4])   # 二叉樹
    cases.append([0, 0, 1, 1, 2, 3, 4, 5, 6])   # 毛毛蟲
    cases.append([0, 0, 0, 1, 3, 4, 5, 2, 7, 8, 9, 10])  # 兩條長鏈
    for _ in range(200):
        n = rnd(200) + 1
        par = [0] * n
        for i in range(1, n):
            par[i] = rnd(i)
        cases.append(par)

    for par in cases:
        n = len(par)
        dep, maxlen, top, down, up, upk = build(n, par)
        ch = [[] for _ in range(n)]
        for i in range(1, n):
            ch[par[i]].append(i)
        # hson 由 build 內部邏輯重算一份做一致性檢查
        hson = [-1] * n
        for i in range(1, n):
            p = par[i]
            if hson[p] < 0 and maxlen[i] + 1 == maxlen[p]:
                hson[p] = i
        for u in range(n):
            # dep 暴力：往根爬到 0 的步數
            v = u
            steps = 0
            while v != 0:
                v = par[v]
                steps += 1
            assert dep[u] == steps, (n, u, "dep")
            assert maxlen[u] == maxlen_brute(u, ch), (n, u, "maxlen")
            assert top[u] == top_brute(u, par, hson), (n, u, "top")
            for k in range(0, dep[u] + 3):
                want = kth_brute(u, k, par, dep)
                got = kth(u, k, dep, maxlen, top, down, up, upk)
                assert got == want, (n, u, k, want, got)

    # 鏈長合計 = n（鏈互不相交），抽樣驗證
    for par in cases[:20]:
        n = len(par)
        dep, maxlen, top, down, up, upk = build(n, par)
        total = sum(len(down[t]) for t in range(n) if top[t] == t)
        assert total == n, (n, total)

    print("all tests passed")


def main():
    data = sys.stdin.buffer.read().split()
    if not data:
        run_tests()
        return
    it = iter(data)
    n = int(next(it))
    q = int(next(it))
    par = [0] * n
    for i in range(1, n):
        par[i] = int(next(it))
    dep, maxlen, top, down, up, upk = build(n, par)
    out = []
    for _ in range(q):
        op = int(next(it))
        u = int(next(it))
        if op == 1:
            k = int(next(it))
            out.append(str(kth(u, k, dep, maxlen, top, down, up, upk)))
        elif op == 2:
            out.append(str(dep[u]))
        elif op == 3:
            out.append(str(maxlen[u]))
        else:
            out.append(str(top[u]))
    sys.stdout.write("\n".join(out) + "\n")


if __name__ == "__main__":
    main()
