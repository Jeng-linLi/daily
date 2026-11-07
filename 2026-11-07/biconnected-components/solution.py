# -*- coding: utf-8 -*-
"""雙連通分量：邊雙連通分量（e-DCC）、點雙連通分量（v-DCC）、橋、割點

題意：給定一張簡單無向圖（n 點 m 邊，無自環、無重邊），求：
        1. 橋（刪掉後圖不連通的邊）；
        2. 割點（刪掉後圖不連通的點）；
        3. 邊雙連通分量（實作上 = 刪去所有橋後的連通分量，
           孤立點各自成一個分量）；
        4. 點雙連通分量（割點可屬於多個 v-DCC，橋的兩端構成大小 2 的
           v-DCC，孤立點自成一個 v-DCC）。

思路（Tarjan，迭代版避免深遞迴）：
  橋：一次 DFS 求 low/disc，樹邊 (p,u) 是橋 iff low[u] > disc[p]。
  邊雙：刪橋後 BFS/DFS 劃分分量（含孤立點），分量按「最小點編號」
        重新編號——與 DFS 順序無關，保證跨語言輸出一致。
  點雙：維護點棧，child c 完成時若 low[c] >= disc[u]，把棧中 c 以上的
        點連同 u 彈出成一個分量；根若無子樹分量則自成單點分量。
  割點：u 是割點 iff（u 是 DFS 根且子樹數 >= 2）或存在子樹 c 使
        low[c] >= disc[u]。

IO 模式（stdin 有內容時）：
  n m
  m 行：u v
  輸出（依序）：
    行 1：邊雙分量數 e
    行 2：n 個整數，每點的邊雙編號（按分量最小點升序編號）
    行 3：點雙分量數 v
    行 4..：v 行，每行「大小 + 升序點列」（分量按字典序升序）
    行：橋數 b；接著 b 行「u v」（端點升序、橋之間字典序升序）
    行：割點數 a；接著 1 行割點升序點列（a=0 時為空行）
  無 stdin 時跑內建斷言（與指數暴力對拍）並輸出 all tests passed。
"""
import sys


def _build_adj(n, edges):
    adj = [[] for _ in range(n)]
    for i, (u, v) in enumerate(edges):
        adj[u].append((v, i))
        adj[v].append((u, i))
    return adj


def tarjan_bridges(n, edges):
    """回傳橋的邊編號集合（迭代 Tarjan）。"""
    adj = _build_adj(n, edges)
    disc = [-1] * n
    low = [0] * n
    ptr = [0] * n
    timer = 0
    bridges = set()
    for s in range(n):
        if disc[s] != -1:
            continue
        disc[s] = low[s] = timer
        timer += 1
        stack = [(s, -1)]
        while stack:
            u, pe = stack[-1]
            advanced = False
            while ptr[u] < len(adj[u]):
                v, eid = adj[u][ptr[u]]
                ptr[u] += 1
                if eid == pe:
                    continue
                if disc[v] == -1:
                    disc[v] = low[v] = timer
                    timer += 1
                    stack.append((v, eid))
                    advanced = True
                    break
                if disc[v] < low[u]:
                    low[u] = disc[v]
            if not advanced:
                stack.pop()
                if stack:
                    p = stack[-1][0]
                    if low[u] < low[p]:
                        low[p] = low[u]
                    if low[u] > disc[p]:
                        bridges.add(pe)
    return bridges


def edge_dcc(n, edges):
    """邊雙：刪橋後連通分量 + 孤立點；按最小點編號重編號。"""
    bridges = tarjan_bridges(n, edges)
    adj = _build_adj(n, edges)
    comp = [-1] * n
    nc = 0
    for s in range(n):
        if comp[s] != -1:
            continue
        comp[s] = nc
        st = [s]
        while st:
            u = st.pop()
            for v, eid in adj[u]:
                if eid not in bridges and comp[v] == -1:
                    comp[v] = nc
                    st.append(v)
        nc += 1
    mins = {}
    for v in range(n):
        c = comp[v]
        if c not in mins or v < mins[c]:
            mins[c] = v
    order = sorted(mins, key=lambda c: mins[c])
    remap = {c: i for i, c in enumerate(order)}
    return [remap[comp[v]] for v in range(n)]


def vertex_dcc(n, edges):
    """點雙（點棧 Tarjan）：回傳分量列表（分量內升序、分量間字典序升序）。"""
    adj = _build_adj(n, edges)
    disc = [-1] * n
    low = [0] * n
    ptr = [0] * n
    timer = 0
    vstk = []
    comps = []
    for s in range(n):
        if disc[s] != -1:
            continue
        created = len(comps)
        disc[s] = low[s] = timer
        timer += 1
        vstk.append(s)
        stack = [(s, -1)]
        while stack:
            u, pe = stack[-1]
            advanced = False
            while ptr[u] < len(adj[u]):
                v, eid = adj[u][ptr[u]]
                ptr[u] += 1
                if eid == pe:
                    continue
                if disc[v] == -1:
                    disc[v] = low[v] = timer
                    timer += 1
                    vstk.append(v)
                    stack.append((v, eid))
                    advanced = True
                    break
                if disc[v] < low[u]:
                    low[u] = disc[v]
            if not advanced:
                stack.pop()
                if stack:
                    p = stack[-1][0]
                    if low[u] < low[p]:
                        low[p] = low[u]
                    if low[u] >= disc[p]:
                        # p 是割點（或根），彈出 u 及其以上 + p 成一個分量
                        comp = []
                        while True:
                            w = vstk.pop()
                            comp.append(w)
                            if w == u:
                                break
                        comp.append(p)
                        comps.append(sorted(comp))
        if len(comps) == created:
            # 根孤立（整棵 DFS 沒產生任何分量）：自成單點分量
            comps.append([vstk.pop()])
        else:
            # 根已含在子分量中，棧底剩餘的根是重複，丟棄
            vstk.pop()
    return sorted(comps)


def articulation_points(n, edges):
    """割點：根看子樹數；其餘看是否存在子樹 c 使 low[c] >= disc[u]。"""
    adj = _build_adj(n, edges)
    disc = [-1] * n
    low = [0] * n
    ptr = [0] * n
    timer = 0
    is_cut = [False] * n
    for s in range(n):
        if disc[s] != -1:
            continue
        root_children = 0
        disc[s] = low[s] = timer
        timer += 1
        stack = [(s, -1)]
        while stack:
            u, pe = stack[-1]
            advanced = False
            while ptr[u] < len(adj[u]):
                v, eid = adj[u][ptr[u]]
                ptr[u] += 1
                if eid == pe:
                    continue
                if disc[v] == -1:
                    disc[v] = low[v] = timer
                    timer += 1
                    stack.append((v, eid))
                    advanced = True
                    if u == s:
                        root_children += 1
                    break
                if disc[v] < low[u]:
                    low[u] = disc[v]
            if not advanced:
                stack.pop()
                if stack:
                    p = stack[-1][0]
                    if low[u] < low[p]:
                        low[p] = low[u]
                    if low[u] >= disc[p] and p != s:
                        is_cut[p] = True
        if root_children >= 2:
            is_cut[s] = True
    return [v for v in range(n) if is_cut[v]]


# ---------- 暴力對拍（僅用於小規模） ----------

def _connected(vertices, adj):
    """vertices 內部連通性。"""
    if not vertices:
        return True
    vs = set(vertices)
    seen = {vertices[0]}
    st = [vertices[0]]
    while st:
        u = st.pop()
        for v, _ in adj[u]:
            if v in vs and v not in seen:
                seen.add(v)
                st.append(v)
    return len(seen) == len(vs)


def brute_bridges(n, edges):
    """刪邊看端點是否仍連通。"""
    res = []
    for i in range(len(edges)):
        rest = [e for j, e in enumerate(edges) if j != i]
        radj = _build_adj(n, rest)
        u, v = edges[i]
        seen = {u}
        st = [u]
        while st:
            x = st.pop()
            for y, _ in radj[x]:
                if y not in seen:
                    seen.add(y)
                    st.append(y)
        if v not in seen:
            res.append((min(u, v), max(u, v)))
    return sorted(res)


def brute_articulation(n, edges):
    """標準定義：刪去 x 後剩餘點的連通分量數比原圖增加，則 x 是割點。
    （孤立點刪去後分量數反而減一，不算割點。）"""
    adj = _build_adj(n, edges)

    def comp_count(vs):
        vs = set(vs)
        seen = set()
        cnt = 0
        for s in vs:
            if s in seen:
                continue
            cnt += 1
            seen.add(s)
            st = [s]
            while st:
                u = st.pop()
                for v, _ in adj[u]:
                    if v in vs and v not in seen:
                        seen.add(v)
                        st.append(v)
        return cnt

    full = comp_count(range(n))
    res = []
    for x in range(n):
        rest = [v for v in range(n) if v != x]
        # x 孤立時刪去後分量數減一，閾值也要減一
        if comp_count(rest) > full - (0 if adj[x] else 1):
            res.append(x)
    return sorted(res)


def brute_edge_dcc(n, edges):
    """刪所有橋後的連通分量 + 孤立點，按最小點編號重編號。"""
    bridges = set(brute_bridges(n, edges))
    adj = _build_adj(n, edges)
    comp = [-1] * n
    nc = 0
    for s in range(n):
        if comp[s] != -1:
            continue
        comp[s] = nc
        st = [s]
        while st:
            u = st.pop()
            for v, eid in adj[u]:
                key = (min(u, v), max(u, v))
                if key not in bridges and comp[v] == -1:
                    comp[v] = nc
                    st.append(v)
        nc += 1
    mins = {}
    for v in range(n):
        c = comp[v]
        if c not in mins or v < mins[c]:
            mins[c] = v
    order = sorted(mins, key=lambda c: mins[c])
    remap = {c: i for i, c in enumerate(order)}
    return [remap[comp[v]] for v in range(n)]


def brute_vertex_dcc(n, edges):
    """定義暴力：v-DCC = 極大的「2-點連通 / 橋對 / 孤立點」子集。
    n <= 8 時枚舉全部子集取極大者。"""
    eset = set(tuple(sorted(e)) for e in edges)
    deg = [0] * n
    for u, v in edges:
        deg[u] += 1
        deg[v] += 1

    def valid(sub):
        k = len(sub)
        if k == 1:
            return deg[sub[0]] == 0
        if k == 2:
            return tuple(sorted(sub)) in eset
        for skip in sub:
            rest = [v for v in sub if v != skip]
            if not _connected(rest, _build_adj(n, edges)):
                return False
        return True

    vsets = []
    for mask in range(1, 1 << n):
        sub = [i for i in range(n) if mask >> i & 1]
        if valid(sub):
            vsets.append(sub)
    maximal = []
    for a in vsets:
        sa = set(a)
        if not any(set(b) > sa for b in vsets):
            maximal.append(sorted(a))
    return sorted(maximal)


def run_tests():
    seed = 20261107

    def rnd(mv):
        nonlocal seed
        seed = (seed * 1103515245 + 12345) & 0x7FFFFFFF
        return seed % mv

    fixed = [
        (6, [(0, 1), (1, 2), (2, 0), (3, 4), (4, 5), (5, 3), (2, 3)]),  # 兩環一橋
        (4, [(0, 1), (1, 2), (2, 3)]),            # 鏈：全橋
        (4, [(0, 1), (1, 2), (2, 3), (3, 0)]),    # 四環
        (5, [(0, 1), (0, 2), (0, 3), (0, 4)]),    # 星：中心割點
        (5, [(0, 1), (1, 2), (2, 0), (2, 3), (3, 4), (4, 2)]),  # 8 字
        (3, []),                                   # 全孤立
        (5, [(0, 1), (1, 2), (2, 3), (3, 4), (4, 0), (0, 2), (0, 3)]),
    ]
    for n, es in fixed:
        assert brute_vertex_dcc(n, es) == vertex_dcc(n, es), (n, es)
        assert brute_edge_dcc(n, es) == edge_dcc(n, es), (n, es)
        assert brute_articulation(n, es) == articulation_points(n, es), (n, es)
        bs = set(brute_bridges(n, es))
        bset = tarjan_bridges(n, es)
        want = set((min(es[i][0], es[i][1]), max(es[i][0], es[i][1]))
                   for i in bset)
        assert bs == want, (n, es, bs, want)

    # 隨機小圖（n <= 7）與全暴力對拍
    for _ in range(250):
        n = rnd(7) + 1
        pairs = [(a, b) for a in range(n) for b in range(a + 1, n)]
        m = rnd(len(pairs) + 1)
        idx = list(range(len(pairs)))
        for i in range(len(idx) - 1, 0, -1):
            j = rnd(i + 1)
            idx[i], idx[j] = idx[j], idx[i]
        es = [pairs[i] for i in idx[:m]]
        assert brute_vertex_dcc(n, es) == vertex_dcc(n, es), (n, es)
        assert brute_edge_dcc(n, es) == edge_dcc(n, es), (n, es)
        assert brute_articulation(n, es) == articulation_points(n, es), (n, es)
        bs = set(brute_bridges(n, es))
        bset = tarjan_bridges(n, es)
        want = set((min(es[i][0], es[i][1]), max(es[i][0], es[i][1]))
                   for i in bset)
        assert bs == want, (n, es, bs, want)

    # 中等隨機圖冒煙：邊雙分量內部連通、點雙覆蓋所有點
    for _ in range(30):
        n = rnd(60) + 10
        m = rnd(2 * n) + n // 2
        es = []
        seen = set()
        guard = 0
        while len(es) < m and guard < m * 20:
            guard += 1
            a, b = rnd(n), rnd(n)
            if a == b:
                continue
            key = (min(a, b), max(a, b))
            if key in seen:
                continue
            seen.add(key)
            es.append(key)
        ecomp = edge_dcc(n, es)
        vcomps = vertex_dcc(n, es)
        mins = {}
        for v in range(n):
            mins.setdefault(ecomp[v], []).append(v)
        for c, vs in mins.items():
            if len(vs) >= 2:
                assert _connected(vs, _build_adj(n, es)), (n, es, c, vs)
        cov = set()
        for c in vcomps:
            cov |= set(c)
        assert cov == set(range(n)), (n, es)

    print("all tests passed")


def main():
    data = sys.stdin.buffer.read().split()
    if not data:
        run_tests()
        return
    it = iter(data)
    n = int(next(it))
    m = int(next(it))
    edges = []
    for _ in range(m):
        u = int(next(it))
        v = int(next(it))
        edges.append((u, v))
    ecomp = edge_dcc(n, edges)
    vcomps = vertex_dcc(n, edges)
    bset = tarjan_bridges(n, edges)
    cuts = articulation_points(n, edges)
    lines = []
    lines.append(str(len(set(ecomp))))
    lines.append(" ".join(map(str, ecomp)))
    lines.append(str(len(vcomps)))
    for c in vcomps:
        lines.append(str(len(c)) + " " + " ".join(map(str, c)))
    bl = sorted((min(edges[i][0], edges[i][1]), max(edges[i][0], edges[i][1]))
                for i in bset)
    lines.append(str(len(bl)))
    for p in bl:
        lines.append("%d %d" % p)
    lines.append(str(len(cuts)))
    lines.append(" ".join(map(str, cuts)))
    sys.stdout.write("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
