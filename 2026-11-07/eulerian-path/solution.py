# -*- coding: utf-8 -*-
"""歐拉路徑 / 歐拉迴路（Hierholzer 演算法）

題意：給定一張圖（可帶自環與重邊），求一條經過**每一條邊恰好一次**的
      走訪（歐拉路徑）；若起點與終點相同則為歐拉迴路。不存在輸出 -1。

判定條件：
  無向圖：度數為奇數的點恰為 0 個（迴路）或 2 個（路徑），
          且所有度數 > 0 的點在同一連通分量內。
  有向圖：每點 out-in ∈ {-1, 0, 1}，且 +1 與 -1 的點各至多一個
          （全 0 為迴路；有 +1 必有 -1 為路徑），
          且所有 out+in > 0 的點弱連通。

Hierholzer（迭代版）：沿未走過的邊不斷深入，走不動就把當前點彈進答案棧；
最後把棧倒序即為歐拉路徑。複雜度 O(n + m)。

IO 模式（stdin 有內容時）：
  T           （1 = 無向圖，2 = 有向圖）
  n m
  m 行：u v    （0-based，允許 u == v 自環與重邊）
  存在則輸出一行 m+1 個點（空格分隔，m == 0 時輸出單點 0），
  否則輸出 -1。
無 stdin 時跑內建斷言（與指數回溯暴力對拍存在性）並輸出 all tests passed。
"""
import sys


def solve_undirected(n, edges):
    """回傳歐拉路徑（點序列）或 None。"""
    m = len(edges)
    if m == 0:
        return [0] if n >= 1 else None
    deg = [0] * n
    adj = [[] for _ in range(n)]
    for i, (u, v) in enumerate(edges):
        adj[u].append((v, 2 * i))
        adj[v].append((u, 2 * i + 1))
        deg[u] += 1
        deg[v] += 1
    odd = [v for v in range(n) if deg[v] & 1]
    if len(odd) not in (0, 2):
        return None
    act = [v for v in range(n) if deg[v] > 0]
    seen = {act[0]}
    st = [act[0]]
    while st:
        u = st.pop()
        for v, _ in adj[u]:
            if v not in seen:
                seen.add(v)
                st.append(v)
    if len(seen) != len(act):
        return None
    start = odd[0] if odd else act[0]
    used = [False] * m
    it = [0] * n
    stack = [start]
    path = []
    while stack:
        u = stack[-1]
        advanced = False
        while it[u] < len(adj[u]):
            v, arc = adj[u][it[u]]
            it[u] += 1
            if not used[arc >> 1]:
                used[arc >> 1] = True
                stack.append(v)
                advanced = True
                break
        if not advanced:
            path.append(stack.pop())
    path.reverse()
    if len(path) != m + 1:
        return None
    return path


def solve_directed(n, edges):
    """回傳歐拉路徑（點序列）或 None。"""
    m = len(edges)
    if m == 0:
        return [0] if n >= 1 else None
    outdeg = [0] * n
    indeg = [0] * n
    adj = [[] for _ in range(n)]
    radj = [[] for _ in range(n)]
    for i, (u, v) in enumerate(edges):
        adj[u].append((v, i))
        radj[v].append(u)
        outdeg[u] += 1
        indeg[v] += 1
    bal = [outdeg[i] - indeg[i] for i in range(n)]
    plus = sum(1 for b in bal if b == 1)
    minus = sum(1 for b in bal if b == -1)
    if any(b < -1 or b > 1 for b in bal):
        return None
    if plus != minus:
        return None
    act = [v for v in range(n) if outdeg[v] + indeg[v] > 0]
    seen = {act[0]}
    st = [act[0]]
    while st:
        u = st.pop()
        for v, _ in adj[u]:
            if v not in seen:
                seen.add(v)
                st.append(v)
        for v in radj[u]:
            if v not in seen:
                seen.add(v)
                st.append(v)
    if len(seen) != len(act):
        return None
    start = -1
    for v in range(n):
        if bal[v] == 1:
            start = v
            break
    if start < 0:
        start = act[0]
    used = [False] * m
    it = [0] * n
    stack = [start]
    path = []
    while stack:
        u = stack[-1]
        advanced = False
        while it[u] < len(adj[u]):
            v, eid = adj[u][it[u]]
            it[u] += 1
            if not used[eid]:
                used[eid] = True
                stack.append(v)
                advanced = True
                break
        if not advanced:
            path.append(stack.pop())
    path.reverse()
    if len(path) != m + 1:
        return None
    return path


# ---------- 暴力對拍（指數回溯，僅用於小規模） ----------

def brute_undirected(n, edges):
    m = len(edges)
    if m == 0:
        return [0] if n >= 1 else None
    deg = [0] * n
    for u, v in edges:
        deg[u] += 1
        deg[v] += 1
    odd = [v for v in range(n) if deg[v] & 1]
    if len(odd) not in (0, 2):
        return None
    starts = odd if odd else [v for v in range(n) if deg[v] > 0]
    used = [False] * m

    def dfs(u, cnt):
        if cnt == m:
            return [u]
        for i, (a, b) in enumerate(edges):
            if used[i]:
                continue
            nxt = None
            if a == u:
                nxt = b
            elif b == u:
                nxt = a
            if nxt is None:
                continue
            used[i] = True
            r = dfs(nxt, cnt + 1)
            used[i] = False
            if r is not None:
                return [u] + r
        return None

    for s in starts:
        r = dfs(s, 0)
        if r is not None:
            return r
    return None


def brute_directed(n, edges):
    m = len(edges)
    if m == 0:
        return [0] if n >= 1 else None
    outdeg = [0] * n
    indeg = [0] * n
    for u, v in edges:
        outdeg[u] += 1
        indeg[v] += 1
    bal = [outdeg[i] - indeg[i] for i in range(n)]
    plus = sum(1 for b in bal if b == 1)
    minus = sum(1 for b in bal if b == -1)
    if any(b < -1 or b > 1 for b in bal):
        return None
    if plus != minus:
        return None
    starts = [v for v in range(n) if bal[v] == 1] or \
             [v for v in range(n) if outdeg[v] > 0]
    used = [False] * m

    def dfs(u, cnt):
        if cnt == m:
            return [u]
        for i, (a, b) in enumerate(edges):
            if not used[i] and a == u:
                used[i] = True
                r = dfs(b, cnt + 1)
                used[i] = False
                if r is not None:
                    return [u] + r
        return None

    for s in starts:
        r = dfs(s, 0)
        if r is not None:
            return r
    return None


def validate_path(n, edges, path, directed):
    """驗證一條宣稱的歐拉路徑確實每條邊恰好一次、每步合法。"""
    m = len(edges)
    if path is None:
        return True
    if len(path) != m + 1:
        return False
    cnt = {}
    for i in range(m):
        key = edges[i] if directed else tuple(sorted(edges[i]))
        cnt[key] = cnt.get(key, 0) + 1
    for i in range(m):
        key = (path[i], path[i + 1]) if directed else \
              tuple(sorted((path[i], path[i + 1])))
        if cnt.get(key, 0) == 0:
            return False
        cnt[key] -= 1
    return all(v == 0 for v in cnt.values())


def run_tests():
    seed = 20261107

    def rnd(mv):
        nonlocal seed
        seed = (seed * 1103515245 + 12345) & 0x7FFFFFFF
        return seed % mv

    # 固定用例（存在性 + 路徑合法性）
    fixed = [
        (1, 3, [(0, 1), (1, 2), (2, 0)]),           # 三角形 -> 迴路
        (1, 3, [(0, 1), (1, 2)]),                   # 鏈 -> 路徑
        (1, 4, [(0, 1), (0, 2), (0, 3), (0, 4)][:3]),  # 截短星形
        (1, 5, [(0, 1), (1, 2), (2, 0), (2, 3), (3, 4), (4, 2)]),  # 8 字 -> 迴路
        (1, 4, [(0, 1), (0, 2), (0, 3), (3, 1)]),   # 四奇點 -> 無解
        (1, 2, [(0, 0), (0, 1)]),                   # 自環 + 邊
        (1, 2, [(0, 1), (0, 1), (0, 1)]),           # 重邊
        (2, 3, [(0, 1), (1, 2), (2, 0)]),           # 有向三角形 -> 迴路
        (2, 3, [(0, 1), (1, 2), (0, 2)]),           # 0->1->2 路徑可排
        (2, 4, [(0, 1), (1, 2), (3, 2)]),           # bal 檢查失敗 -> 無解
        (2, 3, [(0, 1), (1, 0), (1, 2)]),           # 有向路徑
    ]
    for t, n, es in fixed:
        got = solve_undirected(n, es) if t == 1 else solve_directed(n, es)
        want = brute_undirected(n, es) if t == 1 else brute_directed(n, es)
        assert (got is None) == (want is None), (t, n, es, got, want)
        assert validate_path(n, es, got, t == 2), (t, n, es, got)

    # 隨機用例：n, m <= 8，與回溯暴力對拍存在性
    for _ in range(300):
        t = rnd(2) + 1
        n = rnd(7) + 1
        m = rnd(9)
        es = [(rnd(n), rnd(n)) for _ in range(m)]
        got = solve_undirected(n, es) if t == 1 else solve_directed(n, es)
        want = brute_undirected(n, es) if t == 1 else brute_directed(n, es)
        assert (got is None) == (want is None), (t, n, es, got, want)
        assert validate_path(n, es, got, t == 2), (t, n, es, got)

    # 大圖冒煙：隨機歐拉迴路（把環上點連成環 + 額外孤立點）
    for _ in range(20):
        k = rnd(500) + 2
        cyc = [(i, (i + 1) % k) for i in range(k)]
        p = solve_undirected(k, cyc)
        assert p is not None and len(p) == k + 1
        assert validate_path(k, cyc, p, False)

    print("all tests passed")


def main():
    data = sys.stdin.buffer.read().split()
    if not data:
        run_tests()
        return
    it = iter(data)
    t = int(next(it))
    n = int(next(it))
    m = int(next(it))
    edges = []
    for _ in range(m):
        u = int(next(it))
        v = int(next(it))
        edges.append((u, v))
    if t == 1:
        path = solve_undirected(n, edges)
    else:
        path = solve_directed(n, edges)
    if path is None:
        sys.stdout.write("-1\n")
    else:
        sys.stdout.write(" ".join(map(str, path)) + "\n")


if __name__ == "__main__":
    main()
