"""Bellman-Ford 单源最短路（支持负权边 + 负环检测 + 路径还原）

题意：
    给定一张 n 个点、m 条边的**带权有向图**（边权可以为负），以及源点 s。要求：
      1) 判断从 s 出发**是否能到达某个负权环**（环上边权之和 < 0）；
      2) 若不存在这样的负环，输出 s 到每个点的最短距离，不可达输出 INF。

思路：
    Dijkstra 依赖「已经出队的点距离不再变小」这个贪心性质，一旦有负权边就不成立，
    所以负权图要换一套思路。Bellman-Ford 的出发点是一个朴素事实：

        **一条最短路最多经过 n-1 条边**（再多就一定绕了环；正环/零环可以删掉，
         负环则根本不存在「最短路」）。

    于是算法就非常简单：把所有边**整体松弛** n-1 轮。第 k 轮结束后，
    所有「不超过 k 条边」的最短路径都已经求出来了。

    再做第 n 轮：如果还有边能被松弛，就说明存在一条含 n 条边还能不断变短的路径，
    也就是绕了负环 —— 这就是**负环检测**，也是 Bellman-Ford 相比 Dijkstra 的核心能力。

    几个要点：
      1) 只松弛 `dist[u]` 有限的边，否则 `INF + w` 会污染结果；
      2) 某一轮一条边都没松弛就可以提前退出（`changed` 标记），实际运行远快于 n 轮；
      3) 记录前驱 `pre[v]` 就能还原最短路；无负环时前驱构成一棵最短路树；
      4) **队列优化版（SPFA）**：只有上一轮被更新过的点，它的出边才可能继续松弛，
         用队列维护这些点即可。期望快很多，最坏仍是 O(n·m)。
         负环判定换成「某个点入队次数 > n」。

    注意：负环必须是**从 s 可达**的才会被检测到；图里另一头的负环与 s 无关，
    不影响 s 到其他点的最短路（本文件有专门的用例验证这一点）。

输入格式（stdin，所有数字按空白分隔即可）：
    n m s
    u v w               （共 m 行；n = 0 时这两部分都省略）
输出格式（stdout）：
    第 1 行：1 表示存在从 s 可达的负环，0 表示不存在
    第 2 行（仅当第 1 行为 0）：dist[0] ... dist[n-1]，空格分隔，不可达输出 INF
无 stdin 输入时运行内置断言测试并输出 `all tests passed`。
"""

import sys
from collections import deque
from typing import Iterator, List, Tuple

import random

INF = float("inf")
INF_STR = "INF"


def bellman_ford(n: int, edges: List[Tuple[int, int, int]], s: int):
    """Bellman-Ford。返回 (has_neg_cycle, dist, pre)。

    dist[v] 为 s 到 v 的最短距离（不可达为 INF）；pre[v] 为最短路上的前驱（-1 表示无）。
    时间 O(n·m)，空间 O(n)。
    """
    dist = [INF] * n
    pre = [-1] * n
    if 0 <= s < n:
        dist[s] = 0
    has_neg = False
    for it in range(n):
        changed = False
        for u, v, w in edges:
            # 只从已可达的点往外松弛，避免 INF + w 污染结果
            if dist[u] != INF and dist[u] + w < dist[v]:
                dist[v] = dist[u] + w
                pre[v] = u
                changed = True
        if not changed:
            break                      # 这一轮没人被更新，后面也不可能再变
        if it == n - 1:
            has_neg = True             # 第 n 轮还能松弛 → 绕了负环
    return has_neg, dist, pre


def spfa(n: int, edges: List[Tuple[int, int, int]], s: int):
    """Bellman-Ford 的队列优化版（SPFA）。返回 (has_neg_cycle, dist)。

    只有距离被更新过的点才需要再次松弛其出边。期望远快于朴素版，最坏仍 O(n·m)。
    """
    adj: List[List[Tuple[int, int]]] = [[] for _ in range(n)]
    for u, v, w in edges:
        adj[u].append((v, w))

    dist = [INF] * n
    inq = [False] * n
    cnt = [0] * n
    if 0 <= s < n:
        dist[s] = 0
        inq[s] = True
        cnt[s] = 1

    q = deque([s]) if 0 <= s < n else deque()
    while q:
        u = q.popleft()
        inq[u] = False
        if dist[u] == INF:
            continue
        for v, w in adj[u]:
            if dist[u] + w < dist[v]:
                dist[v] = dist[u] + w
                if not inq[v]:
                    inq[v] = True
                    cnt[v] += 1
                    if cnt[v] > n:     # 入队超过 n 次 → 有负环
                        return True, dist
                    q.append(v)
    return False, dist


def build_path(pre: List[int], s: int, t: int) -> List[int]:
    """沿前驱数组还原 s → t 的最短路；不可达返回空列表。时间 O(路径长度)。"""
    if not pre or not (0 <= t < len(pre)):
        return []
    path = []
    cur = t
    while cur != -1:
        path.append(cur)
        if cur == s:
            return path[::-1]
        cur = pre[cur]
    return []


def path_weight(edges: List[Tuple[int, int, int]], path: List[int]):
    """校验路径：返回 (是否每条边都存在, 路径总权重)。"""
    best = {}
    for u, v, w in edges:
        if (u, v) not in best or w < best[(u, v)]:
            best[(u, v)] = w
    total = 0
    for i in range(len(path) - 1):
        e = (path[i], path[i + 1])
        if e not in best:
            return False, 0
        total += best[e]
    return True, total


# ---------------- 对照用的其他算法 ----------------

def dijkstra(n: int, edges: List[Tuple[int, int, int]], s: int) -> List[float]:
    """O(n^2) 版 Dijkstra，仅用于**非负权**图的对拍。"""
    adj: List[List[Tuple[int, int]]] = [[] for _ in range(n)]
    for u, v, w in edges:
        adj[u].append((v, w))
    dist = [INF] * n
    used = [False] * n
    if 0 <= s < n:
        dist[s] = 0
    for _ in range(n):
        u = -1
        for i in range(n):
            if not used[i] and dist[i] != INF and (u == -1 or dist[i] < dist[u]):
                u = i
        if u == -1:
            break
        used[u] = True
        for v, w in adj[u]:
            if dist[u] + w < dist[v]:
                dist[v] = dist[u] + w
    return dist


def floyd_warshall(n: int, edges: List[Tuple[int, int, int]]):
    """Floyd-Warshall 全源最短路。返回 (图中是否存在任意负环, 距离矩阵)。时间 O(n^3)。"""
    d = [[INF] * n for _ in range(n)]
    for i in range(n):
        d[i][i] = 0
    for u, v, w in edges:
        if 0 <= u < n and 0 <= v < n and w < d[u][v]:
            d[u][v] = w
    for k in range(n):
        for i in range(n):
            if d[i][k] == INF:
                continue
            for j in range(n):
                if d[k][j] != INF and d[i][k] + d[k][j] < d[i][j]:
                    d[i][j] = d[i][k] + d[k][j]
    has_neg = any(d[i][i] < 0 for i in range(n))
    return has_neg, d


# ---------------- IO ----------------

def _next_int(it: Iterator[str], default: int = 0) -> int:
    """取下一个整数；输入被截断时用默认值兜底，避免直接抛异常。"""
    try:
        return int(next(it))
    except (StopIteration, ValueError):
        return default


def run_io(data: str) -> None:
    """按统一输入输出格式处理 stdin 数据。"""
    it = iter(data.split())
    n = _next_int(it)
    m = _next_int(it)
    s = _next_int(it)
    edges = []
    for _ in range(m):
        u = _next_int(it)
        v = _next_int(it)
        w = _next_int(it)
        if 0 <= u < n and 0 <= v < n:
            edges.append((u, v, w))

    has_neg, dist, _ = bellman_ford(n, edges, s)
    print(1 if has_neg else 0)
    if not has_neg:
        print(" ".join(INF_STR if d == INF else str(d) for d in dist))   # n = 0 时输出空行


def run_tests() -> None:
    # 示例一：含负权边但无负环的图
    #   0 -(4)-> 1 -(2)-> 3 -(2)-> 4 -(1)-> 1（正环，不影响）
    #   0 -(2)-> 2 -(-3)-> 1，所以 0→1 走 0→2→1 只需 -1
    e1 = [(0, 1, 4), (0, 2, 2), (2, 1, -3), (1, 3, 2), (2, 3, 5), (3, 4, 2), (4, 1, 1)]
    neg, dist, pre = bellman_ford(5, e1, 0)
    assert neg is False
    assert dist == [0, -1, 2, 1, 3]
    assert spfa(5, e1, 0) == (False, dist)
    ok, w = path_weight(e1, build_path(pre, 0, 1))
    assert ok and w == -1 and build_path(pre, 0, 1) == [0, 2, 1]

    # 示例二：负环 1 -> 2 -> 1，权值和 -2，且从 0 可达
    e2 = [(0, 1, 1), (1, 2, -3), (2, 1, 1)]
    neg2, dist2, _ = bellman_ford(3, e2, 0)
    assert neg2 is True
    assert spfa(3, e2, 0)[0] is True

    # 负环存在但**从 s 不可达**：不影响 s 的最短路，检测也应当报告「无」
    e3 = [(0, 1, 1), (2, 3, -5), (3, 2, 2)]
    neg3, dist3, _ = bellman_ford(4, e3, 0)
    assert neg3 is False
    assert dist3[0] == 0 and dist3[1] == 1
    assert dist3[2] == INF and dist3[3] == INF
    assert spfa(4, e3, 0)[0] is False

    # 自环负权 → 负环
    neg4, _, _ = bellman_ford(1, [(0, 0, -1)], 0)
    assert neg4 is True
    # 自环零权 / 正权 → 不是负环
    assert bellman_ford(1, [(0, 0, 0)], 0)[0] is False
    assert bellman_ford(1, [(0, 0, 5)], 0)[0] is False

    # 边界：单点无边 / 无边图 / 源点即终点
    assert bellman_ford(1, [], 0) == (False, [0], [-1])
    assert bellman_ford(3, [], 0)[1] == [0, INF, INF]
    assert bellman_ford(0, [], 0) == (False, [], [])

    # 重边取最小：两条 0→1，权值 7 和 3
    neg5, dist5, _ = bellman_ford(2, [(0, 1, 7), (0, 1, 3)], 0)
    assert neg5 is False and dist5 == [0, 3]

    # 路径还原：链状图 0→1→2→3
    e6 = [(0, 1, 2), (1, 2, 3), (2, 3, 4)]
    neg6, dist6, pre6 = bellman_ford(4, e6, 0)
    assert dist6 == [0, 2, 5, 9]
    assert build_path(pre6, 0, 3) == [0, 1, 2, 3]
    assert build_path(pre6, 0, 0) == [0]
    # 不可达的点还原不出路径
    neg7, dist7, pre7 = bellman_ford(4, [(2, 3, 1)], 0)
    assert build_path(pre7, 0, 3) == []

    random.seed(20260928)

    # 随机对拍一：**非负权**图，与 O(n^2) Dijkstra 比对
    for _ in range(400):
        n = random.randint(1, 8)
        m = random.randint(0, n * 2)
        edges = [(random.randrange(n), random.randrange(n), random.randint(0, 9))
                 for _ in range(m)]
        s = random.randrange(n)
        neg, dist, _ = bellman_ford(n, edges, s)
        assert neg is False                     # 非负权不可能有负环
        assert dist == dijkstra(n, edges, s)

    # 随机对拍二：**允许负权**，与 SPFA 互相印证，并用 Floyd 与最优性条件兜底
    for _ in range(400):
        n = random.randint(1, 7)
        m = random.randint(0, n * 2)
        edges = [(random.randrange(n), random.randrange(n), random.randint(-6, 9))
                 for _ in range(m)]
        s = random.randrange(n)

        neg, dist, pre = bellman_ford(n, edges, s)
        neg_s, dist_s = spfa(n, edges, s)
        assert neg == neg_s                     # 两种判据必须一致

        if neg:
            continue                            # 有负环时距离无意义，只校验判据一致

        assert dist == dist_s                   # 两版距离必须一致
        assert dist[s] == 0

        # 最优性条件：无负环时不存在还能被松弛的边（三角不等式成立）
        for u, v, w in edges:
            if dist[u] != INF:
                assert dist[u] + w >= dist[v]

        # 每个可达点的距离都必须对应一条真实存在的、权重相等的路径
        for t in range(n):
            if dist[t] == INF:
                assert build_path(pre, s, t) == []
                continue
            p = build_path(pre, s, t)
            assert p and p[0] == s and p[-1] == t
            ok, w = path_weight(edges, p)
            assert ok and w == dist[t]

        # 与 Floyd-Warshall 交叉验证（只在整张图都没有负环时才可比）
        has_neg_any, mat = floyd_warshall(n, edges)
        if not has_neg_any:
            for t in range(n):
                assert dist[t] == mat[s][t]

    print("all tests passed")


if __name__ == "__main__":
    raw = sys.stdin.read() if not sys.stdin.isatty() else ""
    if raw.strip():
        run_io(raw)
    else:
        run_tests()
