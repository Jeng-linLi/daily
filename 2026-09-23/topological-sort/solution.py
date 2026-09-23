"""拓扑排序（Topological Sort，Kahn 算法 + DFS 逆后序）

题意：给定一张 n 个点、m 条边的有向图（点编号 0..n-1），求一个拓扑序：
    一个把所有点排成一列的顺序，使得对每条有向边 u -> v，u 都排在 v 之前。
    若图中存在环，则不存在拓扑序，报告无解。

思路（两种等价的实现）：

    1. Kahn 算法（BFS / 剥洋葱）：
       不断把「入度为 0 的点」拿掉，拿掉一个点就把它指向的边的入度减一，
       于是新的入度为 0 的点会出现。若最终拿掉了全部 n 个点，拿掉的先后顺序
       就是一个拓扑序；若中途再也找不到入度为 0 的点却还剩点没拿，说明
       剩下的点互相卡住——它们必定在一个环里。
       用「最小堆」而不是普通队列来挑选入度为 0 的点，就得到**字典序最小**
       的拓扑序；用普通队列则得到某个（取决于入边顺序的）合法拓扑序。

    2. DFS 逆后序：
       对图做深度优先搜索，在「一个点的所有后继都访问完之后」把该点压入栈，
       最后把栈倒过来输出。正确性来自：DFS 的后序天然保证「后继先于前驱完成」，
       反过来就是「前驱排在 successors 之前」。
       检测环靠三色标记：递归栈上的点是灰色，若 DFS 走到灰色点说明有回边，即有环。

    两者的一致性与差异：Kahn 是迭代的（不会有递归深度问题），还能顺手判环；
    DFS 版更短，但需要注意递归深度，且得到的序一般不是字典序最小的。
    本文件用同一个固定邻接表（升序排列）驱动两种实现，保证结果可复现。

输入格式（stdin）：
    第一行：n m
    接下来 m 行：u v（一条有向边 u -> v）
输出格式（stdout）：
    一行：若存在拓扑序，输出 n 个点编号（空格分隔，字典序最小的那个）；
          若存在环，输出 -1
无 stdin 输入时运行内置断言测试。
"""

import heapq
import sys
from typing import List, Optional, Tuple


def kahn_topological_sort(n: int, edges: List[Tuple[int, int]]) -> Optional[List[int]]:
    """Kahn 算法 + 最小堆，返回字典序最小的拓扑序；有环返回 None。时间 O((n+m) log n)。"""
    adj: List[List[int]] = [[] for _ in range(n)]
    indeg = [0] * n
    for u, v in edges:
        adj[u].append(v)
        indeg[v] += 1

    # 最小堆：每步取编号最小的入度 0 点 -> 结果字典序最小
    heap = [i for i in range(n) if indeg[i] == 0]
    heapq.heapify(heap)

    order: List[int] = []
    while heap:
        u = heapq.heappop(heap)
        order.append(u)
        for v in adj[u]:
            indeg[v] -= 1
            if indeg[v] == 0:
                heapq.heappush(heap, v)

    # 拿掉的点不足 n 个 -> 剩下的点都还在环里
    return order if len(order) == n else None


def dfs_topological_sort(n: int, edges: List[Tuple[int, int]]) -> Optional[List[int]]:
    """DFS 逆后序，返回一个合法拓扑序；有环返回 None。时间 O(n + m)。"""
    adj: List[List[int]] = [[] for _ in range(n)]
    for u, v in edges:
        adj[u].append(v)
    for lst in adj:
        lst.sort()  # 邻接表升序，保证两版实现结果可复现

    WHITE, GRAY, BLACK = 0, 1, 2
    color = [WHITE] * n
    post: List[int] = []

    # 迭代版 DFS：避免深图把 Python 递归栈打爆
    for start in range(n):
        if color[start] != WHITE:
            continue
        stack: List[Tuple[int, int]] = [(start, 0)]  # (当前点, 下一条要走的边下标)
        color[start] = GRAY
        while stack:
            u, idx = stack[-1]
            if idx < len(adj[u]):
                stack[-1] = (u, idx + 1)
                v = adj[u][idx]
                if color[v] == GRAY:
                    return None          # 回边 -> 有环
                if color[v] == WHITE:
                    color[v] = GRAY
                    stack.append((v, 0))
            else:
                color[u] = BLACK
                post.append(u)           # 后继都已完成，本点才算完成
                stack.pop()

    post.reverse()                        # 逆后序即拓扑序
    return post


def is_valid_topological_order(
    n: int, edges: List[Tuple[int, int]], order: Optional[List[int]]
) -> bool:
    """校验：是 n 个点的一个排列，且每条边 u -> v 都满足 u 在 v 之前。"""
    if order is None or len(order) != n or sorted(order) != list(range(n)):
        return False
    pos = {v: i for i, v in enumerate(order)}
    return all(pos[u] < pos[v] for u, v in edges)


def brute_lex_topological_order(n: int, edges: List[Tuple[int, int]]) -> Optional[List[int]]:
    """对照用的全排列枚举，返回字典序最小的拓扑序；无解返回 None。仅用于 n 很小的测试。"""
    from itertools import permutations

    best: Optional[Tuple[int, ...]] = None
    for perm in permutations(range(n)):
        pos = {v: i for i, v in enumerate(perm)}
        if all(pos[u] < pos[v] for u, v in edges):
            if best is None or perm < best:
                best = perm
    return list(best) if best is not None else None


def run_io(data: str) -> None:
    """按统一输入输出格式处理 stdin 数据。"""
    tokens = data.split()
    if not tokens:
        return
    n, m = int(tokens[0]), int(tokens[1])
    edges: List[Tuple[int, int]] = []
    idx = 2
    for _ in range(m):
        edges.append((int(tokens[idx]), int(tokens[idx + 1])))
        idx += 2
    order = kahn_topological_sort(n, edges)
    print(" ".join(str(x) for x in order) if order is not None else -1)


def run_tests() -> None:
    # README 中的示例：0 1 2 3 5 4
    edges = [(0, 1), (0, 2), (1, 3), (2, 3), (2, 4), (5, 4)]
    assert kahn_topological_sort(6, edges) == [0, 1, 2, 3, 5, 4]
    assert brute_lex_topological_order(6, edges) == [0, 1, 2, 3, 5, 4]
    assert is_valid_topological_order(6, edges, dfs_topological_sort(6, edges))

    # 有环：0 -> 1 -> 2 -> 0
    cyc = [(0, 1), (1, 2), (2, 0)]
    assert kahn_topological_sort(3, cyc) is None
    assert dfs_topological_sort(3, cyc) is None
    assert brute_lex_topological_order(3, cyc) is None

    # 自环
    assert kahn_topological_sort(1, [(0, 0)]) is None
    assert dfs_topological_sort(1, [(0, 0)]) is None

    # 空图：任意排列都是拓扑序。Kahn 用最小堆，得到字典序最小的 0 1 2 3；
    # DFS 版按 start = 0,1,2,... 依次完成，逆后序恰好是 3 2 1 0，同样是合法拓扑序
    assert kahn_topological_sort(0, []) == []
    assert kahn_topological_sort(4, []) == [0, 1, 2, 3]
    assert dfs_topological_sort(4, []) == [3, 2, 1, 0]
    assert is_valid_topological_order(4, [], dfs_topological_sort(4, []))

    # 单点无边
    assert kahn_topological_sort(1, []) == [0]

    # 重边不应导致入度被多减（入度按「边的条数」计数，去重才不会出错）
    assert kahn_topological_sort(2, [(0, 1), (0, 1)]) == [0, 1]
    assert dfs_topological_sort(2, [(0, 1), (0, 1)]) == [0, 1]

    # 链：1 -> 2 -> 3，孤立点 0 排最前
    chain = [(1, 2), (2, 3)]
    assert kahn_topological_sort(4, chain) == [0, 1, 2, 3]

    # 环 + 无环部分混合：只要有一个环就整体无解
    assert kahn_topological_sort(4, [(0, 1), (1, 2), (2, 1), (0, 3)]) is None

    # 与全排列暴力解随机对拍：校验 Kahn 的字典序最小性、DFS 结果的合法性、
    # 以及两者「有解 / 无解」的判断必须一致
    import random

    random.seed(20260923)
    for _ in range(200):
        n = random.randint(1, 6)
        edges = []
        for u in range(n):
            for v in range(n):
                if u != v and random.random() < 0.25:
                    edges.append((u, v))
        k = kahn_topological_sort(n, edges)
        d = dfs_topological_sort(n, edges)
        b = brute_lex_topological_order(n, edges)
        assert k == b                                        # Kahn = 字典序最小
        assert (k is None) == (d is None) == (b is None)     # 有解/无解判断一致
        if k is not None:
            assert is_valid_topological_order(n, edges, k)
            assert is_valid_topological_order(n, edges, d)   # DFS 结果也是合法拓扑序

    # 较大规模的链状图，验证迭代版 DFS 不会爆递归栈
    big_n = 20000
    big_edges = [(i, i + 1) for i in range(big_n - 1)]
    assert kahn_topological_sort(big_n, big_edges) == list(range(big_n))
    assert dfs_topological_sort(big_n, big_edges) == list(range(big_n))

    print("all tests passed")


if __name__ == "__main__":
    raw = sys.stdin.read() if not sys.stdin.isatty() else ""
    if raw.strip():
        run_io(raw)
    else:
        run_tests()
