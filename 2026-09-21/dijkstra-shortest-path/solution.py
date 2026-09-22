"""Dijkstra 单源最短路（非负权图，二叉堆实现）

返回从起点 start 到所有点的最短距离；不可达的点距离为 inf。
要求所有边权非负——若存在负权边，该算法的贪心正确性不成立，应改用 Bellman-Ford / SPFA。
"""

import heapq
from typing import List, Tuple, Optional

INF = float("inf")


def dijkstra(n: int, graph: List[List[Tuple[int, float]]], start: int) -> List[float]:
    """graph[u] = [(v, w), ...] 为邻接表；返回长度 n 的距离数组。"""
    dist = [INF] * n
    dist[start] = 0
    # 小根堆：元素为 (当前距离, 节点)。允许同一节点多次入堆，靠 visited 跳过过期条目
    heap: List[Tuple[float, int]] = [(0.0, start)]
    visited = [False] * n

    while heap:
        d, u = heapq.heappop(heap)
        if visited[u]:
            continue          # 该节点已被更短的距离定型过，堆里这条是过期记录
        visited[u] = True
        for v, w in graph[u]:
            if visited[v]:
                continue
            nd = d + w
            if nd < dist[v]:  # 松弛
                dist[v] = nd
                heapq.heappush(heap, (nd, v))
    return dist


def shortest_path(n: int, edges: List[Tuple[int, int, float]],
                  start: int, target: int) -> Optional[float]:
    """便捷接口：传入无向边列表 (u, v, w)，返回 start→target 的最短距离，不可达返回 None。"""
    graph: List[List[Tuple[int, float]]] = [[] for _ in range(n)]
    for u, v, w in edges:
        graph[u].append((v, w))
        graph[v].append((u, w))
    d = dijkstra(n, graph, start)
    return None if d[target] == INF else d[target]


if __name__ == "__main__":
    n = 5
    edges = [
        (0, 1, 4), (0, 2, 1),
        (2, 1, 2), (2, 3, 5),
        (1, 3, 1), (3, 4, 3),
    ]
    assert shortest_path(n, edges, 0, 4) == 7   # 0→2→1→3→4 = 1+2+1+3
    assert shortest_path(n, edges, 0, 1) == 3   # 0→2→1 = 1+2，比直连 4 更短
    assert shortest_path(n, edges, 4, 4) == 0

    graph = [[] for _ in range(n)]
    for u, v, w in edges:
        graph[u].append((v, w))
        graph[v].append((u, w))
    dist = dijkstra(n, graph, 0)
    assert dist == [0, 3, 1, 4, 7]

    # 不连通的图：节点 2 不可达
    g2 = [[(1, 1.0)], []]
    assert dijkstra(2, g2, 0)[1] == 1.0
    assert dijkstra(2, [[], []], 0)[1] == INF
    print("all tests passed")
