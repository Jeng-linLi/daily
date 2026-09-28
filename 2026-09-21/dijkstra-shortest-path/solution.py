"""Dijkstra 單源最短路（非負權圖，二叉堆實現）

返回從起點 start 到所有點的最短距離；不可達的點距離爲 inf。
要求所有邊權非負——若存在負權邊，該算法的貪心正確性不成立，應改用 Bellman-Ford / SPFA。
"""

import heapq
from typing import List, Tuple, Optional

INF = float("inf")


def dijkstra(n: int, graph: List[List[Tuple[int, float]]], start: int) -> List[float]:
    """graph[u] = [(v, w), ...] 爲鄰接表；返回長度 n 的距離數組。"""
    dist = [INF] * n
    dist[start] = 0
    # 小根堆：元素爲 (當前距離, 節點)。允許同一節點多次入堆，靠 visited 跳過過期條目
    heap: List[Tuple[float, int]] = [(0.0, start)]
    visited = [False] * n

    while heap:
        d, u = heapq.heappop(heap)
        if visited[u]:
            continue          # 該節點已被更短的距離定型過，堆裏這條是過期記錄
        visited[u] = True
        for v, w in graph[u]:
            if visited[v]:
                continue
            nd = d + w
            if nd < dist[v]:  # 鬆弛
                dist[v] = nd
                heapq.heappush(heap, (nd, v))
    return dist


def shortest_path(n: int, edges: List[Tuple[int, int, float]],
                  start: int, target: int) -> Optional[float]:
    """便捷接口：傳入無向邊列表 (u, v, w)，返回 start→target 的最短距離，不可達返回 None。"""
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
    assert shortest_path(n, edges, 0, 1) == 3   # 0→2→1 = 1+2，比直連 4 更短
    assert shortest_path(n, edges, 4, 4) == 0

    graph = [[] for _ in range(n)]
    for u, v, w in edges:
        graph[u].append((v, w))
        graph[v].append((u, w))
    dist = dijkstra(n, graph, 0)
    assert dist == [0, 3, 1, 4, 7]

    # 不連通的圖：節點 2 不可達
    g2 = [[(1, 1.0)], []]
    assert dijkstra(2, g2, 0)[1] == 1.0
    assert dijkstra(2, [[], []], 0)[1] == INF
    print("all tests passed")
