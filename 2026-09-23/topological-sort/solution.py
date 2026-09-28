"""拓撲排序（Topological Sort，Kahn 算法 + DFS 逆後序）

題意：給定一張 n 個點、m 條邊的有向圖（點編號 0..n-1），求一個拓撲序：
    一個把所有點排成一列的順序，使得對每條有向邊 u -> v，u 都排在 v 之前。
    若圖中存在環，則不存在拓撲序，報告無解。

思路（兩種等價的實現）：

    1. Kahn 算法（BFS / 剝洋蔥）：
       不斷把「入度爲 0 的點」拿掉，拿掉一個點就把它指向的邊的入度減一，
       於是新的入度爲 0 的點會出現。若最終拿掉了全部 n 個點，拿掉的先後順序
       就是一個拓撲序；若中途再也找不到入度爲 0 的點卻還剩點沒拿，說明
       剩下的點互相卡住——它們必定在一個環裏。
       用「最小堆」而不是普通隊列來挑選入度爲 0 的點，就得到**字典序最小**
       的拓撲序；用普通隊列則得到某個（取決於入邊順序的）合法拓撲序。

    2. DFS 逆後序：
       對圖做深度優先搜索，在「一個點的所有後繼都訪問完之後」把該點壓入棧，
       最後把棧倒過來輸出。正確性來自：DFS 的後序天然保證「後繼先於前驅完成」，
       反過來就是「前驅排在 successors 之前」。
       檢測環靠三色標記：遞歸棧上的點是灰色，若 DFS 走到灰色點說明有回邊，即有環。

    兩者的一致性與差異：Kahn 是迭代的（不會有遞歸深度問題），還能順手判環；
    DFS 版更短，但需要注意遞歸深度，且得到的序一般不是字典序最小的。
    本文件用同一個固定鄰接表（升序排列）驅動兩種實現，保證結果可復現。

輸入格式（stdin）：
    第一行：n m
    接下來 m 行：u v（一條有向邊 u -> v）
輸出格式（stdout）：
    一行：若存在拓撲序，輸出 n 個點編號（空格分隔，字典序最小的那個）；
          若存在環，輸出 -1
無 stdin 輸入時運行內置斷言測試。
"""

import heapq
import sys
from typing import List, Optional, Tuple


def kahn_topological_sort(n: int, edges: List[Tuple[int, int]]) -> Optional[List[int]]:
    """Kahn 算法 + 最小堆，返回字典序最小的拓撲序；有環返回 None。時間 O((n+m) log n)。"""
    adj: List[List[int]] = [[] for _ in range(n)]
    indeg = [0] * n
    for u, v in edges:
        adj[u].append(v)
        indeg[v] += 1

    # 最小堆：每步取編號最小的入度 0 點 -> 結果字典序最小
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

    # 拿掉的點不足 n 個 -> 剩下的點都還在環裏
    return order if len(order) == n else None


def dfs_topological_sort(n: int, edges: List[Tuple[int, int]]) -> Optional[List[int]]:
    """DFS 逆後序，返回一個合法拓撲序；有環返回 None。時間 O(n + m)。"""
    adj: List[List[int]] = [[] for _ in range(n)]
    for u, v in edges:
        adj[u].append(v)
    for lst in adj:
        lst.sort()  # 鄰接表升序，保證兩版實現結果可復現

    WHITE, GRAY, BLACK = 0, 1, 2
    color = [WHITE] * n
    post: List[int] = []

    # 迭代版 DFS：避免深圖把 Python 遞歸棧打爆
    for start in range(n):
        if color[start] != WHITE:
            continue
        stack: List[Tuple[int, int]] = [(start, 0)]  # (當前點, 下一條要走的邊下標)
        color[start] = GRAY
        while stack:
            u, idx = stack[-1]
            if idx < len(adj[u]):
                stack[-1] = (u, idx + 1)
                v = adj[u][idx]
                if color[v] == GRAY:
                    return None          # 回邊 -> 有環
                if color[v] == WHITE:
                    color[v] = GRAY
                    stack.append((v, 0))
            else:
                color[u] = BLACK
                post.append(u)           # 後繼都已完成，本點才算完成
                stack.pop()

    post.reverse()                        # 逆後序即拓撲序
    return post


def is_valid_topological_order(
    n: int, edges: List[Tuple[int, int]], order: Optional[List[int]]
) -> bool:
    """校驗：是 n 個點的一個排列，且每條邊 u -> v 都滿足 u 在 v 之前。"""
    if order is None or len(order) != n or sorted(order) != list(range(n)):
        return False
    pos = {v: i for i, v in enumerate(order)}
    return all(pos[u] < pos[v] for u, v in edges)


def brute_lex_topological_order(n: int, edges: List[Tuple[int, int]]) -> Optional[List[int]]:
    """對照用的全排列枚舉，返回字典序最小的拓撲序；無解返回 None。僅用於 n 很小的測試。"""
    from itertools import permutations

    best: Optional[Tuple[int, ...]] = None
    for perm in permutations(range(n)):
        pos = {v: i for i, v in enumerate(perm)}
        if all(pos[u] < pos[v] for u, v in edges):
            if best is None or perm < best:
                best = perm
    return list(best) if best is not None else None


def run_io(data: str) -> None:
    """按統一輸入輸出格式處理 stdin 數據。"""
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

    # 有環：0 -> 1 -> 2 -> 0
    cyc = [(0, 1), (1, 2), (2, 0)]
    assert kahn_topological_sort(3, cyc) is None
    assert dfs_topological_sort(3, cyc) is None
    assert brute_lex_topological_order(3, cyc) is None

    # 自環
    assert kahn_topological_sort(1, [(0, 0)]) is None
    assert dfs_topological_sort(1, [(0, 0)]) is None

    # 空圖：任意排列都是拓撲序。Kahn 用最小堆，得到字典序最小的 0 1 2 3；
    # DFS 版按 start = 0,1,2,... 依次完成，逆後序恰好是 3 2 1 0，同樣是合法拓撲序
    assert kahn_topological_sort(0, []) == []
    assert kahn_topological_sort(4, []) == [0, 1, 2, 3]
    assert dfs_topological_sort(4, []) == [3, 2, 1, 0]
    assert is_valid_topological_order(4, [], dfs_topological_sort(4, []))

    # 單點無邊
    assert kahn_topological_sort(1, []) == [0]

    # 重邊不應導致入度被多減（入度按「邊的條數」計數，去重才不會出錯）
    assert kahn_topological_sort(2, [(0, 1), (0, 1)]) == [0, 1]
    assert dfs_topological_sort(2, [(0, 1), (0, 1)]) == [0, 1]

    # 鏈：1 -> 2 -> 3，孤立點 0 排最前
    chain = [(1, 2), (2, 3)]
    assert kahn_topological_sort(4, chain) == [0, 1, 2, 3]

    # 環 + 無環部分混合：只要有一個環就整體無解
    assert kahn_topological_sort(4, [(0, 1), (1, 2), (2, 1), (0, 3)]) is None

    # 與全排列暴力解隨機對拍：校驗 Kahn 的字典序最小性、DFS 結果的合法性、
    # 以及兩者「有解 / 無解」的判斷必須一致
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
        assert (k is None) == (d is None) == (b is None)     # 有解/無解判斷一致
        if k is not None:
            assert is_valid_topological_order(n, edges, k)
            assert is_valid_topological_order(n, edges, d)   # DFS 結果也是合法拓撲序

    # 較大規模的鏈狀圖，驗證迭代版 DFS 不會爆遞歸棧
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
