"""Bellman-Ford 單源最短路（支持負權邊 + 負環檢測 + 路徑還原）

題意：
    給定一張 n 個點、m 條邊的**帶權有向圖**（邊權可以爲負），以及源點 s。要求：
      1) 判斷從 s 出發**是否能到達某個負權環**（環上邊權之和 < 0）；
      2) 若不存在這樣的負環，輸出 s 到每個點的最短距離，不可達輸出 INF。

思路：
    Dijkstra 依賴「已經出隊的點距離不再變小」這個貪心性質，一旦有負權邊就不成立，
    所以負權圖要換一套思路。Bellman-Ford 的出發點是一個樸素事實：

        **一條最短路最多經過 n-1 條邊**（再多就一定繞了環；正環/零環可以刪掉，
         負環則根本不存在「最短路」）。

    於是算法就非常簡單：把所有邊**整體鬆弛** n-1 輪。第 k 輪結束後，
    所有「不超過 k 條邊」的最短路徑都已經求出來了。

    再做第 n 輪：如果還有邊能被鬆弛，就說明存在一條含 n 條邊還能不斷變短的路徑，
    也就是繞了負環 —— 這就是**負環檢測**，也是 Bellman-Ford 相比 Dijkstra 的核心能力。

    幾個要點：
      1) 只鬆弛 `dist[u]` 有限的邊，否則 `INF + w` 會污染結果；
      2) 某一輪一條邊都沒鬆弛就可以提前退出（`changed` 標記），實際運行遠快於 n 輪；
      3) 記錄前驅 `pre[v]` 就能還原最短路；無負環時前驅構成一棵最短路樹；
      4) **隊列優化版（SPFA）**：只有上一輪被更新過的點，它的出邊才可能繼續鬆弛，
         用隊列維護這些點即可。期望快很多，最壞仍是 O(n·m)。
         負環判定換成「某個點入隊次數 > n」。

    注意：負環必須是**從 s 可達**的才會被檢測到；圖裏另一頭的負環與 s 無關，
    不影響 s 到其他點的最短路（本文件有專門的用例驗證這一點）。

輸入格式（stdin，所有數字按空白分隔即可）：
    n m s
    u v w               （共 m 行；n = 0 時這兩部分都省略）
輸出格式（stdout）：
    第 1 行：1 表示存在從 s 可達的負環，0 表示不存在
    第 2 行（僅當第 1 行爲 0）：dist[0] ... dist[n-1]，空格分隔，不可達輸出 INF
無 stdin 輸入時運行內置斷言測試並輸出 `all tests passed`。
"""

import sys
from collections import deque
from typing import Iterator, List, Tuple

import random

INF = float("inf")
INF_STR = "INF"


def bellman_ford(n: int, edges: List[Tuple[int, int, int]], s: int):
    """Bellman-Ford。返回 (has_neg_cycle, dist, pre)。

    dist[v] 爲 s 到 v 的最短距離（不可達爲 INF）；pre[v] 爲最短路上的前驅（-1 表示無）。
    時間 O(n·m)，空間 O(n)。
    """
    dist = [INF] * n
    pre = [-1] * n
    if 0 <= s < n:
        dist[s] = 0
    has_neg = False
    for it in range(n):
        changed = False
        for u, v, w in edges:
            # 只從已可達的點往外鬆弛，避免 INF + w 污染結果
            if dist[u] != INF and dist[u] + w < dist[v]:
                dist[v] = dist[u] + w
                pre[v] = u
                changed = True
        if not changed:
            break                      # 這一輪沒人被更新，後面也不可能再變
        if it == n - 1:
            has_neg = True             # 第 n 輪還能鬆弛 → 繞了負環
    return has_neg, dist, pre


def spfa(n: int, edges: List[Tuple[int, int, int]], s: int):
    """Bellman-Ford 的隊列優化版（SPFA）。返回 (has_neg_cycle, dist)。

    只有距離被更新過的點才需要再次鬆弛其出邊。期望遠快於樸素版，最壞仍 O(n·m)。
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
                    if cnt[v] > n:     # 入隊超過 n 次 → 有負環
                        return True, dist
                    q.append(v)
    return False, dist


def build_path(pre: List[int], s: int, t: int) -> List[int]:
    """沿前驅數組還原 s → t 的最短路；不可達返回空列表。時間 O(路徑長度)。"""
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
    """校驗路徑：返回 (是否每條邊都存在, 路徑總權重)。"""
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


# ---------------- 對照用的其他算法 ----------------

def dijkstra(n: int, edges: List[Tuple[int, int, int]], s: int) -> List[float]:
    """O(n^2) 版 Dijkstra，僅用於**非負權**圖的對拍。"""
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
    """Floyd-Warshall 全源最短路。返回 (圖中是否存在任意負環, 距離矩陣)。時間 O(n^3)。"""
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
    """取下一個整數；輸入被截斷時用默認值兜底，避免直接拋異常。"""
    try:
        return int(next(it))
    except (StopIteration, ValueError):
        return default


def run_io(data: str) -> None:
    """按統一輸入輸出格式處理 stdin 數據。"""
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
        print(" ".join(INF_STR if d == INF else str(d) for d in dist))   # n = 0 時輸出空行


def run_tests() -> None:
    # 示例一：含負權邊但無負環的圖
    #   0 -(4)-> 1 -(2)-> 3 -(2)-> 4 -(1)-> 1（正環，不影響）
    #   0 -(2)-> 2 -(-3)-> 1，所以 0→1 走 0→2→1 只需 -1
    e1 = [(0, 1, 4), (0, 2, 2), (2, 1, -3), (1, 3, 2), (2, 3, 5), (3, 4, 2), (4, 1, 1)]
    neg, dist, pre = bellman_ford(5, e1, 0)
    assert neg is False
    assert dist == [0, -1, 2, 1, 3]
    assert spfa(5, e1, 0) == (False, dist)
    ok, w = path_weight(e1, build_path(pre, 0, 1))
    assert ok and w == -1 and build_path(pre, 0, 1) == [0, 2, 1]

    # 示例二：負環 1 -> 2 -> 1，權值和 -2，且從 0 可達
    e2 = [(0, 1, 1), (1, 2, -3), (2, 1, 1)]
    neg2, dist2, _ = bellman_ford(3, e2, 0)
    assert neg2 is True
    assert spfa(3, e2, 0)[0] is True

    # 負環存在但**從 s 不可達**：不影響 s 的最短路，檢測也應當報告「無」
    e3 = [(0, 1, 1), (2, 3, -5), (3, 2, 2)]
    neg3, dist3, _ = bellman_ford(4, e3, 0)
    assert neg3 is False
    assert dist3[0] == 0 and dist3[1] == 1
    assert dist3[2] == INF and dist3[3] == INF
    assert spfa(4, e3, 0)[0] is False

    # 自環負權 → 負環
    neg4, _, _ = bellman_ford(1, [(0, 0, -1)], 0)
    assert neg4 is True
    # 自環零權 / 正權 → 不是負環
    assert bellman_ford(1, [(0, 0, 0)], 0)[0] is False
    assert bellman_ford(1, [(0, 0, 5)], 0)[0] is False

    # 邊界：單點無邊 / 無邊圖 / 源點即終點
    assert bellman_ford(1, [], 0) == (False, [0], [-1])
    assert bellman_ford(3, [], 0)[1] == [0, INF, INF]
    assert bellman_ford(0, [], 0) == (False, [], [])

    # 重邊取最小：兩條 0→1，權值 7 和 3
    neg5, dist5, _ = bellman_ford(2, [(0, 1, 7), (0, 1, 3)], 0)
    assert neg5 is False and dist5 == [0, 3]

    # 路徑還原：鏈狀圖 0→1→2→3
    e6 = [(0, 1, 2), (1, 2, 3), (2, 3, 4)]
    neg6, dist6, pre6 = bellman_ford(4, e6, 0)
    assert dist6 == [0, 2, 5, 9]
    assert build_path(pre6, 0, 3) == [0, 1, 2, 3]
    assert build_path(pre6, 0, 0) == [0]
    # 不可達的點還原不出路徑
    neg7, dist7, pre7 = bellman_ford(4, [(2, 3, 1)], 0)
    assert build_path(pre7, 0, 3) == []

    random.seed(20260928)

    # 隨機對拍一：**非負權**圖，與 O(n^2) Dijkstra 比對
    for _ in range(400):
        n = random.randint(1, 8)
        m = random.randint(0, n * 2)
        edges = [(random.randrange(n), random.randrange(n), random.randint(0, 9))
                 for _ in range(m)]
        s = random.randrange(n)
        neg, dist, _ = bellman_ford(n, edges, s)
        assert neg is False                     # 非負權不可能有負環
        assert dist == dijkstra(n, edges, s)

    # 隨機對拍二：**允許負權**，與 SPFA 互相印證，並用 Floyd 與最優性條件兜底
    for _ in range(400):
        n = random.randint(1, 7)
        m = random.randint(0, n * 2)
        edges = [(random.randrange(n), random.randrange(n), random.randint(-6, 9))
                 for _ in range(m)]
        s = random.randrange(n)

        neg, dist, pre = bellman_ford(n, edges, s)
        neg_s, dist_s = spfa(n, edges, s)
        assert neg == neg_s                     # 兩種判據必須一致

        if neg:
            continue                            # 有負環時距離無意義，只校驗判據一致

        assert dist == dist_s                   # 兩版距離必須一致
        assert dist[s] == 0

        # 最優性條件：無負環時不存在還能被鬆弛的邊（三角不等式成立）
        for u, v, w in edges:
            if dist[u] != INF:
                assert dist[u] + w >= dist[v]

        # 每個可達點的距離都必須對應一條真實存在的、權重相等的路徑
        for t in range(n):
            if dist[t] == INF:
                assert build_path(pre, s, t) == []
                continue
            p = build_path(pre, s, t)
            assert p and p[0] == s and p[-1] == t
            ok, w = path_weight(edges, p)
            assert ok and w == dist[t]

        # 與 Floyd-Warshall 交叉驗證（只在整張圖都沒有負環時才可比）
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
