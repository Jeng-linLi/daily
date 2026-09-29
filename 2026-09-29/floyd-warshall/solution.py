"""Floyd-Warshall 全源最短路（含負環檢測、路徑還原、傳遞閉包）

題意：
    給定一張 n 個點（編號 0..n-1）、m 條帶權有向邊的圖（邊權可以爲負），求：
      1. 圖中是否存在負權環（從任一點出發沿有向邊能繞到的總權和爲負的環）；
      2. 任意兩點間的最短路距離矩陣（不可達記作 INF，自己到自己爲 0）；
      3. 所有有限距離中的最大值，即圖的「直徑」（沒有任何有限距離時輸出 INF）。

思路：
    Floyd-Warshall 本質是一個**區間 / 集合上的動態規劃**，三個字母的狀態定義是：
        dist[k][i][j] = 從 i 到 j、中間點只允許取自 {0,1,...,k-1} 的最短路長度
    轉移只有兩種選擇：
        - 不經過點 k：dist[k][i][j] = dist[k-1][i][j]
        - 經過點 k：dist[k-1][i][k] + dist[k-1][k][j]
    取 min 即可。而「第 k 層只依賴第 k-1 層」這個性質允許把第一維**滾動掉**，
    直接在原矩陣上原地更新（這正是它只要 O(n^2) 額外空間的原因），
    代價是三重循環的順序必須寫成 `for k: for i: for j:`，順序寫錯結果就錯。

    負環檢測：算法跑完後若存在 dist[i][i] < 0，說明從 i 出發繞一圈回到自己
    還能更短 —— 這就是負環。注意 Floyd 檢測的是**全圖任意位置**的負環，
    不像 Bellman-Ford 只檢測從單一源點可達的那些。

    路徑還原：維護 next[i][j] = 從 i 走向 j 的第一步落在哪個點。
    dist[i][j] 被 i→k→j 更新時，next[i][j] 就跟著改成 next[i][k]。
    還原時從 i 沿 next 一路跳到 j 即可；跳不動（next 爲 -1）說明不可達。

    複雜度 O(n^3) 時間、O(n^2) 空間。n 在幾百以內很合適（n = 1000 就上億次了），
    優勢是寫起來極短，而且一次算出全部點對 —— 這是跑 n 次 Dijkstra 拿不到的簡潔性。

輸入格式（stdin，數字按空白分隔即可）：
    n m
    u1 v1 w1
    ...
    um vm wm          （m = 0 時沒有這些行）
輸出格式（stdout）：
    第 1 行：1 表示存在負環，0 表示不存在
    第 2 .. n+1 行：距離矩陣，每行 n 個數，空格分隔，不可達輸出 INF
    第 n+2 行：圖的直徑（所有有限距離的最大值；沒有則輸出 INF）
無 stdin 輸入時運行內置斷言測試並輸出 `all tests passed`。
"""

import random
import sys
from typing import Iterator, List, Optional, Tuple

INF = float("inf")          # Python 版用 inf 表示不可達，輸出時統一打印成 INF
INF_STR = "INF"


# ---------------- 核心算法 ----------------

def floyd_warshall(n: int, edges: List[Tuple[int, int, int]]):
    """返回 (dist, nxt, has_negative_cycle)。時間 O(n^3)，空間 O(n^2)。"""
    dist = [[INF] * n for _ in range(n)]
    nxt = [[-1] * n for _ in range(n)]          # next[i][j]：從 i 到 j 的第一步

    for i in range(n):
        dist[i][i] = 0                          # 自己到自己距離 0（負自環會把它改小）
        nxt[i][i] = i                           # 保證「dist 有限 ⟹ next 有效」這條不變量
    for u, v, w in edges:
        if not (0 <= u < n and 0 <= v < n):     # 越界邊直接忽略，和 C++ 版保持一致
            continue
        if w < dist[u][v]:
            dist[u][v] = w
            nxt[u][v] = v

    # 三重循環順序必須是 k → i → j：k 是「允許使用的中間點」，必須放在最外層
    for k in range(n):
        for i in range(n):
            dik = dist[i][k]
            if dik == INF:                      # 小剪枝：i 根本到不了 k
                continue
            row_i, row_k = dist[i], dist[k]
            nxt_i = nxt[i]
            for j in range(n):
                if row_k[j] == INF:
                    continue
                nd = dik + row_k[j]
                if nd < row_i[j]:
                    row_i[j] = nd
                    # 走 i → k → j，所以從 i 出發的第一步就是「i 走向 k」的第一步
                    nxt_i[j] = nxt_i[k]

    has_neg = any(dist[i][i] < 0 for i in range(n))
    return dist, nxt, has_neg


def restore_path(nxt: List[List[int]], u: int, v: int) -> List[int]:
    """用 next 矩陣還原 u → v 的一條最短路；不可達返回空列表。

    注意：圖中存在負環時，路徑可能無意義，不要在負環情形下依賴它。
    """
    n = len(nxt)
    if not (0 <= u < n and 0 <= v < n):
        return []
    if nxt[u][v] == -1:                         # 不可達
        return []
    path = [u]
    cur = u
    while cur != v:
        cur = nxt[cur][v]
        if cur == -1:
            return []
        path.append(cur)
        if len(path) > n + 1:                   # 負環時 next 指針可能繞圈，防禦一下
            return []
    return path


def diameter(dist: List[List[float]]) -> float:
    """所有有限距離中的最大值；沒有則返回 INF。"""
    best: float = INF
    for row in dist:
        for d in row:
            if d != INF and (best == INF or d > best):
                best = d
    return best


def warshall_closure(n: int, edges: List[Tuple[int, int, int]]) -> List[List[bool]]:
    """Warshall 傳遞閉包：reach[i][j] 表示 i 能否到達 j。複雜度同爲 O(n^3)。

    它就是 Floyd 把 min/+ 換成 or/and 的版本（同一個三重循環骨架）。
    """
    reach = [[i == j for j in range(n)] for i in range(n)]
    for u, v, _w in edges:
        if 0 <= u < n and 0 <= v < n:
            reach[u][v] = True
    for k in range(n):
        for i in range(n):
            if not reach[i][k]:
                continue
            for j in range(n):
                if reach[k][j]:
                    reach[i][j] = True
    return reach


# ---------------- 對照實現（Bellman-Ford，用於對拍） ----------------

def bellman_ford_from(n: int, edges: List[Tuple[int, int, int]], s: int) -> Tuple[List[float], bool]:
    """從單點 s 出發的 Bellman-Ford，返回 (dist, 是否檢測到從 s 可達的負環)。"""
    dist = [INF] * n
    dist[s] = 0
    for it in range(n):
        changed = False
        for u, v, w in edges:
            if not (0 <= u < n and 0 <= v < n):
                continue
            if dist[u] != INF and dist[u] + w < dist[v]:
                dist[v] = dist[u] + w
                changed = True
        if not changed:
            break
        if it == n - 1:
            return dist, True                   # 第 n 輪還能鬆弛 → 負環
    return dist, False


def brute_all_pairs(n: int, edges: List[Tuple[int, int, int]]):
    """跑 n 次 Bellman-Ford 得到全源最短路，用於對拍 Floyd 的結果。"""
    all_dist, has_neg = [], False
    for s in range(n):
        d, neg = bellman_ford_from(n, edges, s)
        all_dist.append(d)
        has_neg = has_neg or neg
    return all_dist, has_neg


# ---------------- IO ----------------

def _next_int(it: Iterator[str], default: int = 0) -> int:
    """取下一個整數；輸入被截斷時用默認值兜底，避免直接拋異常。"""
    try:
        return int(next(it))
    except (StopIteration, ValueError):
        return default


def _fmt(x) -> str:
    """距離的統一打印格式：不可達打印 INF。"""
    return INF_STR if x == INF else str(x)


def run_io(data: str) -> None:
    """按統一輸入輸出格式處理 stdin 數據。"""
    it = iter(data.split())
    n = _next_int(it)
    m = _next_int(it)
    edges: List[Tuple[int, int, int]] = []
    for _ in range(max(m, 0)):
        u = _next_int(it)
        v = _next_int(it)
        w = _next_int(it)
        edges.append((u, v, w))

    dist, _nxt, has_neg = floyd_warshall(n, edges)
    print(1 if has_neg else 0)
    for i in range(n):
        print(" ".join(_fmt(dist[i][j]) for j in range(n)))
    print(_fmt(diameter(dist)))


def run_tests() -> None:
    # README 示例：4 個點、5 條邊（有負權邊但沒有負環）
    n, edges = 4, [(0, 1, 3), (0, 2, 7), (1, 2, -2), (1, 3, 5), (2, 3, 1)]
    dist, nxt, has_neg = floyd_warshall(n, edges)
    assert has_neg is False
    expected = [
        [0, 3, 1, 2],
        [INF, 0, -2, -1],
        [INF, INF, 0, 1],
        [INF, INF, INF, 0],
    ]
    for i in range(n):
        for j in range(n):
            assert dist[i][j] == expected[i][j], (i, j, dist[i][j], expected[i][j])
    assert diameter(dist) == 3

    # 路徑還原：0 → 3 的最短路是 0 → 1 → 2 → 3
    p = restore_path(nxt, 0, 3)
    assert p == [0, 1, 2, 3]
    # 還原出來的路徑長度必須等於 dist
    total = 0
    wmap = {(u, v): w for u, v, w in edges}
    for a, b in zip(p, p[1:]):
        total += wmap[(a, b)]
    assert total == dist[0][3]

    # 傳遞閉包與距離矩陣必須一致：有限距離 ⟺ 可達
    reach = warshall_closure(n, edges)
    for i in range(n):
        for j in range(n):
            assert reach[i][j] == (dist[i][j] != INF)

    # 負環：1 → 2 → 3 → 1 權和 -1
    dist2, _nxt2, neg2 = floyd_warshall(4, [(0, 1, 1), (1, 2, 1), (2, 3, -3), (3, 1, 1)])
    assert neg2 is True
    assert any(dist2[i][i] < 0 for i in range(4))

    # 負自環也算負環
    _, _, neg_self = floyd_warshall(2, [(0, 0, -1)])
    assert neg_self is True

    # 零環 / 正環不算負環
    _, _, neg_zero = floyd_warshall(3, [(0, 1, 1), (1, 2, 1), (2, 0, -2)])
    assert neg_zero is False
    _, _, neg_pos = floyd_warshall(3, [(0, 1, 1), (1, 2, 1), (2, 0, 1)])
    assert neg_pos is False

    # 空圖 / 單點 / 沒有邊
    d0, _, neg0 = floyd_warshall(0, [])
    assert d0 == [] and neg0 is False and diameter(d0) == INF
    d1, nx1, neg1 = floyd_warshall(1, [])
    assert d1 == [[0]] and neg1 is False and diameter(d1) == 0
    assert restore_path(nx1, 0, 0) == [0]
    d2, _, neg2b = floyd_warshall(3, [])
    assert neg2b is False
    assert d2 == [[0, INF, INF], [INF, 0, INF], [INF, INF, 0]]
    assert diameter(d2) == 0                      # 只有自己到自己這些 0

    # 不可達：單向邊
    d3, nx3, _ = floyd_warshall(3, [(0, 1, 5)])
    assert d3[0][1] == 5 and d3[1][0] == INF
    assert restore_path(nx3, 1, 0) == []
    assert diameter(d3) == 5

    # 重複邊取最小的那條
    d4, _, _ = floyd_warshall(2, [(0, 1, 9), (0, 1, 4)])
    assert d4[0][1] == 4

    random.seed(20260929)

    # 隨機對拍：與跑 n 次 Bellman-Ford 的結果逐格比對
    for _ in range(400):
        n = random.randint(1, 6)
        m = random.randint(0, 10)
        edges = [(random.randrange(n), random.randrange(n), random.randint(-6, 12))
                 for _ in range(m)]
        dist, nxt, has_neg = floyd_warshall(n, edges)
        ref, ref_neg = brute_all_pairs(n, edges)
        assert has_neg == ref_neg
        if not has_neg:
            for i in range(n):
                for j in range(n):
                    assert dist[i][j] == ref[i][j], (n, edges, i, j, dist[i][j], ref[i][j])
            # 無負環時，三角不等式必須成立：dist[i][j] <= dist[i][k] + dist[k][j]
            for i in range(n):
                for j in range(n):
                    for k in range(n):
                        if dist[i][k] != INF and dist[k][j] != INF:
                            assert dist[i][j] <= dist[i][k] + dist[k][j]
            # 傳遞閉包一致性
            reach = warshall_closure(n, edges)
            for i in range(n):
                for j in range(n):
                    assert reach[i][j] == (dist[i][j] != INF)
            # 還原出來的路徑確實是一條合法路徑，且總長等於 dist
            for i in range(n):
                for j in range(n):
                    path = restore_path(nxt, i, j)
                    if dist[i][j] == INF:
                        assert path == []
                    elif i == j:
                        assert path == [i]
                    else:
                        assert path[0] == i and path[-1] == j
                        cost = 0
                        ok = True
                        for a, b in zip(path, path[1:]):
                            if (a, b) not in wmap_of(edges):
                                ok = False
                                break
                            cost += min(w for (_u, _v, w) in edges if (_u, _v) == (a, b))
                        assert ok and cost == dist[i][j], (path, cost, dist[i][j])

    print("all tests passed")


def wmap_of(edges):
    """測試輔助：把邊集轉成 (u, v) → 最小權 的映射（重複邊取最小）。"""
    d = {}
    for u, v, w in edges:
        if (u, v) not in d or w < d[(u, v)]:
            d[(u, v)] = w
    return d


if __name__ == "__main__":
    raw = sys.stdin.read() if not sys.stdin.isatty() else ""
    if raw.strip():
        run_io(raw)
    else:
        run_tests()
