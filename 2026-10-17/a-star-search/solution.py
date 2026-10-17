"""A* 與 IDA*（啟發式搜索）

題意：
    三小題，一次輸入全部做完，每一題都同時給出「啟發式搜索」與「無資訊的最優對照組」並對拍。

    1. 網格最短路：R×C 格子，'#' 是障礙，四方向移動、每步代價 1。
       A*（曼哈頓距離啟發）vs BFS，比較「最短路長度」與「展開節點數」。
    2. 帶權圖最短路：n 點 m 條有向邊（邊權非負），外加每個點的啟發值 h[v]
       （h 為到終點的估計距離）。A* vs Dijkstra，同樣比較距離與展開節點數。
    3. 八碼難題（8-puzzle）：3×3 棋盤，給定初始局面，用 IDA* 求最少步數；
       與「從目標局面 BFS 得到的精確距離」對照（BFS 只做到深度 12，查不到就標記未知）。

思路：
    ### A*：把「已經花的代價 g」和「預估還要花的代價 h」加起來排序
    優先佇列依照 f = g + h 由小到大取出節點。
      - h ≡ 0 時，A* 就退化成 Dijkstra（本題的對照組）。
      - h 若**可採納**（h(v) ≤ 真實距離 d(v, t)）且**一致**
        （h(u) ≤ w(u,v) + h(v)），A* 保證第一次取出終點就是最短路。
        網格四方向 + 曼哈頓距離正好滿足這兩個條件（障礙只會讓真實距離變長，
        所以 h 仍是可採納的下界）。
      - A* 的價值在於**展開的節點更少**：h 越準，搜索越像「直接朝終點走」。
    為了讓 Python 與 C++ 的展開數完全一致，佇列的排序鍵固定為 (f, h, 節點編號)，
    鄰居的枚舉順序也固定為 上 / 下 / 左 / 右。

    ### IDA*：迭代加深 + 記憶「這輪超標多少」當下一輪門檻
    A* 要存 open/closed 表，八碼難題的狀態太多會爆炸。IDA* 用 DFS 換空間：
      - 設門檻 bound，DFS 只展開 f = g + h ≤ bound 的節點；
      - 這輪找不到就記錄「最小的超標值」當下一輪的 bound，從頭再搜一次；
      - 因為每輪門檻只增不減，且 h 可採納，**第一次找到的解一定是最優解**。
    空間複雜度只有 O(最長路徑)，代價是重複搜索（但每輪的分支都被 h 砍掉大半）。
    本題的 h 用「曼哈頓距離 + 線性衝突」：
      曼哈頓距離給出每塊磁磚至少要走幾步；
      同一行（列）裡有兩塊磁磚的目標位置也在這一行（列）但左右（上下）順序相反時，
      它們必定要「繞出去」一次，至少多 2 步 —— 這個額外資訊能砍掉大量分支。

    ### 八碼難題的可解性
    3×3 棋盤寬度為奇數，所以「逆序數（不含空格）為偶數」↔ 可解。
    不可解時直接回報 -1，不必展開任何節點。

應用場景：
    - A*：遊戲 NPC 尋路（網格 / 導航網）、地圖導航（真實路網 + 直線距離或
      landmark 啟發）、機械手臂與機器人路徑規劃、拼圖 / 魔方的狀態空間搜索
    - 權衡：h 越強（越接近真實距離）展開越少，但計算 h 本身也要時間；
      h = 0 是 Dijkstra（最穩但最慢），h 過強且不可採納則會得到次優解
    - IDA*：狀態空間大、記憶體吃緊時的首選（15-puzzle、旅行推銷員的分支限界）、
      也是「迭代加深 + 啟發剪枝」這個通用框架的典型代表

複雜度：
    網格 A* / BFS       最壞 O(R·C) 時間 / O(R·C) 空間
    帶權圖 A* / Dijkstra O((n + m) log n) 時間 / O(n) 空間
    IDA*                時間 O(b^d)（b 分支數、d 解答深度），空間 O(d)
"""

import heapq
import random
import re
import sys
from collections import deque
from typing import Dict, List, Optional, Tuple

INT_RE = re.compile(r"^[+-]?[0-9]+$")
BOARD_RE = re.compile(r"^[0-9]{9}$")

MAXR = 20            # IO 模式下網格最大邊長
MAXN = 40            # 帶權圖最大點數
MAXM = 200           # 最大邊數
WMAX = 1000          # 邊權上限
HMAX = 10 ** 6       # 啟發值上限
IDA_LIMIT = 500000   # IDA* 節點預算（超標就放棄，回報 -2）
ORACLE_DEPTH = 12    # 對照用 BFS 只做到這個深度

GOAL = (1, 2, 3, 4, 5, 6, 7, 8, 0)
DIRS4 = ((-1, 0), (1, 0), (0, -1), (0, 1))     # 上 / 下 / 左 / 右，順序固定

# ---------------------------------------------------------------- 八碼難題：前置表

# blank 在某格時可以移動到哪些格（順序固定為 上/下/左/右）
MOVE_TABLE: List[List[int]] = []
for p in range(9):
    r, c = divmod(p, 3)
    lst = []
    for dr, dc in DIRS4:
        nr, nc = r + dr, c + dc
        if 0 <= nr < 3 and 0 <= nc < 3:
            lst.append(nr * 3 + nc)
    MOVE_TABLE.append(lst)


def manhattan(state: Tuple[int, ...]) -> int:
    """每塊磁磚到目標位置的曼哈頓距離總和（不含空格）。"""
    total = 0
    for i, tile in enumerate(state):
        if tile == 0:
            continue
        gr, gc = divmod(tile - 1, 3)
        r, c = divmod(i, 3)
        total += abs(gr - r) + abs(gc - c)
    return total


def linear_conflict(state: Tuple[int, ...]) -> int:
    """同一行/列中「目標也在同一行/列、但相對順序相反」的磁磚對數 × 2。"""
    extra = 0
    for line in range(3):
        rows = [state[line * 3 + c] for c in range(3)]
        rows = [t for t in rows if t != 0 and (t - 1) // 3 == line]
        for i in range(len(rows)):
            for j in range(i + 1, len(rows)):
                if (rows[i] - 1) % 3 > (rows[j] - 1) % 3:
                    extra += 2
        cols = [state[r * 3 + line] for r in range(3)]
        cols = [t for t in cols if t != 0 and (t - 1) % 3 == line]
        for i in range(len(cols)):
            for j in range(i + 1, len(cols)):
                if (cols[i] - 1) // 3 > (cols[j] - 1) // 3:
                    extra += 2
    return extra


def heur(state: Tuple[int, ...]) -> int:
    return manhattan(state) + linear_conflict(state)


def solvable(state: Tuple[int, ...]) -> bool:
    """3×3 寬度為奇數：逆序數（不含空格）為偶數 ⟺ 可解。"""
    tiles = [t for t in state if t != 0]
    inv = 0
    for i in range(len(tiles)):
        for j in range(i + 1, len(tiles)):
            if tiles[i] > tiles[j]:
                inv += 1
    return inv % 2 == 0


# ---------------------------------------------------------------- 八碼難題：IDA*

INF = 10 ** 9


class Budget(Exception):
    """節點預算用盡。"""


def ida_star(start: Tuple[int, ...], limit: int = IDA_LIMIT) -> Tuple[int, int]:
    """回傳 (最少步數, 展開節點數)；無解回傳 (-1, nodes)，超預算回傳 (-2, nodes)。"""
    if not solvable(start):
        return -1, 0
    if start == GOAL:
        return 0, 1

    counter = [0]

    def dfs(state: Tuple[int, ...], g: int, bound: int, prev_blank: int) -> Tuple[int, int]:
        """回傳 (結果, 下一個門檻)；結果 ≥ 0 表示找到解，-1 表示這輪沒找到。"""
        counter[0] += 1
        if counter[0] > limit:
            raise Budget()
        f = g + heur(state)
        if f > bound:
            return -1, f
        if state == GOAL:
            return g, INF
        blank = state.index(0)
        nxt_bound = INF
        for nb in MOVE_TABLE[blank]:
            if nb == prev_blank:
                continue
            lst = list(state)
            lst[blank], lst[nb] = lst[nb], lst[blank]
            child = tuple(lst)
            res, nb2 = dfs(child, g + 1, bound, blank)
            if res >= 0:
                return res, INF
            if nb2 < nxt_bound:
                nxt_bound = nb2
        return -1, nxt_bound

    bound = heur(start)
    try:
        while bound <= 40:
            res, nxt = dfs(start, 0, bound, -1)
            if res >= 0:
                return res, counter[0]
            if nxt == INF:
                return -1, counter[0]
            bound = nxt
    except Budget:
        return -2, counter[0]
    return -2, counter[0]


# ---------------------------------------------------------------- 八碼難題：BFS 對照組

_ORACLE: Optional[Dict[Tuple[int, ...], int]] = None


def oracle_dist(state: Tuple[int, ...]) -> int:
    """從目標局面 BFS 出來的精確距離；深度超過 ORACLE_DEPTH 就回傳 -1（未知）。"""
    global _ORACLE
    if _ORACLE is None:
        dist = {GOAL: 0}
        frontier = [GOAL]
        d = 0
        while d < ORACLE_DEPTH:
            d += 1
            nxt_frontier = []
            for st in frontier:
                blank = st.index(0)
                for nb in MOVE_TABLE[blank]:
                    lst = list(st)
                    lst[blank], lst[nb] = lst[nb], lst[blank]
                    child = tuple(lst)
                    if child not in dist:
                        dist[child] = d
                        nxt_frontier.append(child)
            frontier = nxt_frontier
            if not frontier:
                break
        _ORACLE = dist
    return _ORACLE.get(state, -1)


# ---------------------------------------------------------------- 網格：A* 與 BFS


def astar_grid(grid: List[str], sr: int, sc: int, tr: int, tc: int) -> Tuple[int, int]:
    """四方向單位權網格最短路。回傳 (距離, 展開節點數)，不可達回傳 (-1, ...)。"""
    R, C = len(grid), len(grid[0]) if grid else 0
    if R == 0 or C == 0:
        return -1, 0

    def h(r: int, c: int) -> int:
        return abs(r - tr) + abs(c - tc)

    INF_D = 1 << 60
    dist = [[INF_D] * C for _ in range(R)]
    closed = [[False] * C for _ in range(R)]
    dist[sr][sc] = 0
    heap = [(h(sr, sc), h(sr, sc), sr * C + sc)]
    expanded = 0
    while heap:
        _f, _hv, code = heapq.heappop(heap)
        r, c = divmod(code, C)
        if closed[r][c]:
            continue
        closed[r][c] = True
        expanded += 1
        if r == tr and c == tc:
            break
        for dr, dc in DIRS4:
            nr, nc = r + dr, c + dc
            if 0 <= nr < R and 0 <= nc < C and grid[nr][nc] != "#":
                if not closed[nr][nc] and dist[r][c] + 1 < dist[nr][nc]:
                    dist[nr][nc] = dist[r][c] + 1
                    hv = h(nr, nc)
                    heapq.heappush(heap, (dist[nr][nc] + hv, hv, nr * C + nc))
    d = dist[tr][tc]
    return (d if d < INF_D else -1), expanded


def bfs_grid(grid: List[str], sr: int, sc: int, tr: int, tc: int) -> Tuple[int, int]:
    R, C = len(grid), len(grid[0]) if grid else 0
    if R == 0 or C == 0:
        return -1, 0
    INF_D = 1 << 60
    dist = [[INF_D] * C for _ in range(R)]
    dist[sr][sc] = 0
    q = deque([(sr, sc)])
    expanded = 0
    while q:
        r, c = q.popleft()
        if dist[r][c] == INF_D:
            continue
        expanded += 1
        if r == tr and c == tc:
            break
        for dr, dc in DIRS4:
            nr, nc = r + dr, c + dc
            if 0 <= nr < R and 0 <= nc < C and grid[nr][nc] != "#" and dist[nr][nc] == INF_D:
                dist[nr][nc] = dist[r][c] + 1
                q.append((nr, nc))
    d = dist[tr][tc]
    return (d if d < INF_D else -1), expanded


# ---------------------------------------------------------------- 帶權圖：A* 與 Dijkstra


def dijkstra(n: int, adj: List[List[Tuple[int, int]]], s: int, t: int) -> Tuple[int, int]:
    INF_D = 1 << 60
    dist = [INF_D] * n
    closed = [False] * n
    dist[s] = 0
    heap = [(0, s)]
    expanded = 0
    while heap:
        d, v = heapq.heappop(heap)
        if closed[v]:
            continue
        closed[v] = True
        expanded += 1
        if v == t:
            break
        for u, w in adj[v]:
            if not closed[u] and d + w < dist[u]:
                dist[u] = d + w
                heapq.heappush(heap, (dist[u], u))
    return (dist[t] if dist[t] < INF_D else -1), expanded


def astar_weighted(n: int, adj: List[List[Tuple[int, int]]], h: List[int],
                   s: int, t: int) -> Tuple[int, int]:
    INF_D = 1 << 60
    dist = [INF_D] * n
    closed = [False] * n
    dist[s] = 0
    heap = [(h[s], h[s], s)]
    expanded = 0
    while heap:
        _f, _hv, v = heapq.heappop(heap)
        if closed[v]:
            continue
        closed[v] = True
        expanded += 1
        if v == t:
            break
        for u, w in adj[v]:
            if not closed[u] and dist[v] + w < dist[u]:
                dist[u] = dist[v] + w
                hv = h[u]
                heapq.heappush(heap, (dist[u] + hv, hv, u))
    return (dist[t] if dist[t] < INF_D else -1), expanded


# ---------------------------------------------------------------- 輸出


def run(grid: List[str], sr: int, sc: int, tr: int, tc: int,
        n: int, edges: List[Tuple[int, int, int]], gs: int, gt: int,
        h: List[int], board: Tuple[int, ...]) -> str:
    adj: List[List[Tuple[int, int]]] = [[] for _ in range(n)]
    for u, v, w in edges:
        adj[u].append((v, w))

    a_dist, a_nodes = astar_grid(grid, sr, sc, tr, tc)
    b_dist, b_nodes = bfs_grid(grid, sr, sc, tr, tc)
    if n == 0:
        w_dist, w_nodes, d_dist, d_nodes = -1, 0, -1, 0
    else:
        w_dist, w_nodes = astar_weighted(n, adj, h, gs, gt)
        d_dist, d_nodes = dijkstra(n, adj, gs, gt)
    i_dist, i_nodes = ida_star(board)
    o_dist = oracle_dist(board)

    ok = 1
    if a_dist != b_dist:
        ok = 0
    if w_dist != d_dist:
        ok = 0                       # h 可採納時兩者必須相等
    if o_dist >= 0 and i_dist != o_dist:
        ok = 0

    return "\n".join([
        "%d %d %d %d" % (a_dist, a_nodes, b_dist, b_nodes),
        "%d %d %d %d" % (w_dist, w_nodes, d_dist, d_nodes),
        "%d %d %d" % (i_dist, i_nodes, o_dist),
        str(ok),
    ])


# ---------------------------------------------------------------- 內置測試


def gen_maze(rng: random.Random, R: int, C: int) -> List[str]:
    return ["".join("#" if rng.random() < 0.28 else "." for _ in range(C)) for _ in range(R)]


def gen_graph(rng: random.Random) -> Tuple[int, List[Tuple[int, int, int]], int, int, List[int]]:
    """產生帶權圖 + 一組可採納且一致的啟發值。"""
    n = rng.randint(1, 12)
    m = rng.randint(0, min(30, n * n))
    edges = []
    for _ in range(m):
        u = rng.randrange(0, n)
        v = rng.randrange(0, n)
        w = rng.randint(0, 9)
        edges.append((u, v, w))
    s = rng.randrange(0, n)
    t = rng.randrange(0, n)
    # 先用 Dijkstra 從 t 反推真實距離（反向圖）
    rev: List[List[Tuple[int, int]]] = [[] for _ in range(n)]
    for u, v, w in edges:
        rev[v].append((u, w))
    INF_D = 1 << 60
    dist = [INF_D] * n
    dist[t] = 0
    heap = [(0, t)]
    while heap:
        d, v = heapq.heappop(heap)
        if d > dist[v]:
            continue
        for u, w in rev[v]:
            if d + w < dist[u]:
                dist[u] = d + w
                heapq.heappush(heap, (dist[u], u))
    num, den = rng.choice([(0, 1), (1, 2), (3, 4), (1, 1)])   # 整數運算，避免浮點差異
    h = []
    for i in range(n):
        true_d = dist[i]
        if true_d >= INF_D:
            h.append(0)              # 不可達：給 0（對不可達點而言 h 不影響最優性）
        else:
            h.append(true_d * num // den)
    return n, edges, s, t, h


def run_tests() -> None:
    # ---- 網格固定用例 ----
    g0 = ["...", "...", "..."]
    assert astar_grid(g0, 0, 0, 2, 2)[0] == 4
    a, an = astar_grid(g0, 0, 0, 2, 2)
    b, bn = bfs_grid(g0, 0, 0, 2, 2)
    assert a == b
    assert an <= bn                       # 曼哈頓啟發不會展開得比 BFS 多

    g1 = ["....", ".##.", "....", "...."]
    a, an = astar_grid(g1, 0, 0, 3, 3)
    b, bn = bfs_grid(g1, 0, 0, 3, 3)
    assert a == b and a > 0
    assert an <= bn

    # 目標被完全圍住 → 不可達
    g2 = [".....", ".###.", ".#.#.", ".###.", "....."]
    assert astar_grid(g2, 0, 0, 2, 2)[0] == -1
    assert bfs_grid(g2, 0, 0, 2, 2)[0] == -1

    # 退化輸入
    assert astar_grid([""], 0, 0, 0, 0)[0] == -1
    assert astar_grid([".#"], 0, 0, 0, 1)[0] == -1
    assert astar_grid([], 0, 0, 0, 0)[0] == -1

    # ---- 帶權圖固定用例 ----
    n = 5
    edges = [(0, 1, 4), (0, 2, 1), (2, 1, 1), (1, 3, 2), (2, 3, 5), (3, 4, 3)]
    adj: List[List[Tuple[int, int]]] = [[] for _ in range(n)]
    for u, v, w in edges:
        adj[u].append((v, w))
    h0 = [0] * n
    assert dijkstra(n, adj, 0, 4) == astar_weighted(n, adj, h0, 0, 4)
    assert dijkstra(n, adj, 0, 4)[0] == 7
    # 真實距離當啟發 → A* 幾乎不展開多餘節點
    htrue = [7, 5, 6, 3, 0]
    assert astar_weighted(n, adj, htrue, 0, 4)[0] == 7
    assert astar_weighted(n, adj, htrue, 0, 4)[1] <= dijkstra(n, adj, 0, 4)[1]
    assert dijkstra(n, adj, 4, 0)[0] == -1

    # ---- 八碼難題固定用例 ----
    assert manhattan(GOAL) == 0 and heur(GOAL) == 0
    assert ida_star(GOAL) == (0, 1)
    one = (1, 2, 3, 4, 5, 0, 7, 8, 6)      # 空格上移一步
    assert oracle_dist(one) == 1
    assert ida_star(one)[0] == 1
    assert solvable(one)
    bad = (2, 1, 3, 4, 5, 6, 7, 8, 0)      # 逆序數 1 → 不可解
    assert not solvable(bad)
    assert ida_star(bad) == (-1, 0)
    assert oracle_dist(bad) == -1

    # 線性衝突：第一行 3 1 2 的目標順序是 1 2 3，有兩對順序反了 → +4
    conflict = (3, 1, 2, 4, 5, 6, 7, 8, 0)
    assert manhattan(conflict) == 4
    assert linear_conflict(conflict) == 4
    assert heur(conflict) == manhattan(conflict) + linear_conflict(conflict)

    # ---- 隨機對拍 ----
    rng = random.Random(20261017)
    for _ in range(300):
        R = rng.randint(1, 8)
        C = rng.randint(1, 8)
        grid = gen_maze(rng, R, C)
        sr, sc = rng.randrange(0, R), rng.randrange(0, C)
        tr, tc = rng.randrange(0, R), rng.randrange(0, C)
        n, edges, gs, gt, h = gen_graph(rng)
        board = list(GOAL)
        pos = 8
        for _ in range(rng.randint(0, 9)):
            cand = MOVE_TABLE[pos]
            np_ = cand[rng.randrange(0, len(cand))]
            board[pos], board[np_] = board[np_], board[pos]
            pos = np_
        out = run(grid, sr, sc, tr, tc, n, edges, gs, gt, h, tuple(board))
        assert out.split("\n")[-1] == "1", (grid, sr, sc, tr, tc, n, edges, gs, gt, h, board, out)

    # ---- 深一點的八碼局面（只驗證 IDA* 與 BFS 對照組一致）----
    for depth in (10, 12):
        rng2 = random.Random(depth)
        for _ in range(3):
            board = list(GOAL)
            pos = 8
            for _ in range(depth):
                cand = MOVE_TABLE[pos]
                np_ = cand[rng2.randrange(0, len(cand))]
                board[pos], board[np_] = board[np_], board[pos]
                pos = np_
            st = tuple(board)
            assert ida_star(st)[0] == oracle_dist(st), (st, ida_star(st), oracle_dist(st))


# ---------------------------------------------------------------- IO


def emit(text: str) -> None:
    sys.stdout.buffer.write(text.encode("latin-1"))
    sys.stdout.buffer.flush()


def run_io(raw: str) -> None:
    toks = raw.split()
    idx = 0

    def next_int(default: int) -> int:
        nonlocal idx
        if idx >= len(toks):
            return default
        t = toks[idx]
        idx += 1
        if not INT_RE.match(t):
            return default
        v = int(t)
        return max(-(10 ** 15), min(10 ** 15, v))

    def next_str() -> str:
        nonlocal idx
        if idx >= len(toks):
            return ""
        t = toks[idx]
        idx += 1
        return t

    R = max(0, min(next_int(0), MAXR))
    C = max(0, min(next_int(0), MAXR))
    grid = []
    for r in range(R):
        row = next_str()
        row = "".join(ch if ch in ".#" else "#" for ch in row[:C])
        grid.append(row + "." * (C - len(row)))
    if R == 0:
        C = 0

    def cell(default: int, hi: int) -> int:
        v = next_int(default)
        return max(0, min(v, max(0, hi - 1)))

    sr = cell(0, R)
    sc = cell(0, C)
    tr = cell(0, R)
    tc = cell(0, C)

    n = max(0, min(next_int(0), MAXN))
    m = max(0, min(next_int(0), MAXM))
    edges = []
    for _ in range(m):
        u = next_int(0)
        v = next_int(0)
        w = next_int(0)
        w = max(0, min(w, WMAX))
        if 0 <= u < n and 0 <= v < n:
            edges.append((u, v, w))

    gs = max(0, min(next_int(0), max(0, n - 1)))
    gt = max(0, min(next_int(0), max(0, n - 1)))
    h = []
    for _ in range(n):
        hv = next_int(0)
        h.append(max(0, min(hv, HMAX)))

    bt = next_str()
    if not BOARD_RE.match(bt):
        board = GOAL
    else:
        cand = tuple(int(ch) for ch in bt)
        board = cand if sorted(cand) == list(range(9)) else GOAL

    emit(run(grid, sr, sc, tr, tc, n, edges, gs, gt, h, board) + "\n")


def main() -> None:
    sys.setrecursionlimit(1000000)
    raw = sys.stdin.buffer.read().decode("latin-1")
    if raw.strip():
        run_io(raw)
    else:
        run_tests()
        emit("all tests passed\n")


if __name__ == "__main__":
    main()
