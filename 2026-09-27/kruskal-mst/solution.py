"""最小生成樹（Kruskal 算法 + 併查集）

題意：
    給定一個 n 個點、m 條邊的**無向帶權圖**（頂點編號 0..n-1，邊權爲整數，可負），
    求一棵最小生成樹：選出若干條邊把全部點連通，且邊權之和最小。
    若圖本身不連通（不存在生成樹），輸出 disconnected。

思路：
    Kruskal 是**貪心 + 併查集**：把所有邊按權值從小到大排序，依次嘗試加入；
    如果這條邊的兩個端點當前還不連通，就把它選進生成樹並合併兩端點所在集合；
    否則這條邊會和已選的邊成環，直接丟棄。選夠 n-1 條邊就結束。

    爲什麼這樣是對的（切分性質）：掃描到某條邊 (u, v, w) 且 u、v 尚未連通時，
    把「u 所在的連通塊」看作一個切分的一側，那麼所有跨這個切分的邊裏，
    (u, v, w) 是當前未處理的最小者 —— 因爲更小的邊全都已經被考慮過了，
    它們要麼沒能跨這個切分，要麼會把點合併進來（那樣 u、v 就已經連通了）。
    跨切分的最小邊必然屬於某棵最小生成樹，所以選它不會錯。

    併查集負責「u、v 是否已連通」的判定與合併，用路徑壓縮 + 按大小合併，
    均攤接近 O(1)。注意排序要**穩定**（同權邊保持輸入順序），
    這樣 Python 與 C++ 兩版輸出才能逐字節一致。

輸入格式（stdin，所有數字按空白分隔）：
    n m
    u1 v1 w1
    u2 v2 w2
    ...（共 m 行，頂點 0-based）
輸出格式（stdout）：
    連通時：第一行爲最小生成樹的總權重，隨後每行一條被選中的邊 `u v w`
            （按 Kruskal 的選中順序輸出，即權值升序、同權值按輸入順序）
    不連通時：只輸出一行 disconnected
無 stdin 輸入時運行內置斷言測試並輸出 `all tests passed`。
"""

import sys
from typing import Iterator, List, Optional, Tuple

Edge = Tuple[int, int, int]


class DSU:
    """併查集：路徑壓縮 + 按大小合併，同時維護連通塊個數。"""

    def __init__(self, n: int) -> None:
        self.parent = list(range(n))
        self.size = [1] * n
        self.components = n

    def find(self, x: int) -> int:
        while self.parent[x] != x:
            self.parent[x] = self.parent[self.parent[x]]   # 路徑壓縮（折半）
            x = self.parent[x]
        return x

    def union(self, a: int, b: int) -> bool:
        """合併 a、b 所在集合；返回是否真的合併了（原本不連通才返回 True）。"""
        ra, rb = self.find(a), self.find(b)
        if ra == rb:
            return False
        if self.size[ra] < self.size[rb]:     # 按大小合併，小樹掛到大樹下
            ra, rb = rb, ra
        self.parent[rb] = ra
        self.size[ra] += self.size[rb]
        self.components -= 1
        return True


def kruskal(n: int, edges: List[Edge]) -> Tuple[int, List[Edge], bool]:
    """返回 (總權重, 被選中的邊, 是否連通)。

    時間 O(m log m)（瓶頸在排序），空間 O(n + m)。
    """
    dsu = DSU(n)
    chosen: List[Edge] = []
    total = 0
    # 穩定排序：只按權值排序，同權值保持輸入先後順序
    for idx in sorted(range(len(edges)), key=lambda i: edges[i][2]):
        u, v, w = edges[idx]
        if dsu.union(u, v):
            chosen.append((u, v, w))
            total += w
            if len(chosen) == n - 1:
                break                          # 已經是一棵生成樹，可以提前結束
    connected = (n <= 1) or (len(chosen) == n - 1)
    return total, chosen, connected


def prim(n: int, edges: List[Edge]) -> Tuple[int, bool]:
    """對照用的 Prim（鄰接表 + O(n^2) 選最小），用於與 Kruskal 交叉驗證總權重。"""
    if n == 0:
        return 0, True
    adj: List[List[Tuple[int, int]]] = [[] for _ in range(n)]
    for u, v, w in edges:
        if u == v:
            continue                           # 自環對 MST 無意義
        adj[u].append((v, w))
        adj[v].append((u, w))

    INF = float("inf")
    dist = [INF] * n
    used = [False] * n
    dist[0] = 0
    total = 0
    picked = 0
    for _ in range(n):
        best = -1
        for i in range(n):
            if not used[i] and (best == -1 or dist[i] < dist[best]):
                best = i
        if best == -1 or dist[best] == INF:
            return total, False                # 剩下的點都夠不着，圖不連通
        used[best] = True
        total += dist[best]
        picked += 1
        for v, w in adj[best]:
            if not used[v] and w < dist[v]:
                dist[v] = w
    return total, picked == n


def mst_brute(n: int, edges: List[Edge]) -> Tuple[Optional[int], bool]:
    """對照用的指數級枚舉：枚舉所有邊子集，挑出權和最小的生成樹。僅用於極小規模測試。"""
    if n <= 1:
        return 0, True
    m = len(edges)
    best: Optional[int] = None
    for mask in range(1 << m):
        if bin(mask).count("1") != n - 1:      # 生成樹恰好 n-1 條邊
            continue
        dsu = DSU(n)
        ok = True
        total = 0
        for i in range(m):
            if (mask >> i) & 1:
                u, v, w = edges[i]
                if dsu.union(u, v):
                    total += w
                else:
                    ok = False                 # 成環，不是樹
                    break
        if ok and dsu.components == 1:
            if best is None or total < best:
                best = total
    return best, best is not None


def total_weight(edges: List[Edge]) -> int:
    return sum(w for _, _, w in edges)


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
    edges: List[Edge] = []
    for _ in range(m):
        u = _next_int(it)
        v = _next_int(it)
        w = _next_int(it)
        edges.append((u, v, w))

    total, chosen, connected = kruskal(n, edges)
    if not connected:
        print("disconnected")
        return
    print(total)
    for u, v, w in chosen:
        print(f"{u} {v} {w}")


def run_tests() -> None:
    # README 示例：4 點 5 邊，MST = (0,1,1) + (1,2,2) + (2,3,3) = 6
    edges = [(0, 1, 1), (0, 2, 4), (1, 2, 2), (1, 3, 5), (2, 3, 3)]
    total, chosen, ok = kruskal(4, edges)
    assert ok is True
    assert total == 6
    assert chosen == [(0, 1, 1), (1, 2, 2), (2, 3, 3)]
    assert prim(4, edges) == (6, True)
    assert mst_brute(4, edges) == (6, True)

    # 不連通：兩個點之間只有一條邊，第三個點孤立
    total2, chosen2, ok2 = kruskal(3, [(0, 1, 5)])
    assert ok2 is False
    assert prim(3, [(0, 1, 5)])[1] is False
    assert mst_brute(3, [(0, 1, 5)])[1] is False

    # 退化情形
    assert kruskal(0, []) == (0, [], True)            # 空圖視爲連通，權重 0
    assert kruskal(1, []) == (0, [], True)            # 單點，不需要邊
    assert kruskal(2, [(0, 1, 7)]) == (7, [(0, 1, 7)], True)
    assert kruskal(2, []) == (0, [], False)           # 兩點無邊，不連通

    # 自環與重邊：自環必被丟棄，重邊只留一條
    assert kruskal(2, [(0, 0, 1), (0, 1, 3), (0, 1, 3)]) == (3, [(0, 1, 3)], True)
    assert kruskal(3, [(0, 1, 2), (1, 2, 2), (0, 2, 2)])[0] == 4

    # 負權邊同樣成立
    assert kruskal(3, [(0, 1, -5), (1, 2, -1), (0, 2, 10)]) == (-6, [(0, 1, -5), (1, 2, -1)], True)

    # 同權邊按輸入順序穩定選中（保證兩語言輸出一致）
    same = [(2, 3, 1), (0, 1, 1), (1, 2, 1)]
    assert kruskal(4, same)[1] == [(2, 3, 1), (0, 1, 1), (1, 2, 1)]

    import random

    random.seed(20260927)

    # 隨機對拍一：小規模圖上 Kruskal、Prim、指數級枚舉三者結果一致
    for _ in range(150):
        n = random.randint(1, 6)
        m = random.randint(0, 8)
        cand = [(u, v) for u in range(n) for v in range(u + 1, n)]
        random.shuffle(cand)
        edges = [(u, v, random.randint(1, 20)) for u, v in cand[:m]]
        total, chosen, ok = kruskal(n, edges)
        pt, pok = prim(n, edges)
        bt, bok = mst_brute(n, edges)
        assert ok == pok == bok
        if ok:
            assert total == pt == bt, (n, edges, total, pt, bt)
            assert len(chosen) == n - 1
            assert total_weight(chosen) == total
            # 選出來的邊確實構成生成樹：n-1 條且把所有點連成一塊
            dsu = DSU(n)
            for u, v, _w in chosen:
                assert dsu.union(u, v) is True
            assert dsu.components == 1
            # 權值和不超過隨便一棵生成樹：與 Prim 再比一次即可
        else:
            assert total == 0 or len(chosen) < n - 1

    # 隨機對拍二：保證連通的隨機圖（先造一條鏈再補隨機邊），Kruskal 與 Prim 對拍
    for _ in range(150):
        n = random.randint(2, 9)
        edges: List[Edge] = []
        for i in range(1, n):                  # 先連成鏈，確保一定連通
            j = random.randint(0, i - 1)
            edges.append((j, i, random.randint(1, 30)))
        for _ in range(random.randint(0, 6)):  # 再補一些隨機邊
            u, v = random.sample(range(n), 2)
            edges.append((u, v, random.randint(1, 30)))
        total, chosen, ok = kruskal(n, edges)
        pt, pok = prim(n, edges)
        assert ok is True and pok is True
        assert total == pt
        assert len(chosen) == n - 1
        dsu = DSU(n)
        for u, v, _w in chosen:
            assert dsu.union(u, v) is True
        assert dsu.components == 1

    print("all tests passed")


if __name__ == "__main__":
    raw = sys.stdin.read() if not sys.stdin.isatty() else ""
    if raw.strip():
        run_io(raw)
    else:
        run_tests()
