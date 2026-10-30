"""樹上啟發式合併（DSU on Tree / Sack）：O(n log n) 做子樹顏色統計

題意：
    給定一棵有 n 個節點的根樹（根為 1），節點 u 有顏色 c[u]（1 <= c[u] <= n）。
    有 q 個詢問，每個詢問給出一對 (u, f)，要求回答 u 的**子樹**中：
        distinct —— 出現過的顏色數
        eq       —— 出現次數**恰好等於 f** 的顏色數
        shared   —— 出現過、而且在 u 的子樹**之外**也出現過的顏色數
    限制：1 <= n, q <= 100000，0 <= f <= n。

思路：
    ### 1. 為什麼需要這招
    樸素做法是對每個節點 u 單獨遍歷一次子樹，O(n^2)。
    「樹上啟發式合併」（又稱 DSU on Tree / Sack / 同步 bags）
    藉由**重輕劃分**把總複雜度降到 O(n log n)。

    ### 2. 重輕劃分（heavy-light by subtree size）
    對每個節點 u，取子樹最大的那個兒子當「重兒子」（heavy child），
    其餘都是「輕兒子」（light child）。
    關鍵性質：從根往下走，**輕邊最多只有 O(log n) 條** ——
    因為走過一條輕邊 (u,v) 時 size[v] <= size[u]/2，
    子樹大小每走一條輕邊就至少減半，所以輕邊數是 O(log n)。

    ### 3. Sack 的流程
    維護一個「全域狀態」（bag），記錄目前被加入的節點的顏色統計。
    對節點 u 的處理分成三段：
        (a) 先處理所有**輕兒子**，處理完立刻清空它們的貢獻
            （keep = false）
        (b) 再處理**重兒子**，但**保留**它的貢獻（keep = true）
        (c) 把所有輕兒子的子樹重新加回來，再加上 u 自己，
            此刻 bag 恰好就是 u 的子樹 —— 據此回答 u 的詢問
        (d) 若 keep = false，把 u 的整個子樹移出 bag
    輕兒子的貢獻被加入又移除了兩次、重兒子只加入一次，
    而每個節點只會是它祖先的 O(log n) 個輕兒子之一，
    所以總加入次數是 O(n log n)。

    不變式：**每次進入 dfs(u, ·) 時 bag 都是空的**，
    而 dfs(u, true) 返回時 bag 恰好是 u 的子樹。
    所以 (c) 段之後 bag = u 的子樹，可以直接回答 u 的詢問。

    ### 4. 維護哪些統計量
    cnt[c]   = 顏色 c 在 bag 中的出現次數
    freq[f]  = 出現次數恰好為 f 的顏色數（f >= 0）
    distinct = 出現過的顏色數
    shared   = 出現過、且在子樹之外也出現過的顏色數
    total[c] = 顏色 c 在整棵樹中的總出現次數（一次預處理）
    palette  = 整棵樹用到的不同顏色數

    加入一個顏色 c（cnt 從 a 變成 a+1）：
        - freq[a] -= 1；freq[a+1] += 1
        - distinct：a == 0 時 +1
        - shared：判斷 (0 < cnt[c] < total[c]) 在變化前後是否成立，
          成立狀態改變就相應 +1 / -1
    移出一個顏色 c 完全對稱。於是每次加入 / 移出都是 O(1)。

    顏色的「宇宙」取整棵樹用到的 palette 種顏色，
    因此 freq[0] 初始值 = palette，eq(0) = freq[0] = palette - distinct
    （也就是「在整棵樹出現過、但沒出現在這個子樹裡」的顏色數），
    與 shared 的語意保持一致。

    ### 5. 為什麼不問「最多的顏色」
    常見的「子樹內最多的顏色」需要維護 max frequency，
    但在 Sack 裡移除是「整批移除子樹」，max 會難以 O(1) 修正。
    本題改問「恰好 f 次的顏色數」，用 freq[] 陣列即可 O(1) 回答，
    避開了 max 追蹤的陷阱。

複雜度：
    預處理        O(n)
    單次詢問      O(1)（在 DFS 到 u 的當下直接算）
    總計          O(n log n + q)，空間 O(n)

應用場景：
    - 專案目錄的檔案副檔名統計（每個子目錄有多少種副檔名、
      有多少副檔名只出現一次、哪些副檔名與其他目錄共用）
    - 基因組 / 蛋白質序列的子樹特徵統計
    - 程式碼倉庫的語言分布分析（每個資料夾用了幾種語言）
    - 線上統計「某個資料夾內有多少種唯一的檔案類型」
"""

import sys
from typing import Dict, List, Optional, Tuple

NMAX = 100000
QMAX = 100000


class Solver:
    """樹上啟發式合併：一次 DFS 回答所有節點的詢問。"""

    def __init__(self, n: int, color: List[int], parent: List[int],
                 queries: Dict[int, List[Tuple[int, int]]]) -> None:
        self.n = n
        self.color = color
        self.parent = parent
        self.queries = queries  # u -> [(f, query_index), ...]

        # ---- 建立子節點清單（依編號遞增）----
        self.children: List[List[int]] = [[] for _ in range(n + 1)]
        for u in range(2, n + 1):
            self.children[parent[u]].append(u)
        self.light: List[List[int]] = [[] for _ in range(n + 1)]

        # ---- 迭代 DFS：取得先序序列與 tin ----
        # 先序序列中，u 的子樹恰好佔據連續區間 [tin[u], tin[u]+size[u])
        self.order: List[int] = []
        self.tin: List[int] = [0] * (n + 1)
        stack = [1]
        while stack:
            u = stack.pop()
            self.tin[u] = len(self.order)
            self.order.append(u)
            # 反序推入，使同層兄弟依編號遞增被拜訪
            for c in reversed(self.children[u]):
                stack.append(c)

        # ---- 子樹大小、重兒子 ----
        self.size: List[int] = [1] * (n + 1)
        for u in reversed(self.order):
            for c in self.children[u]:
                self.size[u] += self.size[c]
        self.heavy: List[int] = [0] * (n + 1)
        for u in range(2, n + 1):
            p = parent[u]
            if self.heavy[p] == 0 or self.size[u] > self.size[self.heavy[p]]:
                self.heavy[p] = u
        for u in range(1, n + 1):
            h = self.heavy[u]
            for c in self.children[u]:
                if c != h:
                    self.light[u].append(c)

        # ---- 整棵樹的顏色總數 ----
        self.total: List[int] = [0] * (n + 1)
        for u in range(1, n + 1):
            self.total[color[u]] += 1
        self.palette = sum(1 for c in range(1, n + 1) if self.total[c] > 0)

        # ---- 全域狀態 ----
        self.cnt: List[int] = [0] * (n + 1)
        # 顏色的宇宙 = palette 種，所以 freq[0] 初始為 palette
        self.freq: List[int] = [0] * (n + 2)
        self.freq[0] = self.palette
        self.distinct = 0
        self.shared = 0

        # ---- 答案 ----
        nq = 0
        for lst in queries.values():
            nq += len(lst)
        self.answers: List[Tuple[int, int, int]] = [(0, 0, 0)] * nq

    # ------------------------------------------------------------ 狀態操作

    def _add(self, c: int) -> None:
        a = self.cnt[c]
        b = a + 1
        # shared：變化前後的 (0 < cnt < total) 指標
        before = 0 < a < self.total[c]
        after = 0 < b < self.total[c]
        if after and not before:
            self.shared += 1
        elif before and not after:
            self.shared -= 1
        self.cnt[c] = b
        self.freq[a] -= 1
        self.freq[b] += 1
        if a == 0:
            self.distinct += 1

    def _remove(self, c: int) -> None:
        a = self.cnt[c]
        b = a - 1
        before = 0 < a < self.total[c]
        after = 0 < b < self.total[c]
        if after and not before:
            self.shared += 1
        elif before and not after:
            self.shared -= 1
        self.cnt[c] = b
        self.freq[a] -= 1
        self.freq[b] += 1
        if b == 0:
            self.distinct -= 1

    # ------------------------------------------------------------ 主流程

    def run(self) -> List[Tuple[int, int, int]]:
        answers = self.answers
        color = self.color
        order = self.order
        tin = self.tin
        size = self.size
        freq = self.freq

        # 迭代版 Sack：frame = (u, keep, stage)
        #   stage 0 = 進入（安排子節點的處理順序）
        #   stage 1 = 子節點都處理完，輪到自己
        stack: List[Tuple[int, bool, int]] = [(1, True, 0)]
        while stack:
            u, keep, stage = stack.pop()
            if stage == 0:
                # 堆疊後進先出：push 順序要反過來，
                # pop 出來的順序才是「輕兒子 -> 重兒子 -> 自己收尾」
                stack.append((u, keep, 1))
                h = self.heavy[u]
                if h:
                    stack.append((h, True, 0))
                for c in reversed(self.light[u]):
                    stack.append((c, False, 0))
            else:
                # (c) 把所有輕兒子的子樹加回來
                for c in self.light[u]:
                    start = tin[c]
                    for i in range(start, start + size[c]):
                        self._add(color[order[i]])
                # 再加上 u 自己 —— 此刻 bag 恰好是 u 的子樹
                self._add(color[u])
                # 回答 u 的詢問
                for f, qi in self.queries.get(u, ()):
                    answers[qi] = (self.distinct, freq[f], self.shared)
                # (d) 若不需要保留，把 u 的整個子樹移出 bag
                if not keep:
                    start = tin[u]
                    for i in range(start, start + size[u]):
                        self._remove(color[order[i]])
        return answers


def solve_brute(n: int, color: List[int], parent: List[int],
                queries: List[Tuple[int, int]]) -> List[Tuple[int, int, int]]:
    """暴力參考實作：對每個詢問單獨遍歷該節點的子樹。"""
    children: List[List[int]] = [[] for _ in range(n + 1)]
    for u in range(2, n + 1):
        children[parent[u]].append(u)
    total: List[int] = [0] * (n + 1)
    for u in range(1, n + 1):
        total[color[u]] += 1
    palette = sum(1 for c in range(1, n + 1) if total[c] > 0)

    out: List[Tuple[int, int, int]] = []
    for u, f in queries:
        cnt: Dict[int, int] = {}
        stack = [u]
        while stack:
            x = stack.pop()
            cnt[color[x]] = cnt.get(color[x], 0) + 1
            for c in children[x]:
                stack.append(c)
        distinct = len(cnt)
        if f == 0:
            eq = palette - distinct
        else:
            eq = sum(1 for v in cnt.values() if v == f)
        shared = sum(1 for c, v in cnt.items() if v < total[c])
        out.append((distinct, eq, shared))
    return out


def run_tests() -> None:
    # ---- 固定用例 ----
    # 樹：1-2, 1-3, 2-4, 3-5, 3-6（根 = 1）
    # 顏色：c1=1, c2=2, c3=1, c4=3, c5=1, c6=1
    #   total[1] = 4（節點 1,3,5,6）、total[2] = 1、total[3] = 1，palette = 3
    n = 6
    color = [0, 1, 2, 1, 3, 1, 1]
    parent = [0, 0, 1, 1, 2, 3, 3]
    q1 = [(1, 0), (1, 1), (1, 2), (1, 4), (1, 3), (2, 1), (3, 1), (4, 1),
          (5, 1), (2, 0), (3, 0), (1, 6), (6, 1)]
    queries: Dict[int, List[Tuple[int, int]]] = {}
    for i, (u, f) in enumerate(q1):
        queries.setdefault(u, []).append((f, i))
    s = Solver(n, color, parent, queries)
    got = s.run()
    want = solve_brute(n, color, parent, q1)
    assert got == want, list(zip(q1, got, want))
    # 手算驗證關鍵值
    # 整棵樹：顏色 1 出現 4 次、2 出現 1 次、3 出現 1 次
    #   distinct = 3；f=1 的顏色數 = 2；f=4 的顏色數 = 1；f=2、f=3 的顏色數 = 0
    #   shared = 0（每種顏色在子樹外都沒有出現）
    assert got[0] == (3, 0, 0), got[0]   # (u=1, f=0)：palette - distinct = 3-3 = 0
    assert got[1] == (3, 2, 0), got[1]   # (u=1, f=1)：顏色 2、3 各出現 1 次
    assert got[2] == (3, 0, 0), got[2]   # (u=1, f=2)：沒有顏色恰好出現 2 次
    assert got[3] == (3, 1, 0), got[3]   # (u=1, f=4)：顏色 1 出現 4 次
    assert got[4] == (3, 0, 0), got[4]   # (u=1, f=3)：沒有顏色恰好出現 3 次
    # 節點 2 的子樹 = {2,4}：顏色 2、3 各一次；兩者在整樹都只出現一次 → shared = 0
    assert got[5] == (2, 2, 0), got[5]   # (u=2, f=1)
    # 節點 3 的子樹 = {3,5,6}：顏色 1 出現 3 次
    #   distinct = 1；f=1 的顏色數 = 0；顏色 1 全樹 4 次、子樹內 3 次 → shared = 1
    assert got[6] == (1, 0, 1), got[6]   # (u=3, f=1)
    # 葉節點 4（顏色 3）：顏色 3 全樹只出現 1 次 → shared = 0
    assert got[7] == (1, 1, 0), got[7]   # (u=4, f=1)
    # 葉節點 5（顏色 1）：顏色 1 全樹 4 次、子樹內 1 次 → shared = 1
    assert got[8] == (1, 1, 1), got[8]   # (u=5, f=1)
    # f = 0 的退化定義
    assert got[9] == (2, 1, 0), got[9]   # (u=2, f=0)：palette - distinct = 3-2 = 1
    assert got[10] == (1, 2, 1), got[10]  # (u=3, f=0)：3-1 = 2
    # f 超過最大可能次數
    assert got[11] == (3, 0, 0), got[11]  # (u=1, f=6)
    # 葉節點 6（顏色 1）：shared = 1
    assert got[12] == (1, 1, 1), got[12]  # (u=6, f=1)

    # 單節點樹
    color2 = [0, 1]
    parent2 = [0, 0]
    s2 = Solver(1, color2, parent2, {1: [(0, 0), (1, 1), (2, 2)]})
    got2 = s2.run()
    assert got2 == solve_brute(1, color2, parent2, [(1, 0), (1, 1), (1, 2)]), got2

    # 所有節點同色
    n3 = 5
    color3 = [0] + [5] * 5
    parent3 = [0, 0, 1, 2, 3, 4]
    q3 = [(u, f) for u in range(1, 6) for f in (0, 1, 5)]
    queries3: Dict[int, List[Tuple[int, int]]] = {}
    for i, (u, f) in enumerate(q3):
        queries3.setdefault(u, []).append((f, i))
    s3 = Solver(n3, color3, parent3, queries3)
    assert s3.run() == solve_brute(n3, color3, parent3, q3)

    # ---- 隨機對拍 ----
    import random

    rng = random.Random(20261030)
    for _ in range(120):
        n = rng.randint(1, 12)
        color = [0] + [rng.randint(1, n) for _ in range(n)]
        parent = [0] * (n + 1)
        for u in range(2, n + 1):
            parent[u] = rng.randint(1, u - 1)
        q = rng.randint(1, 15)
        qlist = [(rng.randint(1, n), rng.randint(0, n)) for _ in range(q)]
        queries = {}
        for i, (u, f) in enumerate(qlist):
            queries.setdefault(u, []).append((f, i))
        s = Solver(n, color, parent, queries)
        assert s.run() == solve_brute(n, color, parent, qlist), (n, color, qlist)

    # 稍大規模：隨機樹 + 每個節點都問
    rng = random.Random(13572468)
    for _ in range(10):
        n = rng.randint(50, 250)
        color = [0] + [rng.randint(1, n) for _ in range(n)]
        parent = [0] * (n + 1)
        for u in range(2, n + 1):
            parent[u] = rng.randint(1, u - 1)
        qlist = [(u, rng.randint(0, n)) for u in range(1, n + 1)]
        queries = {}
        for i, (u, f) in enumerate(qlist):
            queries.setdefault(u, []).append((f, i))
        s = Solver(n, color, parent, queries)
        assert s.run() == solve_brute(n, color, parent, qlist)

    # 深鏈：驗證迭代實作不爆堆
    n4 = 4000
    color4 = [0] + [(u % 3) + 1 for u in range(1, n4 + 1)]
    parent4 = [0] * (n4 + 1)
    for u in range(2, n4 + 1):
        parent4[u] = u - 1
    qlist4 = [(u, 1) for u in range(1, n4 + 1)]
    queries4 = {}
    for i, (u, f) in enumerate(qlist4):
        queries4.setdefault(u, []).append((f, i))
    s4 = Solver(n4, color4, parent4, queries4)
    assert s4.run() == solve_brute(n4, color4, parent4, qlist4)


# ---------------------------------------------------------------- IO 模式

def _parse_int(tok: str) -> Optional[int]:
    """整數字串解析；非純整數或超出 ±10^18 回 None。"""
    if not tok:
        return None
    i = 0
    if tok[0] in "+-":
        i = 1
    if i == len(tok):
        return None
    for ch in tok[i:]:
        if ch < "0" or ch > "9":
            return None
    val = int(tok)
    if val > 10 ** 18:
        return None
    return val


def run_io(raw: str) -> None:
    toks = raw.split()
    pos = 0
    out: List[str] = []

    def nxt() -> Optional[str]:
        nonlocal pos
        if pos >= len(toks):
            return None
        t = toks[pos]
        pos += 1
        return t

    def bad() -> None:
        emit("error\n")

    n = _parse_int(nxt() or "")
    if n is None or not (1 <= n <= NMAX):
        bad()
        return
    color = [0] * (n + 1)
    for u in range(1, n + 1):
        v = _parse_int(nxt() or "")
        if v is None or not (1 <= v <= n):
            bad()
            return
        color[u] = v
    parent = [0] * (n + 1)
    for u in range(2, n + 1):
        p = _parse_int(nxt() or "")
        if p is None or not (1 <= p < u):
            bad()
            return
        parent[u] = p
    q = _parse_int(nxt() or "")
    if q is None or not (1 <= q <= QMAX):
        bad()
        return
    queries: Dict[int, List[Tuple[int, int]]] = {}
    for i in range(q):
        u = _parse_int(nxt() or "")
        f = _parse_int(nxt() or "")
        if (u is None or f is None or not (1 <= u <= n) or not (0 <= f <= n)):
            bad()
            return
        queries.setdefault(u, []).append((f, i))
    solver = Solver(n, color, parent, queries)
    answers = solver.run()
    for d, e, sh in answers:
        out.append(f"{d} {e} {sh}")
    for line in out:
        emit(line + "\n")


def emit(s: str) -> None:
    sys.stdout.buffer.write(s.encode("latin-1"))


def main() -> None:
    raw = sys.stdin.buffer.read().decode("latin-1")
    if raw.strip():
        run_io(raw)
    else:
        run_tests()
        emit("all tests passed\n")


if __name__ == "__main__":
    main()
