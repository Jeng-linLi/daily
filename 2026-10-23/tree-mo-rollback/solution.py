"""樹上莫隊（Mo's algorithm on tree）與迴滾莫隊（rollback Mo's algorithm）

題意：
    給一棵 n 個點的有根樹（根為 1 號點），每個點帶一個整數權值。
    q 個詢問，每次給 (u, v)，問 u 到 v 的**路徑上**有幾個不同的權值
    （路徑含端點；u = v 時路徑就是單一個點）。

    IO 格式（有 stdin 時）：
        第 1 行    n q
        第 2 行    a_1 a_2 ... a_n         權值
        第 3 行    p_2 p_3 ... p_n         點 i（2..n）的父節點
        接下來 q 行  u v
    輸出 q 行，每行一個整數。

思路：
    ### 1. 把「路徑查詢」拍平成一維「區間查詢」
    對樹做一次 DFS，記錄每個點的**進入時間** tin 與**離開時間** tout，
    得到一條長度 2n 的歐拉序（每個點恰好出現兩次：進入、離開）。
    對詢問 (u, v)（不妨設 tin[u] < tin[v]），路徑 u→v 在歐拉序上對應：

        - 若 u 是 v 的祖先：一段連續區間 [tin[u], tin[v]]，
          路徑上每個點在區間裡恰好出現一次，其餘點出現 0 或 2 次；
        - 否則：區間 [tin[u], tin[v]] 裡，u 出現兩次（進+出）、
          LCA(u, v) 出現零次，其餘路徑點各一次。
          於是「區間 + 額外翻轉 u 與 LCA」就正好得到整條路徑。

    翻轉用 vis[] 標記實現：區間裡出現兩次的點來回翻轉兩次等於沒翻，
    出現一次的點留在「啟用」狀態——這正是莫隊維護「區間內點集」的標準技巧。

    ### 2. 迴滾莫隊（rollback Mo）
    普通莫隊在區間上移動左右端點，需要「加」與「減」兩種操作；
    迴滾莫隊只保留「加」，左端點向左移產生的修改全部記在一棵**棧**上，
    回答完當前詢問後統一彈棧撤銷。適用於「刪不掉」的資料結構
    （例如只支援 insert 的平衡樹、二元組統計等）。

    具體結構（塊大小 B ≈ √m，m 為歐拉序長度）：
        - 詢問按 (l 所在塊, r) 排序，同一塊內 r 單調遞增；
        - 進入新塊時清空全部狀態，左右端點都停在**塊尾**；
        - r 只增不減（只加）；l 從塊尾往左加，記棧；
        - 回答完把棧彈回檢查點，l 復原到塊尾，r 保留；
        - 若某詢問的 r 落在塊內（r < 塊尾，整窗在塊中），
          先整體回滾到空窗，再一次性加入 [l, r]，答完再回滾——
          這類「短詢問」代價 O(B)，不破壞「只加」的不變式。

    本題的「加」= 翻轉一個點（維護每種權值的出現次數 cnt 與不同權值數 d），
    翻轉是自逆操作，棧裡記下被翻轉的點，反序再翻一次即精確撤銷。

    ### 3. LCA 用倍增法
    O(n log n) 預處理、O(log n) 查詢；也可以用歐拉序 + Sparse Table 做到 O(1)，
    但莫隊本身是 O(√n) 級別，倍增常數小、好寫，兩語言行為也完全一致。

    ### 4. 兩語言一致性
    - 答案與 DFS 順序無關（路徑點集是樹的固有屬性），歐拉序不同不影響輸出；
    - 莫隊的排序鍵 (l // B, r) 兩邊一致；平手不影響答案；
    - 詢問按輸入順序輸出，與處理順序解耦。

複雜度：
    預處理（DFS + 倍增）    O(n log n)
    單次詢問                O(√n)（均攤，迴滾莫隊標準複雜度）
    q 次詢問合計            O((n + q) √n)
    空間                    O(n log n)

應用場景：
    - 遊戲 / 社交網路：「共同好友」「幫派勢力範圍」等路徑聚合統計
    - 版本控制系統：兩次 commit 之間修改過的檔案種類數（diff 統計）
    - 資料庫查詢最佳化：樹形目錄下「某兩節點路徑上的標籤數」
    - 離線演算法競賽經典題型：靜態樹路徑第 k 小 / 種類數 / 權值和
"""

import random
import re
import sys
from typing import List, Tuple

INT_RE = re.compile(r"^[+-]?[0-9]+$")

NMAX = 20000          # IO 模式的節點數上限
QMAX = 20000          # IO 模式的詢問數上限
VMAX = 10 ** 9        # 權值絕對值上限
CLAMP_ABS = 10 ** 18  # 通用整數夾取範圍


# ---------------------------------------------------------------- 樹與 LCA


class TreeMoSolver:
    """樹上莫隊（迴滾版）的前處理：歐拉序 + 深度 + 倍增 LCA。"""

    def __init__(self, n: int, values: List[int], parent: List[int]) -> None:
        """values / parent 均為 0 起算；parent[root] 必須是 -1。"""
        self.n = n
        self.values = values
        self.parent = parent
        if n == 0:
            self.root = 0
            self.first = []
            self.tout = []
            self.tour = []
            self.depth = []
            self.up = [[]]
            self.log = 1
            return

        children: List[List[int]] = [[] for _ in range(n)]
        root = 0
        for v in range(n):
            p = parent[v]
            if p < 0:
                root = v
            elif 0 <= p < n and p != v:
                children[p].append(v)
        # 自環一律接回根，保證任何輸入都能終止
        for v in range(n):
            if parent[v] == v:
                children[root].append(v)
        self.root = root

        # ---- 迭代 DFS 求歐拉序（進入 / 離開各記一次）----
        first = [0] * n
        tout = [0] * n
        tour: List[int] = []
        stack: List[Tuple[int, bool]] = [(root, True)]
        while stack:
            v, enter = stack.pop()
            if enter:
                first[v] = len(tour)
                tour.append(v)
                stack.append((v, False))
                # 反序入棧，讓編號小的子樹先遍歷（決定性，方便對拍）
                for c in reversed(children[v]):
                    stack.append((c, True))
            else:
                tout[v] = len(tour)
                tour.append(v)
        self.first = first
        self.tout = tout
        self.tour = tour

        # ---- BFS 求深度 ----
        depth = [0] * n
        order = [root]
        seen = [False] * n
        seen[root] = True
        head = 0
        while head < len(order):
            v = order[head]
            head += 1
            for c in children[v]:
                if not seen[c]:
                    seen[c] = True
                    depth[c] = depth[v] + 1
                    order.append(c)
        self.depth = depth

        # ---- 倍增 ----
        log = max(1, n.bit_length())
        up = [[0] * n for _ in range(log)]
        up[0] = [p if 0 <= p < n and p != v else root
                 for v, p in enumerate(parent)]
        up[0][root] = root
        for k in range(1, log):
            upk = up[k]
            upk1 = up[k - 1]
            for v in range(n):
                upk[v] = upk1[upk1[v]]
        self.up = up
        self.log = log

    def lca(self, u: int, v: int) -> int:
        if self.depth[u] < self.depth[v]:
            u, v = v, u
        diff = self.depth[u] - self.depth[v]
        k = 0
        while diff:
            if diff & 1:
                u = self.up[k][u]
            diff >>= 1
            k += 1
        if u == v:
            return u
        for k in range(self.log - 1, -1, -1):
            if self.up[k][u] != self.up[k][v]:
                u = self.up[k][u]
                v = self.up[k][v]
        return self.parent[u] if self.parent[u] >= 0 else u

    def to_window(self, u: int, v: int) -> Tuple[int, int, Tuple[int, ...]]:
        """詢問 → (歐拉序區間 [l, r], 需要額外翻轉的點)。"""
        i, j = self.first[u], self.first[v]
        if i > j:
            u, v = v, u
            i, j = j, i
        l = self.lca(u, v)
        if l == u:
            return (i, j, ())
        return (i, j, (u, l))


# ---------------------------------------------------------------- 迴滾莫隊


class RollbackMo:
    """只支援「加」的莫隊：左端點的修改記棧，答完即撤銷。"""

    def __init__(self, tour: List[int], values: List[int],
                 nvals: int, n: int) -> None:
        self.tour = tour
        self.values = values
        self.m = len(tour)
        self.B = max(1, int(self.m ** 0.5))
        self.n = n
        self.cnt = [0] * nvals
        self.vis = [False] * n
        self.d = 0
        self.stack: List[int] = []
        self.curL = 0
        self.curR = -1

    def _toggle(self, x: int) -> None:
        """翻轉點 x 的啟用狀態（自逆操作，不記棧）。"""
        c = self.values[x]
        if self.vis[x]:
            self.cnt[c] -= 1
            if self.cnt[c] == 0:
                self.d -= 1
            self.vis[x] = False
        else:
            if self.cnt[c] == 0:
                self.d += 1
            self.cnt[c] += 1
            self.vis[x] = True

    def _add_node(self, x: int) -> None:
        self._toggle(x)
        self.stack.append(x)

    def _rollback(self, snap_size: int, cur_l: int, cur_r: int) -> None:
        while len(self.stack) > snap_size:
            x = self.stack.pop()
            self._toggle(x)
        self.curL = cur_l
        self.curR = cur_r

    def solve(self, windows: List[Tuple[int, int, Tuple[int, ...]]]) -> List[int]:
        """windows[i] = (l, r, extras)；回傳與輸入同序的答案。"""
        q = len(windows)
        ans = [0] * q
        order = sorted(range(q),
                       key=lambda idx: (windows[idx][0] // self.B, windows[idx][1]))
        cur_blk = -1
        block_r = 0
        for idx in order:
            l, r, extras = windows[idx]
            blk = l // self.B
            if blk != cur_blk:
                # 新塊：清空全部狀態，端點停在塊尾
                self.stack.clear()
                self.cnt = [0] * len(self.cnt)
                self.vis = [False] * self.n
                self.d = 0
                block_r = min((blk + 1) * self.B, self.m) - 1
                self.curL = block_r
                self.curR = block_r - 1
                cur_blk = blk
            if r < block_r:
                # 短詢問：整個窗口在塊內 → 回滾到空窗再一次加滿 [l, r]
                self._rollback(0, block_r, block_r - 1)
                self.curR = r
                self.curL = r + 1
                while self.curL > l:
                    self.curL -= 1
                    self._add_node(self.tour[self.curL])
                for x in extras:
                    self._add_node(x)
                ans[idx] = self.d
                self._rollback(0, block_r, block_r - 1)
            else:
                # 長詢問：r 只增（只加），l 從塊尾往左加（記棧，答完撤銷）
                while self.curR < r:
                    self.curR += 1
                    self._add_node(self.tour[self.curR])
                save = len(self.stack)
                while self.curL > l:
                    self.curL -= 1
                    self._add_node(self.tour[self.curL])
                for x in extras:
                    self._add_node(x)
                ans[idx] = self.d
                # 只撤銷 l 的部分；r 的部分保留給同一塊的後續詢問
                self._rollback(save, block_r, self.curR)
        return ans


# ---------------------------------------------------------------- 暴力對拍


def bf_distinct(n: int, values: List[int], parent: List[int],
                u: int, v: int) -> int:
    """沿父指標爬行收集路徑點集，統計不同權值（完全獨立的參考實作）。"""
    anc = set()
    x, steps = u, 0
    while x >= 0 and steps <= n:
        anc.add(x)
        x = parent[x]
        steps += 1
    y, steps = v, 0
    while y >= 0 and steps <= n and y not in anc:
        y = parent[y]
        steps += 1
    l = y if y >= 0 else u
    seen = set()
    x, steps = u, 0
    while x != l and steps <= n:
        seen.add(values[x])
        x = parent[x] if parent[x] >= 0 else l
        steps += 1
    seen.add(values[l])
    x, steps = v, 0
    while x != l and steps <= n:
        seen.add(values[x])
        x = parent[x] if parent[x] >= 0 else l
        steps += 1
    return len(seen)


# ---------------------------------------------------------------- 主流程


def run(n: int, values: List[int], parent: List[int],
        queries: List[Tuple[int, int]]) -> str:
    solver = TreeMoSolver(n, values, parent)
    # 權值離散化
    uniq = sorted(set(values))
    comp = {v: i for i, v in enumerate(uniq)}
    cvals = [comp[v] for v in values]
    windows = [solver.to_window(u, v) for (u, v) in queries]
    mo = RollbackMo(solver.tour, cvals, len(uniq), n)
    ans = mo.solve(windows)
    out = [str(a) for a in ans]
    return "\n".join(out) + ("\n" if out else "")


# ---------------------------------------------------------------- 內置測試


def run_tests() -> None:
    # ---- 固定用例：樹 1-2, 1-3, 2-4, 2-5（根 1）----
    # 權值（0 起算）：[10, 20, 30, 20, 10]
    n = 5
    values = [10, 20, 30, 20, 10]
    parent = [-1, 0, 0, 1, 1]
    qs = [(3, 4), (2, 3), (3, 3), (0, 4), (1, 2), (4, 2)]
    # (3,4): 路徑 3-1-4      → {20,10}       → 2
    # (2,3): 路徑 2-0-1-3    → {30,10,20}    → 3
    # (3,3): 單點            → {20}          → 1
    # (0,4): 路徑 0-1-4      → {10,20}       → 2
    # (1,2): 路徑 1-0-2      → {20,10,30}    → 3
    # (4,2): 路徑 4-1-0-2    → {10,20,30}    → 3
    want = [2, 3, 1, 2, 3, 3]
    got = run(n, values, parent, qs)
    assert got.rstrip("\n").split("\n") == [str(w) for w in want], got
    # 與暴力對拍一致
    for (u, v), w in zip(qs, want):
        assert bf_distinct(n, values, parent, u, v) == w

    # ---- 單點樹 ----
    assert run(1, [7], [-1], [(0, 0)]) == "1\n"

    # ---- 沒有詢問 → 空輸出 ----
    assert run(3, [1, 2, 3], [-1, 0, 0], []) == ""

    # ---- 全部權值相同 → 每個答案都是 1 ----
    assert run(4, [5, 5, 5, 5], [-1, 0, 0, 0],
               [(0, 1), (1, 3), (2, 3)]) == "1\n1\n1\n"

    # ---- 鏈狀樹（最深路徑，壓力測試塊邏輯）----
    n = 1500
    values = [((i * 37) % 9) for i in range(n)]
    parent = list(range(-1, n - 1))
    rng = random.Random(777)
    queries = [(rng.randrange(n), rng.randrange(n)) for _ in range(150)]
    got = [int(x) for x in run(n, values, parent, queries).rstrip("\n").split("\n")]
    want2 = [bf_distinct(n, values, parent, u, v) for (u, v) in queries]
    assert got == want2

    # ---- 隨機對拍：四種樹形 × 兩種值域 ----
    rng = random.Random(20261023)
    for t in range(400):
        n = rng.randint(1, 60)
        shape = t % 4
        if shape == 0:                      # 隨機樹
            parent = [-1] + [rng.randrange(v) for v in range(1, n)]
        elif shape == 1:                    # 鏈
            parent = list(range(-1, n - 1))
        elif shape == 2:                    # 星（根帶全部子樹）
            parent = [-1] + [0] * (n - 1)
        else:                               # 雙層：一半掛根、一半掛 1 號長子
            parent = [-1] + [0] * (n - 1)
            for v in range(2, n):
                if rng.random() < 0.5:
                    parent[v] = 1
        small = (t % 2 == 0)
        if small:
            values = [rng.randint(0, 3) for _ in range(n)]     # 強制大量重複
        else:
            values = [rng.randint(-10 ** 9, 10 ** 9) for _ in range(n)]
        q = rng.randint(0, 20)
        queries = [(rng.randrange(n), rng.randrange(n)) for _ in range(q)]
        out = run(n, values, parent, queries)
        got = [int(x) for x in out.rstrip("\n").split("\n")] if q else []
        want2 = [bf_distinct(n, values, parent, u, v) for (u, v) in queries]
        assert got == want2, (t, shape, small, n, values, parent, queries)

    # ---- 隨機大樹：2000 點、400 詢問 ----
    rng2 = random.Random(777)
    n = 2000
    values = [rng2.randint(0, 50) for _ in range(n)]
    parent = [-1] + [rng2.randrange(v) for v in range(1, n)]
    queries = [(rng2.randrange(n), rng2.randrange(n)) for _ in range(400)]
    got = [int(x) for x in run(n, values, parent, queries).rstrip("\n").split("\n")]
    want2 = [bf_distinct(n, values, parent, u, v) for (u, v) in queries]
    assert got == want2


# ---------------------------------------------------------------- IO


def emit(text: str) -> None:
    """二進位輸出，避免 Windows 把 \n 換成 \r\n。"""
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
        try:
            v = int(t)
        except ValueError:
            return default
        return max(-CLAMP_ABS, min(CLAMP_ABS, v))

    def cl(v: int, lo: int, hi: int) -> int:
        return max(lo, min(hi, v))

    n = cl(next_int(0), 0, NMAX)
    q = cl(next_int(0), 0, QMAX)
    values = [cl(next_int(0), -VMAX, VMAX) for _ in range(n)]
    parent = [-1] * n
    for v in range(1, n):
        p = cl(next_int(1), 1, n)
        # 自環一律當成指向根，保證後續爬行終止
        parent[v] = p - 1 if p - 1 != v else 0
    queries: List[Tuple[int, int]] = []
    for _ in range(q):
        if n == 0:
            break
        u = cl(next_int(1), 1, n) - 1
        v = cl(next_int(1), 1, n) - 1
        queries.append((u, v))
    emit(run(n, values, parent, queries))


def main() -> None:
    raw = sys.stdin.buffer.read().decode("latin-1")
    if raw.strip():
        run_io(raw)
    else:
        run_tests()
        emit("all tests passed\n")


if __name__ == "__main__":
    main()
