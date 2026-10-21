"""Link-Cut Tree（LCT）：用 splay 維護「動態森林」，支援加邊 / 刪邊 / 路徑查詢 / 路徑修改

題意：
    有 n 個節點（編號 1..n），每個節點有一個整數權值。一開始整張圖**沒有邊**。
    依序處理 q 個操作，操作會動態地連邊、斷邊（整個過程圖始終保持是森林），
    並在任意時刻查詢／修改「u 到 v 這條路徑」上的資訊。

    操作碼（每行一個操作）：
        1  u v        連邊（若 u、v 已連通則忽略，避免產生環）
        2  u v        斷邊（若邊 (u,v) 不存在則忽略）
        3  u v        路徑查詢，輸出「權值和 最大值 最小值 邊數」；不連通輸出 -1
        4  u v w      路徑修改：u 到 v 路徑上每個點權都加 w（不連通則忽略）
        5  u v        連通性判定，輸出 1 / 0
        6  u x        單點改值：把 u 的權值設成 x
        7  u v        輸出路徑上的點權（由 u 走到 v 的順序）；不連通輸出 -1

思路：
    ### 1. 為什麼需要 LCT
    靜態樹上做路徑查詢，樹鏈剖分（heavy-light）就夠了；但它**不能加邊刪邊**。
    LCT 用「**實鏈剖分**」解決這件事：把每個點的兒子中**至多一個**標成「實邊」，
    其餘是「虛邊」。實邊串成的鏈用一棵 **splay** 維護（鍵是中序對應的樹上深度），
    虛邊只記錄「這棵 splay 的根的父親是誰」（path-parent，單向）。

    於是「森林」被拆成一棵棵 splay。兩個關鍵操作：

    - **`access(x)`**：把根到 x 的路徑整條變成**實鏈**。
      反覆「splay(x) → 把 x 的右子樹換成上一輪處理過的鏈 → 往上跳 path-parent」。
      做完 x 所在的那棵 splay 就恰好含根到 x 的全部節點，且 x 是 splay 根。
    - **`makeroot(x)`**：`access(x)` 之後翻轉整棵 splay（懶標記 rev），
      把 x 變成整棵樹的**代表根**。這正是 LCT 能做到「無根樹任意路徑」的原因。

    有了 `split(u,v) = makeroot(u); access(v)`，v 的那棵 splay 就**恰好**是 u..v 的路徑，
    路徑查詢／修改全部退化成「對一個 splay 根做一次 O(1) 讀取 / O(log n) 懶標記」。

    ### 2. 各操作的複雜度
    全部都是 **均攤 O(log n)**（splay 的攤還分析）。這是目前已知動態樹的主流最優解。

    ### 3. 兩語言必須一致的五個細節
    - `access` 的寫法統一用 `for y in ...` 迴圈，且**不要**額外維護 `fa[last]`：
      `x = fa[y = x]` 這個賦值順序本身就保證了 `fa[y] == x`，多寫一次反而可能寫錯。
    - `isRoot(x)` 的定義：`fa[x] == 0` 或 x 不是 `fa[x]` 的**任何一個**兒子。
      （splay 的父親 vs 虛邊的 path-parent 靠這個區分。）
    - `cut(u,v)` 的成立條件是 `ch[v][0] == u and ch[u][1] == 0`：
      前者表示 u 是 v 的前驅，後者表示 u..v 之間沒有第三個點 —— 兩者合起來才是「直接相鄰」。
    - `rev` 與 `add` 兩個懶標記的**下推順序**固定為「先 rev 再 add」。
    - 路徑輸出用 splay 的**中序遍歷**：`makeroot(u); access(v)` 之後中序就是 u → v。
      遍歷時每一步都要先 `pushdown` 再往下走。

    ### 4. 容易踩的坑
    - 哨兵節點 0 的 `mx[0]` 必須是 -INF、`mn[0]` 必須是 +INF，否則取 max/min 會被 0 污染。
    - `pushadd` 要乘的是 **splay 的子樹大小**（不是鏈長），所以 `sz` 必須維護。
    - `findroot` 要先 `access(x)`，往左下角走到底，中途每一層都要 `pushdown`，
      否則 rev 標記沒解開會走錯方向。

應用場景：
    - 動態連通性：網路拓撲變動（連線／斷線）下的即時連通查詢
    - 線上維護生成樹：邊權會改、邊會加刪的最小生成樹（動態 MST）
    - 交通路網：施工封路（cut）與新路開通（link）後的兩地最長／最短瓶頸路段
    - 版本管理：分支合併圖上查詢兩版本間的修改總量
    - 遊戲伺服器：動態變化的地圖傳送網路上的路徑統計

複雜度（n 個點，q 次操作）：
    link / cut / 路徑查詢 / 路徑修改     均攤 O(log n)
    單點改值                             均攤 O(log n)
    空間                                 O(n)
"""

import random
import re
import sys
from collections import deque
from typing import Deque, List, Optional

INT_RE = re.compile(r"^[+-]?[0-9]+$")

NMAX = 300          # IO 模式下的點數上限
QMAX = 400          # IO 模式下的操作數上限
VMAX = 10 ** 6      # 點權與路徑加值的絕對值上限
INF = 1 << 60


# ------------------------------------------------------------------ LCT


class LCT:
    """Link-Cut Tree：陣列版（0 號是哨兵，節點編號 1..n）。"""

    def __init__(self, vals: List[int]) -> None:
        n = len(vals)
        self.n = n
        self.ch = [[0, 0] for _ in range(n + 1)]
        self.fa = [0] * (n + 1)
        self.rev = [0] * (n + 1)
        self.lazy = [0] * (n + 1)
        self.sz = [0] * (n + 1)
        self.val = [0] * (n + 1)
        self.sm = [0] * (n + 1)
        self.mx = [0] * (n + 1)
        self.mn = [0] * (n + 1)
        self.mx[0] = -INF
        self.mn[0] = INF
        for i in range(1, n + 1):
            v = vals[i - 1]
            self.val[i] = v
            self.sm[i] = v
            self.mx[i] = v
            self.mn[i] = v
            self.sz[i] = 1

    # -- 維護 --

    def _pushup(self, x: int) -> None:
        l = self.ch[x][0]
        r = self.ch[x][1]
        self.sz[x] = 1 + self.sz[l] + self.sz[r]
        self.sm[x] = self.val[x] + self.sm[l] + self.sm[r]
        m = self.val[x]
        if self.mx[l] > m:
            m = self.mx[l]
        if self.mx[r] > m:
            m = self.mx[r]
        self.mx[x] = m
        m = self.val[x]
        if self.mn[l] < m:
            m = self.mn[l]
        if self.mn[r] < m:
            m = self.mn[r]
        self.mn[x] = m

    def _pushrev(self, x: int) -> None:
        if x == 0:
            return
        c = self.ch[x]
        c[0], c[1] = c[1], c[0]
        self.rev[x] ^= 1

    def _pushadd(self, x: int, w: int) -> None:
        if x == 0:
            return
        self.val[x] += w
        self.sm[x] += w * self.sz[x]
        self.mx[x] += w
        self.mn[x] += w
        self.lazy[x] += w

    def _pushdown(self, x: int) -> None:
        if self.rev[x]:
            self._pushrev(self.ch[x][0])
            self._pushrev(self.ch[x][1])
            self.rev[x] = 0
        if self.lazy[x] != 0:
            self._pushadd(self.ch[x][0], self.lazy[x])
            self._pushadd(self.ch[x][1], self.lazy[x])
            self.lazy[x] = 0

    def _is_root(self, x: int) -> bool:
        f = self.fa[x]
        if f == 0:
            return True
        return self.ch[f][0] != x and self.ch[f][1] != x

    def _rotate(self, x: int) -> None:
        y = self.fa[x]
        z = self.fa[y]
        k = 1 if self.ch[y][1] == x else 0
        if not self._is_root(y):
            if self.ch[z][1] == y:
                self.ch[z][1] = x
            else:
                self.ch[z][0] = x
        self.fa[x] = z
        self.ch[y][k] = self.ch[x][k ^ 1]
        if self.ch[x][k ^ 1]:
            self.fa[self.ch[x][k ^ 1]] = y
        self.ch[x][k ^ 1] = y
        self.fa[y] = x
        self._pushup(y)
        self._pushup(x)

    def _splay(self, x: int) -> None:
        st: List[int] = []
        y = x
        st.append(y)
        while not self._is_root(y):
            y = self.fa[y]
            st.append(y)
        for i in range(len(st) - 1, -1, -1):
            self._pushdown(st[i])
        while not self._is_root(x):
            y = self.fa[x]
            z = self.fa[y]
            if not self._is_root(y):
                if (self.ch[y][0] == x) != (self.ch[z][0] == y):
                    self._rotate(x)
                else:
                    self._rotate(y)
            self._rotate(x)
        self._pushup(x)

    def _access(self, x: int) -> None:
        last = 0
        y = x
        while y != 0:
            self._splay(y)
            self.ch[y][1] = last
            self._pushup(y)
            last = y
            y = self.fa[y]
        self._splay(x)

    def _makeroot(self, x: int) -> None:
        self._access(x)
        self._pushrev(x)

    def _findroot(self, x: int) -> int:
        self._access(x)
        y = x
        while self.ch[y][0] != 0:
            self._pushdown(y)
            y = self.ch[y][0]
        self._splay(y)
        return y

    def _split(self, u: int, v: int) -> None:
        self._makeroot(u)
        self._access(v)

    # -- 對外介面 --

    def connected(self, u: int, v: int) -> bool:
        if u == v:
            return True
        self._makeroot(u)
        return self._findroot(v) == u

    def link(self, u: int, v: int) -> None:
        if u == v:
            return
        self._makeroot(u)
        if self._findroot(v) == u:
            return          # 已連通，連了會成環
        self.fa[u] = v

    def cut(self, u: int, v: int) -> None:
        if u == v:
            return
        self._makeroot(u)
        self._access(v)
        # v 已是 splay 根；u 必須緊貼在 v 左邊、且 u 右邊沒有別的點
        if self.ch[v][0] == u and self.ch[u][1] == 0:
            self.ch[v][0] = 0
            self.fa[u] = 0
            self._pushup(v)

    def path_query(self, u: int, v: int) -> Optional[List[int]]:
        if not self.connected(u, v):
            return None
        self._split(u, v)
        return [self.sm[v], self.mx[v], self.mn[v], self.sz[v] - 1]

    def path_add(self, u: int, v: int, w: int) -> None:
        if not self.connected(u, v):
            return
        self._split(u, v)
        self._pushadd(v, w)

    def set_value(self, u: int, x: int) -> None:
        self._access(u)          # u 變成 splay 根，直接改權值即可
        self.val[u] = x
        self._pushup(u)

    def path_values(self, u: int, v: int) -> Optional[List[int]]:
        if not self.connected(u, v):
            return None
        self._split(u, v)
        res: List[int] = []
        st: List[int] = []
        cur = v
        while st or cur != 0:
            while cur != 0:
                self._pushdown(cur)
                st.append(cur)
                cur = self.ch[cur][0]
            cur = st.pop()
            res.append(self.val[cur])
            cur = self.ch[cur][1]
        return res


# ------------------------------------------------------- 暴力對照（一般圖森林）


class Brute:
    """O(n) 的參考實作：鄰接表 + BFS，只求正確、不求快。"""

    def __init__(self, vals: List[int]) -> None:
        self.val = [0] + list(vals)
        self.adj: List[set] = [set() for _ in range(len(vals) + 1)]

    def _path(self, u: int, v: int) -> Optional[List[int]]:
        if u == v:
            return [u]
        n = len(self.adj) - 1
        par = [0] * (n + 1)
        par[u] = -1
        q: Deque[int] = deque([u])
        while q:
            x = q.popleft()
            for y in sorted(self.adj[x]):
                if par[y] == 0:
                    par[y] = x
                    if y == v:
                        path = [v]
                        while path[-1] != u:
                            path.append(par[path[-1]])
                        path.reverse()
                        return path
                    q.append(y)
        return None

    def connected(self, u: int, v: int) -> bool:
        return self._path(u, v) is not None

    def link(self, u: int, v: int) -> None:
        if self.connected(u, v):
            return
        self.adj[u].add(v)
        self.adj[v].add(u)

    def cut(self, u: int, v: int) -> None:
        self.adj[u].discard(v)
        self.adj[v].discard(u)

    def path_query(self, u: int, v: int) -> Optional[List[int]]:
        p = self._path(u, v)
        if p is None:
            return None
        vs = [self.val[i] for i in p]
        return [sum(vs), max(vs), min(vs), len(p) - 1]

    def path_add(self, u: int, v: int, w: int) -> None:
        p = self._path(u, v)
        if p is None:
            return
        for i in p:
            self.val[i] += w

    def set_value(self, u: int, x: int) -> None:
        self.val[u] = x

    def path_values(self, u: int, v: int) -> Optional[List[int]]:
        p = self._path(u, v)
        if p is None:
            return None
        return [self.val[i] for i in p]


# ------------------------------------------------------------------ 驅動


def run(vals: List[int], ops: List[List[int]]) -> str:
    n = len(vals)
    lct = LCT(vals)
    bf = Brute(vals)
    out: List[str] = []

    def ok(u: int) -> bool:
        return 1 <= u <= n

    for op in ops:
        if not op:
            continue
        code = op[0]
        if code == 1:
            u, v = op[1], op[2]
            if ok(u) and ok(v):
                lct.link(u, v)
                bf.link(u, v)
        elif code == 2:
            u, v = op[1], op[2]
            if ok(u) and ok(v):
                lct.cut(u, v)
                bf.cut(u, v)
        elif code == 3:
            u, v = op[1], op[2]
            if ok(u) and ok(v):
                r = lct.path_query(u, v)
                e = bf.path_query(u, v)
                assert r == e, ("path_query", u, v, r, e)
                out.append("-1" if r is None else "%d %d %d %d" % (r[0], r[1], r[2], r[3]))
        elif code == 4:
            u, v, w = op[1], op[2], op[3]
            if ok(u) and ok(v):
                lct.path_add(u, v, w)
                bf.path_add(u, v, w)
        elif code == 5:
            u, v = op[1], op[2]
            if ok(u) and ok(v):
                r = lct.connected(u, v)
                e = bf.connected(u, v)
                assert r == e, ("connected", u, v, r, e)
                out.append("1" if r else "0")
        elif code == 6:
            u, x = op[1], op[2]
            if ok(u):
                lct.set_value(u, x)
                bf.set_value(u, x)
        elif code == 7:
            u, v = op[1], op[2]
            if ok(u) and ok(v):
                r = lct.path_values(u, v)
                e = bf.path_values(u, v)
                assert r == e, ("path_values", u, v, r, e)
                if r is None:
                    out.append("-1")
                else:
                    out.append("%d %s" % (len(r), " ".join(str(x) for x in r)))
    return "".join(s + "\n" for s in out)


def _check_case(vals: List[int], ops: List[List[int]]) -> None:
    run(vals, ops)


def run_tests() -> None:
    # 固定用例 1：鏈 1-2-3-4，權值 1,2,3,4
    _check_case([1, 2, 3, 4], [
        [1, 1, 2], [1, 2, 3], [1, 3, 4],
        [3, 1, 4], [7, 1, 4], [5, 1, 4], [5, 1, 3],
        [4, 2, 3, 10], [3, 1, 4],
        [2, 2, 3], [5, 1, 4], [3, 1, 4],
    ])

    # 固定用例 2：兩棵樹，反覆 link / cut
    _check_case([5, -3, 7, 0, 2], [
        [1, 1, 2], [1, 3, 4],
        [5, 1, 3], [3, 3, 4],
        [1, 2, 3], [5, 1, 5], [3, 1, 5],
        [1, 4, 5], [3, 1, 5], [7, 5, 1],
        [6, 5, 100], [3, 1, 5],
        [2, 3, 4], [5, 3, 4],
    ])

    # 固定用例 3：單點 / 非法編號 / 重複連邊
    _check_case([42], [
        [3, 1, 1], [7, 1, 1], [5, 1, 1],
        [1, 1, 1], [5, 1, 1],
        [2, 1, 1], [3, 1, 2], [5, 0, 1], [3, 1, 0],
    ])

    # 固定用例 4：n = 0
    _check_case([], [[3, 1, 1], [5, 1, 2]])

    # 隨機對拍：隨機 link / cut / 查詢
    rnd = random.Random(20261021)
    for _ in range(220):
        n = rnd.randint(1, 16)
        vals = [rnd.randint(-40, 40) for _ in range(n)]
        ops: List[List[int]] = []
        for _ in range(rnd.randint(1, 40)):
            code = rnd.choice([1, 1, 2, 3, 3, 4, 5, 6, 7])
            u = rnd.randint(1, n)
            v = rnd.randint(1, n)
            if code == 4:
                ops.append([4, u, v, rnd.randint(-20, 20)])
            elif code == 6:
                ops.append([6, u, rnd.randint(-50, 50)])
            else:
                ops.append([code, u, v])
        _check_case(vals, ops)

    # 隨機對拍：只連邊（保證森林）、長路徑
    for _ in range(40):
        n = rnd.randint(2, 30)
        vals = [rnd.randint(-100, 100) for _ in range(n)]
        ops: List[List[int]] = []
        for _ in range(rnd.randint(1, 60)):
            code = rnd.choice([1, 3, 4, 5, 7, 2])
            u = rnd.randint(1, n)
            v = rnd.randint(1, n)
            if code == 4:
                ops.append([4, u, v, rnd.randint(-30, 30)])
            else:
                ops.append([code, u, v])
        _check_case(vals, ops)


def emit(text: str) -> None:
    sys.stdout.buffer.write(text.encode("latin-1"))


def run_io(raw: str) -> None:
    toks = raw.split()
    pos = 0

    def nxt(default: int = 0) -> int:
        nonlocal pos
        if pos >= len(toks):
            return default
        t = toks[pos]
        pos += 1
        if not INT_RE.match(t):
            return default
        try:
            v = int(t)
        except ValueError:
            return default
        return max(-(10 ** 18), min(10 ** 18, v))

    def cl(v: int) -> int:
        return max(-VMAX, min(VMAX, v))

    n = nxt(0)
    if n < 0:
        n = 0
    if n > NMAX:
        n = NMAX
    vals = [cl(nxt(0)) for _ in range(n)]
    q = nxt(0)
    if q < 0:
        q = 0
    if q > QMAX:
        q = QMAX
    ops: List[List[int]] = []
    for _ in range(q):
        code = nxt(0)
        if code == 1 or code == 2 or code == 3 or code == 5 or code == 7:
            ops.append([code, nxt(1), nxt(1)])
        elif code == 4:
            u, v = nxt(1), nxt(1)
            ops.append([4, u, v, cl(nxt(0))])
        elif code == 6:
            ops.append([6, nxt(1), cl(nxt(0))])
        else:
            continue
    emit(run(vals, ops))


def main() -> None:
    raw = sys.stdin.buffer.read().decode("latin-1")
    if raw.strip():
        run_io(raw)
    else:
        run_tests()
        emit("all tests passed\n")


if __name__ == "__main__":
    main()
