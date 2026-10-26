"""帶權併查集：關係推斷（奇偶性 / 食物鏈）與擴展域建模

題意：
    給定 n 個元素，元素之間的關係是一個模 d 的「循環群」：
    d = 2 時 0 = 同類、1 = 不同類；d = 3 時 0 = 同類、1 = x 吃 y、2 = x 被 y 吃。
    關係可傳遞：已知 rel(a,b) 與 rel(b,c) 就能推出
        rel(a,c) = (rel(a,b) + rel(b,c)) mod d
    依序給出 m 條陳述 (x, y, r)，表示「rel(x,y) = r」。
    與已知資訊矛盾的陳述稱為**假話**；假話不會被記錄，不影響後續推斷。
    統計假話的數量。

    限制：d ∈ {2, 3}，1 <= n, m <= 100000，1 <= x, y <= n，0 <= r < d。

思路：
    ### 1. 關係的代數結構
    把每個元素看成一個「種類值」k[i] ∈ Z_d，定義
        rel(x, y) = (k[x] - k[y]) mod d
    於是 rel(a,c) = rel(a,b) + rel(b,c) 自動成立 —— 傳遞性免費得到。
    一條陳述 (x, y, r) 本質上是對兩個元素的種類差給出一個約束；
    矛盾 ⟺ 兩個約束給出了不同的差。

    ### 2. 帶權併查集（主實作）
    併查集的每個節點附一個權重 rel[x] = rel(x, parent[x])。
    查找時做路徑壓縮，順手把權重累加成「相對於根」的
        rel(x, root) = rel(x, parent) + rel(parent, root)
    合併兩棵樹時，根據要維護的 rel(x, y) = r 解出新邊的權重：
        rel(ry, rx) = (rel[x] - rel[y] - r) mod d
    查詢時若 x, y 已連通，直接比較 (rel[x] - rel[y]) mod d 與 r。
    兩次操作都接近 O(α(n))。

    ### 3. 擴展域併查集（第二實作，測試對拍用）
    另一個經典建模：把每個元素複製成 d 個「域節點」(i, t)
    表示「元素 i 的種類是 t」。陳述 (x, y, r) 等價於
        對所有 u：union(x + u·n, y + ((u - r) mod d)·n)
    假話的判定也變得很直接：rel(x,y) = t 已知
        ⟺ find(x + t·n) == find(y)
    所以只要存在 t ≠ r 使得 find(x + t·n) == find(y)，就是假話。
    x == y 且 r ≠ 0 的假話也會被正確抓到（此時 t = 0 就會命中）。
    域節點之間永不矛盾，正是「種類唯一」的編碼。

    ### 4. 三套實作互相對拍
    測試裡同時跑 帶權併查集 / 擴展域併查集 / 樸素傳遞閉包（O(n^3)），
    隨機數據三方結果必須一致 —— 關係取模的邊界（負數、x==y、
    跨合併才暴露的矛盾）極易寫錯，多一套實作就多一道保險。

複雜度（α 為阿克曼函數的反函數）：
    帶權併查集        每次操作接近 O(α(n))
    擴展域併查集      每次操作 O(d·α(n))，此處 d <= 3
    樸素閉包（僅測試）O(n^3) 預處理 + O(n^2) 每條陳述

應用場景：
    - 食物鏈推斷（NOI2001 經典題）：根據觀察記錄判定哪些記錄是假的
    - 奇偶性 / 二分圖判定：d = 2 就是「同色 / 異色」，離線判斷圖是否
      二分圖、找出破壞二分性的矛盾邊
    - 種類合一（unification）：Prolog / 類型推導裡「變數相等或不等」
      的約束求解，關係群從 Z_d 推廣到更一般的置換群
    - 版本相容性：模組 A 與 B 相容、B 與 C 不相容……傳遞推斷衝突組合
    - 密碼學份額：Shamir 份額之間「同組 / 不同組」關係的批量校驗
"""

import random
import sys
from typing import List, Optional, Tuple

NMAX = 100000
MMAX = 100000


class WeightedDSU:
    """帶權併查集：rel[x] = x 相對於 parent[x] 的關係（mod d）。"""

    def __init__(self, n: int, d: int) -> None:
        self.n = n
        self.d = d
        self.parent = list(range(n + 1))
        self.size = [1] * (n + 1)
        self.rel = [0] * (n + 1)

    def find(self, x: int) -> int:
        # 第一遍：找根
        root = x
        while self.parent[root] != root:
            root = self.parent[root]
        # 第二遍：自頂向下壓縮，權重累加為相對於根
        path: List[int] = []
        cur = x
        while self.parent[cur] != root:
            path.append(cur)
            cur = self.parent[cur]
        for node in reversed(path):
            p = self.parent[node]
            self.rel[node] = (self.rel[node] + self.rel[p]) % self.d
            self.parent[node] = root
        return root

    def union(self, x: int, y: int, r: int) -> None:
        """記錄一條確定為真的陳述 rel(x,y) = r。"""
        rx = self.find(x)
        ry = self.find(y)
        if rx == ry:
            return
        dx = self.rel[x]  # rel(x, rx)
        dy = self.rel[y]  # rel(y, ry)
        if self.size[rx] < self.size[ry]:
            # 把 rx 掛到 ry 下面：rel(rx, ry) = (r - dx + dy) mod d
            self.parent[rx] = ry
            self.rel[rx] = (r - dx + dy) % self.d
            self.size[ry] += self.size[rx]
        else:
            # 把 ry 掛到 rx 下面：rel(ry, rx) = (dx - dy - r) mod d
            self.parent[ry] = rx
            self.rel[ry] = (dx - dy - r) % self.d
            self.size[rx] += self.size[ry]

    def relation(self, x: int, y: int) -> Optional[int]:
        """x, y 已連通時回傳 rel(x,y)，否則回 None。"""
        rx = self.find(x)
        ry = self.find(y)
        if rx != ry:
            return None
        return (self.rel[x] - self.rel[y]) % self.d


class ExtendedDSU:
    """擴展域併查集：每個元素複製成 d 個域節點 (i, t)，i 從 1 開始。

    節點編號：i + t·n（i ∈ [1, n]，t ∈ [0, d)）。
    約定：rel(x, y) = t 已知  ⟺  find(x + t·n) == find(y)。
    """

    def __init__(self, n: int, d: int) -> None:
        self.n = n
        self.d = d
        self.parent = list(range(n * d + 1))
        self.size = [1] * (n * d + 1)

    def find(self, x: int) -> int:
        root = x
        while self.parent[root] != root:
            root = self.parent[root]
        while self.parent[x] != root:
            self.parent[x], x = root, self.parent[x]
        return root

    def union(self, a: int, b: int) -> None:
        ra = self.find(a)
        rb = self.find(b)
        if ra == rb:
            return
        if self.size[ra] < self.size[rb]:
            ra, rb = rb, ra
        self.parent[rb] = ra
        self.size[ra] += self.size[rb]

    def add_statement(self, x: int, y: int, r: int) -> bool:
        """記錄一條陳述；回傳 True 表示為真並已記錄，False 表示是假話。"""
        for t in range(self.d):
            if t != r and self.find(x + t * self.n) == self.find(y):
                return False  # 假話：已知關係是 t 而非 r
        for u in range(self.d):
            self.union(x + u * self.n, y + ((u - r) % self.d) * self.n)
        return True


def solve_weighted(n: int, d: int, stmts: List[Tuple[int, int, int]]) -> int:
    dsu = WeightedDSU(n, d)
    false_cnt = 0
    for x, y, r in stmts:
        known = dsu.relation(x, y)
        if known is None:
            dsu.union(x, y, r)
        elif known != r:
            false_cnt += 1
    return false_cnt


def solve_extended(n: int, d: int, stmts: List[Tuple[int, int, int]]) -> int:
    dsu = ExtendedDSU(n, d)
    false_cnt = 0
    for x, y, r in stmts:
        if not dsu.add_statement(x, y, r):
            false_cnt += 1
    return false_cnt


def solve_brute(n: int, d: int, stmts: List[Tuple[int, int, int]]) -> int:
    """樸素傳遞閉包：R[a][b] = 已知的 rel(a,b) 或 None。"""
    R: List[List[Optional[int]]] = [[None] * (n + 1) for _ in range(n + 1)]
    for i in range(1, n + 1):
        R[i][i] = 0
    false_cnt = 0
    for x, y, r in stmts:
        if R[x][y] is not None:
            if R[x][y] != r:
                false_cnt += 1
            continue
        # 試探性地寫入並取閉包，有矛盾就是假話
        R[x][y] = r
        R[y][x] = (d - r) % d
        if not _close(R, n, d):
            false_cnt += 1
            R[x][y] = None
            R[y][x] = None
    return false_cnt


def _close(R: List[List[Optional[int]]], n: int, d: int) -> bool:
    """對關係矩陣取傳遞閉包；發現矛盾回 False。"""
    changed = True
    while changed:
        changed = False
        for a in range(1, n + 1):
            for b in range(1, n + 1):
                if R[a][b] is None:
                    continue
                for c in range(1, n + 1):
                    if R[b][c] is None:
                        continue
                    v = (R[a][b] + R[b][c]) % d
                    if R[a][c] is None:
                        R[a][c] = v
                        changed = True
                    elif R[a][c] != v:
                        return False
    return True


def solve(n: int, d: int, stmts: List[Tuple[int, int, int]]) -> int:
    return solve_weighted(n, d, stmts)


def run_tests() -> None:
    # ---- 固定用例（手工驗算）----
    # d=3 食物鏈：第 3、5 條是假話
    stmts = [(1, 2, 1), (2, 3, 1), (3, 1, 2), (1, 1, 0), (1, 1, 1), (1, 3, 2)]
    assert solve(3, 3, stmts) == 2
    assert solve_extended(3, 3, stmts) == 2
    assert solve_brute(3, 3, stmts) == 2
    # d=2 奇偶：1-2 異、2-3 異 ⇒ 1-3 同；1-4 異、4-2 異 ⇒ 2-4 同，
    # 但「2 4 異」與之矛盾；最後一條 3-4 異一致
    stmts2 = [(1, 2, 1), (2, 3, 1), (1, 3, 0), (1, 4, 1), (4, 2, 1), (3, 4, 1)]
    assert solve(4, 2, stmts2) == 1
    assert solve_extended(4, 2, stmts2) == 1
    assert solve_brute(4, 2, stmts2) == 1
    # d=2 奇環：三條「異」圍成奇數環，第三條矛盾
    stmts3 = [(1, 2, 1), (2, 3, 1), (1, 3, 1)]
    assert solve(3, 2, stmts3) == 1
    assert solve_extended(3, 2, stmts3) == 1
    assert solve_brute(3, 2, stmts3) == 1
    # x == x 的各種情形
    assert solve(1, 2, [(1, 1, 0), (1, 1, 1)]) == 1
    assert solve(1, 3, [(1, 1, 0), (1, 1, 1), (1, 1, 2)]) == 2
    assert solve_extended(1, 3, [(1, 1, 0), (1, 1, 1), (1, 1, 2)]) == 2
    # 空語句（m = 0 在 IO 模式不合法，此處直接測函數）
    assert solve(5, 3, []) == 0
    # 跨合併才暴露的矛盾：A 吃 B、C 吃 D，再知 B 與 C 同類
    #   ⇒ A 吃 D；此時「A 與 D 同類」為假
    stmts4 = [(1, 2, 1), (3, 4, 1), (2, 3, 0), (1, 4, 0)]
    assert solve(4, 3, stmts4) == 1
    assert solve_extended(4, 3, stmts4) == 1
    assert solve_brute(4, 3, stmts4) == 1

    # ---- 隨機對拍：三套實作結果必須一致 ----
    rng = random.Random(20261026)
    for _ in range(300):
        d = rng.choice([2, 3])
        n = rng.randint(1, 10)
        m = rng.randint(1, 40)
        stmts = [(rng.randint(1, n), rng.randint(1, n), rng.randrange(d))
                 for _ in range(m)]
        a = solve_weighted(n, d, stmts)
        b = solve_extended(n, d, stmts)
        c = solve_brute(n, d, stmts)
        assert a == b == c, (d, n, stmts, a, b, c)

    # ---- 稍大規模：帶權 vs 擴展域（樸素閉包太慢，不參與）----
    rng = random.Random(13572468)
    n = 3000
    m = 6000
    stmts = [(rng.randint(1, n), rng.randint(1, n), rng.randrange(3))
             for _ in range(m)]
    assert solve_weighted(n, 3, stmts) == solve_extended(n, 3, stmts)


# ---------------------------------------------------------------- IO 模式

def _parse_int(tok: str) -> Optional[int]:
    """整數字串解析；越界（超出 ±10^18）或非純整數回 None。"""
    if not tok:
        return None
    s = tok
    i = 0
    if s[0] in "+-":
        i = 1
    if i == len(s):
        return None
    for ch in s[i:]:
        if ch < "0" or ch > "9":
            return None
    val = int(s)
    if val > 10 ** 18:
        return None
    return val


def run_io(raw: str) -> None:
    toks = raw.split()
    pos = 0

    def nxt() -> Optional[str]:
        nonlocal pos
        if pos >= len(toks):
            return None
        t = toks[pos]
        pos += 1
        return t

    def bad() -> None:
        emit("error\n")

    d = _parse_int(nxt() or "")
    if d not in (2, 3):
        bad()
        return
    n = _parse_int(nxt() or "")
    if n is None or not (1 <= n <= NMAX):
        bad()
        return
    m = _parse_int(nxt() or "")
    if m is None or not (1 <= m <= MMAX):
        bad()
        return
    stmts: List[Tuple[int, int, int]] = []
    for _ in range(m):
        x = _parse_int(nxt() or "")
        y = _parse_int(nxt() or "")
        r = _parse_int(nxt() or "")
        if (x is None or y is None or r is None
                or not (1 <= x <= n) or not (1 <= y <= n)
                or not (0 <= r < d)):
            bad()
            return
        stmts.append((x, y, r))
    emit(str(solve(n, d, stmts)) + "\n")


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
