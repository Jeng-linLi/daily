"""伸展樹（Splay Tree）：自調整二元搜尋樹

題意：
    維護一個可重複元素的有序集合，支援：
        insert(x)      插入 x
        erase(x)       刪除一個 x（不存在則無動作）
        contains(x)    是否存在
        kth(k)         第 k 小（1-indexed）
        rank(x)        嚴格小於 x 的元素個數
        predecessor(x) 嚴格小於 x 的最大值
        successor(x)   嚴格大於 x 的最小值
    所有操作都靠**伸展（splay）**把被訪問的節點一路旋轉到根。
    另外用一個「有序陣列」參考實作（bisect）逐項對拍，驗證 splay 的結果正確。

思路：
    ### 核心想法：把剛訪問過的節點搬到根
    伸展樹不維護任何平衡資訊（沒有 AVL 的高度、沒有紅黑樹的顏色），
    它只遵守一條規則：**每次訪問一個節點，就把它旋轉到根**。
    靠「均攤分析」（access lemma / 勢能法）可證：
    m 次操作的總代價是 O((n + m) log n)，也就是**均攤 O(log n)**。
    代價是單次操作可能退化成 O(n)，好處是常數小、程式短，
    而且自帶「快取效應」：最近被訪問的元素永遠在根附近，
    對「有局部性（locality）的存取序列」特別快——這是它真正的價值。

    ### 三種旋轉情況
    設 x 是要伸展的節點，f = fa[x]，g = fa[f]：
      - **zig**（g 為空）：只轉一次 x。
      - **zig-zig**（x 與 f 同側，都是左兒子或都是右兒子）：**先轉 f 再轉 x**。
      - **zig-zag**（x 與 f 異側）：**轉兩次 x**。
    zig-zig 與 zig-zag 的差別是伸展樹的精髓：
    zig-zig 先轉父親會把整條路徑「折半」，把樹壓扁；
    如果這裡也寫成「轉兩次 x」，均攤複雜度就會退化（變成 Ω(n) 的競爭比）。

    ### 刪除：伸展後靠前驅頂替
    先把目標 splay 到根：
      - 計數 > 1：直接減一（本題用 cnt 支援重複元素）；
      - 沒有兒子 / 只有一個兒子：讓那個兒子當根；
      - 兩個兒子都在：取左子樹的最大值（前驅），把它在左子樹裡 splay 到根，
        再把右子樹接成它的右兒子。這跟 BST 的「中序後繼頂替」是同一套思路，
        差別在於 splay 版順便把前驅送到根，讓下次存取更快。

    ### 與 Treap / 紅黑樹的比較
      - Treap 靠隨機優先級保證期望平衡，單次操作嚴格 O(log n)；
      - 紅黑樹靠顏色規則保證最壞 O(log n)，但常數大、程式長；
      - Splay 不保證單次操作，但均攤 O(log n) 且有局部性加速，
        是「序列維護」（區間翻轉 / 區間移動）與 LCT 的基礎結構。

應用場景：
    - 需要平衡樹但存取有局部性：編譯器的符號表、垃圾回收的物件快取
    - 區間維護（文藝平衡樹）：隱式 Splay 用中序位置當下標
    - Link-Cut Tree 的輔助樹就是 Splay
    - 競賽：普通平衡樹（插入 / 刪除 / 排名 / 前驅後繼）、序列翻轉

複雜度（均攤）：插入 / 刪除 / 查詢 / 第 k 小 / 排名 O(log n)；
                單次操作最壞 O(n)；空間 O(n)。
"""

import bisect
import re
import sys
from typing import List, Optional

NUM_RE = re.compile(r"^[+-]?[0-9]+$")

MAXN = 3000             # n 與 q 的截斷上限
VMAX = 10 ** 12         # 數值截斷上限

NULL = 0                # 0 號是空節點哨兵，真實節點從 1 開始編號


class Splay:
    """陣列池實作的伸展樹（支援重複元素，cnt 為該 key 的個數）。"""

    __slots__ = ("ch", "fa", "key", "cnt", "sz", "root", "tot")

    def __init__(self, capacity: int) -> None:
        cap = capacity + 5
        self.ch = [[0, 0] for _ in range(cap)]
        self.fa = [0] * cap
        self.key = [0] * cap
        self.cnt = [0] * cap
        self.sz = [0] * cap
        self.root = 0
        self.tot = 0

    # ---- 基礎工具 ----

    def _new(self, key: int) -> int:
        self.tot += 1
        i = self.tot
        self.ch[i][0] = self.ch[i][1] = 0
        self.fa[i] = 0
        self.key[i] = key
        self.cnt[i] = 1
        self.sz[i] = 1
        return i

    def _pushup(self, x: int) -> None:
        self.sz[x] = self.sz[self.ch[x][0]] + self.sz[self.ch[x][1]] + self.cnt[x]

    def _dir(self, f: int, x: int) -> int:
        return 1 if self.ch[f][1] == x else 0

    def _rotate(self, x: int) -> None:
        f = self.fa[x]
        g = self.fa[f]
        d = self._dir(f, x)
        b = self.ch[x][1 - d]
        self.ch[f][d] = b
        if b:
            self.fa[b] = f
        self.ch[x][1 - d] = f
        self.fa[f] = x
        self.fa[x] = g
        if g:
            self.ch[g][1 if self.ch[g][1] == f else 0] = x
        else:
            self.root = x
        self._pushup(f)
        self._pushup(x)

    def _splay(self, x: int) -> None:
        """把 x 旋轉到根（三種情況：zig / zig-zig / zig-zag）。"""
        while self.fa[x]:
            f = self.fa[x]
            g = self.fa[f]
            if g == 0:
                self._rotate(x)                      # zig
            elif self._dir(g, f) == self._dir(f, x):
                self._rotate(f)                      # zig-zig：先轉父親
                self._rotate(x)
            else:
                self._rotate(x)                      # zig-zag：轉兩次自己
                self._rotate(x)
        self.root = x

    # ---- 對外操作 ----

    def size(self) -> int:
        return self.sz[self.root]

    def insert(self, key: int) -> None:
        if self.root == 0:
            self.root = self._new(key)
            return
        x = self.root
        while True:
            if key == self.key[x]:
                self.cnt[x] += 1
                self._pushup(x)
                self._splay(x)
                return
            d = 1 if key > self.key[x] else 0
            if self.ch[x][d] == 0:
                i = self._new(key)
                self.ch[x][d] = i
                self.fa[i] = x
                self._pushup(x)
                self._splay(i)
                return
            x = self.ch[x][d]

    def contains(self, key: int) -> int:
        x = self.root
        last = 0
        while x:
            last = x
            if key == self.key[x]:
                self._splay(x)
                return 1
            x = self.ch[x][1] if key > self.key[x] else self.ch[x][0]
        if last:
            self._splay(last)
        return 0

    def erase(self, key: int) -> None:
        x = self.root
        last = 0
        while x:
            last = x
            if key == self.key[x]:
                break
            x = self.ch[x][1] if key > self.key[x] else self.ch[x][0]
        if x == 0:
            if last:
                self._splay(last)
            return
        self._splay(x)
        if self.cnt[x] > 1:
            self.cnt[x] -= 1
            self._pushup(x)
            return
        l, r = self.ch[x][0], self.ch[x][1]
        if l == 0:
            self.root = r
            self.fa[r] = 0
        elif r == 0:
            self.root = l
            self.fa[l] = 0
        else:
            # 取左子樹的最大值（前驅），先在左子樹裡把它 splay 到根
            self.fa[l] = 0
            self.fa[r] = 0
            y = l
            while self.ch[y][1]:
                y = self.ch[y][1]
            self._splay(y)
            self.ch[y][1] = r
            self.fa[r] = y
            self._pushup(y)
            self.root = y

    def kth(self, k: int) -> Optional[int]:
        if k < 1 or k > self.sz[self.root]:
            return None
        x = self.root
        last = 0
        while x:
            last = x
            ls = self.sz[self.ch[x][0]]
            if k <= ls:
                x = self.ch[x][0]
            elif k > ls + self.cnt[x]:
                k -= ls + self.cnt[x]
                x = self.ch[x][1]
            else:
                self._splay(x)
                return self.key[x]
        if last:
            self._splay(last)
        return None

    def rank(self, key: int) -> int:
        """嚴格小於 key 的元素個數。"""
        x = self.root
        last = 0
        ans = 0
        while x:
            last = x
            if key < self.key[x]:
                x = self.ch[x][0]
            elif key > self.key[x]:
                ans += self.sz[self.ch[x][0]] + self.cnt[x]
                x = self.ch[x][1]
            else:
                ans += self.sz[self.ch[x][0]]
                self._splay(x)
                return ans
        if last:
            self._splay(last)
        return ans

    def predecessor(self, key: int) -> Optional[int]:
        x = self.root
        last = 0
        ans = None
        while x:
            last = x
            if self.key[x] < key:
                ans = self.key[x]
                x = self.ch[x][1]
            else:
                x = self.ch[x][0]
        if last:
            self._splay(last)
        return ans

    def successor(self, key: int) -> Optional[int]:
        x = self.root
        last = 0
        ans = None
        while x:
            last = x
            if self.key[x] > key:
                ans = self.key[x]
                x = self.ch[x][0]
            else:
                x = self.ch[x][1]
        if last:
            self._splay(last)
        return ans

    def inorder(self) -> List[int]:
        """中序遍歷（含重複元素），迭代版避免深度過大。"""
        res: List[int] = []
        st: List[int] = []
        x = self.root
        while x or st:
            while x:
                st.append(x)
                x = self.ch[x][0]
            x = st.pop()
            res.extend([self.key[x]] * self.cnt[x])
            x = self.ch[x][1]
        return res


# ---------------------------------------------------------------- 參考實作


class RefSet:
    """有序陣列版的有序集合（bisect），用來對拍 splay 的結果。"""

    __slots__ = ("a",)

    def __init__(self) -> None:
        self.a: List[int] = []

    def insert(self, key: int) -> None:
        bisect.insort(self.a, key)

    def erase(self, key: int) -> None:
        i = bisect.bisect_left(self.a, key)
        if i < len(self.a) and self.a[i] == key:
            del self.a[i]

    def contains(self, key: int) -> int:
        i = bisect.bisect_left(self.a, key)
        return 1 if i < len(self.a) and self.a[i] == key else 0

    def kth(self, k: int) -> Optional[int]:
        return self.a[k - 1] if 1 <= k <= len(self.a) else None

    def rank(self, key: int) -> int:
        return bisect.bisect_left(self.a, key)

    def predecessor(self, key: int) -> Optional[int]:
        i = bisect.bisect_left(self.a, key)
        return self.a[i - 1] if i > 0 else None

    def successor(self, key: int) -> Optional[int]:
        i = bisect.bisect_right(self.a, key)
        return self.a[i] if i < len(self.a) else None


# ---------------------------------------------------------------- 驅動


def fmt(v: Optional[int]) -> str:
    return "-1" if v is None else str(v)


def run(init: List[int], ops: List[tuple]) -> str:
    cap = len(init) + len(ops) + 10
    tree = Splay(cap)
    ref = RefSet()
    for v in init:
        tree.insert(v)
        ref.insert(v)

    out: List[str] = []
    ok = 1
    for op, a in ops:
        if op == 1:
            tree.insert(a)
            ref.insert(a)
        elif op == 2:
            tree.erase(a)
            ref.erase(a)
        elif op == 3:
            v = tree.contains(a)
            out.append(str(v))
            if v != ref.contains(a):
                ok = 0
        elif op == 4:
            v = tree.kth(a)
            out.append(fmt(v))
            if v != ref.kth(a):
                ok = 0
        elif op == 5:
            v = tree.rank(a)
            out.append(str(v))
            if v != ref.rank(a):
                ok = 0
        elif op == 6:
            v = tree.predecessor(a)
            out.append(fmt(v))
            if v != ref.predecessor(a):
                ok = 0
        elif op == 7:
            v = tree.successor(a)
            out.append(fmt(v))
            if v != ref.successor(a):
                ok = 0

    order = tree.inorder()
    if tree.size() != len(ref.a) or order != ref.a:
        ok = 0
    # BST 性質檢查：中序遍歷必須非遞減（同一 key 的重複元素相鄰）
    for i in range(1, len(order)):
        if order[i] < order[i - 1]:
            ok = 0
            break

    return "\n".join([
        str(tree.size()),
        " ".join(str(v) for v in order),
        " ".join(out),
        str(ok),
    ])


# ---------------------------------------------------------------- 內置測試


def run_tests() -> None:
    import random

    # ---- 固定用例 ----
    assert run([], []) == "0\n\n\n1"
    assert run([5, 3, 7], []) == "3\n3 5 7\n\n1"
    assert run([], [(1, 5), (1, 3), (1, 7)]) == "3\n3 5 7\n\n1"
    assert run([1, 2, 3], [(3, 2), (3, 9), (2, 2), (3, 2)]) == "2\n1 3\n1 0 0\n1"
    assert run([1, 2, 3], [(4, 1), (4, 3), (4, 4), (4, 0)]) == "3\n1 2 3\n1 3 -1 -1\n1"
    assert run([1, 2, 3], [(5, 1), (5, 2), (5, 3), (5, 0), (5, 9)]) == "3\n1 2 3\n0 1 2 0 3\n1"
    assert run([10, 20, 30], [(6, 25), (7, 25), (6, 10), (7, 30)]) == "3\n10 20 30\n20 30 -1 -1\n1"
    # 重複元素：insert 三次 5，size 應該是 3，中序有兩個 5
    assert run([], [(1, 5), (1, 5), (1, 5), (2, 5)]) == "2\n5 5\n\n1"

    # ---- 隨機對拍 ----
    rng = random.Random(20261015)
    for _ in range(400):
        n = rng.randint(0, 40)
        init = [rng.randint(-20, 20) for _ in range(n)]
        q = rng.randint(0, 80)
        ops = [(rng.randint(1, 7), rng.randint(-25, 25)) for _ in range(q)]
        assert run(init, ops).endswith("\n1")

    # ---- 大序列壓測（只跑 splay，不跑參考實作的完整對拍輸出比對）----
    big = Splay(2 * MAXN + 20)
    vals = [rng.randint(-10 ** 9, 10 ** 9) for _ in range(MAXN)]
    for v in vals:
        big.insert(v)
    assert big.size() == MAXN
    assert big.inorder() == sorted(vals)
    for v in vals:
        assert big.contains(v) == 1
    for v in vals:
        big.erase(v)
    assert big.size() == 0


# ---------------------------------------------------------------- IO


def emit(text: str) -> None:
    """二進位輸出（latin-1），避免 Windows 把 \\n 換成 \\r\\n。"""
    sys.stdout.buffer.write(text.encode("latin-1"))
    sys.stdout.buffer.flush()


def run_io(raw: str) -> None:
    """stdin 解析：嚴格只吃 `^[+-]?[0-9]+$` 的 token，非數字一律當作輸入結束。"""
    tokens = raw.split()
    pos = 0

    def read_int() -> Optional[int]:
        nonlocal pos
        if pos < len(tokens) and NUM_RE.match(tokens[pos]):
            v = int(tokens[pos])
            pos += 1
            return v
        return None

    def clamp(v: Optional[int], lo: int, hi: int, default: int = 0) -> int:
        if v is None:
            return default
        return max(lo, min(hi, v))

    n = clamp(read_int(), 0, MAXN, 0)
    init: List[int] = []
    for _ in range(n):
        v = read_int()
        if v is None:
            break
        init.append(clamp(v, -VMAX, VMAX))

    q = clamp(read_int(), 0, MAXN, 0)
    ops: List[tuple] = []
    for _ in range(q):
        o = read_int()
        a = read_int()
        if o is None or a is None:
            break
        ops.append((clamp(o, 0, 7, 0), clamp(a, -VMAX, VMAX)))

    emit(run(init, ops) + "\n")


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
