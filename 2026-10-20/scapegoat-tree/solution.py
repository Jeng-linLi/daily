"""替罪羊樹（Scapegoat Tree）：不旋轉，靠「整棵子樹重建」維持平衡的權重平衡樹

題意：
    維護一個**有序集合**（元素唯一），支援插入 / 刪除 / 查找 / 排名 / 第 k 小 /
    前驅 / 後繼，並與「有序陣列 + 二分」的參考實作逐項對拍。

    操作碼（每行 `op x`）：
        1  insert x    插入 x（已存在則忽略）
        2  erase  x    刪除 x（不存在則忽略）
        3  count  x    x 是否在集合中（1 / 0）
        4  rank   x    嚴格小於 x 的元素個數
        5  kth    k    第 k 小的元素（k 從 1 起算），越界回傳 -1
        6  pred   x    嚴格小於 x 的最大元素，不存在回傳 -1
        7  succ   x    嚴格大於 x 的最小元素，不存在回傳 -1

思路：
    ### 1. 什麼是「權重平衡」
    紅黑樹 / AVL / Treap / Splay 都靠**旋轉**（rotation）局部調整形狀。
    替罪羊樹完全不旋轉，改成一個更粗暴的條件：

        對每個節點 u，其左右子樹大小都必須 ≤ α · size(u)   （本題 α = 7/10）

    一旦某個節點違反這個條件，就把它的**整棵子樹**拉直（中序遍歷成有序陣列），
    再從中點遞迴建成一棵完美平衡樹。這個動作叫 **rebuild（重建）**。

    ### 2. 為什麼這樣還能是 O(log n)
    - α 權重平衡 ⟹ 從根往下每走一層，子樹大小至少縮為 α 倍
      ⟹ 樹高 ≤ log_{1/α}(n) = ln n / ln(10/7) ≈ 2.81 · log₂ n。
    - 重建的攤還代價：一個節點被重建，代表它底下剛累積了 Θ(size) 次插入，
      於是把 O(size) 的重建成本分攤到這 Θ(size) 次插入上，每次插入攤還 O(log n)。

    ### 3. 插入後要重建「哪一個」節點（本題的關鍵設計）
    插入後從根沿著 key 走一趟，找**最淺**（離根最近）的那個失衡節點 v 重建。
    重建完之後，subtree(v) 內部全部完美平衡，而 v 的祖先們大小完全沒變
    （重建不改變節點數），所以「下一次失衡的位置」只能嚴格地往根的方向移動。
    於是把這個動作寫成迴圈：

        while 存在失衡節點 v（取最淺的）:  rebuild(v)

    每一輪 v 的深度嚴格遞減，最多走 depth 輪就結束；一旦 v 就是根，
    重建後整棵樹完美平衡，必然收斂。這比「找最深的 / 只重建一次」更穩健——
    因為插入會讓祖先的 size 與子樹 size 同時 +1，而 α < 1，
    `child+1 ≤ α·(size+1)` 是有可能被打破的（1 > α），只修一次不保證全域平衡。

    ### 4. 刪除：靠「縮水太多就整棵重建」
    刪除走一般 BST 刪除（兩個小孩時用**中序前驅**頂替），沿途把 size -1。
    維護 `max_size`（歷史最大節點數）：當 `size · 10 < max_size · 7` 時，
    表示節點已經掉到歷史高峰的 70% 以下，整棵樹重建並把 max_size 設回 size。
    這樣攤還起來每次刪除也是 O(log n)。

    ### 5. 兩語言一致性的兩個陷阱
    - 平衡條件用**整數**比較 `max(l, r) · 10 > size · 7`，不要用浮點 α。
    - 建平衡樹取中點用 `mid = (l + r) // 2`（下標非負，Python 的 floor 與
      C++ 的 truncate 結果相同），不要寫 `(l + r) / 2`。

應用場景：
    - 需要「有序集合」但又想避開旋轉邏輯的場合（程式碼量比紅黑樹少一個數量級）
    - 函數式 / 持久化結構：重建天然產生新節點，不需要可變指標
    - 資料庫索引中「批次重建比逐次旋轉更快」的場景（LSM 的 read-only segment 合併）
    - 競賽中 Treap / Splay 的替代：特別是需要「按子樹大小攤平」的題
    - 不希望旋轉造成迭代器 / 外部指標失效的嵌入式與即時系統

複雜度（n = 元素個數，α = 0.7）：
    查找 / rank / kth / pred / succ   O(log n) 最壞
    插入                             攤還 O(log n)，單次最壞 O(n)（重建）
    刪除                             攤還 O(log n)，單次最壞 O(n)
    樹高                             ≤ log_{1/α}(n) + O(1) ≈ 2.81 log₂ n
    空間                             O(n)
"""

import random
import re
import sys
from bisect import bisect_left, bisect_right, insort
from typing import List, Optional, Tuple

INT_RE = re.compile(r"^[+-]?[0-9]+$")

QMAX = 300          # IO 模式下的操作數上限
VMAX = 10 ** 6      # 元素值上限（夾到 [0, VMAX]，讓 -1 可以當「不存在」哨兵）

ALPHA_NUM = 7       # α = 7/10
ALPHA_DEN = 10


# ---------------------------------------------------------------- 節點與小工具


class Node:
    """替罪羊樹節點：key + 左右小孩 + 子樹大小。"""

    __slots__ = ("key", "left", "right", "size")

    def __init__(self, key: int) -> None:
        self.key = key
        self.left: Optional["Node"] = None
        self.right: Optional["Node"] = None
        self.size = 1


def _sz(t: Optional[Node]) -> int:
    return t.size if t is not None else 0


def _pull(t: Node) -> None:
    t.size = 1 + _sz(t.left) + _sz(t.right)


def _collect(t: Optional[Node], out: List[int]) -> None:
    """中序遍歷，把子樹拉直成有序陣列。"""
    if t is None:
        return
    _collect(t.left, out)
    out.append(t.key)
    _collect(t.right, out)


def _build(keys: List[int], l: int, r: int) -> Optional[Node]:
    """把有序陣列 keys[l:r] 建成一棵高度平衡的樹。"""
    if l >= r:
        return None
    mid = (l + r) // 2
    t = Node(keys[mid])
    t.left = _build(keys, l, mid)
    t.right = _build(keys, mid + 1, r)
    _pull(t)
    return t


def _height(t: Optional[Node]) -> int:
    if t is None:
        return 0
    return 1 + max(_height(t.left), _height(t.right))


def _true_size(t: Optional[Node]) -> int:
    """實際走一趟算出來的子樹大小，用來驗證 size 欄位沒有算錯。"""
    if t is None:
        return 0
    return 1 + _true_size(t.left) + _true_size(t.right)


def _check_sizes(t: Optional[Node]) -> bool:
    if t is None:
        return True
    if t.size != 1 + _sz(t.left) + _sz(t.right):
        return False
    return _check_sizes(t.left) and _check_sizes(t.right)


def _unbalanced(t: Node) -> bool:
    """α 權重平衡條件（α = 7/10）。整數運算，兩語言完全一致。"""
    return max(_sz(t.left), _sz(t.right)) * ALPHA_DEN > t.size * ALPHA_NUM


# ---------------------------------------------------------------- 替罪羊樹


class ScapegoatTree:
    def __init__(self) -> None:
        self.root: Optional[Node] = None
        self.max_size = 0

    # -- 讀取 --

    def size(self) -> int:
        return _sz(self.root)

    def height(self) -> int:
        return _height(self.root)

    def inorder(self) -> List[int]:
        out: List[int] = []
        _collect(self.root, out)
        return out

    def contains(self, key: int) -> bool:
        t = self.root
        while t is not None:
            if key == t.key:
                return True
            t = t.left if key < t.key else t.right
        return False

    def count_less(self, key: int) -> int:
        cnt = 0
        t = self.root
        while t is not None:
            if key <= t.key:
                t = t.left
            else:
                cnt += 1 + _sz(t.left)
                t = t.right
        return cnt

    def kth(self, k: int) -> int:
        """第 k 小（k 從 1 起算）。"""
        if k < 1 or k > _sz(self.root):
            return -1
        t = self.root
        while t is not None:
            ls = _sz(t.left)
            if k <= ls:
                t = t.left
            elif k == ls + 1:
                return t.key
            else:
                k -= ls + 1
                t = t.right
        return -1

    def pred(self, key: int) -> int:
        best = -1
        t = self.root
        while t is not None:
            if t.key < key:
                best = t.key
                t = t.right
            else:
                t = t.left
        return best

    def succ(self, key: int) -> int:
        best = -1
        t = self.root
        while t is not None:
            if t.key > key:
                best = t.key
                t = t.left
            else:
                t = t.right
        return best

    # -- 修改 --

    def insert(self, key: int) -> bool:
        if self.root is None:
            self.root = Node(key)
            self.max_size = 1
            return True
        path: List[Node] = []
        t: Optional[Node] = self.root
        while t is not None:
            if key == t.key:
                return False
            path.append(t)
            t = t.left if key < t.key else t.right
        parent = path[-1]
        if key < parent.key:
            parent.left = Node(key)
        else:
            parent.right = Node(key)
        for u in path:
            u.size += 1
        if self.root.size > self.max_size:
            self.max_size = self.root.size
        self._fix(key)
        return True

    def erase(self, key: int) -> bool:
        path: List[Node] = []
        t: Optional[Node] = self.root
        while t is not None and t.key != key:
            path.append(t)
            t = t.left if key < t.key else t.right
        if t is None:
            return False
        if t.left is not None and t.right is not None:
            # 兩個小孩：用中序前驅（左子樹最大值）頂替
            ppath: List[Node] = [t]
            p: Node = t.left
            while p.right is not None:
                ppath.append(p)
                p = p.right
            t.key = p.key
            parent = ppath[-1]
            if parent.left is p:
                parent.left = p.left
            else:
                parent.right = p.left
            # 注意：t 的祖先（path）也要 -1，ppath 只包含 t 往下到 p 的父節點
            for u in path:
                u.size -= 1
            for u in ppath:
                u.size -= 1
        else:
            child = t.left if t.left is not None else t.right
            if path:
                parent = path[-1]
                if parent.left is t:
                    parent.left = child
                else:
                    parent.right = child
            else:
                self.root = child
            for u in path:
                u.size -= 1
        if self.root is not None and self.root.size * ALPHA_DEN < self.max_size * ALPHA_NUM:
            keys: List[int] = []
            _collect(self.root, keys)
            self.root = _build(keys, 0, len(keys))
            self.max_size = self.root.size
        elif self.root is None:
            self.max_size = 0
        return True

    # -- 內部：插入後反覆重建最淺的失衡節點 --

    def _fix(self, key: int) -> None:
        while True:
            path: List[Node] = []
            t: Optional[Node] = self.root
            while t is not None:
                path.append(t)
                if key == t.key:
                    break
                t = t.left if key < t.key else t.right
            bad = -1
            for i, u in enumerate(path):
                if _unbalanced(u):
                    bad = i
                    break
            if bad < 0:
                return
            v = path[bad]
            keys: List[int] = []
            _collect(v, keys)
            nb = _build(keys, 0, len(keys))
            if bad == 0:
                self.root = nb
            else:
                p = path[bad - 1]
                if p.left is v:
                    p.left = nb
                else:
                    p.right = nb


# ---------------------------------------------------------------- 參考實作對拍


def run(ops: List[Tuple[int, int]]) -> str:
    """跑一串操作：替罪羊樹 vs 有序陣列參考實作，回傳三行文字。"""
    tree = ScapegoatTree()
    ref: List[int] = []
    qres: List[str] = []

    for op, x in ops:
        if op == 1:
            a = 1 if tree.insert(x) else 0
            b = 0
            if all(v != x for v in ref):
                insort(ref, x)
                b = 1
            assert a == b, ("insert", op, x)
        elif op == 2:
            a = 1 if tree.erase(x) else 0
            b = 0
            i = bisect_left(ref, x)
            if i < len(ref) and ref[i] == x:
                ref.pop(i)
                b = 1
            assert a == b, ("erase", op, x)
        elif op == 3:
            qres.append("1" if tree.contains(x) else "0")
            i = bisect_left(ref, x)
            assert tree.contains(x) == (i < len(ref) and ref[i] == x)
        elif op == 4:
            qres.append(str(tree.count_less(x)))
            assert tree.count_less(x) == bisect_left(ref, x)
        elif op == 5:
            qres.append(str(tree.kth(x)))
            want = ref[x - 1] if 1 <= x <= len(ref) else -1
            assert tree.kth(x) == want, ("kth", x, tree.kth(x), want)
        elif op == 6:
            qres.append(str(tree.pred(x)))
            i = bisect_left(ref, x)
            want = ref[i - 1] if i > 0 else -1
            assert tree.pred(x) == want, ("pred", x)
        elif op == 7:
            qres.append(str(tree.succ(x)))
            i = bisect_right(ref, x)          # 嚴格大於：用右插入點
            want = ref[i] if i < len(ref) else -1
            assert tree.succ(x) == want, ("succ", x)

    # size 欄位的一致性檢查：中序內容對但 size 錯，kth / rank 仍會給錯答案
    assert _check_sizes(tree.root), "size field inconsistent"
    assert _true_size(tree.root) == tree.size()
    ok = 1 if tree.inorder() == ref else 0
    line1 = " ".join(qres)
    line2 = " ".join([str(tree.size()), str(tree.height())] + [str(v) for v in tree.inorder()])
    return line1 + "\n" + line2 + "\n" + str(ok) + "\n"


# ---------------------------------------------------------------- 內置測試


def run_tests() -> None:
    # ---- 固定用例：基本插入與查詢 ----
    out = run([(1, 5), (1, 3), (1, 8), (1, 1), (1, 5), (3, 3), (3, 4),
               (4, 5), (5, 2), (6, 5), (7, 5)])
    assert out.split("\n")[0] == "1 0 2 3 3 8", out
    assert out.split("\n")[1] == "4 3 1 3 5 8", out
    assert out.split("\n")[2] == "1", out

    # ---- 空集合 ----
    out = run([(5, 1), (6, 0), (7, 0), (3, 0), (4, 0), (2, 0)])
    assert out.split("\n")[0] == "-1 -1 -1 0 0", out
    assert out.split("\n")[1] == "0 0", out
    assert out.split("\n")[2] == "1", out

    # ---- 樹高：連續插入 1..2000（已排序輸入，對 BST 是最壞情況）----
    tree = ScapegoatTree()
    for i in range(1, 2001):
        tree.insert(i)
    assert tree.size() == 2000
    h = tree.height()
    assert h <= 30, h           # 理論上界 ≈ 2.81·log₂2000 ≈ 31，實測遠低於此
    assert tree.inorder() == list(range(1, 2001))

    # ---- 樹高：隨機插入 2000 個 ----
    rng0 = random.Random(1234)
    tree2 = ScapegoatTree()
    vals = [rng0.randint(0, 10 ** 6) for _ in range(2000)]
    for v in vals:
        tree2.insert(v)
    assert tree2.inorder() == sorted(set(vals))
    assert tree2.height() <= 30, tree2.height()

    # ---- 隨機對拍 ----
    rng = random.Random(20261020)
    worst_h = 0
    for _ in range(700):
        q = rng.randint(0, 40)
        ops = []
        for _ in range(q):
            ops.append((rng.randint(1, 7), rng.randint(0, 12)))
        out = run(ops)
        assert out.split("\n")[2] == "1", (ops, out)
        h = int(out.split("\n")[1].split()[1])
        n = int(out.split("\n")[1].split()[0])
        if n > 0:
            worst_h = max(worst_h, h - 2.81 * max(1.0, n.bit_length()))

    # ---- 大量插入後隨機查詢（順便檢查樹高）----
    for trial in range(60):
        rng2 = random.Random(9000 + trial)
        ops = [(1, rng2.randint(0, 200)) for _ in range(rng2.randint(0, 150))]
        for _ in range(30):
            ops.append((rng2.randint(2, 7), rng2.randint(0, 210)))
        out = run(ops)
        assert out.split("\n")[2] == "1", (trial, out)
        parts = out.split("\n")[1].split()
        n, h = int(parts[0]), int(parts[1])
        if n >= 2:
            assert h <= 3 * n.bit_length() + 6, (trial, n, h)

    # ---- 刪除壓力測試：插入 800 個再全部刪掉 ----
    rng3 = random.Random(555)
    vals3 = [rng3.randint(0, 500) for _ in range(800)]
    ops = [(1, v) for v in vals3] + [(2, v) for v in vals3]
    out = run(ops)
    assert out.split("\n")[1] == "0 0", out
    assert out.split("\n")[2] == "1", out


# ---------------------------------------------------------------- IO


def emit(text: str) -> None:
    """二進位輸出，避免 Windows 把 \\n 換成 \\r\\n。"""
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
        return max(-(10 ** 18), min(10 ** 18, v))

    q = next_int(0)
    if q < 0:
        q = 0
    if q > QMAX:
        q = QMAX
    ops: List[Tuple[int, int]] = []
    for _ in range(q):
        op = next_int(0)
        x = next_int(0)
        if op < 1 or op > 7:
            op = 0
        if x < 0:
            x = 0
        if x > VMAX:
            x = VMAX
        ops.append((op, x))
    emit(run(ops))


def main() -> None:
    raw = sys.stdin.buffer.read().decode("latin-1")
    if raw.strip():
        run_io(raw)
    else:
        run_tests()
        emit("all tests passed\n")


if __name__ == "__main__":
    main()
