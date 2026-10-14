"""隨機 Treap：旋轉版 / FHQ（split-merge）版 / 隱式 Treap 區間翻轉

題意：
    維護一棵以 key 為序的平衡二元搜尋樹，支援插入、刪除、查詢存在性、
    第 k 小、排名（有多少元素嚴格小於 x）、中序遍歷；
    另外用**隱式 Treap** 維護一個序列，支援區間翻轉。
    同時用「旋轉版」與「FHQ 版」兩套實作互相對拍，並驗證堆性質與 BST 性質。

思路：
    ### Treap = Binary Search Tree + Heap
    Treap 的每個節點有兩個權值：
      - `key`  滿足 **BST 性質**：左子樹 < 自己 ≤ 右子樹；
      - `prio` 滿足 **堆性質**（本題用小根堆）：父節點的 prio 小於兩個兒子。
    若 prio 是隨機獨立均勻分布的，則「key 決定的 BST 形狀」恰好等價於
    「按 prio 從小到大插入的順序」，也就是一棵**隨機 BST**，期望深度 O(log n)。
    這就是 Treap 不需要像 AVL / 紅黑樹那樣維護平衡資訊的原因：
    **隨機性本身就是平衡。**

    ### 旋轉版：插入後靠旋轉把 prio 提上去
    先按普通 BST 規則遞迴插到底，回程時若「兒子的 prio 比自己小」就旋轉：
      - 左兒子 prio 更小 → 右旋（zig）；
      - 右兒子 prio 更小 → 左旋（zag）。
    旋轉只改變父子關係、不破壞 BST 的中序序列，因此是合法的平衡化手段。
    刪除則是「把目標節點旋轉到葉子再丟棄」：
      若只有一個兒子，直接用兒子頂替；
      若兩個兒子都在，把 prio 較小的那個兒子旋轉上來，再往另一側繼續刪。

    ### FHQ Treap：只用 split 與 merge
    FHQ（范浩強）Treap 完全不寫旋轉，只靠兩個操作：
      - `merge(a, b)`：合併兩棵樹（要求 a 的所有 key ≤ b 的所有 key），
        prio 小者當根，遞迴合併另一側；
      - `split(t, key)`：按 key 切成 (≤ key, > key) 兩棵樹。
    於是：
        insert(x)  = split(t, x) → merge(A, new) → merge(..., B)
        erase(x)   = split(t, x) → split(A, x-1) → 丟掉中間那一棵
    FHQ 的優勢是**結構持久化友善、區間操作直觀**，而且天然支援隱式 Treap。

    ### 隱式 Treap：用「下標」當 key
    把序列的第 i 個元素放在「中序遍歷第 i 個位置」的節點上，
    不需要真的存 key，只需要維護子樹大小 `sz`。
      - `split_size(t, k)` 切成「前 k 個」與「剩下的」；
      - 區間翻轉 [l, r] = split 出三段 → 給中段打懶標記 → 合併回去。
    懶標記 `rev` 用 XOR 下傳（`push`）：打兩次等於沒打，這正是翻轉的語義。

    ### 為什麼兩個語言能逐位元組一致
    Treap 的形狀完全由 prio 序列決定。本題用 **xorshift64** 確定性亂數
    （種子固定），並且「隱式 Treap / FHQ / 旋轉版」各自用**獨立的亂數流**，
    因此 Python 與 C++ 構造出的樹**完全同構**，輸出必然一致。

應用場景：
    - 需要平衡樹但不想寫紅黑樹：資料庫索引、記憶體 KV 儲存的排序結構
    - 隱式 Treap 是「可持久化序列」與「區間操作」的基礎：文字編輯器的行緩衝、
      區間翻轉 / 區間移動（rope 結構的低成本替代）
    - split / merge 的語義讓 Treap 很容易改成 **持久化 Treap**（函數式資料結構）
    - 競賽：普通平衡樹（插入 / 刪除 / 排名 / 前驅後繼）、序列維護（文藝平衡樹）

複雜度（n 為元素個數，期望值）：
    插入 / 刪除 / 查詢 / 第 k 小 / 排名   O(log n)
    隱式 Treap 區間翻轉                   O(log n)（懶標記 O(1)）
    中序遍歷                              O(n)
    空間                                  O(n)

輸入格式（stdin，全部以空白分隔的整數）：
    n                      初始序列長度（隱式 Treap 的內容）
    v_1 ... v_n            初始序列
    q                      操作數
    接下來 q 組，每組 3 個整數 `op a b`：
        op = 1  a=x  → insert(x)
        op = 2  a=x  → erase(x)（若不存在則無動作）
        op = 3  a=x  → contains(x)，輸出 1 / 0
        op = 4  a=k  → kth(k)（1-indexed），輸出值；越界輸出 -1
        op = 5  a=x  → rank(x)，輸出「嚴格小於 x 的元素個數」
        op = 6  a=l b=r → 隱式 Treap 翻轉區間 [l, r]（1-indexed，閉區間）
    非數字 token 一律視為「輸入到此為止」；讀不滿的那一組直接丟棄。
    n、q 截斷到 [0, 5000]；數值截斷到 [-10^12, 10^12]。
輸出格式（stdout，共 5 行）：
    第 1 行  BST 的元素個數
    第 2 行  BST 的中序遍歷（空格分隔；空則為空行）
    第 3 行  op = 3 / 4 / 5 的查詢結果（依執行順序，空格分隔）
    第 4 行  隱式 Treap 最終序列（空格分隔）
    第 5 行  1 表示「FHQ 版與旋轉版中序一致、且堆性質與 BST 性質都成立」
無 stdin 輸入時運行內置斷言測試並輸出 `all tests passed`。
"""

import random
import re
import sys
from typing import List, Optional, Sequence, Tuple

NUM_RE = re.compile(r"^[+-]?[0-9]+$")

MAXN = 5000
VMAX = 10 ** 12
MASK64 = (1 << 64) - 1

SEED_SEQ = 0x9E3779B97F4A7C15 & MASK64
SEED_FHQ = 0xC2B2AE3D27D4EB4F & MASK64
SEED_ROT = 0x165667B19E3779F9 & MASK64


# ---------------------------------------------------------------- 確定性亂數


class Rng:
    """xorshift64：兩個語言逐位元一致（uint64 語義由 MASK64 保證）。"""

    __slots__ = ("s",)

    def __init__(self, seed: int) -> None:
        self.s = seed & MASK64

    def next_u64(self) -> int:
        x = self.s
        x ^= (x << 13) & MASK64
        x ^= x >> 7
        x ^= (x << 17) & MASK64
        self.s = x & MASK64
        return self.s

    def next_prio(self) -> int:
        """回傳 1 .. 2^31 的優先級（小根堆，越小越靠近根）。"""
        return ((self.next_u64() >> 33) & 0x7FFFFFFF) + 1


# ---------------------------------------------------------------- Treap


class Treap:
    """陣列池實作的 Treap（-1 表示空節點）。"""

    __slots__ = ("key", "prio", "lc", "rc", "sz", "rev")

    def __init__(self) -> None:
        self.key: List[int] = []
        self.prio: List[int] = []
        self.lc: List[int] = []
        self.rc: List[int] = []
        self.sz: List[int] = []
        self.rev: List[int] = []

    # ---- 基礎 ----

    def new_node(self, key: int, prio: int) -> int:
        idx = len(self.key)
        self.key.append(key)
        self.prio.append(prio)
        self.lc.append(-1)
        self.rc.append(-1)
        self.sz.append(1)
        self.rev.append(0)
        return idx

    def size(self, t: int) -> int:
        return 0 if t < 0 else self.sz[t]

    def pull(self, t: int) -> None:
        self.sz[t] = 1 + self.size(self.lc[t]) + self.size(self.rc[t])

    def push(self, t: int) -> None:
        """下傳翻轉懶標記（XOR 語義：打兩次等於沒打）。"""
        if t < 0 or self.rev[t] == 0:
            return
        self.lc[t], self.rc[t] = self.rc[t], self.lc[t]
        if self.lc[t] >= 0:
            self.rev[self.lc[t]] ^= 1
        if self.rc[t] >= 0:
            self.rev[self.rc[t]] ^= 1
        self.rev[t] = 0

    # ---- FHQ：merge / split ----

    def merge(self, a: int, b: int) -> int:
        if a < 0:
            return b
        if b < 0:
            return a
        if self.prio[a] < self.prio[b]:
            self.push(a)
            self.rc[a] = self.merge(self.rc[a], b)
            self.pull(a)
            return a
        self.push(b)
        self.lc[b] = self.merge(a, self.lc[b])
        self.pull(b)
        return b

    def split_key(self, t: int, key: int) -> Tuple[int, int]:
        """切成 (所有 key ≤ key, 所有 key > key)。"""
        if t < 0:
            return (-1, -1)
        self.push(t)
        if self.key[t] <= key:
            a, b = self.split_key(self.rc[t], key)
            self.rc[t] = a
            self.pull(t)
            return (t, b)
        a, b = self.split_key(self.lc[t], key)
        self.lc[t] = b
        self.pull(t)
        return (a, t)

    def split_size(self, t: int, k: int) -> Tuple[int, int]:
        """切成 (中序前 k 個, 剩下的)；k ≤ 0 時左為空，k ≥ size 時右為空。"""
        if t < 0:
            return (-1, -1)
        self.push(t)
        ls = self.size(self.lc[t])
        if k <= ls:
            a, b = self.split_size(self.lc[t], k)
            self.lc[t] = b
            self.pull(t)
            return (a, t)
        a, b = self.split_size(self.rc[t], k - ls - 1)
        self.rc[t] = a
        self.pull(t)
        return (t, b)

    # ---- FHQ 版插入 / 刪除 ----

    def fhq_insert(self, root: int, key: int, prio: int) -> int:
        a, b = self.split_key(root, key)
        return self.merge(self.merge(a, self.new_node(key, prio)), b)

    def fhq_erase(self, root: int, key: int) -> int:
        a, b = self.split_key(root, key)
        c, d = self.split_key(a, key - 1)
        # d 裡面是「恰好等於 key」的那一批；丟掉其中一個
        # （拿 d 的根開刀，把它的左右子樹接回去），剩下的合併回來
        if d >= 0:
            d = self.merge(self.lc[d], self.rc[d])
        return self.merge(self.merge(c, d), b)

    # ---- 旋轉版插入 / 刪除 ----

    def rot_right(self, t: int) -> int:
        l = self.lc[t]
        self.push(t)
        self.push(l)
        self.lc[t] = self.rc[l]
        self.rc[l] = t
        self.pull(t)
        self.pull(l)
        return l

    def rot_left(self, t: int) -> int:
        r = self.rc[t]
        self.push(t)
        self.push(r)
        self.rc[t] = self.lc[r]
        self.lc[r] = t
        self.pull(t)
        self.pull(r)
        return r

    def rot_insert(self, t: int, key: int, prio: int) -> int:
        if t < 0:
            return self.new_node(key, prio)
        if key <= self.key[t]:
            self.lc[t] = self.rot_insert(self.lc[t], key, prio)
            self.pull(t)
            if self.prio[self.lc[t]] < self.prio[t]:
                t = self.rot_right(t)
        else:
            self.rc[t] = self.rot_insert(self.rc[t], key, prio)
            self.pull(t)
            if self.prio[self.rc[t]] < self.prio[t]:
                t = self.rot_left(t)
        return t

    def rot_erase(self, t: int, key: int) -> int:
        """把目標節點一路旋轉到葉子再丟棄。"""
        if t < 0:
            return -1
        if key < self.key[t]:
            self.lc[t] = self.rot_erase(self.lc[t], key)
            self.pull(t)
            return t
        if key > self.key[t]:
            self.rc[t] = self.rot_erase(self.rc[t], key)
            self.pull(t)
            return t
        if self.lc[t] < 0:
            return self.rc[t]
        if self.rc[t] < 0:
            return self.lc[t]
        if self.prio[self.lc[t]] < self.prio[self.rc[t]]:
            t = self.rot_right(t)
            self.rc[t] = self.rot_erase(self.rc[t], key)
            self.pull(t)
            return t
        t = self.rot_left(t)
        self.lc[t] = self.rot_erase(self.lc[t], key)
        self.pull(t)
        return t

    # ---- 查詢 ----

    def contains(self, root: int, key: int) -> int:
        t = root
        while t >= 0:
            self.push(t)
            if key == self.key[t]:
                return 1
            t = self.lc[t] if key < self.key[t] else self.rc[t]
        return 0

    def kth(self, root: int, k: int) -> int:
        """第 k 小（1-indexed）；越界回傳 -1。"""
        if k < 1 or k > self.size(root):
            return -1
        t = root
        while t >= 0:
            self.push(t)
            ls = self.size(self.lc[t])
            if k <= ls:
                t = self.lc[t]
            elif k == ls + 1:
                return self.key[t]
            else:
                k -= ls + 1
                t = self.rc[t]
        return -1

    def rank(self, root: int, key: int) -> int:
        """嚴格小於 key 的元素個數。"""
        res = 0
        t = root
        while t >= 0:
            self.push(t)
            if self.key[t] < key:
                res += self.size(self.lc[t]) + 1
                t = self.rc[t]
            else:
                t = self.lc[t]
        return res

    def to_list(self, t: int) -> List[int]:
        """中序遍歷（疊代版，避免遞迴深度風險）。"""
        out: List[int] = []
        st: List[int] = []
        cur = t
        while cur >= 0 or st:
            while cur >= 0:
                st.append(cur)
                self.push(cur)
                cur = self.lc[cur]
            cur = st.pop()
            out.append(self.key[cur])
            cur = self.rc[cur]
        return out

    # ---- 隱式 Treap 建構與區間翻轉 ----

    def build_sequence(self, vals: Sequence[int], rng: Rng) -> int:
        root = -1
        for v in vals:
            root = self.merge(root, self.new_node(v, rng.next_prio()))
        return root

    def reverse_range(self, root: int, l: int, r: int) -> int:
        """翻轉 1-indexed 閉區間 [l, r]（l、r 需已合法化到 [1, n]）。"""
        n = self.size(root)
        if n == 0 or l > r or l > n:
            return root
        if r > n:
            r = n
        if l < 1:
            l = 1
        left, rest = self.split_size(root, r)
        mid_a, mid = self.split_size(left, l - 1)
        if mid >= 0:
            self.rev[mid] ^= 1
        return self.merge(self.merge(mid_a, mid), rest)

    # ---- 結構性質檢查 ----

    def check_bst(self, t: int, lo: Optional[int], hi: Optional[int]) -> bool:
        if t < 0:
            return True
        self.push(t)
        if lo is not None and self.key[t] < lo:
            return False
        if hi is not None and self.key[t] > hi:
            return False
        if self.size(t) != 1 + self.size(self.lc[t]) + self.size(self.rc[t]):
            return False
        return (self.check_bst(self.lc[t], lo, self.key[t])
                and self.check_bst(self.rc[t], self.key[t], hi))

    def check_heap(self, t: int) -> bool:
        if t < 0:
            return True
        if self.lc[t] >= 0 and self.prio[self.lc[t]] < self.prio[t]:
            return False
        if self.rc[t] >= 0 and self.prio[self.rc[t]] < self.prio[t]:
            return False
        return self.check_heap(self.lc[t]) and self.check_heap(self.rc[t])


# ---------------------------------------------------------------- 主流程


def run(init: Sequence[int], ops: Sequence[Tuple[int, int, int]]) -> str:
    """回傳 5 行輸出（不含末尾換行）。"""
    tr = Treap()

    # 1) 隱式 Treap（獨立亂數流）
    rng_seq = Rng(SEED_SEQ)
    seq_root = tr.build_sequence(init, rng_seq)

    # 2) 兩套 BST（各自獨立亂數流）
    rng_fhq = Rng(SEED_FHQ)
    rng_rot = Rng(SEED_ROT)
    fhq_root = -1
    rot_root = -1

    queries: List[int] = []

    for op, a, b in ops:
        if op == 1:
            p1 = rng_fhq.next_prio()
            p2 = rng_rot.next_prio()
            fhq_root = tr.fhq_insert(fhq_root, a, p1)
            rot_root = tr.rot_insert(rot_root, a, p2)
        elif op == 2:
            fhq_root = tr.fhq_erase(fhq_root, a)
            rot_root = tr.rot_erase(rot_root, a)
        elif op == 3:
            queries.append(tr.contains(fhq_root, a))
        elif op == 4:
            queries.append(tr.kth(fhq_root, a))
        elif op == 5:
            queries.append(tr.rank(fhq_root, a))
        elif op == 6:
            n = tr.size(seq_root)
            if n > 0:
                l = max(1, min(a, n))
                r = max(1, min(b, n))
                if l > r:
                    l, r = r, l
                seq_root = tr.reverse_range(seq_root, l, r)

    inorder_fhq = tr.to_list(fhq_root)
    inorder_rot = tr.to_list(rot_root)
    seq = tr.to_list(seq_root)

    # 隱式 Treap 的 key 是序列值、不滿足 BST 性質，因此只檢查它的堆性質
    ok = 1 if (inorder_fhq == inorder_rot
               and tr.check_bst(fhq_root, None, None)
               and tr.check_heap(fhq_root)
               and tr.check_bst(rot_root, None, None)
               and tr.check_heap(rot_root)
               and tr.check_heap(seq_root)) else 0

    lines = [
        str(tr.size(fhq_root)),
        " ".join(str(v) for v in inorder_fhq),
        " ".join(str(v) for v in queries),
        " ".join(str(v) for v in seq),
        str(ok),
    ]
    return "\n".join(lines)


# ---------------------------------------------------------------- 對拍用暴力


def bf_run(init: Sequence[int], ops: Sequence[Tuple[int, int, int]]) -> str:
    """用 list 直接模擬，驗證 Treap 的語義。"""
    seq = list(init)
    bst: List[int] = []
    queries: List[int] = []

    for op, a, b in ops:
        if op == 1:
            # 插入：保持有序（允許重複）
            lo, hi = 0, len(bst)
            while lo < hi:
                mid = (lo + hi) // 2
                if bst[mid] < a:
                    lo = mid + 1
                else:
                    hi = mid
            bst.insert(lo, a)
        elif op == 2:
            if a in bst:
                bst.remove(a)
        elif op == 3:
            queries.append(1 if a in bst else 0)
        elif op == 4:
            queries.append(bst[a - 1] if 1 <= a <= len(bst) else -1)
        elif op == 5:
            lo, hi = 0, len(bst)
            while lo < hi:
                mid = (lo + hi) // 2
                if bst[mid] < a:
                    lo = mid + 1
                else:
                    hi = mid
            queries.append(lo)
        elif op == 6:
            n = len(seq)
            if n > 0:
                l = max(1, min(a, n))
                r = max(1, min(b, n))
                if l > r:
                    l, r = r, l
                seq[l - 1:r] = list(reversed(seq[l - 1:r]))

    lines = [
        str(len(bst)),
        " ".join(str(v) for v in bst),
        " ".join(str(v) for v in queries),
        " ".join(str(v) for v in seq),
        "1",
    ]
    return "\n".join(lines)


# ---------------------------------------------------------------- 測試


def check(init: Sequence[int], ops: Sequence[Tuple[int, int, int]]) -> None:
    got = run(init, ops)
    want = bf_run(init, ops)
    assert got == want, (list(init), list(ops), got, want)


def run_tests() -> None:
    # ---- 固定用例 ----
    assert run([1, 2, 3, 4, 5], []) == (
        "0\n\n\n1 2 3 4 5\n1")
    assert run([1, 2, 3, 4, 5], [(6, 2, 4)]) == (
        "0\n\n\n1 4 3 2 5\n1")
    assert run([5, 4, 3, 2, 1], [(1, 3, 0), (1, 1, 0), (1, 2, 0),
                                 (3, 2, 0), (4, 2, 0), (5, 3, 0)]) == (
        "3\n1 2 3\n1 2 2\n5 4 3 2 1\n1")
    assert run([1, 2, 3], [(1, 7, 0), (1, 7, 0), (2, 7, 0),
                           (3, 7, 0), (4, 1, 0)]) == (
        "1\n7\n1 7\n1 2 3\n1")

    check([], [])
    check([1], [])
    check([1, 2, 3], [(1, 2, 0)])
    check([1, 2, 3], [(6, 1, 3)])
    check([1, 2, 3], [(6, 3, 1)])
    check([1, 2, 3, 4, 5], [(6, 2, 4), (6, 1, 5)])
    check([1, 2, 3, 4, 5], [(6, 0, 9)])
    check([10, 20, 30], [(1, 20, 0), (1, 20, 0), (1, 10, 0),
                         (3, 20, 0), (4, 1, 0), (4, 3, 0), (5, 20, 0)])
    check([1], [(1, 5, 0), (2, 5, 0), (2, 5, 0), (3, 5, 0)])
    check([1, 2], [(4, 0, 0), (4, 1, 0)])

    # ---- 隨機對拍 ----
    rng = random.Random(20261014)
    for _ in range(220):
        n = rng.randint(0, 14)
        init = [rng.randint(-9, 9) for _ in range(n)]
        q = rng.randint(0, 16)
        ops: List[Tuple[int, int, int]] = []
        for _ in range(q):
            op = rng.randint(1, 6)
            a = rng.randint(-3, 12)
            b = rng.randint(-3, 12)
            ops.append((op, a, b))
        check(init, ops)

    # ---- 中量級：結構性質與效能 ----
    for _ in range(6):
        n = rng.randint(200, 600)
        init = list(range(n))
        ops = []
        for _ in range(400):
            op = rng.randint(1, 6)
            ops.append((op, rng.randint(0, n + 50), rng.randint(0, n + 50)))
        check(init, ops)

    # ---- 只翻轉的大序列（不跑 O(q·n) 暴力，只驗證可逆性）----
    base_seq = list(range(1, 801))
    base = "0\n\n\n" + " ".join(str(v) for v in base_seq) + "\n1"
    assert run(base_seq, [(6, 1, 800), (6, 100, 200)]) != base
    assert run(base_seq, [(6, 100, 200), (6, 1, 800)]) != base
    # 翻轉兩次回到原狀（懶標記 XOR 語義）
    assert run(base_seq, [(6, 1, 800), (6, 1, 800)]) == base


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

    def clamp(v: Optional[int], lo: int, hi: int, default: int) -> int:
        if v is None:
            return default
        return max(lo, min(hi, v))

    def clamp_val(v: Optional[int]) -> int:
        if v is None:
            return 0
        return max(-VMAX, min(VMAX, v))

    n = clamp(read_int(), 0, MAXN, 0)
    init: List[int] = []
    for _ in range(n):
        v = read_int()
        if v is None:
            break
        init.append(clamp_val(v))

    q = clamp(read_int(), 0, MAXN, 0)
    ops: List[Tuple[int, int, int]] = []
    for _ in range(q):
        o = read_int()
        a = read_int()
        b = read_int()
        if o is None or a is None or b is None:
            break
        ops.append((clamp(o, 0, 6, 0), clamp_val(a), clamp_val(b)))

    out = run(init, ops)
    emit(out + "\n")


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
