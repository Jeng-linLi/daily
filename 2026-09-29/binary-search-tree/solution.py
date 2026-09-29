"""二元搜尋樹（BST）：插入 / 搜尋 / 刪除 / 中序遍歷 / 樹高 / 合法性檢查

題意：
    先給出 n 個待插入的鍵值，再給出 m 個待刪除的鍵值，維護一棵二元搜尋樹並依次輸出：
      1. 全部插入後的中序遍歷（非嚴格升序）；
      2. 全部插入後的樹高（空樹爲 0，只有根節點時爲 1）；
      3. 在插入後的樹上依次搜尋 m 個鍵值的結果（找到 1、否則 0），空格分隔；
      4. 依次刪除這 m 個鍵值（每個鍵值只刪掉一個副本，不存在則忽略）後的中序遍歷；
      5. 刪除完成後的樹高。
    重複鍵值的約定：插入時「大於等於」當前節點值的走右子樹，因此中序遍歷是非嚴格升序。

思路：
    BST 的核心不變量：對任意節點 x，**左子樹所有鍵值 < x.val <= 右子樹所有鍵值**。
    有了這個不變量，搜尋就退化成二分：每比較一次就排除掉一棵子樹，
    插入 / 搜尋 / 刪除都只沿著一條從根到葉的路徑走，複雜度 O(h)，h 爲樹高。
    AVL / 紅黑樹之類的平衡 BST 能把 h 壓到 O(log n)；而**不做平衡**的樸素 BST
    在「有序插入」時會退化成一條鏈（h = n），這是樹類題目最常見的性能坑。

    刪除是最麻煩的一步，分三種情形：
      - 葉子節點：直接摘掉；
      - 只有一棵子樹：讓那棵子樹頂替它的位置；
      - 左右子樹都在：不能直接摘，否則兩棵子樹會一起丟失。
        標準做法是找**中序後繼**（右子樹中的最小節點）——它是「比當前值大的最小者」。
        把它的值搬到當前節點（不破壞 BST 不變量），再把那個後繼節點刪掉。
        後繼節點至多只有右孩子，所以第二趟刪除必定落在前兩種簡單情形。
        （用中序前驅、即左子樹的最大節點，同理可行。）

    本實作的插入 / 刪除 / 遍歷 / 求高全部寫成**迭代**版本，
    避免樹退化成鏈時 Python 遞迴深度爆掉。

輸入格式（stdin，數字按空白分隔即可）：
    n m
    a1 a2 ... an        （依次插入；n = 0 時這行省略）
    b1 b2 ... bm        （依次刪除；m = 0 時這行省略）
輸出格式（stdout）：
    第 1 行：插入後的中序遍歷，空格分隔（空樹輸出空行）
    第 2 行：插入後的樹高
    第 3 行：m 個鍵值的搜尋結果（1/0），空格分隔（m = 0 時輸出空行）
    第 4 行：刪除後的中序遍歷，空格分隔
    第 5 行：刪除後的樹高
無 stdin 輸入時運行內置斷言測試並輸出 `all tests passed`。
"""

import math
import random
import sys
from typing import Iterator, List, Optional


class Node:
    """BST 節點。__slots__ 省掉每節點的 __dict__，節點多時省內存。"""

    __slots__ = ("val", "left", "right")

    def __init__(self, val: int) -> None:
        self.val = val
        self.left: Optional["Node"] = None
        self.right: Optional["Node"] = None


# ---------------- 基本操作 ----------------

def insert(root: Optional[Node], val: int) -> Node:
    """插入 val（重複值走右子樹）。迭代版，時間 O(h)，空間 O(1)。"""
    if root is None:
        return Node(val)
    cur = root
    while True:
        if val < cur.val:                      # 嚴格小於走左邊
            if cur.left is None:
                cur.left = Node(val)
                return root
            cur = cur.left
        else:                                  # 大於或等於都走右邊（重複鍵值策略）
            if cur.right is None:
                cur.right = Node(val)
                return root
            cur = cur.right


def search(root: Optional[Node], val: int) -> bool:
    """判斷 val 是否存在。時間 O(h)，空間 O(1)。"""
    cur = root
    while cur is not None:
        if val == cur.val:
            return True
        cur = cur.left if val < cur.val else cur.right
    return False


def _replace_child(parent: Node, old: Node, new: Optional[Node]) -> None:
    """把 parent 的某個孩子 old 換成 new（用身份比較，不用值比較，避免重複值歧義）。"""
    if parent.left is old:
        parent.left = new
    else:
        parent.right = new


def delete(root: Optional[Node], val: int) -> Optional[Node]:
    """刪除一個值爲 val 的節點（不存在則原樣返回）。時間 O(h)，空間 O(1)。

    用一個哨兵節點當「根的父節點」，這樣連刪根節點都不需要特判。
    """
    dummy = Node(0)
    dummy.right = root
    parent, cur = dummy, root
    while cur is not None and cur.val != val:
        parent = cur
        cur = cur.left if val < cur.val else cur.right
    if cur is None:                            # 沒找到，什麼都不做
        return dummy.right

    if cur.left is None or cur.right is None:  # 情形 1 / 2：至多一棵子樹
        child = cur.left if cur.left is not None else cur.right
        _replace_child(parent, cur, child)
    else:                                      # 情形 3：兩棵子樹都在，找中序後繼
        succ_parent, succ = cur, cur.right
        while succ.left is not None:
            succ_parent = succ
            succ = succ.left
        cur.val = succ.val                     # 值搬上來，結構不動
        _replace_child(succ_parent, succ, succ.right)   # 後繼至多只有右孩子
    return dummy.right


def inorder(root: Optional[Node]) -> List[int]:
    """中序遍歷，結果即升序序列。迭代版（顯式棧），時間 O(n)，空間 O(h)。"""
    out: List[int] = []
    stack: List[Node] = []
    cur = root
    while cur is not None or stack:
        while cur is not None:                 # 一路向左到底
            stack.append(cur)
            cur = cur.left
        cur = stack.pop()
        out.append(cur.val)
        cur = cur.right                        # 轉向右子樹
    return out


def height(root: Optional[Node]) -> int:
    """樹高：空樹 0，單節點 1。迭代 DFS，時間 O(n)，空間 O(n)。"""
    if root is None:
        return 0
    best = 0
    stack = [(root, 1)]
    while stack:
        node, depth = stack.pop()
        if depth > best:
            best = depth
        if node.left is not None:
            stack.append((node.left, depth + 1))
        if node.right is not None:
            stack.append((node.right, depth + 1))
    return best


def size(root: Optional[Node]) -> int:
    """節點個數。迭代版，時間 O(n)。"""
    cnt = 0
    stack = [root] if root is not None else []
    while stack:
        node = stack.pop()
        cnt += 1
        if node.left is not None:
            stack.append(node.left)
        if node.right is not None:
            stack.append(node.right)
    return cnt


def is_valid(root: Optional[Node]) -> bool:
    """檢查 BST 不變量是否成立：左子樹嚴格小於，右子樹允許等於（重複鍵值策略）。

    做法是給每個節點帶一組取值區間 (lo, hi) 往下傳，區間邊界是否可取用
    lo_strict / hi_strict 兩個旗標區分。時間 O(n)，空間 O(n)。
    """
    if root is None:
        return True
    neg_inf, pos_inf = -math.inf, math.inf
    stack = [(root, neg_inf, False, pos_inf, False)]
    while stack:
        node, lo, lo_strict, hi, hi_strict = stack.pop()
        if node.val < lo or (lo_strict and node.val == lo):
            return False
        if node.val > hi or (hi_strict and node.val == hi):
            return False
        if node.left is not None:
            # 左子樹：上界是當前值，且必須嚴格小於
            stack.append((node.left, lo, lo_strict, node.val, True))
        if node.right is not None:
            # 右子樹：下界是當前值，可以取等號（重複值往右放）
            stack.append((node.right, node.val, False, hi, hi_strict))
    return True


def build(values: List[int]) -> Optional[Node]:
    """依序插入整個序列，返回根節點。"""
    root: Optional[Node] = None
    for v in values:
        root = insert(root, v)
    return root


# ---------------- 獨立對照實現（用於對拍） ----------------

def height_bfs(root: Optional[Node]) -> int:
    """另一種求高寫法：層序 BFS 數層數。用於交叉驗證 height()。"""
    if root is None:
        return 0
    level = [root]
    h = 0
    while level:
        h += 1
        nxt: List[Node] = []
        for node in level:
            if node.left is not None:
                nxt.append(node.left)
            if node.right is not None:
                nxt.append(node.right)
        level = nxt
    return h


def inorder_rec(root: Optional[Node]) -> List[int]:
    """遞迴版中序遍歷，用於交叉驗證迭代版（測試數據規模很小，不會爆棧）。"""
    if root is None:
        return []
    return inorder_rec(root.left) + [root.val] + inorder_rec(root.right)


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
    inserts = [_next_int(it) for _ in range(n)]
    deletes = [_next_int(it) for _ in range(m)]

    root = build(inserts)
    print(" ".join(str(v) for v in inorder(root)))                    # 插入後的中序
    print(height(root))                                               # 插入後的樹高
    print(" ".join("1" if search(root, v) else "0" for v in deletes))  # 逐個搜尋
    for v in deletes:
        root = delete(root, v)
    print(" ".join(str(v) for v in inorder(root)))                    # 刪除後的中序
    print(height(root))                                               # 刪除後的樹高


def run_tests() -> None:
    # README 示例：插入 [5, 3, 7, 3, 6]，刪除 [3, 7]
    root = build([5, 3, 7, 3, 6])
    assert inorder(root) == [3, 3, 5, 6, 7]
    assert height(root) == 3
    assert size(root) == 5
    assert is_valid(root)
    assert [search(root, v) for v in (3, 7, 9)] == [True, True, False]
    root = delete(root, 3)
    assert inorder(root) == [3, 5, 6, 7]      # 重複鍵值只刪掉一個
    assert is_valid(root)
    root = delete(root, 7)                    # 有兩個孩子的節點：用中序後繼頂替
    assert inorder(root) == [3, 5, 6]
    assert is_valid(root)
    assert height(root) == 2

    # 空樹
    assert inorder(None) == []
    assert height(None) == 0
    assert size(None) == 0
    assert is_valid(None)
    assert delete(None, 1) is None

    # 單節點
    one = build([42])
    assert inorder(one) == [42]
    assert height(one) == 1
    assert delete(one, 42) is None

    # 有序插入 → 退化成鏈（這正是樸素 BST 的軟肋）
    chain = build([1, 2, 3, 4, 5])
    assert inorder(chain) == [1, 2, 3, 4, 5]
    assert height(chain) == 5
    assert is_valid(chain)
    # 隨機插入通常會矮一些
    assert height(build([5, 3, 8, 1, 4])) == 3

    # 刪根節點 / 刪葉子 / 刪不存在的鍵
    t = build([4, 2, 6, 1, 3, 5, 7])
    assert height(t) == 3
    t = delete(t, 4)                          # 刪根：中序後繼是 5
    assert inorder(t) == [1, 2, 3, 5, 6, 7]
    assert is_valid(t)
    t = delete(t, 1)                          # 刪葉子
    assert inorder(t) == [2, 3, 5, 6, 7]
    t = delete(t, 100)                        # 鍵不存在：原樣返回
    assert inorder(t) == [2, 3, 5, 6, 7]
    assert is_valid(t)

    # 全部刪空
    t = build([2, 1, 3])
    for v in (2, 1, 3):
        t = delete(t, v)
        assert is_valid(t)
    assert t is None

    # 重複鍵值：全相同
    dup = build([7, 7, 7])
    assert inorder(dup) == [7, 7, 7]
    assert height(dup) == 3
    assert is_valid(dup)
    dup = delete(dup, 7)
    assert inorder(dup) == [7, 7]
    assert is_valid(dup)

    random.seed(20260929)

    # 隨機對拍：用有序列表模擬多重集合，每一步都與樹對比
    for _ in range(300):
        n = random.randint(0, 60)
        inserts = [random.randint(-8, 8) for _ in range(n)]   # 值域小 → 大量重複
        root = build(inserts)
        model = sorted(inserts)
        assert inorder(root) == model
        assert inorder(root) == inorder_rec(root)
        assert is_valid(root)
        assert size(root) == len(model)
        assert height(root) == height_bfs(root)
        if n > 0:
            # 樹高介於「完全二叉樹的下界」與「退化成鏈的上界」之間
            assert math.ceil(math.log2(n + 1)) <= height(root) <= n

        m = random.randint(0, 20)
        deletes = [random.randint(-8, 8) for _ in range(m)]
        for v in deletes:
            assert search(root, v) == (v in model)
            root = delete(root, v)
            if v in model:
                model.remove(v)               # 只刪一個副本
            assert inorder(root) == model     # 中序必須始終等於有序列表
            assert is_valid(root)
            assert size(root) == len(model)
            assert height(root) == height_bfs(root)

    print("all tests passed")


if __name__ == "__main__":
    raw = sys.stdin.read() if not sys.stdin.isatty() else ""
    if raw.strip():
        run_io(raw)
    else:
        run_tests()
