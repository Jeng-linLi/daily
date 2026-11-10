# -*- coding: utf-8 -*-
"""配對堆 Pairing Heap（two-pass combine / 支援 decrease-key）

題意：維護 10 個（編號 0..9）可合併的小根堆，支援六種操作：
  1 h x     把值 x 插入堆 h
  2 h       彈出堆 h 的最小值並輸出；空堆輸出 -1
  3 h       詢問堆 h 的最小值（不彈出）；空堆輸出 -1
  4 a b     把堆 b 併入堆 a（之後 b 為空；a == b 時視為無操作）
  5 h       輸出堆 h 的元素個數
  6 h x y   把堆 h 中的值 x 改為 y（decrease-key，無輸出）。
            x 不在堆 h、y 已存在或 x == y 時視為無操作。

配對堆要點：
  - 節點用 child（第一個兒子）+ sibling（下一個兄弟）的左兒子-右兄弟表示。
  - meld(a, b)：值大的根作為值小根的第一個兒子，O(1)。
  - delete-min：根的所有兒子先**兩兩配對**合併（第一遍），再把結果
    從右往左**依序合併**（第二遍），均攤 O(log n)。
  - decrease-key：把目標節點從父親的兒子鏈上**摘下**、改值後 meld 回根，
    攤還 O(log n)。這是配對堆相對左偏樹/二叉堆的特色操作。
  - 注意：本題的「6」允許任意改值（含調大）。被改節點的兒子鏈必須
    **摘下 two-pass 重組後再與自身 meld**——調大時兒子 ≥ 舊值卻可能
    < 新值，子樹內部堆序會被破壞；調小時重組無害。

IO 模式（stdin 有內容時）：q 行操作，輸出每個 2/3/5 操作的結果。
  值域 |x| ≤ 10^9；測試保證同一時刻全域值互異（decrease-key 靠值定位節點）。
無 stdin 時跑內建斷言（與有序陣列參考實作對拍 + 堆序不變量實走檢查）
並輸出 all tests passed。
"""
import sys
from bisect import insort


class Node:
    __slots__ = ("val", "child", "sib", "par")

    def __init__(self, val):
        self.val = val
        self.child = None
        self.sib = None
        self.par = None


def meld(a, b):
    """合併兩棵配對堆：值大的根變成值小根的第一個兒子。"""
    if a is None:
        return b
    if b is None:
        return a
    if b.val < a.val:
        a, b = b, a
    # b 掛到 a 的兒子鏈頭；par 指向父親，a 的舊兒子 par 仍是 a，不動
    b.sib = a.child
    b.par = a
    a.child = b
    return a


def _collect(head):
    """把兒子鏈摘成 list，並清掉 par/sib 指針。"""
    out = []
    c = head
    while c is not None:
        nxt = c.sib
        c.sib = None
        c.par = None
        out.append(c)
        c = nxt
    return out


def two_pass_combine(children):
    """兩遍合併：第一遍兩兩配對，第二遍從右往左依序合併。"""
    if not children:
        return None
    # 第一遍：相鄰兩兩 meld
    pairs = []
    for i in range(0, len(children) - 1, 2):
        pairs.append(meld(children[i], children[i + 1]))
    if len(children) % 2 == 1:
        pairs.append(children[-1])
    # 第二遍：從右往左
    cur = pairs[-1]
    for i in range(len(pairs) - 2, -1, -1):
        cur = meld(pairs[i], cur)
    return cur


class PairingHeap:
    def __init__(self):
        self.root = None
        self.size = 0
        self.pos = {}  # 值 -> 節點（同時刻全域值互異）

    def push(self, x):
        node = Node(x)
        self.root = meld(self.root, node)
        self.pos[x] = node
        self.size += 1

    def pop(self):
        if self.root is None:
            return -1
        node = self.root
        v = node.val
        self.root = two_pass_combine(_collect(node.child))
        self.pos.pop(v, None)
        self.size -= 1
        return v

    def top(self):
        return self.root.val if self.root is not None else -1

    def decrease_key(self, x, y):
        """把值 x 的節點改為 y（允許調大或調小）。

        摘下節點後，其兒子鏈也要 two-pass 重組再與自身 meld：
        調大時（y > x）兒子 ≥ 舊值 x 卻可能 < 新值 y，子樹內部
        堆序被破壞，必須重建；調小時重組亦無害。
        """
        if x == y or x not in self.pos or y in self.pos:
            return
        node = self.pos.pop(x)
        if node is not self.root:
            # 從父親的兒子鏈上摘下 node（兄弟的 par 仍指向父親，無需改動）
            par = node.par
            if par.child is node:
                par.child = node.sib
            else:
                prev = par.child
                while prev.sib is not node:
                    prev = prev.sib
                prev.sib = node.sib
        # 摘下兒子鏈重組（清掉 child，防止 meld 掛到陳舊鏈）
        children = _collect(node.child)
        node.child = None
        node.sib = None
        node.par = None
        node.val = y
        tree = meld(two_pass_combine(children), node)
        if node is self.root:
            self.root = tree
        else:
            self.root = meld(self.root, tree)
        self.pos[y] = node

    def meld_into(self, other):
        """把堆 other 併入 self。"""
        self.root = meld(self.root, other.root)
        self.size += other.size
        self.pos.update(other.pos)
        other.root = None
        other.size = 0
        other.pos = {}


def check_invariants(root):
    """實走檢查：堆序（父 ≤ 子）+ 節點計數；回傳節點數。"""
    if root is None:
        return 0
    cnt = 0
    stk = [(root, None)]
    while stk:
        u, lo = stk.pop()
        cnt += 1
        assert lo is None or lo <= u.val, "heap order violated"
        c = u.child
        while c is not None:
            assert c.par is u, "parent pointer broken"
            stk.append((c, u.val))
            c = c.sib
    return cnt


class RefHeap:
    """參考實作：有序陣列（值互異，remove/insort 直接可用）。"""

    def __init__(self):
        self.a = []

    def push(self, x):
        insort(self.a, x)

    def pop(self):
        return self.a.pop(0) if self.a else -1

    def top(self):
        return self.a[0] if self.a else -1

    def meld_into(self, other):
        self.a = sorted(self.a + other.a)
        other.a = []


def run_ops(ops, hl, hr):
    """同一組操作分別跑配對堆與參考實作，回傳雙方輸出。"""
    out_l, out_r = [], []
    for op in ops:
        c = op[0]
        if c == 1:
            h, x = op[1], op[2]
            hl[h].push(x)
            hr[h].push(x)
        elif c == 2:
            h = op[1]
            out_l.append(hl[h].pop())
            out_r.append(hr[h].pop())
        elif c == 3:
            h = op[1]
            out_l.append(hl[h].top())
            out_r.append(hr[h].top())
        elif c == 4:
            a, b = op[1], op[2]
            if a != b:
                hl[a].meld_into(hl[b])
                hr[a].meld_into(hr[b])
        elif c == 5:
            h = op[1]
            out_l.append(hl[h].size)
            out_r.append(len(hr[h].a))
        else:
            h, x, y = op[1], op[2], op[3]
            hl[h].decrease_key(x, y)
            if x in hr[h].a and y not in hr[h].a and x != y:
                hr[h].a.remove(x)
                insort(hr[h].a, y)
    return out_l, out_r


def fresh_heaps():
    return ([PairingHeap() for _ in range(10)], [RefHeap() for _ in range(10)])


def run_tests():
    # 1) 空堆行為
    hl, hr = fresh_heaps()
    assert hl[0].pop() == -1 and hl[0].top() == -1 and hl[0].size == 0
    assert check_invariants(hl[0].root) == 0

    # 2) 固定情境（全域值互異；期望值由實跑取得）
    hl, hr = fresh_heaps()
    ops = [(1, 0, 5), (1, 0, 2), (1, 0, 8), (3, 0), (2, 0), (5, 0),
           (6, 0, 8, 1), (3, 0), (1, 1, 7), (4, 0, 1), (5, 0), (5, 1),
           (2, 0), (2, 0), (2, 0), (2, 0)]
    out_l, out_r = run_ops(ops, hl, hr)
    assert out_l == out_r, "fixed scenario mismatch vs reference"
    # 實跑取值：top=2, pop=2, size=2；decrease 8→1 後 top=1；
    # 合併後 size=3, h1=0；pop 依序 1, 5, 7, 最後空堆 -1
    assert out_l == [2, 2, 2, 1, 3, 0, 1, 5, 7, -1], out_l
    for h in hl:
        assert check_invariants(h.root) == h.size

    # 3) 自製 LCG 隨機對拍（全域值互異）+ 不變量檢查
    seed = 20261109
    mod = 1 << 31

    def rnd(n):
        nonlocal seed
        seed = (seed * 1103515245 + 12345) % mod
        return seed % n

    used = set()

    def fresh_val():
        while True:
            v = rnd(2_000_001) - 1_000_000
            if v not in used:
                used.add(v)
                return v

    for tc in range(120):
        hl, hr = fresh_heaps()
        used.clear()
        q = 60 + rnd(120)
        ops = []
        for _ in range(q):
            c = rnd(100)
            if c < 40:
                ops.append((1, rnd(10), fresh_val()))
            elif c < 62:
                ops.append((2, rnd(10)))
            elif c < 77:
                ops.append((3, rnd(10)))
            elif c < 87:
                ops.append((4, rnd(10), rnd(10)))
            elif c < 92:
                ops.append((5, rnd(10)))
            else:
                # decrease-key：挑堆 h 裡已存在的值改成全新值
                h = rnd(10)
                if hl[h].pos:
                    x = rnd(2_000_001) - 1_000_000
                    while x not in hl[h].pos:
                        x = rnd(2_000_001) - 1_000_000
                    y = fresh_val()
                    ops.append((6, h, x, y))
                else:
                    ops.append((1, h, fresh_val()))
        out_l, out_r = run_ops(ops, hl, hr)
        assert out_l == out_r, "random scenario %d mismatch" % tc
        for h in hl:
            assert check_invariants(h.root) == h.size, "size mismatch"


def main():
    data = sys.stdin.buffer.read().split()
    if not data:
        run_tests()
        print("all tests passed")
        return
    pos = 0

    def read_int():
        # token 耗盡或非數字時回傳 None（與 C++ 的 scanf 失敗語義對齊：
        # 一律當「輸入到此為止」，不拋例外）
        nonlocal pos
        if pos >= len(data):
            return None
        tok = data[pos]
        pos += 1
        try:
            return int(tok)
        except ValueError:
            return None

    q = read_int()
    if q is None:
        return
    heaps = [PairingHeap() for _ in range(10)]
    out = []
    for _ in range(q):
        c = read_int()
        if c is None:
            break
        if c == 1:
            h, x = read_int(), read_int()
            if h is None or x is None:
                break
            if 0 <= h <= 9:
                heaps[h].push(x)
        elif c == 2:
            h = read_int()
            if h is None:
                break
            out.append(heaps[h].pop() if 0 <= h <= 9 else -1)
        elif c == 3:
            h = read_int()
            if h is None:
                break
            out.append(heaps[h].top() if 0 <= h <= 9 else -1)
        elif c == 4:
            a, b = read_int(), read_int()
            if a is None or b is None:
                break
            if 0 <= a <= 9 and 0 <= b <= 9 and a != b:
                heaps[a].meld_into(heaps[b])
        elif c == 5:
            h = read_int()
            if h is None:
                break
            out.append(heaps[h].size if 0 <= h <= 9 else 0)
        elif c == 6:
            h, x, y = read_int(), read_int(), read_int()
            if h is None or x is None or y is None:
                break
            if 0 <= h <= 9:
                heaps[h].decrease_key(x, y)
    if out:
        sys.stdout.write("\n".join(map(str, out)) + "\n")


if __name__ == "__main__":
    main()
