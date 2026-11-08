# -*- coding: utf-8 -*-
"""左偏樹 Leftist Heap（可合併堆 Mergeable Heap）

題意：維護 10 個（編號 0..9）可合併的小根堆，支援五種操作：
  1 h x   把值 x 插入堆 h
  2 h     彈出堆 h 的最小值並輸出；空堆輸出 -1
  3 h     詢問堆 h 的最小值（不彈出）；空堆輸出 -1
  4 a b   把堆 b 併入堆 a（之後 b 為空；a == b 時視為無操作）
  5 h     輸出堆 h 的元素個數

左偏樹要點：
  - 節點帶 dist（外節點距離 / null 路徑長度）：dist(空) = -1、
    葉子 dist = 0、dist(x) = dist(右子) + 1。
  - 左偏性質：左兒子的 dist ≥ 右兒子的 dist ⟹ 樹高 O(log n)。
  - merge：值大的根掛到值小的根的右子樹上遞歸合併，回溯時若違反
    左偏性質則交換左右兒子。插入 = 與單節點 merge；彈出 = merge(左, 右)。
  - push/pop/top/merge 全部 O(log n)。

IO 模式（stdin 有內容時）：
  q
  q 行操作（如上），輸出每個 2/3/5 操作的結果，一行一個。
  值域 |x| ≤ 10^9（隨機測試保證同時刻各堆內值互異，輸出順序才能唯一）。
無 stdin 時跑內建斷言（與「有序陣列 + bisect」參考實作對拍 + 左偏性/堆序
不變量實走檢查）並輸出 all tests passed。
"""
import sys
from bisect import insort


class Node:
    __slots__ = ("val", "dist", "l", "r")

    def __init__(self, val):
        self.val = val
        self.dist = 0
        self.l = None
        self.r = None


def meld(a, b):
    """合併兩棵左偏樹，回傳新根。"""
    if a is None:
        return b
    if b is None:
        return a
    if b.val < a.val:
        a, b = b, a
    a.r = meld(a.r, b)
    if a.l is None or (a.r is not None and a.l.dist < a.r.dist):
        a.l, a.r = a.r, a.l
    a.dist = 0 if a.r is None else a.r.dist + 1
    return a


class LeftistHeap:
    def __init__(self):
        self.root = None
        self.size = 0

    def push(self, x):
        self.root = meld(self.root, Node(x))
        self.size += 1

    def pop(self):
        if self.root is None:
            return -1
        v = self.root.val
        self.root = meld(self.root.l, self.root.r)
        self.size -= 1
        return v

    def top(self):
        return self.root.val if self.root is not None else -1


def check_invariants(root):
    """實走檢查：堆序性 + 左偏性 + dist 定義一致；回傳節點數。"""
    if root is None:
        return 0
    cnt = 0
    stk = [root]
    while stk:
        u = stk.pop()
        cnt += 1
        # dist 慣例：dist(空) = -1，葉子 dist = 0，dist(x) = dist(右子) + 1
        dl = u.l.dist if u.l is not None else -1
        dr = u.r.dist if u.r is not None else -1
        assert u.l is None or u.l.val >= u.val, "heap order violated"
        assert u.r is None or u.r.val >= u.val, "heap order violated"
        assert dl >= dr, "leftist property violated"
        assert u.dist == dr + 1, "dist field inconsistent"
        if u.l is not None:
            stk.append(u.l)
        if u.r is not None:
            stk.append(u.r)
    return cnt


class RefHeap:
    """參考實作：有序陣列（O(n) 插入，測試規模足夠）。"""

    def __init__(self):
        self.a = []

    def push(self, x):
        insort(self.a, x)

    def pop(self):
        if not self.a:
            return -1
        return self.a.pop(0)

    def top(self):
        return self.a[0] if self.a else -1


def run_ops(ops, heaps_l, heaps_r):
    """對同一組操作分別跑左偏樹與參考實作，回傳雙方輸出（用於比對）。"""
    out_l, out_r = [], []
    for op in ops:
        c = op[0]
        if c == 1:
            h, x = op[1], op[2]
            heaps_l[h].push(x)
            heaps_r[h].push(x)
        elif c == 2:
            h = op[1]
            out_l.append(heaps_l[h].pop())
            out_r.append(heaps_r[h].pop())
        elif c == 3:
            h = op[1]
            out_l.append(heaps_l[h].top())
            out_r.append(heaps_r[h].top())
        elif c == 4:
            a, b = op[1], op[2]
            if a != b:
                la, lb = heaps_l[a], heaps_l[b]
                la.root = meld(la.root, lb.root)
                la.size += lb.size
                lb.root, lb.size = None, 0
                ra, rb = heaps_r[a], heaps_r[b]
                ra.a = sorted(ra.a + rb.a)
                rb.a = []
        else:
            h = op[1]
            out_l.append(heaps_l[h].size)
            out_r.append(len(heaps_r[h].a))
    return out_l, out_r


def fresh_heaps():
    return [LeftistHeap() for _ in range(10)], [RefHeap() for _ in range(10)]


def run_tests():
    # 1) 空堆行為
    hl, hr = fresh_heaps()
    assert hl[0].pop() == -1 and hl[0].top() == -1 and hl[0].size == 0

    # 2) 固定情境（值互異）
    hl, hr = fresh_heaps()
    ops = [(1, 0, 5), (1, 0, 2), (1, 0, 8), (3, 0), (2, 0), (5, 0),
           (1, 1, 7), (4, 0, 1), (5, 0), (5, 1), (2, 0), (2, 0), (2, 0)]
    out_l, out_r = run_ops(ops, hl, hr)
    assert out_l == out_r, "fixed scenario mismatch"
    for h in hl:
        if h.root is not None:
            check_invariants(h.root)
    # 實跑取值：top=2、pop=2、size=2 → 合併後 size=3、h1=0 → pop 5,7,8
    assert out_l == [2, 2, 2, 3, 0, 5, 7, 8]
    assert hl[0].size == 0 and hl[0].pop() == -1

    # 3) 自製 LCG 隨機對拍（同時刻全域值互異）+ 不變量檢查
    seed = 20261108
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
            if c < 45:
                ops.append((1, rnd(10), fresh_val()))
            elif c < 70:
                ops.append((2, rnd(10)))
            elif c < 85:
                ops.append((3, rnd(10)))
            elif c < 95:
                a, b = rnd(10), rnd(10)
                ops.append((4, a, b))
            else:
                ops.append((5, rnd(10)))
        out_l, out_r = run_ops(ops, hl, hr)
        assert out_l == out_r, "random scenario %d mismatch" % tc
        for h in hl:
            if h.root is not None:
                assert check_invariants(h.root) == h.size, "size mismatch"


def main():
    data = sys.stdin.buffer.read().split()
    if not data:
        run_tests()
        print("all tests passed")
        return
    pos = 0

    def read_int():
        nonlocal pos
        v = int(data[pos])
        pos += 1
        return v

    q = read_int()
    heaps = [LeftistHeap() for _ in range(10)]
    out = []
    for _ in range(q):
        c = read_int()
        if c == 1:
            h, x = read_int(), read_int()
            if 0 <= h <= 9:
                heaps[h].push(x)
        elif c == 2:
            h = read_int()
            out.append(heaps[h].pop() if 0 <= h <= 9 else -1)
        elif c == 3:
            h = read_int()
            out.append(heaps[h].top() if 0 <= h <= 9 else -1)
        elif c == 4:
            a, b = read_int(), read_int()
            if 0 <= a <= 9 and 0 <= b <= 9 and a != b:
                ha, hb = heaps[a], heaps[b]
                ha.root = meld(ha.root, hb.root)
                ha.size += hb.size
                hb.root, hb.size = None, 0
        elif c == 5:
            h = read_int()
            out.append(heaps[h].size if 0 <= h <= 9 else 0)
    if out:
        sys.stdout.write("\n".join(map(str, out)) + "\n")


if __name__ == "__main__":
    main()
