#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""費波那契堆（Fibonacci heap）：可合併堆，支援攤還 O(1) 的插入與 decrease-key。

命令（IO 模式，從 stdin 讀入）：
  push k      -> 插入鍵值 k，輸出分配到的 id（從 0 遞增）
  popmin      -> 彈出最小鍵值並輸出該鍵值；空堆輸出 -1
  min         -> 輸出當前最小鍵值（不彈出）；空堆輸出 -1
  dec id nk   -> 把 id 的鍵值降為 nk（僅當 nk <= 現值且 id 仍存活，否則忽略）
  size        -> 輸出當前元素個數

stdin 無內容時跑內建斷言測試，輸出 all tests passed。

核心要點：
  * 結構：一組「根串列」樹 + 每節點 children 串列；minRoot 恆指向最小根。
  * 堆序採字典序 (key, nid)：父節點 (key,nid) 嚴格小於孩子。如此全域
    最小節點必為根，且平手時 popmin 彈出的節點完全確定（鍵值相同取
    id 最小者），與參考實作、C++ 版逐位元組一致——dec 操作之後
    id 的死活才不會因平手分岔。
  * push：新節點掛進根串列，必要時更新 minRoot —— O(1)。
  * dec：改鍵後若破壞堆序（比父小）就 cut（剪下丟回根串列），
    父節點失去第二個孩子時級聯 cut（cascading cut）—— 攤還 O(1)。
  * popmin：把最小根的孩子全部升為根，然後 consolidate——
    依 degree 合併同度數的樹（key 小者當根），確保每個度數至多一棵樹，
    使第 k 小的樹至少有 F(k) 個節點 —— 攤還 O(log n)。
  * 本實作以 Python list 當根串列 / children 串列（語義與環狀雙向鏈等價）；
    popmin 的輸出只含鍵值，與平手時具體彈出哪個同鍵節點無關，輸出確定。
  * 應用：Dijkstra / Prim 的 O(E + V log V) 理論最優實作。
"""

import sys

INF = float('inf')


class Node:
    __slots__ = ('key', 'nid', 'deg', 'mark', 'parent', 'children', 'alive')

    def __init__(self, key, nid):
        self.key = key
        self.nid = nid
        self.deg = 0
        self.mark = False
        self.parent = None
        self.children = []
        self.alive = True


class FibonacciHeap:
    def __init__(self):
        self.roots = []
        self.min_root = None
        self.nodes = []   # id -> Node（含已彈出的死節點）
        self.cnt = 0

    def _link(self, y, x):
        """把根 y 掛到根 x 下面（(y.key,y.nid) >= (x.key,x.nid)）。"""
        y.parent = x
        x.children.append(y)
        x.deg += 1
        y.mark = False

    def push(self, key):
        node = Node(key, self.cnt)
        self.cnt += 1
        self.nodes.append(node)
        self.roots.append(node)
        if self.min_root is None or (node.key, node.nid) < (self.min_root.key, self.min_root.nid):
            self.min_root = node
        return node.nid

    def _consolidate(self):
        by_deg = {}
        for w in self.roots:
            x = w
            d = x.deg
            while d in by_deg:
                y = by_deg.pop(d)
                if (y.key, y.nid) < (x.key, x.nid):
                    x, y = y, x
                self._link(y, x)
                d += 1
            by_deg[d] = x
        self.roots = list(by_deg.values())
        self.min_root = None
        for r in self.roots:
            if self.min_root is None or (r.key, r.nid) < (self.min_root.key, self.min_root.nid):
                self.min_root = r

    def popmin(self):
        z = self.min_root
        if z is None:
            return -1
        z.alive = False
        # 最小根的孩子全部升為根（清 mark：根不該帶標記）
        for c in z.children:
            c.parent = None
            c.mark = False
        new_roots = [r for r in self.roots if r is not z]
        new_roots.extend(z.children)
        self.roots = new_roots
        z.children = []
        self._consolidate()
        return z.key

    def _cut(self, x):
        p = x.parent
        p.children.remove(x)
        p.deg -= 1
        x.parent = None
        x.mark = False
        self.roots.append(x)

    def _cascading_cut(self, y):
        z = y.parent
        if z is None:
            return
        if not y.mark:
            y.mark = True
        else:
            self._cut(y)
            self._cascading_cut(z)

    def dec(self, nid, nk):
        if not (0 <= nid < len(self.nodes)):
            return
        node = self.nodes[nid]
        if not node.alive or nk > node.key:
            return
        node.key = nk
        p = node.parent
        # 堆序是字典序 (key, nid)：與父節點「平手但 id 更小」也算破壞堆序
        if p is not None and (node.key, node.nid) < (p.key, p.nid):
            self._cut(node)
            self._cascading_cut(p)
        if self.min_root is None or (node.key, node.nid) < (self.min_root.key, self.min_root.nid):
            self.min_root = node

    def peek_min(self):
        return -1 if self.min_root is None else self.min_root.key

    def size(self):
        return sum(1 for nd in self.nodes if nd.alive)


class RefHeap:
    """線性掃描參考實作（僅供內建測試對拍）。"""

    def __init__(self):
        self.alive = {}
        self.cnt = 0

    def push(self, key):
        self.alive[self.cnt] = key
        self.cnt += 1
        return self.cnt - 1

    def popmin(self):
        if not self.alive:
            return -1
        best_id = min(self.alive, key=lambda i: (self.alive[i], i))
        v = self.alive.pop(best_id)
        return v

    def peek_min(self):
        if not self.alive:
            return -1
        return min(self.alive.values())

    def dec(self, nid, nk):
        if nid in self.alive and nk <= self.alive[nid]:
            self.alive[nid] = nk

    def size(self):
        return len(self.alive)


def run_io(data):
    tokens = data.split()
    pos = 0
    out = []

    def nxt():
        nonlocal pos
        if pos >= len(tokens):
            return None
        t = tokens[pos]
        pos += 1
        return t

    def nxt_int():
        t = nxt()
        if t is None:
            return None
        try:
            return int(t)
        except ValueError:
            return None

    heap = FibonacciHeap()
    while True:
        cmd = nxt()
        if cmd is None:
            break
        if cmd == 'push':
            k = nxt_int()
            if k is None:
                break
            out.append(str(heap.push(k)))
        elif cmd == 'popmin':
            out.append(str(heap.popmin()))
        elif cmd == 'min':
            out.append(str(heap.peek_min()))
        elif cmd == 'dec':
            i = nxt_int()
            k = nxt_int()
            if i is None or k is None:
                break
            heap.dec(i, k)
        elif cmd == 'size':
            out.append(str(heap.size()))
        else:
            break
    return ''.join(line + '\n' for line in out)


def run_tests():
    # 1) 基本序列（期望值以實跑為準，先驗證再寫死）
    h = FibonacciHeap()
    assert h.push(5) == 0
    assert h.push(3) == 1
    assert h.push(8) == 2
    assert h.peek_min() == 3
    assert h.size() == 3
    assert h.popmin() == 3
    assert h.peek_min() == 5
    h.dec(2, 1)
    assert h.peek_min() == 1
    assert h.popmin() == 1
    assert h.popmin() == 5
    assert h.popmin() == -1
    assert h.size() == 0

    # 2) 空堆行為
    h2 = FibonacciHeap()
    assert h2.popmin() == -1
    assert h2.peek_min() == -1
    assert h2.size() == 0
    h2.dec(99, 0)  # 非法 id 忽略
    assert h2.size() == 0

    # 3) 大量重複鍵值 + dec 打通 cut / cascading cut 路徑
    h3 = FibonacciHeap()
    ids = [h3.push(i % 7) for i in range(50)]
    assert h3.size() == 50
    for i in range(0, 50, 5):
        h3.dec(ids[i], -i)   # 負值觸發剪枝
    assert h3.peek_min() == -45
    vals = [h3.popmin() for _ in range(50)]
    expect = [(-i if i % 5 == 0 else i % 7) for i in range(50)]
    assert sorted(vals) == sorted(expect), (vals, expect)
    assert h3.size() == 0

    # 4) 與線性掃描參考實作對拍（隨機操作流，小鍵值域製造平手）
    state = 20261111

    def rnd(mm):
        nonlocal state
        state = (state * 1103515245 + 12345) & 0x7FFFFFFF
        return state % mm

    for _ in range(300):
        h = FibonacciHeap()
        ref = RefHeap()
        ops = rnd(200) + 30
        for _ in range(ops):
            op = rnd(10)
            if op < 4:
                k = rnd(20)
                a = h.push(k)
                b = ref.push(k)
                assert a == b
            elif op < 6:
                assert h.popmin() == ref.popmin()
            elif op < 7:
                assert h.peek_min() == ref.peek_min()
            elif op < 9:
                if h.cnt:
                    nid = rnd(h.cnt)
                    k = rnd(20)
                    h.dec(nid, k)
                    ref.dec(nid, k)
            else:
                assert h.size() == ref.size()
        # 收尾：逐一彈空，序列必須一致
        while ref.alive:
            assert h.popmin() == ref.popmin()
        assert h.popmin() == -1 == ref.popmin()


def main():
    data = sys.stdin.read()
    if data.strip():
        # 二進位輸出：避免 Windows 文字模式把 \n 轉成 \r\n（與 C++ 版不一致）
        sys.stdout.buffer.write(run_io(data).encode('utf-8'))
    else:
        run_tests()
        sys.stdout.buffer.write(b'all tests passed\n')


if __name__ == '__main__':
    main()
