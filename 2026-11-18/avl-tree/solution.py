#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""AVL 樹（自平衡二元搜尋樹）與順序統計。

命令（IO 模式，stdin 第一個 token 為操作數 n，之後 n 個操作）：
  i x  插入 x          -> 輸出插入後元素總數
  e x  刪除一個 x      -> 刪除後元素總數
  f x  是否存在 x      -> 1 / 0
  r x  嚴格小於 x 的個數
  k x  第 x 小（1-based），越界輸出 -1
  p x  嚴格小於 x 的最大者，無則 -1
  s x  嚴格大於 x 的最小者，無則 -1
  o    中序輸出全部元素（含重複），空則輸出 -

這是一個實際可用的有序索引結構：insert/erase/find/rank/kth/pred/succ 全 O(log n)。
stdin 無內容時跑內建對拍（與排序陣列 + bisect 參考實作比較），輸出 all tests passed。
"""

import sys
import re
from bisect import bisect_left, bisect_right, insort

_INT = re.compile(r'^[+-]?[0-9]+$')


class AVL:
    """以陣列模擬節點的 AVL 樹；節點 0 為空 sentinel。"""

    def __init__(self):
        self.key = [0]
        self.cnt = [0]
        self.h = [0]
        self.sz = [0]
        self.lc = [0]
        self.rc = [0]

    def _new(self, k):
        self.key.append(k)
        self.cnt.append(1)
        self.h.append(1)
        self.sz.append(1)
        self.lc.append(0)
        self.rc.append(0)
        return len(self.key) - 1

    def _pull(self, v):
        if v == 0:
            return
        self.h[v] = 1 + max(self.h[self.lc[v]], self.h[self.rc[v]])
        self.sz[v] = self.cnt[v] + self.sz[self.lc[v]] + self.sz[self.rc[v]]

    def _rot_right(self, y):
        x = self.lc[y]
        self.lc[y] = self.rc[x]
        self.rc[x] = y
        self._pull(y)
        self._pull(x)
        return x

    def _rot_left(self, x):
        y = self.rc[x]
        self.rc[x] = self.lc[y]
        self.lc[y] = x
        self._pull(x)
        self._pull(y)
        return y

    def _rebalance(self, v):
        if v == 0:
            return 0
        self._pull(v)
        bf = self.h[self.lc[v]] - self.h[self.rc[v]]
        if bf > 1:
            if self.h[self.lc[self.lc[v]]] < self.h[self.rc[self.lc[v]]]:
                self.lc[v] = self._rot_left(self.lc[v])
            return self._rot_right(v)
        if bf < -1:
            if self.h[self.rc[self.rc[v]]] < self.h[self.lc[self.rc[v]]]:
                self.rc[v] = self._rot_right(self.rc[v])
            return self._rot_left(v)
        return v

    def insert(self, v, k):
        if v == 0:
            return self._new(k)
        if k == self.key[v]:
            self.cnt[v] += 1
        elif k < self.key[v]:
            self.lc[v] = self.insert(self.lc[v], k)
        else:
            self.rc[v] = self.insert(self.rc[v], k)
        return self._rebalance(v)

    def _pop_min(self, v):
        """移除 v 子樹的最小節點，回傳 (新根, key, cnt)。"""
        if self.lc[v] == 0:
            return self.rc[v], self.key[v], self.cnt[v]
        self.lc[v], k, c = self._pop_min(self.lc[v])
        return self._rebalance(v), k, c

    def erase(self, v, k):
        if v == 0:
            return 0
        if k < self.key[v]:
            self.lc[v] = self.erase(self.lc[v], k)
        elif k > self.key[v]:
            self.rc[v] = self.erase(self.rc[v], k)
        else:
            if self.cnt[v] > 1:
                self.cnt[v] -= 1
            elif self.lc[v] == 0 or self.rc[v] == 0:
                v = self.lc[v] or self.rc[v]
            else:
                self.rc[v], self.key[v], self.cnt[v] = self._pop_min(self.rc[v])
        return self._rebalance(v)

    def find(self, v, k):
        while v:
            if k == self.key[v]:
                return True
            v = self.lc[v] if k < self.key[v] else self.rc[v]
        return False

    def rank_less(self, v, k):
        """嚴格小於 k 的元素個數。"""
        ans = 0
        while v:
            if k <= self.key[v]:
                v = self.lc[v]
            else:
                ans += self.sz[self.lc[v]] + self.cnt[v]
                v = self.rc[v]
        return ans

    def kth(self, v, k):
        """第 k 小（1-based），越界回傳 None。"""
        if k < 1 or k > self.sz[v]:
            return None
        while v:
            ls = self.sz[self.lc[v]]
            if k <= ls:
                v = self.lc[v]
            elif k > ls + self.cnt[v]:
                k -= ls + self.cnt[v]
                v = self.rc[v]
            else:
                return self.key[v]
        return None

    def pred(self, v, k):
        """嚴格小於 k 的最大者，無則 None。"""
        ans = None
        while v:
            if self.key[v] < k:
                ans = self.key[v]
                v = self.rc[v]
            else:
                v = self.lc[v]
        return ans

    def succ(self, v, k):
        """嚴格大於 k 的最小者，無則 None。"""
        ans = None
        while v:
            if self.key[v] > k:
                ans = self.key[v]
                v = self.lc[v]
            else:
                v = self.rc[v]
        return ans

    def inorder(self, v, out):
        if v == 0:
            return
        self.inorder(self.lc[v], out)
        out.extend([self.key[v]] * self.cnt[v])
        self.inorder(self.rc[v], out)

    def check(self, v, lo, hi):
        """驗證 BST 順序、高度與 size 不變量，回傳高度。"""
        if v == 0:
            return 0
        assert lo is None or self.key[v] > lo
        assert hi is None or self.key[v] < hi
        lh = self.check(self.lc[v], lo, self.key[v])
        rh = self.check(self.rc[v], self.key[v], hi)
        assert abs(lh - rh) <= 1
        assert self.h[v] == 1 + max(lh, rh)
        assert self.sz[v] == self.cnt[v] + self.sz[self.lc[v]] + self.sz[self.rc[v]]
        return self.h[v]


# ----------------------------- 參考實作（排序陣列） -----------------------------

def ref_op(ref, op, arg):
    if op == 'i':
        insort(ref, arg)
        return str(len(ref))
    if op == 'e':
        j = bisect_left(ref, arg)
        if j < len(ref) and ref[j] == arg:
            ref.pop(j)
        return str(len(ref))
    if op == 'f':
        j = bisect_left(ref, arg)
        return '1' if j < len(ref) and ref[j] == arg else '0'
    if op == 'r':
        return str(bisect_left(ref, arg))
    if op == 'k':
        return str(ref[arg - 1]) if 1 <= arg <= len(ref) else '-1'
    if op == 'p':
        j = bisect_left(ref, arg)
        return str(ref[j - 1]) if j > 0 else '-1'
    if op == 's':
        j = bisect_right(ref, arg)
        return str(ref[j]) if j < len(ref) else '-1'
    return ' '.join(map(str, ref)) if ref else '-'


def fmt(x):
    return '-1' if x is None else str(x)


def run_ops(ops):
    """回傳 AVL 輸出與參考輸出（二者應逐行相同）。"""
    tree = AVL()
    root = 0
    ref = []
    got = []
    want = []
    for op, arg in ops:
        if op == 'i':
            root = tree.insert(root, arg)
            got.append(str(tree.sz[root]))
        elif op == 'e':
            root = tree.erase(root, arg)
            got.append(str(tree.sz[root]))
        elif op == 'f':
            got.append('1' if tree.find(root, arg) else '0')
        elif op == 'r':
            got.append(str(tree.rank_less(root, arg)))
        elif op == 'k':
            got.append(fmt(tree.kth(root, arg)))
        elif op == 'p':
            got.append(fmt(tree.pred(root, arg)))
        elif op == 's':
            got.append(fmt(tree.succ(root, arg)))
        else:
            seq = []
            tree.inorder(root, seq)
            got.append(' '.join(map(str, seq)) if seq else '-')
        want.append(ref_op(ref, op, arg))
    return got, want


# ------------------------------- IO 模式 -------------------------------

def run_io(data):
    toks = data.split()
    pos = 0
    out = []

    def nxt():
        nonlocal pos
        if pos >= len(toks):
            return None
        t = toks[pos]
        pos += 1
        return t

    nt = nxt()
    if nt is None or not _INT.match(nt):
        return ''
    n = int(nt)
    n = min(max(n, 0), 200000)
    tree = AVL()
    root = 0
    for _ in range(n):
        op = nxt()
        if op is None or op not in ('i', 'e', 'f', 'r', 'k', 'p', 's', 'o'):
            break
        arg = 0
        if op != 'o':
            at = nxt()
            if at is None or not _INT.match(at):
                break
            arg = int(at)
            arg = min(max(arg, -10**15), 10**15)
        if op == 'i':
            root = tree.insert(root, arg)
            out.append(str(tree.sz[root]))
        elif op == 'e':
            root = tree.erase(root, arg)
            out.append(str(tree.sz[root]))
        elif op == 'f':
            out.append('1' if tree.find(root, arg) else '0')
        elif op == 'r':
            out.append(str(tree.rank_less(root, arg)))
        elif op == 'k':
            out.append(fmt(tree.kth(root, arg)))
        elif op == 'p':
            out.append(fmt(tree.pred(root, arg)))
        elif op == 's':
            out.append(fmt(tree.succ(root, arg)))
        else:
            seq = []
            tree.inorder(root, seq)
            out.append(' '.join(map(str, seq)) if seq else '-')
    return ''.join(line + '\n' for line in out)


# ------------------------------ 內建測試 ------------------------------

def run_tests():
    # 固定情境：先跑過確認真實輸出
    fixed = [
        ('i', 5), ('i', 3), ('i', 7), ('i', 3), ('o', 0),
        ('f', 3), ('f', 4), ('r', 6), ('k', 2), ('k', 5),
        ('p', 6), ('s', 6), ('e', 3), ('o', 0), ('e', 3), ('o', 0),
    ]
    got, want = run_ops(fixed)
    assert got == want, (got, want)
    assert got == ['1', '2', '3', '4', '3 3 5 7', '1', '0', '3', '3', '-1',
                   '5', '7', '3', '3 5 7', '2', '5 7'], got

    # 隨機對拍：與排序陣列 + bisect 比較，並逐步驗證 AVL 不變量
    state = 20261118

    def rnd(m):
        nonlocal state
        state = (state * 1103515245 + 12345) & 0x7FFFFFFF
        return state % m

    for tc in range(160):
        ops = []
        for _ in range(120):
            op = 'ief rkps o'.replace(' ', '')[rnd(8)]
            arg = rnd(41) - 20
            ops.append((op, arg))
        tree = AVL()
        root = 0
        ref = []
        for op, arg in ops:
            if op == 'i':
                root = tree.insert(root, arg)
                insort(ref, arg)
                assert tree.sz[root] == len(ref)
            elif op == 'e':
                root = tree.erase(root, arg)
                j = bisect_left(ref, arg)
                if j < len(ref) and ref[j] == arg:
                    ref.pop(j)
                assert tree.sz[root] == len(ref)
            elif op in ('f', 'r', 'k', 'p', 's'):
                pass
            tree.check(root, None, None)
        got, want = run_ops(ops)
        assert got == want, (tc, ops, got, want)

    # 單調插入 / 刪除（最考驗旋轉）
    ops = [('i', i) for i in range(400)] + [('e', i) for i in range(400)]
    got, want = run_ops(ops)
    assert got == want


def main():
    data = sys.stdin.read()
    if data.strip():
        sys.stdout.buffer.write(run_io(data).encode('utf-8'))
    else:
        run_tests()
        sys.stdout.buffer.write(b'all tests passed\n')


if __name__ == '__main__':
    main()
