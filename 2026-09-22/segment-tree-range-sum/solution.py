"""线段树（区间和 + 懒标记区间加）

题意：给定数组 a，支持两种操作（下标 0-based，区间为闭区间 [l, r]）：
    query(l, r)      求 a[l..r] 的元素和
    add(l, r, v)     把 a[l..r] 每个元素都加上 v
暴力做法单次 O(n)，m 次操作就是 O(nm)；线段树把两种操作都降到 O(log n)。

思路：
    线段树是一棵描述「区间」的二叉树：
      - 叶子存单个元素，内部节点存左右儿子区间的合并结果（这里是和）。
      - 区间查询：把目标区间拆成 O(log n) 个节点，命中完全覆盖的节点就直接返回。
      - 区间修改：若整段被覆盖，只更新该节点的 sum 并打上懒标记 lazy，不再往下递归；
        下次访问必须下推（push_down）时再把标记分给两个儿子。
      懒标记的本质是「延迟执行」：只要没人问细节，就不必真的改到叶子。
    用数组（堆式存储）实现：节点 p 的左儿子 2p、右儿子 2p+1，空间开 4n。

输入格式（stdin）：
    第一行 n q（数组长度、操作条数）
    第二行 n 个整数（初始数组）
    接下来 q 行：1 l r（查询 sum）/ 2 l r v（区间加 v）
输出格式（stdout）：每个查询操作输出一行区间和
无 stdin 输入时运行内置断言测试。
"""

import sys
from typing import List


class SegmentTree:
    """区间和线段树，支持区间加与区间求和，均为 O(log n)。"""

    def __init__(self, arr: List[int]) -> None:
        self.n = len(arr)
        self.sum = [0] * (4 * self.n + 5)   # sum[p]：节点 p 对应区间的元素和
        self.lazy = [0] * (4 * self.n + 5)  # lazy[p]：待下推给子孙的「整体加多少」
        if self.n > 0:
            self._build(1, 0, self.n - 1, arr)

    def _build(self, p: int, l: int, r: int, arr: List[int]) -> None:
        """自底向上建树，O(n)。"""
        if l == r:
            self.sum[p] = arr[l]
            return
        mid = (l + r) // 2
        self._build(2 * p, l, mid, arr)
        self._build(2 * p + 1, mid + 1, r, arr)
        self.sum[p] = self.sum[2 * p] + self.sum[2 * p + 1]

    def _apply(self, p: int, l: int, r: int, v: int) -> None:
        """给节点 p 整段加 v：更新 sum，并累积懒标记。"""
        self.sum[p] += v * (r - l + 1)
        self.lazy[p] += v

    def _push_down(self, p: int, l: int, r: int) -> None:
        """把 p 的懒标记下推给两个儿子，只在需要深入时使用。"""
        if self.lazy[p] == 0 or l == r:
            return
        mid = (l + r) // 2
        v = self.lazy[p]
        self._apply(2 * p, l, mid, v)
        self._apply(2 * p + 1, mid + 1, r, v)
        self.lazy[p] = 0

    def add(self, ql: int, qr: int, v: int) -> None:
        """区间加：a[ql..qr] += v，O(log n)。"""
        if self.n == 0:
            return
        self._add(1, 0, self.n - 1, ql, qr, v)

    def _add(self, p: int, l: int, r: int, ql: int, qr: int, v: int) -> None:
        if ql <= l and r <= qr:      # 完全覆盖，打标记后直接返回
            self._apply(p, l, r, v)
            return
        self._push_down(p, l, r)
        mid = (l + r) // 2
        if ql <= mid:
            self._add(2 * p, l, mid, ql, qr, v)
        if qr > mid:
            self._add(2 * p + 1, mid + 1, r, ql, qr, v)
        self.sum[p] = self.sum[2 * p] + self.sum[2 * p + 1]

    def query(self, ql: int, qr: int) -> int:
        """区间求和：返回 a[ql..qr] 的和，O(log n)。"""
        if self.n == 0:
            return 0
        return self._query(1, 0, self.n - 1, ql, qr)

    def _query(self, p: int, l: int, r: int, ql: int, qr: int) -> int:
        if ql <= l and r <= qr:      # 完全覆盖，直接返回整段和
            return self.sum[p]
        self._push_down(p, l, r)
        mid = (l + r) // 2
        ans = 0
        if ql <= mid:
            ans += self._query(2 * p, l, mid, ql, qr)
        if qr > mid:
            ans += self._query(2 * p + 1, mid + 1, r, ql, qr)
        return ans


def run_io(data: str) -> None:
    """按统一输入输出格式处理 stdin 数据。"""
    tokens = data.split()
    n = int(tokens[0])
    q = int(tokens[1])
    arr = [int(t) for t in tokens[2:2 + n]]
    st = SegmentTree(arr)
    idx = 2 + n
    out: List[str] = []
    for _ in range(q):
        kind = int(tokens[idx]); idx += 1
        l = int(tokens[idx]); idx += 1
        r = int(tokens[idx]); idx += 1
        if kind == 1:
            out.append(str(st.query(l, r)))
        else:
            v = int(tokens[idx]); idx += 1
            st.add(l, r, v)
    print("\n".join(out))


def run_tests() -> None:
    arr = [1, 3, 5, 7, 9, 11]
    st = SegmentTree(arr)
    assert st.query(0, 5) == 36
    assert st.query(1, 3) == 15      # 3 + 5 + 7
    assert st.query(2, 2) == 5       # 单点查询
    st.add(1, 3, 10)                 # [1, 13, 15, 17, 9, 11]
    assert st.query(1, 3) == 45
    assert st.query(0, 5) == 66      # 36 + 3 * 10
    assert st.query(0, 0) == 1       # 区间外的点不受影响
    assert st.query(4, 5) == 20

    # 懒标记叠加与部分覆盖交叉验证
    st = SegmentTree([0] * 5)
    # 结果数组为 [2, 5, 5, 5, 2]
    st.add(0, 4, 2)
    st.add(1, 3, 3)
    assert st.query(0, 4) == 19       # 2*5 + 3*3
    assert st.query(0, 0) == 2
    assert st.query(1, 3) == 15       # 5 + 5 + 5
    assert st.query(4, 4) == 2

    # 单元素与空数组边界
    st1 = SegmentTree([42])
    assert st1.query(0, 0) == 42
    st1.add(0, 0, -40)
    assert st1.query(0, 0) == 2
    st0 = SegmentTree([])
    assert st0.query(0, 0) == 0

    # 与朴素数组对拍
    import random
    base = [random.randint(-20, 20) for _ in range(40)]
    st2 = SegmentTree(base)
    for _ in range(300):
        l = random.randint(0, 39)
        r = random.randint(l, 39)
        if random.random() < 0.5:
            assert st2.query(l, r) == sum(base[l:r + 1])
        else:
            v = random.randint(-10, 10)
            st2.add(l, r, v)
            for k in range(l, r + 1):
                base[k] += v
    print("all tests passed")


if __name__ == "__main__":
    raw = sys.stdin.read() if not sys.stdin.isatty() else ""
    if raw.strip():
        run_io(raw)
    else:
        run_tests()
