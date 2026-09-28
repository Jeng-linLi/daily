"""線段樹（區間和 + 懶標記區間加）

題意：給定數組 a，支持兩種操作（下標 0-based，區間爲閉區間 [l, r]）：
    query(l, r)      求 a[l..r] 的元素和
    add(l, r, v)     把 a[l..r] 每個元素都加上 v
暴力做法單次 O(n)，m 次操作就是 O(nm)；線段樹把兩種操作都降到 O(log n)。

思路：
    線段樹是一棵描述「區間」的二叉樹：
      - 葉子存單個元素，內部節點存左右兒子區間的合併結果（這裡是和）。
      - 區間查詢：把目標區間拆成 O(log n) 個節點，命中完全覆蓋的節點就直接返回。
      - 區間修改：若整段被覆蓋，只更新該節點的 sum 並打上懶標記 lazy，不再往下遞歸；
        下次訪問必須下推（push_down）時再把標記分給兩個兒子。
      懶標記的本質是「延遲執行」：只要沒人問細節，就不必真的改到葉子。
    用數組（堆式存儲）實現：節點 p 的左兒子 2p、右兒子 2p+1，空間開 4n。

輸入格式（stdin）：
    第一行 n q（數組長度、操作條數）
    第二行 n 個整數（初始數組）
    接下來 q 行：1 l r（查詢 sum）/ 2 l r v（區間加 v）
輸出格式（stdout）：每個查詢操作輸出一行區間和
無 stdin 輸入時運行內置斷言測試。
"""

import sys
from typing import List


class SegmentTree:
    """區間和線段樹，支持區間加與區間求和，均爲 O(log n)。"""

    def __init__(self, arr: List[int]) -> None:
        self.n = len(arr)
        self.sum = [0] * (4 * self.n + 5)   # sum[p]：節點 p 對應區間的元素和
        self.lazy = [0] * (4 * self.n + 5)  # lazy[p]：待下推給子孫的「整體加多少」
        if self.n > 0:
            self._build(1, 0, self.n - 1, arr)

    def _build(self, p: int, l: int, r: int, arr: List[int]) -> None:
        """自底向上建樹，O(n)。"""
        if l == r:
            self.sum[p] = arr[l]
            return
        mid = (l + r) // 2
        self._build(2 * p, l, mid, arr)
        self._build(2 * p + 1, mid + 1, r, arr)
        self.sum[p] = self.sum[2 * p] + self.sum[2 * p + 1]

    def _apply(self, p: int, l: int, r: int, v: int) -> None:
        """給節點 p 整段加 v：更新 sum，並累積懶標記。"""
        self.sum[p] += v * (r - l + 1)
        self.lazy[p] += v

    def _push_down(self, p: int, l: int, r: int) -> None:
        """把 p 的懶標記下推給兩個兒子，只在需要深入時使用。"""
        if self.lazy[p] == 0 or l == r:
            return
        mid = (l + r) // 2
        v = self.lazy[p]
        self._apply(2 * p, l, mid, v)
        self._apply(2 * p + 1, mid + 1, r, v)
        self.lazy[p] = 0

    def add(self, ql: int, qr: int, v: int) -> None:
        """區間加：a[ql..qr] += v，O(log n)。"""
        if self.n == 0:
            return
        self._add(1, 0, self.n - 1, ql, qr, v)

    def _add(self, p: int, l: int, r: int, ql: int, qr: int, v: int) -> None:
        if ql <= l and r <= qr:      # 完全覆蓋，打標記後直接返回
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
        """區間求和：返回 a[ql..qr] 的和，O(log n)。"""
        if self.n == 0:
            return 0
        return self._query(1, 0, self.n - 1, ql, qr)

    def _query(self, p: int, l: int, r: int, ql: int, qr: int) -> int:
        if ql <= l and r <= qr:      # 完全覆蓋，直接返回整段和
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
    """按統一輸入輸出格式處理 stdin 數據。"""
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
    assert st.query(2, 2) == 5       # 單點查詢
    st.add(1, 3, 10)                 # [1, 13, 15, 17, 9, 11]
    assert st.query(1, 3) == 45
    assert st.query(0, 5) == 66      # 36 + 3 * 10
    assert st.query(0, 0) == 1       # 區間外的點不受影響
    assert st.query(4, 5) == 20

    # 懶標記疊加與部分覆蓋交叉驗證
    st = SegmentTree([0] * 5)
    # 結果數組爲 [2, 5, 5, 5, 2]
    st.add(0, 4, 2)
    st.add(1, 3, 3)
    assert st.query(0, 4) == 19       # 2*5 + 3*3
    assert st.query(0, 0) == 2
    assert st.query(1, 3) == 15       # 5 + 5 + 5
    assert st.query(4, 4) == 2

    # 單元素與空數組邊界
    st1 = SegmentTree([42])
    assert st1.query(0, 0) == 42
    st1.add(0, 0, -40)
    assert st1.query(0, 0) == 2
    st0 = SegmentTree([])
    assert st0.query(0, 0) == 0

    # 與樸素數組對拍
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
