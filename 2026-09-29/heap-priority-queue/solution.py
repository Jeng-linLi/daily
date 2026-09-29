"""二元堆與優先佇列（堆排序 / Top-K / 合併 K 個有序陣列 / 雙堆維護中位數）

題意：
    輸入一個無序陣列 a（長度 n）、一個整數 k、以及 m 個**已排序**的陣列，輸出：
      1. a 的堆排序結果（升序）；
      2. a 中最小的 k 個元素（升序；k <= 0 輸出空行，k >= n 輸出全部）；
      3. 合併 m 個有序陣列後的結果（升序）。

思路：
    二元堆（binary heap）是一棵用**陣列**存的完全二元樹：下標 i 的父節點是 (i-1)//2，
    左右孩子是 2i+1、2i+2。完全二元樹的高度是 floor(log2 n)，所以「上浮 / 下沉」
    最多走 log n 層 —— 這就是 push / pop 都是 O(log n) 的來源，而 top 是 O(1)。
    堆只維護**偏序**（父節點優於兩個孩子），不維護全序，所以它比有序陣列便宜得多：
    只取最大值/最小值時，沒有必要把整個集合排好。

    關鍵操作：
      - sift_up（上浮）：新元素放在末尾，若比父節點更優就交換，直到恢復堆性質；
      - sift_down（下沉）：彈出堆頂後把末尾元素搬到堆頂，再與「更優的那個孩子」交換；
      - heapify（自底向上建堆）：從最後一個非葉節點 (n-2)//2 開始往前逐個 sift_down。
        **注意它是 O(n) 而不是 O(n log n)**：大部分節點深度很淺，級數求和後收斂到 2n。

    三個應用：
      1. 堆排序：建堆後反覆彈出堆頂。本題用最小堆，彈出順序天然升序；
         另外給出原地版本（建最大堆，把堆頂換到末尾再下沉），額外空間 O(1)（不穩定排序）。
      2. Top-K 最小：維護一個**大小不超過 k 的最大堆**——堆頂是當前 k 個候選裏最大的，
         來了新元素只要比堆頂小就替換它。時間 O(n log k)，比排序的 O(n log n) 更划算，
         n 很大而 k 很小時優勢明顯（這也是搜索引擎取 Top-K 結果的標準做法）。
      3. 合併 K 個有序陣列：把每個陣列的**當前頭元素**放進最小堆，彈出最小的那個，
         再把該陣列的下一個元素補進去。時間 O(N log K)，N 爲總元素數。
         擴展一下就是作業系統的「多路歸併」與定時器 / 任務調度的優先佇列。

    彩蛋：`MedianMaintainer` 用一個最大堆 + 一個最小堆維護數據流的中位數
    （大堆存較小的一半、小堆存較大的一半，並保持兩堆大小差 <= 1）。
    插入 O(log n)，查詢 O(1)，是「即時中位數」類題目的標準解法。

輸入格式（stdin，數字按空白分隔即可）：
    n k
    a1 a2 ... an            （n = 0 時這行省略）
    m                       （接下來有多少個有序陣列）
    L1 x1 x2 ... x_L1       （第 1 個有序陣列；重複 m 行；m = 0 時全部省略）
    L2 ...
輸出格式（stdout）：
    第 1 行：a 的堆排序結果（升序，空格分隔）
    第 2 行：a 中最小的 k 個元素（升序）
    第 3 行：合併後的結果（升序）
無 stdin 輸入時運行內置斷言測試並輸出 `all tests passed`。
"""

import heapq
import random
import sys
from typing import Iterator, List, Optional


class BinaryHeap:
    """陣列實現的二元堆。max_heap=True 時爲最大堆（堆頂最大），否則爲最小堆。"""

    def __init__(self, values: Optional[List[int]] = None, max_heap: bool = False) -> None:
        self.data: List[int] = []
        self.max_heap = max_heap
        if values:
            self.heapify(values)

    def _better(self, a: int, b: int) -> bool:
        """a 是否應該排在 b 上面（更靠近堆頂）。"""
        return a > b if self.max_heap else a < b

    def __len__(self) -> int:
        return len(self.data)

    def top(self) -> int:
        """堆頂元素，O(1)。空堆時拋 IndexError。"""
        return self.data[0]

    def push(self, x: int) -> None:
        """插入，O(log n)。"""
        self.data.append(x)
        self._sift_up(len(self.data) - 1)

    def pop(self) -> int:
        """彈出堆頂，O(log n)。"""
        if not self.data:
            raise IndexError("pop from empty heap")
        top = self.data[0]
        last = self.data.pop()
        if self.data:                      # 把末尾元素搬到堆頂再下沉
            self.data[0] = last
            self._sift_down(0)
        return top

    def replace_top(self, x: int) -> int:
        """彈出堆頂並插入 x（只做一次下沉，比 pop+push 快一倍）。"""
        top = self.data[0]
        self.data[0] = x
        self._sift_down(0)
        return top

    def heapify(self, values: List[int]) -> None:
        """自底向上建堆，O(n)（不是 O(n log n)）。"""
        self.data = list(values)
        for i in range((len(self.data) - 2) // 2, -1, -1):   # 從最後一個非葉節點往前
            self._sift_down(i)

    def _sift_up(self, i: int) -> None:
        """上浮：與父節點比較，更優則交換。"""
        data = self.data
        x = data[i]
        while i > 0:
            parent = (i - 1) // 2
            if not self._better(x, data[parent]):
                break
            data[i] = data[parent]
            i = parent
        data[i] = x

    def _sift_down(self, i: int) -> None:
        """下沉：與更優的那個孩子交換。"""
        data = self.data
        n = len(data)
        x = data[i]
        while True:
            left = 2 * i + 1
            if left >= n:
                break
            right = left + 1
            child = left
            if right < n and self._better(data[right], data[left]):
                child = right
            if not self._better(data[child], x):
                break
            data[i] = data[child]
            i = child
        data[i] = x


def is_valid_heap(h: BinaryHeap) -> bool:
    """檢查堆性質：每個父節點都不劣於它的孩子。"""
    d = h.data
    for i in range(1, len(d)):
        parent = (i - 1) // 2
        if h._better(d[i], d[parent]):
            return False
    return True


# ---------------- 應用 1：堆排序 ----------------

def heap_sort(values: List[int]) -> List[int]:
    """堆排序（用最小堆逐個彈出），升序。時間 O(n log n)，空間 O(n)。"""
    h = BinaryHeap(values, max_heap=False)
    return [h.pop() for _ in range(len(h))]


def heap_sort_in_place(values: List[int]) -> List[int]:
    """原地堆排序：建最大堆後把堆頂換到末尾。時間 O(n log n)，額外空間 O(1)。

    注意：堆排序**不是穩定排序**（相等元素的相對順序可能改變）。
    """
    a = list(values)
    n = len(a)

    def sift_down(i: int, size: int) -> None:
        x = a[i]
        while True:
            left = 2 * i + 1
            if left >= size:
                break
            right = left + 1
            child = left
            if right < size and a[right] > a[left]:
                child = right
            if a[child] <= x:
                break
            a[i] = a[child]
            i = child
        a[i] = x

    for i in range((n - 2) // 2, -1, -1):      # 1) 建最大堆
        sift_down(i, n)
    for end in range(n - 1, 0, -1):            # 2) 堆頂（最大）換到末尾，縮小堆
        a[0], a[end] = a[end], a[0]
        sift_down(0, end)
    return a


# ---------------- 應用 2：Top-K 最小 ----------------

def top_k_smallest(values: List[int], k: int) -> List[int]:
    """最小的 k 個元素（升序）。維護大小爲 k 的最大堆，時間 O(n log k)，空間 O(k)。"""
    if k <= 0:
        return []
    h = BinaryHeap(max_heap=True)
    for v in values:
        if len(h) < k:
            h.push(v)
        elif v < h.top():
            h.replace_top(v)                   # 比當前第 k 小還小 → 擠掉堆頂
    return sorted(h.data)


# ---------------- 應用 3：合併 K 個有序陣列 ----------------

def merge_k_sorted(arrays: List[List[int]]) -> List[int]:
    """合併若干個已排序陣列。時間 O(N log K)，空間 O(K)。

    堆裏放 (值, 陣列編號, 元素下標)，後兩項保證值相等時順序也確定（與 C++ 版一致）。
    """
    heap: List[tuple] = []
    for i, arr in enumerate(arrays):
        if arr:                                # 空陣列直接跳過
            heapq.heappush(heap, (arr[0], i, 0))
    out: List[int] = []
    while heap:
        val, i, j = heapq.heappop(heap)
        out.append(val)
        if j + 1 < len(arrays[i]):
            heapq.heappush(heap, (arrays[i][j + 1], i, j + 1))
    return out


# ---------------- 彩蛋：雙堆維護數據流中位數 ----------------

class MedianMaintainer:
    """一個最大堆存較小的一半，一個最小堆存較大的一半，兩堆大小差 <= 1。

    插入 O(log n)，查詢中位數 O(1)。
    """

    def __init__(self) -> None:
        self.lo: List[int] = []                # 最大堆：存較小的一半（用負值模擬）
        self.hi: List[int] = []                # 最小堆：存較大的一半

    def add(self, x: int) -> None:
        if not self.lo or x <= -self.lo[0]:
            heapq.heappush(self.lo, -x)
        else:
            heapq.heappush(self.hi, x)
        # 重新平衡：lo 允許比 hi 多一個，反過來不行
        if len(self.lo) > len(self.hi) + 1:
            heapq.heappush(self.hi, -heapq.heappop(self.lo))
        elif len(self.hi) > len(self.lo):
            heapq.heappush(self.lo, -heapq.heappop(self.hi))

    def median(self) -> Optional[float]:
        if not self.lo and not self.hi:
            return None
        if len(self.lo) == len(self.hi):
            return (-self.lo[0] + self.hi[0]) / 2.0
        return float(-self.lo[0])


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
    k = _next_int(it)
    a = [_next_int(it) for _ in range(max(n, 0))]
    m = _next_int(it)
    arrays: List[List[int]] = []
    for _ in range(max(m, 0)):
        length = _next_int(it)
        arrays.append([_next_int(it) for _ in range(max(length, 0))])

    print(" ".join(str(v) for v in heap_sort(a)))            # n = 0 時輸出空行
    print(" ".join(str(v) for v in top_k_smallest(a, k)))
    print(" ".join(str(v) for v in merge_k_sorted(arrays)))


def run_tests() -> None:
    # README 示例：a = [5, 3, 8, 1, 4]，k = 3，兩個有序陣列 [1,4,7] 與 [2,3,9]
    a = [5, 3, 8, 1, 4]
    assert heap_sort(a) == [1, 3, 4, 5, 8]
    assert heap_sort_in_place(a) == [1, 3, 4, 5, 8]
    assert top_k_smallest(a, 3) == [1, 3, 4]
    assert merge_k_sorted([[1, 4, 7], [2, 3, 9]]) == [1, 2, 3, 4, 7, 9]

    # 堆的基本性質
    h = BinaryHeap([3, 1, 6, 5, 2, 4], max_heap=False)
    assert is_valid_heap(h)
    assert h.top() == 1
    assert [h.pop() for _ in range(6)] == [1, 2, 3, 4, 5, 6]
    hm = BinaryHeap([3, 1, 6, 5, 2, 4], max_heap=True)
    assert is_valid_heap(hm)
    assert hm.top() == 6
    assert [hm.pop() for _ in range(6)] == [6, 5, 4, 3, 2, 1]

    # 邊界：空堆 / 單元素 / 重複值 / 負數
    empty = BinaryHeap()
    assert len(empty) == 0
    assert heap_sort([]) == []
    assert heap_sort_in_place([]) == []
    assert heap_sort([7]) == [7]
    assert heap_sort([2, 2, 2]) == [2, 2, 2]
    assert is_valid_heap(BinaryHeap([2, 2, 2]))
    assert heap_sort([-3, 0, -1, 5]) == [-3, -1, 0, 5]
    assert heap_sort([1, 2, 3, 4, 5]) == [1, 2, 3, 4, 5]      # 已升序
    assert heap_sort([5, 4, 3, 2, 1]) == [1, 2, 3, 4, 5]      # 逆序

    # Top-K 邊界
    assert top_k_smallest([5, 3, 8, 1, 4], 0) == []
    assert top_k_smallest([5, 3, 8, 1, 4], -2) == []          # k 爲負：空
    assert top_k_smallest([5, 3, 8, 1, 4], 100) == [1, 3, 4, 5, 8]   # k 超過 n
    assert top_k_smallest([], 3) == []
    assert top_k_smallest([4, 4, 1, 4], 2) == [1, 4]          # 有重複值

    # 合併：空陣列 / 全空 / 單個陣列
    assert merge_k_sorted([]) == []
    assert merge_k_sorted([[], []]) == []
    assert merge_k_sorted([[1, 2, 3]]) == [1, 2, 3]
    assert merge_k_sorted([[], [1, 2], []]) == [1, 2]
    assert merge_k_sorted([[1, 1], [1, 1]]) == [1, 1, 1, 1]   # 值全相等：按下標順序出

    # 中位數維護器
    mm = MedianMaintainer()
    assert mm.median() is None
    for x in [5, 2, 9, 1, 7]:
        mm.add(x)
    assert abs(mm.median() - 5.0) < 1e-9                       # [1,2,5,7,9] → 5
    mm2 = MedianMaintainer()
    for x in [1, 2, 3, 4]:
        mm2.add(x)
    assert abs(mm2.median() - 2.5) < 1e-9                      # 偶數個取中間兩個的平均

    random.seed(20260929)

    # 隨機對拍：全部與 sorted() 的結果比對
    for _ in range(400):
        n = random.randint(0, 40)
        a = [random.randint(-10, 10) for _ in range(n)]        # 值域小 → 大量重複
        assert heap_sort(a) == sorted(a)
        assert heap_sort_in_place(a) == sorted(a)
        k = random.randint(-1, n + 3)
        assert top_k_smallest(a, k) == sorted(a)[:max(k, 0)]

        # 堆操作與排序結果一致
        h = BinaryHeap(a)
        assert is_valid_heap(h)
        assert [h.pop() for _ in range(len(h))] == sorted(a)

        # 合併隨機個有序陣列
        m = random.randint(0, 5)
        arrays = []
        for _ in range(m):
            size = random.randint(0, 6)
            arrays.append(sorted(random.randint(-10, 10) for _ in range(size)))
        merged = merge_k_sorted(arrays)
        flat = [x for arr in arrays for x in arr]
        assert merged == sorted(flat)

        # 中位數維護器與「排序後取中位」逐前綴比對
        mmr = MedianMaintainer()
        seen: List[int] = []
        for x in a:
            mmr.add(x)
            seen.append(x)
            seen.sort()
            size = len(seen)
            want = (seen[size // 2] if size % 2 == 1
                    else (seen[size // 2 - 1] + seen[size // 2]) / 2.0)
            assert abs(mmr.median() - want) < 1e-9

    print("all tests passed")


if __name__ == "__main__":
    raw = sys.stdin.read() if not sys.stdin.isatty() else ""
    if raw.strip():
        run_io(raw)
    else:
        run_tests()
