"""快速排序與快速選擇（隨機化 pivot + 三路分區 + Lomuto 分區）

題意：
    給定 n 個整數（允許重複、允許負數、允許爲空），要求：
      1) 輸出升序排序結果；
      2) 輸出第 k 小元素（k 爲 1-based，越界時輸出 -1）。

思路：
    快排是「分治 + 原地分區」：每輪挑一個 pivot，把區間切成
    「< pivot」「== pivot」「> pivot」三段，然後只對左右兩段遞歸。

    1) **Lomuto 分區**：維護指針 i，把 < pivot 的元素換到左邊。
       寫法最直觀，但遇到大量與 pivot 相等的元素時，兩段會嚴重不平衡，
       所有元素都相等時會退化成 O(n^2)。

    2) **三路分區（Dijkstra / 荷蘭國旗）**：維護 lt / i / gt 三個指針，
       a[lo:lt] < pivot，a[lt:i] == pivot，a[i:gt+1] 待定，a[gt+1:hi+1] > pivot。
       相等元素一次整段歸位，**全相同元素時是 O(n)**，這是它最大的價值。

    3) **隨機化 pivot**：固定挑首/尾元素時，已排序輸入會退化成 O(n^2)。
       隨機挑 pivot 讓任何輸入分布的期望時間都是 O(n log n)。
       實踐中還會「小區間（長度 < 16）改用插入排序」——常數更小。

    4) **遞歸深度控制**：先遞歸短的那半、用 while 循環處理長的那半，
       棧深度被壓到 O(log n)，不會因最壞情況爆棧。

    5) **快速選擇（Quickselect）**：排序只爲了拿第 k 個太浪費。
       分區後看 k 落在哪一段：落在 == pivot 段就直接返回，否則只遞歸
       包含 k 的那一段。平均 O(n)，最壞 O(n^2)（隨機化後幾乎不會碰到）。

    注意：本文件裏的排序是**原地**的，會修改傳入的列表。

輸入格式（stdin，所有數字按空白分隔即可）：
    n
    a1 a2 ... an        （n = 0 時這一行直接省略）
    k
輸出格式（stdout）：
    第 1 行：升序排序結果，空格分隔（n = 0 時輸出空行）
    第 2 行：第 k 小元素的值（1-based；k 越界或 n = 0 時輸出 -1）
無 stdin 輸入時運行內置斷言測試並輸出 `all tests passed`。
"""

import random
import sys
from typing import Iterator, List

SMALL = 16  # 小區間閾值：低於這個長度改用插入排序，常數更小


def insertion_sort(a: List[int], lo: int = 0, hi: int = -1) -> None:
    """對 a[lo:hi+1] 做插入排序。時間 O(len^2)，但小區間裡常數極小。"""
    if hi < 0:
        hi = len(a) - 1
    for i in range(lo + 1, hi + 1):
        x = a[i]
        j = i - 1
        while j >= lo and a[j] > x:
            a[j + 1] = a[j]
            j -= 1
        a[j + 1] = x


def lomuto_partition(a: List[int], lo: int, hi: int) -> int:
    """Lomuto 分區：以 a[hi] 爲 pivot，把 < pivot 的換到左邊，返回 pivot 最終下標。

    時間 O(hi-lo+1)，空間 O(1)。重複元素多時兩段會不平衡。
    """
    pivot = a[hi]
    i = lo
    for j in range(lo, hi):
        if a[j] < pivot:
            a[i], a[j] = a[j], a[i]
            i += 1
    a[i], a[hi] = a[hi], a[i]
    return i


def partition3(a: List[int], lo: int, hi: int, rng=None) -> tuple:
    """三路分區：把 a[lo:hi+1] 切成 < == > 三段，返回 (lt, gt)。

    a[lo:lt] < pivot，a[lt:gt+1] == pivot，a[gt+1:hi+1] > pivot。
    時間 O(hi-lo+1)，空間 O(1)。pivot 隨機挑選。
    """
    idx = rng.randint(lo, hi) if rng is not None else random.randint(lo, hi)
    pivot = a[idx]
    lt, i, gt = lo, lo, hi
    while i <= gt:
        if a[i] < pivot:
            a[lt], a[i] = a[i], a[lt]
            lt += 1
            i += 1
        elif a[i] > pivot:
            a[i], a[gt] = a[gt], a[i]
            gt -= 1          # 換過來的元素還沒看過，i 不動
        else:
            i += 1
    return lt, gt


def _quick_sort3(a: List[int], lo: int, hi: int, rng=None) -> None:
    """三路快排主過程：先遞歸短半邊，長半邊用循環，棧深度 O(log n)。"""
    while lo < hi:
        if hi - lo + 1 <= SMALL:
            insertion_sort(a, lo, hi)
            return
        lt, gt = partition3(a, lo, hi, rng)
        if lt - lo < hi - gt:          # 左段更短：遞歸左段，循環處理右段
            _quick_sort3(a, lo, lt - 1, rng)
            lo = gt + 1
        else:                          # 右段更短：遞歸右段，循環處理左段
            _quick_sort3(a, gt + 1, hi, rng)
            hi = lt - 1


def quick_sort(a: List[int]) -> List[int]:
    """三路快排（原地），返回 a 本身。平均 O(n log n)，最壞 O(n^2)，棧空間 O(log n)。"""
    _quick_sort3(a, 0, len(a) - 1)
    return a


def _quick_sort_lomuto(a: List[int], lo: int, hi: int) -> None:
    """Lomuto 版快排（原地），作爲對照實現保留。"""
    while lo < hi:
        if hi - lo + 1 <= SMALL:
            insertion_sort(a, lo, hi)
            return
        p = lomuto_partition(a, lo, hi)
        if p - lo < hi - p:
            _quick_sort_lomuto(a, lo, p - 1)
            lo = p + 1
        else:
            _quick_sort_lomuto(a, p + 1, hi)
            hi = p - 1


def quick_sort_lomuto(a: List[int]) -> List[int]:
    """Lomuto 版快排（原地），返回 a 本身。"""
    _quick_sort_lomuto(a, 0, len(a) - 1)
    return a


def quick_select(a: List[int], k: int, rng=None) -> int:
    """返回第 k 小元素（k 爲 0-based），會原地修改 a。

    平均 O(n)（每輪期望砍掉一半），最壞 O(n^2)，空間 O(1)。
    """
    lo, hi = 0, len(a) - 1
    while lo <= hi:
        if hi - lo + 1 <= SMALL:
            insertion_sort(a, lo, hi)
            return a[k]
        lt, gt = partition3(a, lo, hi, rng)
        if k < lt:
            hi = lt - 1                 # k 落在「小於」段
        elif k > gt:
            lo = gt + 1                 # k 落在「大於」段
        else:
            return a[k]                 # k 落在「等於」段，pivot 就是答案
    raise IndexError("k out of range")


def kth_smallest(a: List[int], k: int) -> int:
    """第 k 小元素（k 爲 1-based）；k 越界返回 -1。內部拷貝一份以免破壞原數組。"""
    if k < 1 or k > len(a):
        return -1
    return quick_select(list(a), k - 1)


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
    a = [_next_int(it) for _ in range(n)]
    k = _next_int(it)

    quick_sort(a)
    print(" ".join(str(v) for v in a))          # n = 0 時輸出空行
    print(kth_smallest(a, k) if a else -1)      # 排好序後取第 k 個，這裡直接取也等價


def run_tests() -> None:
    # README 示例
    a = [5, 3, 8, 3, 1, 9, 3]
    assert quick_sort(list(a)) == [1, 3, 3, 3, 5, 8, 9]
    assert quick_sort_lomuto(list(a)) == [1, 3, 3, 3, 5, 8, 9]
    assert kth_smallest(a, 1) == 1
    assert kth_smallest(a, 3) == 3
    assert kth_smallest(a, 7) == 9
    assert kth_smallest(a, 0) == -1             # k 越界
    assert kth_smallest(a, 8) == -1

    # 空數組 / 單元素
    assert quick_sort([]) == []
    assert quick_sort_lomuto([]) == []
    assert quick_sort([42]) == [42]
    assert kth_smallest([], 1) == -1
    assert kth_smallest([42], 1) == 42

    # 已排序 / 逆序 / 全相同 —— 三種最容易觸發退化的輸入
    for arr in ([1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
                [10, 9, 8, 7, 6, 5, 4, 3, 2, 1],
                [7] * 20,
                [-3, -1, -2, -5, -4],
                [0, 0, -1, 1, 0]):
        exp = sorted(arr)
        assert quick_sort(list(arr)) == exp
        assert quick_sort_lomuto(list(arr)) == exp
        for k in range(1, len(arr) + 1):
            assert kth_smallest(arr, k) == exp[k - 1]

    # Lomuto 分區的返回值必須自洽：pivot 左邊的都 < pivot，右邊的都 >= pivot
    b = [4, 2, 7, 2, 9, 1]
    p = lomuto_partition(b, 0, len(b) - 1)
    assert all(v < b[p] for v in b[:p])
    assert all(v >= b[p] for v in b[p + 1:])

    # 三路分區：切出來的三段必須真的有序
    c = [5, 1, 5, 3, 5, 2, 5]
    lt, gt = partition3(c, 0, len(c) - 1, random.Random(1))
    pivot_val = c[lt]
    assert all(v < pivot_val for v in c[:lt])
    assert all(v == pivot_val for v in c[lt:gt + 1])
    assert all(v > pivot_val for v in c[gt + 1:])

    # 小區間閾值附近（<= 16 走插入排序分支）
    for n in range(0, 40):
        arr = [random.randint(-5, 5) for _ in range(n)]
        assert quick_sort(list(arr)) == sorted(arr)
        assert quick_sort_lomuto(list(arr)) == sorted(arr)

    random.seed(20260928)
    rng = random.Random(20260928)

    # 隨機對拍：排序結果與 sorted() 比對，quickselect 與排序結果比對
    for _ in range(500):
        n = random.randint(0, 60)
        # 值域故意取小，製造大量重複元素，考驗三路分區
        arr = [random.randint(-9, 9) for _ in range(n)]
        exp = sorted(arr)
        assert quick_sort(list(arr)) == exp
        assert quick_sort_lomuto(list(arr)) == exp
        for k in range(1, n + 1):
            assert kth_smallest(arr, k) == exp[k - 1]
        assert kth_smallest(arr, 0) == -1
        assert kth_smallest(arr, n + 1) == -1

    # 大數組：確認沒有爆棧、也沒有退化到不可接受
    big = [random.randint(-1000, 1000) for _ in range(20000)]
    assert quick_sort(list(big)) == sorted(big)
    assert kth_smallest(big, 12345) == sorted(big)[12344]

    print("all tests passed")


if __name__ == "__main__":
    raw = sys.stdin.read() if not sys.stdin.isatty() else ""
    if raw.strip():
        run_io(raw)
    else:
        run_tests()
