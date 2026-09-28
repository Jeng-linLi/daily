"""快速排序与快速选择（随机化 pivot + 三路分区 + Lomuto 分区）

题意：
    给定 n 个整数（允许重复、允许负数、允许为空），要求：
      1) 输出升序排序结果；
      2) 输出第 k 小元素（k 为 1-based，越界时输出 -1）。

思路：
    快排是「分治 + 原地分区」：每轮挑一个 pivot，把区间切成
    「< pivot」「== pivot」「> pivot」三段，然后只对左右两段递归。

    1) **Lomuto 分区**：维护指针 i，把 < pivot 的元素换到左边。
       写法最直观，但遇到大量与 pivot 相等的元素时，两段会严重不平衡，
       所有元素都相等时会退化成 O(n^2)。

    2) **三路分区（Dijkstra / 荷兰国旗）**：维护 lt / i / gt 三个指针，
       a[lo:lt] < pivot，a[lt:i] == pivot，a[i:gt+1] 待定，a[gt+1:hi+1] > pivot。
       相等元素一次整段归位，**全相同元素时是 O(n)**，这是它最大的价值。

    3) **随机化 pivot**：固定挑首/尾元素时，已排序输入会退化成 O(n^2)。
       随机挑 pivot 让任何输入分布的期望时间都是 O(n log n)。
       实践中还会「小区间（长度 < 16）改用插入排序」——常数更小。

    4) **递归深度控制**：先递归短的那半、用 while 循环处理长的那半，
       栈深度被压到 O(log n)，不会因最坏情况爆栈。

    5) **快速选择（Quickselect）**：排序只为了拿第 k 个太浪费。
       分区后看 k 落在哪一段：落在 == pivot 段就直接返回，否则只递归
       包含 k 的那一段。平均 O(n)，最坏 O(n^2)（随机化后几乎不会碰到）。

    注意：本文件里的排序是**原地**的，会修改传入的列表。

输入格式（stdin，所有数字按空白分隔即可）：
    n
    a1 a2 ... an        （n = 0 时这一行直接省略）
    k
输出格式（stdout）：
    第 1 行：升序排序结果，空格分隔（n = 0 时输出空行）
    第 2 行：第 k 小元素的值（1-based；k 越界或 n = 0 时输出 -1）
无 stdin 输入时运行内置断言测试并输出 `all tests passed`。
"""

import random
import sys
from typing import Iterator, List

SMALL = 16  # 小区间阈值：低于这个长度改用插入排序，常数更小


def insertion_sort(a: List[int], lo: int = 0, hi: int = -1) -> None:
    """对 a[lo:hi+1] 做插入排序。时间 O(len^2)，但小区间里常数极小。"""
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
    """Lomuto 分区：以 a[hi] 为 pivot，把 < pivot 的换到左边，返回 pivot 最终下标。

    时间 O(hi-lo+1)，空间 O(1)。重复元素多时两段会不平衡。
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
    """三路分区：把 a[lo:hi+1] 切成 < == > 三段，返回 (lt, gt)。

    a[lo:lt] < pivot，a[lt:gt+1] == pivot，a[gt+1:hi+1] > pivot。
    时间 O(hi-lo+1)，空间 O(1)。pivot 随机挑选。
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
            gt -= 1          # 换过来的元素还没看过，i 不动
        else:
            i += 1
    return lt, gt


def _quick_sort3(a: List[int], lo: int, hi: int, rng=None) -> None:
    """三路快排主过程：先递归短半边，长半边用循环，栈深度 O(log n)。"""
    while lo < hi:
        if hi - lo + 1 <= SMALL:
            insertion_sort(a, lo, hi)
            return
        lt, gt = partition3(a, lo, hi, rng)
        if lt - lo < hi - gt:          # 左段更短：递归左段，循环处理右段
            _quick_sort3(a, lo, lt - 1, rng)
            lo = gt + 1
        else:                          # 右段更短：递归右段，循环处理左段
            _quick_sort3(a, gt + 1, hi, rng)
            hi = lt - 1


def quick_sort(a: List[int]) -> List[int]:
    """三路快排（原地），返回 a 本身。平均 O(n log n)，最坏 O(n^2)，栈空间 O(log n)。"""
    _quick_sort3(a, 0, len(a) - 1)
    return a


def _quick_sort_lomuto(a: List[int], lo: int, hi: int) -> None:
    """Lomuto 版快排（原地），作为对照实现保留。"""
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
    """返回第 k 小元素（k 为 0-based），会原地修改 a。

    平均 O(n)（每轮期望砍掉一半），最坏 O(n^2)，空间 O(1)。
    """
    lo, hi = 0, len(a) - 1
    while lo <= hi:
        if hi - lo + 1 <= SMALL:
            insertion_sort(a, lo, hi)
            return a[k]
        lt, gt = partition3(a, lo, hi, rng)
        if k < lt:
            hi = lt - 1                 # k 落在「小于」段
        elif k > gt:
            lo = gt + 1                 # k 落在「大于」段
        else:
            return a[k]                 # k 落在「等于」段，pivot 就是答案
    raise IndexError("k out of range")


def kth_smallest(a: List[int], k: int) -> int:
    """第 k 小元素（k 为 1-based）；k 越界返回 -1。内部拷贝一份以免破坏原数组。"""
    if k < 1 or k > len(a):
        return -1
    return quick_select(list(a), k - 1)


# ---------------- IO ----------------

def _next_int(it: Iterator[str], default: int = 0) -> int:
    """取下一个整数；输入被截断时用默认值兜底，避免直接抛异常。"""
    try:
        return int(next(it))
    except (StopIteration, ValueError):
        return default


def run_io(data: str) -> None:
    """按统一输入输出格式处理 stdin 数据。"""
    it = iter(data.split())
    n = _next_int(it)
    a = [_next_int(it) for _ in range(n)]
    k = _next_int(it)

    quick_sort(a)
    print(" ".join(str(v) for v in a))          # n = 0 时输出空行
    print(kth_smallest(a, k) if a else -1)      # 排好序后取第 k 个，这里直接取也等价


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

    # 空数组 / 单元素
    assert quick_sort([]) == []
    assert quick_sort_lomuto([]) == []
    assert quick_sort([42]) == [42]
    assert kth_smallest([], 1) == -1
    assert kth_smallest([42], 1) == 42

    # 已排序 / 逆序 / 全相同 —— 三种最容易触发退化的输入
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

    # Lomuto 分区的返回值必须自洽：pivot 左边的都 < pivot，右边的都 >= pivot
    b = [4, 2, 7, 2, 9, 1]
    p = lomuto_partition(b, 0, len(b) - 1)
    assert all(v < b[p] for v in b[:p])
    assert all(v >= b[p] for v in b[p + 1:])

    # 三路分区：切出来的三段必须真的有序
    c = [5, 1, 5, 3, 5, 2, 5]
    lt, gt = partition3(c, 0, len(c) - 1, random.Random(1))
    pivot_val = c[lt]
    assert all(v < pivot_val for v in c[:lt])
    assert all(v == pivot_val for v in c[lt:gt + 1])
    assert all(v > pivot_val for v in c[gt + 1:])

    # 小区间阈值附近（<= 16 走插入排序分支）
    for n in range(0, 40):
        arr = [random.randint(-5, 5) for _ in range(n)]
        assert quick_sort(list(arr)) == sorted(arr)
        assert quick_sort_lomuto(list(arr)) == sorted(arr)

    random.seed(20260928)
    rng = random.Random(20260928)

    # 随机对拍：排序结果与 sorted() 比对，quickselect 与排序结果比对
    for _ in range(500):
        n = random.randint(0, 60)
        # 值域故意取小，制造大量重复元素，考验三路分区
        arr = [random.randint(-9, 9) for _ in range(n)]
        exp = sorted(arr)
        assert quick_sort(list(arr)) == exp
        assert quick_sort_lomuto(list(arr)) == exp
        for k in range(1, n + 1):
            assert kth_smallest(arr, k) == exp[k - 1]
        assert kth_smallest(arr, 0) == -1
        assert kth_smallest(arr, n + 1) == -1

    # 大数组：确认没有爆栈、也没有退化到不可接受
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
