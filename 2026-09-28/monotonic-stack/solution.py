"""单调栈（下一个更大元素 / 每日温度 / 柱状图最大矩形）

题意：
    给定一个长度为 n 的非负整数数组 a（允许重复、允许为空），同一个数组上回答三个
    经典的单调栈问题：
      1) 下一个更大元素：对每个位置 i，找右边第一个 **严格大于** a[i] 的元素下标，没有则 -1；
      2) 每日温度：与 (1) 同一件事，但输出的是「距离」（下标差），没有则 0；
      3) 柱状图最大矩形：把 a 视为柱状图的柱高，求其中最大的矩形面积（LeetCode 84）。

思路：
    单调栈解决的是一类「为每个元素找左/右边第一个满足某种大小关系的元素」的问题。
    核心只有一句话：**栈里始终保持一个单调序列，新元素入栈前先把被它"破坏"单调性的
    元素弹出去，而那些被弹出的元素，答案恰好就是当前这个新元素。**

    1) 下一个更大元素：栈里存下标，对应的值**单调递减**（从栈底到栈顶）。
       扫描到 a[i] 时，所有栈里比 a[i] 小的元素都被 a[i] "挡住"了——
       a[i] 就是它们右边第一个更大元素，弹出来把答案记为 i。
       每个下标入栈一次、出栈一次，所以总共 O(n)，而不是 O(n^2)。

    2) 每日温度：与 (1) 完全同源，差别只在写进答案的是 `i - j`（距离）而不是下标 i。
       用严格大于（`<`）比较，所以相等的元素不会被互相弹掉。

    3) 柱状图最大矩形：反过来维护**单调递增**栈。扫描时在末尾补一个高度 0 的哨兵，
       强制把栈清空。当遇到一个比栈顶矮的柱子 h[i] 时，栈顶那根柱子能向右延伸的
       边界就被确定了 —— 就是 i - 1；向左的边界是弹出后的新栈顶 + 1。
       于是以这根柱子为高的最大矩形宽 = `i - 新栈顶 - 1`。
       注意这里用的是 `>`（严格大于才弹），相等高度的柱子留在栈里，
       这样最左边那根相等柱子会负责算出跨越整段的最大矩形，不会漏解。

    三个问题共用同一个套路，区别只有两处：栈是递增还是递减、弹栈条件是 `<` 还是 `>`。

输入格式（stdin，所有数字按空白分隔即可）：
    n
    a1 a2 ... an        （n = 0 时这一行直接省略）
输出格式（stdout）：
    第 1 行：下一个更大元素的下标，空格分隔（-1 表示没有；n = 0 时输出空行）
    第 2 行：每日温度（距离下一个更大元素的下标差，0 表示没有）
    第 3 行：柱状图最大矩形的面积（n = 0 时为 0）
无 stdin 输入时运行内置断言测试并输出 `all tests passed`。
"""

import sys
from typing import Iterator, List

import random


def next_greater_index(a: List[int]) -> List[int]:
    """每个位置右边第一个严格大于它的元素下标；没有则 -1。

    栈中下标对应的值单调递减。时间 O(n)，空间 O(n)。
    """
    n = len(a)
    res = [-1] * n
    stack: List[int] = []          # 存下标，值单调递减
    for i, v in enumerate(a):
        while stack and a[stack[-1]] < v:
            res[stack.pop()] = i   # a[i] 就是这些元素右边第一个更大的
        stack.append(i)
    return res


def daily_temperatures(a: List[int]) -> List[int]:
    """与 next_greater_index 同源，但输出「距离」而非下标；没有则 0。

    时间 O(n)，空间 O(n)。
    """
    n = len(a)
    res = [0] * n
    stack: List[int] = []
    for i, v in enumerate(a):
        while stack and a[stack[-1]] < v:
            j = stack.pop()
            res[j] = i - j         # 隔了多少天才等到更暖的一天
        stack.append(i)
    return res


def largest_rectangle(h: List[int]) -> int:
    """柱状图最大矩形面积（LeetCode 84）。时间 O(n)，空间 O(n)。

    维护单调递增栈，末尾补一个高度 0 的哨兵把栈清空。
    """
    n = len(h)
    stack: List[int] = []          # 存下标，值单调递增
    best = 0
    for i in range(n + 1):
        cur = h[i] if i < n else 0          # 哨兵：高度 0 会弹出所有柱子
        while stack and h[stack[-1]] > cur:  # 严格大于才弹：相等高度留在栈里，避免漏解
            top = stack.pop()
            left = stack[-1] if stack else -1
            width = i - left - 1             # 能向右延伸到 i-1，向左到 left+1
            area = h[top] * width
            if area > best:
                best = area
        stack.append(i)
    return best


# ---------------- 对照用的 O(n^2) 暴力实现 ----------------

def next_greater_index_brute(a: List[int]) -> List[int]:
    """向右线性扫描第一个更大的元素，用于对拍。"""
    n = len(a)
    res = []
    for i in range(n):
        j = -1
        for k in range(i + 1, n):
            if a[k] > a[i]:
                j = k
                break
        res.append(j)
    return res


def daily_temperatures_brute(a: List[int]) -> List[int]:
    """暴力版每日温度，用于对拍。"""
    n = len(a)
    res = []
    for i in range(n):
        d = 0
        for k in range(i + 1, n):
            if a[k] > a[i]:
                d = k - i
                break
        res.append(d)
    return res


def largest_rectangle_brute(h: List[int]) -> int:
    """枚举每根柱子并向两侧扩张，用于对拍。时间 O(n^2)。"""
    n = len(h)
    best = 0
    for i in range(n):
        left = i
        while left - 1 >= 0 and h[left - 1] >= h[i]:
            left -= 1
        right = i
        while right + 1 < n and h[right + 1] >= h[i]:
            right += 1
        area = h[i] * (right - left + 1)
        if area > best:
            best = area
    return best


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

    print(" ".join(str(v) for v in next_greater_index(a)))   # n = 0 时输出空行
    print(" ".join(str(v) for v in daily_temperatures(a)))
    print(largest_rectangle(a))


def run_tests() -> None:
    # README 示例：a = [2, 1, 2, 4, 3]
    a = [2, 1, 2, 4, 3]
    assert next_greater_index(a) == [3, 2, 3, -1, -1]
    assert daily_temperatures(a) == [3, 1, 1, 0, 0]
    assert largest_rectangle(a) == 6      # 高 2 宽 3：[2,2,4,3] 中的 2x3

    # 经典用例
    assert next_greater_index([73, 74, 75, 71, 69, 72, 76, 73]) == [1, 2, 6, 5, 5, 6, -1, -1]
    assert daily_temperatures([73, 74, 75, 71, 69, 72, 76, 73]) == [1, 1, 4, 2, 1, 1, 0, 0]
    assert largest_rectangle([2, 1, 5, 6, 2, 3]) == 10      # LeetCode 84 官方用例
    assert largest_rectangle([2, 4]) == 4
    assert largest_rectangle([1, 1, 1, 1]) == 4
    assert largest_rectangle([5]) == 5
    assert largest_rectangle([0]) == 0
    assert largest_rectangle([0, 0, 0]) == 0

    # 递减 / 递增 / 全相同：三种极端形态
    assert next_greater_index([5, 4, 3, 2, 1]) == [-1, -1, -1, -1, -1]
    assert daily_temperatures([5, 4, 3, 2, 1]) == [0, 0, 0, 0, 0]
    assert largest_rectangle([5, 4, 3, 2, 1]) == 9          # 高 3 宽 3：(3,2,1)
    assert next_greater_index([1, 2, 3, 4, 5]) == [1, 2, 3, 4, -1]
    assert daily_temperatures([1, 2, 3, 4, 5]) == [1, 1, 1, 1, 0]
    assert largest_rectangle([1, 2, 3, 4, 5]) == 9          # 高 3 宽 3：(3,4,5)
    assert next_greater_index([3, 3, 3]) == [-1, -1, -1]    # 严格大于，相等不算
    assert daily_temperatures([3, 3, 3]) == [0, 0, 0]
    assert largest_rectangle([3, 3, 3]) == 9

    # 空数组
    assert next_greater_index([]) == []
    assert daily_temperatures([]) == []
    assert largest_rectangle([]) == 0

    # 答案自洽：next_greater 与 daily_temperatures 必须指向同一个位置
    for arr in ([2, 1, 2, 4, 3], [1], [4, 2, 9, 1, 7], [0, 0, 5, 0]):
        ng = next_greater_index(arr)
        dt = daily_temperatures(arr)
        for i in range(len(arr)):
            if ng[i] == -1:
                assert dt[i] == 0
            else:
                assert dt[i] == ng[i] - i
                assert arr[ng[i]] > arr[i]
                assert all(arr[k] <= arr[i] for k in range(i + 1, ng[i]))  # 中间没有更大的

    random.seed(20260928)

    # 随机对拍：三个函数全部与 O(n^2) 暴力解比对
    for _ in range(500):
        n = random.randint(0, 40)
        arr = [random.randint(0, 12) for _ in range(n)]     # 值域小 → 大量重复
        assert next_greater_index(arr) == next_greater_index_brute(arr)
        assert daily_temperatures(arr) == daily_temperatures_brute(arr)
        assert largest_rectangle(arr) == largest_rectangle_brute(arr)
        # 最大矩形面积不会超过「最大高度 x n」，也不会小于最大高度
        if n > 0:
            assert max(arr) <= largest_rectangle(arr) <= max(arr) * n

    print("all tests passed")


if __name__ == "__main__":
    raw = sys.stdin.read() if not sys.stdin.isatty() else ""
    if raw.strip():
        run_io(raw)
    else:
        run_tests()
