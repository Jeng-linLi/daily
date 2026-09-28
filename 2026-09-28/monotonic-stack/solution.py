"""單調棧（下一個更大元素 / 每日溫度 / 柱狀圖最大矩形）

題意：
    給定一個長度爲 n 的非負整數數組 a（允許重複、允許爲空），同一個數組上回答三個
    經典的單調棧問題：
      1) 下一個更大元素：對每個位置 i，找右邊第一個 **嚴格大於** a[i] 的元素下標，沒有則 -1；
      2) 每日溫度：與 (1) 同一件事，但輸出的是「距離」（下標差），沒有則 0；
      3) 柱狀圖最大矩形：把 a 視爲柱狀圖的柱高，求其中最大的矩形面積（LeetCode 84）。

思路：
    單調棧解決的是一類「爲每個元素找左/右邊第一個滿足某種大小關係的元素」的問題。
    核心只有一句話：**棧裏始終保持一個單調序列，新元素入棧前先把被它"破壞"單調性的
    元素彈出去，而那些被彈出的元素，答案恰好就是當前這個新元素。**

    1) 下一個更大元素：棧裏存下標，對應的值**單調遞減**（從棧底到棧頂）。
       掃描到 a[i] 時，所有棧裏比 a[i] 小的元素都被 a[i] "擋住"了——
       a[i] 就是它們右邊第一個更大元素，彈出來把答案記爲 i。
       每個下標入棧一次、出棧一次，所以總共 O(n)，而不是 O(n^2)。

    2) 每日溫度：與 (1) 完全同源，差別只在寫進答案的是 `i - j`（距離）而不是下標 i。
       用嚴格大於（`<`）比較，所以相等的元素不會被互相彈掉。

    3) 柱狀圖最大矩形：反過來維護**單調遞增**棧。掃描時在末尾補一個高度 0 的哨兵，
       強制把棧清空。當遇到一個比棧頂矮的柱子 h[i] 時，棧頂那根柱子能向右延伸的
       邊界就被確定了 —— 就是 i - 1；向左的邊界是彈出後的新棧頂 + 1。
       於是以這根柱子爲高的最大矩形寬 = `i - 新棧頂 - 1`。
       注意這裡用的是 `>`（嚴格大於才彈），相等高度的柱子留在棧裏，
       這樣最左邊那根相等柱子會負責算出跨越整段的最大矩形，不會漏解。

    三個問題共用同一個套路，區別只有兩處：棧是遞增還是遞減、彈棧條件是 `<` 還是 `>`。

輸入格式（stdin，所有數字按空白分隔即可）：
    n
    a1 a2 ... an        （n = 0 時這一行直接省略）
輸出格式（stdout）：
    第 1 行：下一個更大元素的下標，空格分隔（-1 表示沒有；n = 0 時輸出空行）
    第 2 行：每日溫度（距離下一個更大元素的下標差，0 表示沒有）
    第 3 行：柱狀圖最大矩形的面積（n = 0 時爲 0）
無 stdin 輸入時運行內置斷言測試並輸出 `all tests passed`。
"""

import sys
from typing import Iterator, List

import random


def next_greater_index(a: List[int]) -> List[int]:
    """每個位置右邊第一個嚴格大於它的元素下標；沒有則 -1。

    棧中下標對應的值單調遞減。時間 O(n)，空間 O(n)。
    """
    n = len(a)
    res = [-1] * n
    stack: List[int] = []          # 存下標，值單調遞減
    for i, v in enumerate(a):
        while stack and a[stack[-1]] < v:
            res[stack.pop()] = i   # a[i] 就是這些元素右邊第一個更大的
        stack.append(i)
    return res


def daily_temperatures(a: List[int]) -> List[int]:
    """與 next_greater_index 同源，但輸出「距離」而非下標；沒有則 0。

    時間 O(n)，空間 O(n)。
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
    """柱狀圖最大矩形面積（LeetCode 84）。時間 O(n)，空間 O(n)。

    維護單調遞增棧，末尾補一個高度 0 的哨兵把棧清空。
    """
    n = len(h)
    stack: List[int] = []          # 存下標，值單調遞增
    best = 0
    for i in range(n + 1):
        cur = h[i] if i < n else 0          # 哨兵：高度 0 會彈出所有柱子
        while stack and h[stack[-1]] > cur:  # 嚴格大於才彈：相等高度留在棧裏，避免漏解
            top = stack.pop()
            left = stack[-1] if stack else -1
            width = i - left - 1             # 能向右延伸到 i-1，向左到 left+1
            area = h[top] * width
            if area > best:
                best = area
        stack.append(i)
    return best


# ---------------- 對照用的 O(n^2) 暴力實現 ----------------

def next_greater_index_brute(a: List[int]) -> List[int]:
    """向右線性掃描第一個更大的元素，用於對拍。"""
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
    """暴力版每日溫度，用於對拍。"""
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
    """枚舉每根柱子並向兩側擴張，用於對拍。時間 O(n^2)。"""
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

    print(" ".join(str(v) for v in next_greater_index(a)))   # n = 0 時輸出空行
    print(" ".join(str(v) for v in daily_temperatures(a)))
    print(largest_rectangle(a))


def run_tests() -> None:
    # README 示例：a = [2, 1, 2, 4, 3]
    a = [2, 1, 2, 4, 3]
    assert next_greater_index(a) == [3, 2, 3, -1, -1]
    assert daily_temperatures(a) == [3, 1, 1, 0, 0]
    assert largest_rectangle(a) == 6      # 高 2 寬 3：[2,2,4,3] 中的 2x3

    # 經典用例
    assert next_greater_index([73, 74, 75, 71, 69, 72, 76, 73]) == [1, 2, 6, 5, 5, 6, -1, -1]
    assert daily_temperatures([73, 74, 75, 71, 69, 72, 76, 73]) == [1, 1, 4, 2, 1, 1, 0, 0]
    assert largest_rectangle([2, 1, 5, 6, 2, 3]) == 10      # LeetCode 84 官方用例
    assert largest_rectangle([2, 4]) == 4
    assert largest_rectangle([1, 1, 1, 1]) == 4
    assert largest_rectangle([5]) == 5
    assert largest_rectangle([0]) == 0
    assert largest_rectangle([0, 0, 0]) == 0

    # 遞減 / 遞增 / 全相同：三種極端形態
    assert next_greater_index([5, 4, 3, 2, 1]) == [-1, -1, -1, -1, -1]
    assert daily_temperatures([5, 4, 3, 2, 1]) == [0, 0, 0, 0, 0]
    assert largest_rectangle([5, 4, 3, 2, 1]) == 9          # 高 3 寬 3：(3,2,1)
    assert next_greater_index([1, 2, 3, 4, 5]) == [1, 2, 3, 4, -1]
    assert daily_temperatures([1, 2, 3, 4, 5]) == [1, 1, 1, 1, 0]
    assert largest_rectangle([1, 2, 3, 4, 5]) == 9          # 高 3 寬 3：(3,4,5)
    assert next_greater_index([3, 3, 3]) == [-1, -1, -1]    # 嚴格大於，相等不算
    assert daily_temperatures([3, 3, 3]) == [0, 0, 0]
    assert largest_rectangle([3, 3, 3]) == 9

    # 空數組
    assert next_greater_index([]) == []
    assert daily_temperatures([]) == []
    assert largest_rectangle([]) == 0

    # 答案自洽：next_greater 與 daily_temperatures 必須指向同一個位置
    for arr in ([2, 1, 2, 4, 3], [1], [4, 2, 9, 1, 7], [0, 0, 5, 0]):
        ng = next_greater_index(arr)
        dt = daily_temperatures(arr)
        for i in range(len(arr)):
            if ng[i] == -1:
                assert dt[i] == 0
            else:
                assert dt[i] == ng[i] - i
                assert arr[ng[i]] > arr[i]
                assert all(arr[k] <= arr[i] for k in range(i + 1, ng[i]))  # 中間沒有更大的

    random.seed(20260928)

    # 隨機對拍：三個函數全部與 O(n^2) 暴力解比對
    for _ in range(500):
        n = random.randint(0, 40)
        arr = [random.randint(0, 12) for _ in range(n)]     # 值域小 → 大量重複
        assert next_greater_index(arr) == next_greater_index_brute(arr)
        assert daily_temperatures(arr) == daily_temperatures_brute(arr)
        assert largest_rectangle(arr) == largest_rectangle_brute(arr)
        # 最大矩形面積不會超過「最大高度 x n」，也不會小於最大高度
        if n > 0:
            assert max(arr) <= largest_rectangle(arr) <= max(arr) * n

    print("all tests passed")


if __name__ == "__main__":
    raw = sys.stdin.read() if not sys.stdin.isatty() else ""
    if raw.strip():
        run_io(raw)
    else:
        run_tests()
