"""Manacher 算法：最長回文子串 / 回文子串計數（O(n)）

題意：
    給定一個字串 s，要求：
      1. 最長回文子串的長度；
      2. 最長回文子串本身（若有多個同長度的，取最靠左的那個）；
      3. 回文子串的總數（不同位置的同一內容算不同子串）；
      4. 以每個位置為中心的奇回文半徑陣列 d1；
      5. 以每個間隙為中心的偶回文半徑陣列 d2。
    全部要在 O(n) 內完成（n = len(s)）。

思路：
    **中心擴展**是最直觀的做法：枚舉 n 個「奇中心」和 n−1 個「偶中心」，從中心往兩邊擴，
    最壞 O(n^2)（例如 "aaaa…a"）。Manacher 的洞見是：**已經算過的回文可以被複用**。

    維護當前「右端點最遠」的回文區間 `[l, r]`。處理中心 i 時：
      - 若 i 在 (l, r) 內，設 j 為 i 關於回文中心 `(l+r)/2` 的**鏡像位置**，
        則 i 處的半徑至少是 `min(已算出的 d[j], r − i + 1)`（超出右界的部分不能保證，必須老實再擴）；
      - 若 i ≥ r，從最小半徑開始擴。
    由於 `[l, r]` 的右端只會單調往右移動，整個過程的 while 迴圈總執行次數是 O(n)，
    因此整體 O(n)、O(n) 空間。

    兩個半徑陣列的定義（cp-algorithms 的寫法）：
      - `d1[i]`：以 i 為中心的**奇**回文半徑（含中心），最長奇回文長度 = `2*d1[i] − 1`；
        同時 `d1[i]` 也正好等於「以 i 為中心的奇回文個數」。
      - `d2[i]`：以 i−1 與 i 之間的間隙為中心的**偶**回文半徑，最長偶回文長度 = `2*d2[i]`；
        同樣地，`d2[i]` 等於「以該間隙為中心的偶回文個數」。

    因此**回文子串總數 = sum(d1) + sum(d2)** —— 這是 Manacher 一個很漂亮的副產品。

    最長回文子串：掃一遍 d1 / d2 取最大值即可；長度相同時要取最靠左的，
    所以更新條件必須是嚴格 `>`（先掃到的先佔住）。

輸入格式（stdin）：
    s          （第一行；可含空格以外的任意字元。空字串請給一個只有換行的輸入）
輸出格式（stdout）：
    第 1 行：最長回文子串長度
    第 2 行：最長回文子串（長度為 0 時輸出空行）
    第 3 行：回文子串總數
    第 4 行：d1 陣列，空格分隔（n = 0 時輸出空行）
    第 5 行：d2 陣列，空格分隔（n = 0 時輸出空行）
無 stdin 輸入時運行內置斷言測試並輸出 `all tests passed`。
"""

import random
import sys
from typing import List, Tuple


# ---------------------------------------------------------------- Manacher
def manacher_odd(s: str) -> List[int]:
    """d1[i] = 以 i 為中心的奇回文半徑（含中心）；最長奇回文長度 = 2*d1[i] − 1。"""
    n = len(s)
    d1 = [0] * n
    l, r = 0, -1                  # 目前右端最遠的回文區間 [l, r]
    for i in range(n):
        # i 落在已知回文內時，先繼承鏡像位置的半徑，但不能越過右界
        k = 1 if i > r else min(d1[l + r - i], r - i + 1)
        while i - k >= 0 and i + k < n and s[i - k] == s[i + k]:
            k += 1
        d1[i] = k
        if i + k - 1 > r:         # 右端推進了，更新當前最右回文
            l = i - k + 1
            r = i + k - 1
    return d1


def manacher_even(s: str) -> List[int]:
    """d2[i] = 以 i−1 與 i 之間為中心的偶回文半徑；最長偶回文長度 = 2*d2[i]。"""
    n = len(s)
    d2 = [0] * n
    l, r = 0, -1
    for i in range(n):
        # 鏡像關係與奇回文差 1，注意索引是 l + r − i + 1
        k = 0 if i > r else min(d2[l + r - i + 1], r - i + 1)
        while i - k - 1 >= 0 and i + k < n and s[i - k - 1] == s[i + k]:
            k += 1
        d2[i] = k
        if i + k - 1 > r:
            l = i - k
            r = i + k - 1
    return d2


# ---------------------------------------------------------------- 中心擴展（暴力基準）
def center_expansion_odd(s: str) -> List[int]:
    """O(n^2) 求奇回文半徑，只在測試裡當 Manacher 的對拍基準。"""
    n = len(s)
    d1 = [0] * n
    for i in range(n):
        k = 1
        while i - k >= 0 and i + k < n and s[i - k] == s[i + k]:
            k += 1
        d1[i] = k
    return d1


def center_expansion_even(s: str) -> List[int]:
    """O(n^2) 求偶回文半徑，只在測試裡當對拍基準。"""
    n = len(s)
    d2 = [0] * n
    for i in range(n):
        k = 0
        while i - k - 1 >= 0 and i + k < n and s[i - k - 1] == s[i + k]:
            k += 1
        d2[i] = k
    return d2


# ---------------------------------------------------------------- 由半徑推出答案
def longest_palindrome(s: str, d1: List[int], d2: List[int]) -> Tuple[int, str]:
    """最長回文子串的 (長度, 內容)；同長度取最靠左的（嚴格 > 更新）。"""
    n = len(s)
    best = 0
    start = 0
    for i in range(n):
        if 2 * d1[i] - 1 > best:
            best = 2 * d1[i] - 1
            start = i - d1[i] + 1
        if 2 * d2[i] > best:
            best = 2 * d2[i]
            start = i - d2[i]
    return best, s[start:start + best]


def count_palindromes(d1: List[int], d2: List[int]) -> int:
    """回文子串總數 = sum(d1) + sum(d2)。"""
    return sum(d1) + sum(d2)


def longest_palindrome_naive(s: str) -> Tuple[int, str]:
    """O(n^3) 枚舉所有子串求最長回文，只在測試裡當基準。"""
    n = len(s)
    best = 0
    ans = ""
    for i in range(n):
        for j in range(i, n):
            sub = s[i:j + 1]
            if len(sub) > best and sub == sub[::-1]:
                best = len(sub)
                ans = sub
    return best, ans


def count_palindromes_naive(s: str) -> int:
    """O(n^3) 統計回文子串個數，只在測試裡當基準。"""
    n = len(s)
    cnt = 0
    for i in range(n):
        for j in range(i, n):
            if s[i:j + 1] == s[i:j + 1][::-1]:
                cnt += 1
    return cnt


# ---------------------------------------------------------------- IO 模式
def run_io(data: str) -> None:
    s = data.split("\n")[0].rstrip("\r")
    d1 = manacher_odd(s)
    d2 = manacher_even(s)
    length, sub = longest_palindrome(s, d1, d2)
    print(length)
    print(sub)
    print(count_palindromes(d1, d2))
    print(" ".join(map(str, d1)))
    print(" ".join(map(str, d2)))


# ---------------------------------------------------------------- 測試
def run_tests() -> None:
    # 空字串
    assert manacher_odd("") == []
    assert manacher_even("") == []
    assert longest_palindrome("", [], []) == (0, "")
    assert count_palindromes([], []) == 0

    # 單字元
    assert manacher_odd("a") == [1]
    assert manacher_even("a") == [0]
    assert longest_palindrome("a", [1], [0]) == (1, "a")
    assert count_palindromes([1], [0]) == 1

    # README 示例："cbbd" —— 最長回文是偶回文 "bb"
    s = "cbbd"
    d1 = manacher_odd(s)
    d2 = manacher_even(s)
    assert d1 == [1, 1, 1, 1]
    assert d2 == [0, 0, 1, 0]
    assert d1 == center_expansion_odd(s)
    assert d2 == center_expansion_even(s)
    assert longest_palindrome(s, d1, d2) == (2, "bb")
    assert count_palindromes(d1, d2) == count_palindromes_naive(s) == 5

    # 同長度取最左："abacdc" 的 "aba" 與 "cdc" 都是 3，取 "aba"
    s = "abacdc"
    d1 = manacher_odd(s)
    d2 = manacher_even(s)
    assert d1 == [1, 2, 1, 1, 2, 1]
    assert d2 == [0, 0, 0, 0, 0, 0]
    assert d1 == center_expansion_odd(s)
    assert d2 == center_expansion_even(s)
    assert longest_palindrome(s, d1, d2) == (3, "aba")     # 同長度取最左：aba 先於 cdc
    assert count_palindromes(d1, d2) == 8
    assert count_palindromes(d1, d2) == count_palindromes_naive(s)

    # 全同字元："aaaa" 的回文子串數 = 4*5/2 = 10
    s = "aaaa"
    d1 = manacher_odd(s)
    d2 = manacher_even(s)
    assert d1 == [1, 2, 2, 1]
    assert d2 == [0, 1, 2, 1]                              # d2[3] 只到 1（右邊沒有第 5 個 a）
    assert d1 == center_expansion_odd(s)
    assert d2 == center_expansion_even(s)
    assert longest_palindrome(s, d1, d2) == (4, "aaaa")
    assert count_palindromes(d1, d2) == 10
    assert count_palindromes(d1, d2) == count_palindromes_naive(s)

    # 完全無回文（長度 > 1）："abcde"
    s = "abcde"
    d1 = manacher_odd(s)
    d2 = manacher_even(s)
    assert longest_palindrome(s, d1, d2) == (1, "a")
    assert count_palindromes(d1, d2) == 5

    # 偶回文最長："abba"
    s = "abba"
    d1 = manacher_odd(s)
    d2 = manacher_even(s)
    assert longest_palindrome(s, d1, d2) == (4, "abba")
    assert count_palindromes(d1, d2) == count_palindromes_naive(s) == 6

    # 隨機對拍：Manacher vs 中心擴展 vs 暴力枚舉
    random.seed(20261002)
    alphabet = "ab"
    for _ in range(400):
        n = random.randint(0, 12)
        s = "".join(random.choice(alphabet) for _ in range(n))
        d1 = manacher_odd(s)
        d2 = manacher_even(s)

        assert d1 == center_expansion_odd(s)
        assert d2 == center_expansion_even(s)

        # 每個半徑都要真的對應一個回文，且再往外一格就不回文
        for i in range(n):
            assert d1[i] >= 1
            lo, hi = i - d1[i] + 1, i + d1[i] - 1
            assert s[lo:hi + 1] == s[lo:hi + 1][::-1]
            if lo - 1 >= 0 and hi + 1 < n:
                assert s[lo - 1] != s[hi + 1]
            lo, hi = i - d2[i], i + d2[i] - 1
            if d2[i] > 0:
                assert s[lo:hi + 1] == s[lo:hi + 1][::-1]
            if i - d2[i] - 1 >= 0 and i + d2[i] < n:
                assert s[i - d2[i] - 1] != s[i + d2[i]]

        assert longest_palindrome(s, d1, d2) == longest_palindrome_naive(s)
        assert count_palindromes(d1, d2) == count_palindromes_naive(s)

    # 三字元字母表再跑一輪，增加「同長度取最左」的覆蓋
    for _ in range(200):
        n = random.randint(0, 14)
        s = "".join(random.choice("abc") for _ in range(n))
        d1 = manacher_odd(s)
        d2 = manacher_even(s)
        assert d1 == center_expansion_odd(s)
        assert d2 == center_expansion_even(s)
        assert longest_palindrome(s, d1, d2) == longest_palindrome_naive(s)
        assert count_palindromes(d1, d2) == count_palindromes_naive(s)

    print("all tests passed")


if __name__ == "__main__":
    raw = sys.stdin.read() if not sys.stdin.isatty() else ""
    if raw:
        run_io(raw)
    else:
        run_tests()
