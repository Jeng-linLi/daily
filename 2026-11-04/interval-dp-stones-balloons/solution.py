"""區間 DP：合併石子與戳氣球

題意：
    兩個經典區間 DP 模型，輸入第一行為任務編號：
    任務 1（合併石子 / Stone Merge）：
        n 堆石子排成一排，每次只能合併**相鄰**兩堆，
        代價 = 兩堆石子數之和，求把所有堆合成一堆的最小總代價。
    任務 2（戳氣球 / Burst Balloons）：
        n 個氣球排成一排，戳破第 i 個獲得
        nums[left-1] * nums[i] * nums[right+1] 枚金幣
        （left / right 為與 i 相鄰且未被戳破的氣球；越界視為 1），
        求能獲得的最大金幣數。

    輸入：第一行 T（1 或 2）；第二行 n；第三行 n 個整數。
    輸出：一個整數。

思路：
    ### 1. 合併石子：從「最後一步」倒推
    設 dp[i][j] = 把第 i..j 堆合成一堆的最小代價。
    最後一步一定是把「i..k 合成的堆」與「k+1..j 合成的堆」合併，
    這一步的代價固定為 sum(i..j)（最終一堆的總石子數）：
        dp[i][j] = min_k { dp[i][k] + dp[k+1][j] } + sum(i..j)
    區間長度從小到大枚舉。前綴和 O(1) 求區間和。

    ### 2. 戳氣球：把「先戳」改成「最後戳」
    直接模擬戳的順序，戳掉一個會改變左右鄰居，狀態難維護。
    反向思考：定義 dp[l][r] = 「開區間 (l, r) 內的氣球全部戳破」
    的最大收益，並規定**最後一個**戳破的是 k（l < k < r）：
    戳破 k 時，區間 (l, r) 內其他氣球已經沒了，k 的左右鄰居
    恰好是 l 和 r（注意：邊界用補的虛擬氣球 1，不需真的戳破）：
        dp[l][r] = max_k { dp[l][k] + dp[k][r] + a[l]*a[k]*a[r] }
    這樣每個子問題的邊界固定，不依賴「外面」被戳掉的歷史。

    ### 3. 對拍
    合併石子用「模擬所有合併順序」的指數級暴力（n ≤ 10）；
    戳氣球用「枚舉戳破排列」的指數級暴力（n ≤ 8）。逐一比對。

複雜度：
    兩個任務都是狀態 O(n²)、轉移 O(n)，總 O(n³)
    空間 O(n²)

應用場景（Applications）：
    - **合併類成本**：文件 / 區塊合併的最小 IO 代價、
      資料庫歸併排序成本、供應鏈中相鄰倉庫合併調度
    - **戳氣球型收益**：任何「刪除一個元素，收益與當前左右鄰有關」
      的問題（廣告位拆除、任務調度中的窗口收益最大化）
    - 矩陣鏈乘法、多邊形三角剖分、最優三角剖分本質上同構
"""

import sys
from typing import List, Optional, Sequence

NMAX = 300
VAL_LIM = 10 ** 18


def stone_merge_min(piles: Sequence[int]) -> int:
    """相鄰合併最小總代價（O(n^3) 區間 DP + 前綴和）。"""
    n = len(piles)
    if n <= 1:
        return 0
    pre = [0] * (n + 1)
    for i, v in enumerate(piles):
        pre[i + 1] = pre[i] + v
    dp = [[0] * n for _ in range(n)]
    for length in range(2, n + 1):          # 區間長度
        for i in range(0, n - length + 1):
            j = i + length - 1
            total = pre[j + 1] - pre[i]
            best = min(dp[i][k] + dp[k + 1][j] for k in range(i, j))
            dp[i][j] = best + total
    return dp[0][n - 1]


def burst_balloons_max(nums: Sequence[int]) -> int:
    """戳氣球最大收益（開區間 DP，兩端補虛擬氣球 1）。"""
    n = len(nums)
    a = [1] + list(nums) + [1]
    m = n + 2
    dp = [[0] * m for _ in range(m)]
    # len = r - l，從小到大
    for length in range(2, m):
        for l in range(0, m - length):
            r = l + length
            best = 0
            al, ar = a[l], a[r]
            for k in range(l + 1, r):
                cur = dp[l][k] + dp[k][r] + al * a[k] * ar
                if cur > best:
                    best = cur
            dp[l][r] = best
    return dp[0][m - 1]


# ---------------- 暴力（對拍用） ----------------

def stone_brute(piles: List[int]) -> int:
    """模擬所有合併順序（指數級，n 很小時才用）。"""
    if len(piles) <= 1:
        return 0
    best = VAL_LIM
    for i in range(len(piles) - 1):
        merged = piles[i] + piles[i + 1]
        rest = piles[:i] + [merged] + piles[i + 2:]
        cost = merged + stone_brute(rest)
        if cost < best:
            best = cost
    return best


def balloons_brute(nums: List[int]) -> int:
    """枚舉戳破順序（指數級）。"""
    if not nums:
        return 0
    best = 0
    for i in range(len(nums)):
        left = nums[i - 1] if i > 0 else 1
        right = nums[i + 1] if i + 1 < len(nums) else 1
        rest = nums[:i] + nums[i + 1:]
        gain = left * nums[i] * right + balloons_brute(rest)
        if gain > best:
            best = gain
    return best


# 自製 LCG（兩語言一致的隨機序列）
def run_tests() -> None:
    # ---- 合併石子固定用例 ----
    assert stone_merge_min([]) == 0
    assert stone_merge_min([7]) == 0
    assert stone_merge_min([4, 5]) == 9
    assert stone_merge_min([1, 3, 5, 2]) == 22   # 經典例
    assert stone_merge_min([3, 1, 2, 4]) == 19   # (1+2)=3 → (3+3)=6 → 10
    # ---- 戳氣球固定用例 ----
    assert burst_balloons_max([]) == 0
    assert burst_balloons_max([5]) == 5
    assert burst_balloons_max([3, 1, 5, 8]) == 167   # LeetCode 148 經典例
    assert burst_balloons_max([1, 5]) == 10          # 1*1*5 + 1*5*1

    # ---- 對拍（自製 LCG，兩語言同一序列）----
    state = 20261104

    def nxt(mod: int) -> int:
        nonlocal state
        state = (state * 1103515245 + 12345) & 0x7FFFFFFF
        return state % mod

    for _ in range(150):
        n = 1 + nxt(10)                      # 1..10
        piles = [1 + nxt(20) for _ in range(n)]
        assert stone_merge_min(piles) == stone_brute(list(piles))
    for _ in range(120):
        n = 1 + nxt(8)                       # 1..8
        nums = [nxt(10) for _ in range(n)]   # 允許 0
        assert burst_balloons_max(nums) == balloons_brute(list(nums))


# ---------------------------------------------------------------- IO 模式

def _parse_int(tok: str) -> Optional[int]:
    """整數字串解析；非純整數或超出 ±10^18 回 None。"""
    if not tok:
        return None
    i = 0
    if tok[0] in "+-":
        i = 1
    if i == len(tok):
        return None
    for ch in tok[i:]:
        if ch < "0" or ch > "9":
            return None
    val = int(tok)
    if val > VAL_LIM:
        return None
    return val


def emit(s: str) -> None:
    sys.stdout.buffer.write(s.encode("latin-1"))


def run_io(raw: str) -> None:
    toks = raw.split()
    pos = 0

    def nxt() -> Optional[str]:
        nonlocal pos
        if pos >= len(toks):
            return None
        t = toks[pos]
        pos += 1
        return t

    def bad() -> None:
        emit("error\n")

    t = _parse_int(nxt() or "")
    if t not in (1, 2):
        bad()
        return
    n = _parse_int(nxt() or "")
    if n is None or not (0 <= n <= NMAX):
        bad()
        return
    arr: List[int] = []
    for _ in range(n):
        x = _parse_int(nxt() or "")
        if x is None or x < 0 or x > 10 ** 9:
            bad()
            return
        arr.append(x)
    if t == 1:
        emit(str(stone_merge_min(arr)) + "\n")
    else:
        emit(str(burst_balloons_max(arr)) + "\n")


def main() -> None:
    raw = sys.stdin.buffer.read().decode("latin-1")
    if raw.strip():
        run_io(raw)
    else:
        run_tests()
        emit("all tests passed\n")


if __name__ == "__main__":
    main()
