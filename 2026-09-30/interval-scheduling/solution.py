"""區間問題與貪心：最多不重疊區間 / 最少移除 / 最少會議室 / 合併區間 / 引爆氣球

題意：
    給定 n 個區間，依次輸出五個經典貪心問題的答案：
      1. 最多能選出多少個**互不重疊**的區間（區間調度，LeetCode 435 系列）；
      2. 最少移除多少個區間才能使剩下區間互不重疊 = n − 第 1 問；
      3. 同時安排所有區間最少需要幾個「會議室」（LeetCode 253）；
      4. 合併所有重疊區間後的區間個數與合併結果（LeetCode 56）；
      5. 引爆所有「氣球」最少需要幾支箭（LeetCode 452）。

思路：
    這五問的共同前提是**排序**：區間題 90% 的第一步都是「按左端點或右端點排序」，
    排序之後區間之間的關係就變成單調的，可以一次線性掃描解決。

    1. **最多不重疊區間**——按**右端點升序**貪心。
       直覺：每次都選「結束最早」的那個，給後面留下的空間最大。
       這是經典的交換論證：若最優解的第一個區間不是結束最早的，
       用結束最早的那個替換它，不會與任何後續區間衝突，解的規模不變。
    2. **最少移除** = n − 最多保留，同一個貪心的另一種問法。
    3. **最少會議室**——按左端點升序掃描，用小根堆維護「目前正在用的會議室的結束時間」。
       遇到新區間時，先把所有「結束時間 ≤ 新區間左端點」的會議室釋放（堆頂最小，彈完即止），
       再把新區間的結束時間壓入堆；堆的最大容量就是答案。
       也可以用掃描線（左端 +1、右端 −1，同一座標先處理 −1）得到同樣的結果。
    4. **合併區間**——按左端點升序掃描，若下一個區間的左端點 ≤ 當前區間的右端點就合併
       （端點相接也要合併），否則開一個新段。
    5. **引爆氣球**——仍然按右端點升序貪心，但判定條件更嚴格：
       箭射在區間右端點 x 處時，所有滿足 `l ≤ x ≤ r` 的氣球都會爆（閉區間）；
       所以只有當下一個氣球的 `l > x` 時才需要補一支新箭。
       注意它與第 1 問的區別：第 1 問把「端點相接」視為不衝突（`l ≥ last_end` 即可），
       而氣球問題是閉區間，箭落在端點上也能擊中，因此判定是嚴格的 `l > last_x`。

    「端點相接是否算重疊」是區間題最容易踩的坑，本題統一約定寫在下面：
      - 不重疊 / 會議室：半開語義 `[l, r)`，端點相接**不算**衝突，條件用 `≥` / `≤`；
      - 合併區間 / 引爆氣球：閉區間 `[l, r]`，端點相接**算**重疊，條件用 `≤` / `>`。

輸入格式（stdin，數字按空白分隔即可）：
    n
    l1 r1
    l2 r2
    …… （共 n 行；n = 0 時無後續行）
輸出格式（stdout）：
    第 1 行：最多互不重疊區間數
    第 2 行：最少移除區間數
    第 3 行：最少會議室數
    第 4 行：合併後的區間數
    第 5 行：合併後的區間，按 `l1 r1 l2 r2 …` 空格分隔（無區間時輸出空行）
    第 6 行：引爆所有氣球最少需要的箭數
無 stdin 輸入時運行內置斷言測試並輸出 `all tests passed`。
"""

import heapq
import itertools
import random
import sys
from typing import List, Tuple

Interval = Tuple[int, int]


# ---------------------------------------------------------------- 各問解法


def max_non_overlap(intervals: List[Interval]) -> int:
    """最多互不重疊區間數：按右端點升序貪心，端點相接不算衝突（l >= last_end）。"""
    cnt = 0
    last_end = None
    for l, r in sorted(intervals, key=lambda x: (x[1], x[0])):
        if last_end is None or l >= last_end:
            cnt += 1
            last_end = r
    return cnt


def max_non_overlap_dp(intervals: List[Interval]) -> int:
    """同上，但用 O(n^2) 動態規劃暴力求解，只在測試裡當貪心的對拍基準。

    dp[i] = 只看前 i 個（已按右端點排序）能選出的最多不重疊區間數。
    轉移：不選第 i 個 → dp[i-1]；選第 i 個 → 1 + dp[p]，
    其中 p = 前 i−1 個中「右端點 ≤ l_i」的個數（這些都與第 i 個不衝突）。
    """
    v = sorted(intervals, key=lambda x: (x[1], x[0]))
    n = len(v)
    dp = [0] * (n + 1)
    for i in range(1, n + 1):
        l, _ = v[i - 1]
        p = sum(1 for j in range(i - 1) if v[j][1] <= l)
        dp[i] = max(dp[i - 1], 1 + dp[p])
    return dp[n]


def min_removed(intervals: List[Interval]) -> int:
    """最少移除多少個區間才能使剩下的互不重疊。"""
    return len(intervals) - max_non_overlap(intervals)


def min_meeting_rooms_heap(intervals: List[Interval]) -> int:
    """最少會議室數：按左端點掃描 + 小根堆維護最早空出來的房間。"""
    heap: List[int] = []
    busy = 0
    for l, r in sorted(intervals, key=lambda x: (x[0], x[1])):
        if l == r:
            continue                     # 半開語義下 [l, l) 是空區間，不佔用房間
        while heap and heap[0] <= l:     # 半開語義：r <= l 表示房間已空出
            heapq.heappop(heap)
        heapq.heappush(heap, r)
        busy = max(busy, len(heap))
    return busy


def min_meeting_rooms_sweep(intervals: List[Interval]) -> int:
    """最少會議室數：掃描線。同一座標先處理結束（-1）再處理開始（+1）。"""
    events = []
    for l, r in intervals:
        if l == r:
            continue
        events.append((l, 1))
        events.append((r, -1))
    events.sort(key=lambda e: (e[0], e[1]))   # -1 < 1，同座標先退房再入住
    cur = best = 0
    for _, delta in events:
        cur += delta
        best = max(best, cur)
    return best


def merge_intervals(intervals: List[Interval]) -> List[Interval]:
    """合併所有重疊（含端點相接）的區間，回傳按左端點升序的不重疊區間列表。"""
    merged: List[Interval] = []
    for l, r in sorted(intervals, key=lambda x: (x[0], x[1])):
        if merged and l <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], r))
        else:
            merged.append((l, r))
    return merged


def min_arrows(intervals: List[Interval]) -> int:
    """引爆所有氣球的最少箭數（閉區間）：按右端點升序，箭射在右端點，l > x 才補新箭。"""
    arrows = 0
    last_x = None
    for l, r in sorted(intervals, key=lambda x: (x[1], x[0])):
        if last_x is None or l > last_x:
            arrows += 1
            last_x = r
    return arrows


def min_arrows_bruteforce(intervals: List[Interval]) -> int:
    """同上，但枚舉所有候選位置的子集求最優解，只用於測試對拍（座標範圍必須很小）。"""
    if not intervals:
        return 0
    # 最優箭的位置一定可以取在某個區間的端點上（交換論證），候選集取全部端點即可
    cand = sorted({p for l, r in intervals for p in (l, r)})
    for size in range(1, len(cand) + 1):
        for combo in itertools.combinations(cand, size):
            if all(any(l <= x <= r for x in combo) for l, r in intervals):
                return size
    return len(intervals)


# ---------------------------------------------------------------- 輸入輸出


def run_io(data: str) -> None:
    """按題目格式解析 stdin 並輸出結果。"""
    tokens = list(map(int, data.split()))
    if not tokens:
        n = 0
        rest: List[int] = []
    else:
        n, rest = tokens[0], tokens[1:]

    intervals: List[Interval] = []
    for i in range(min(n, len(rest) // 2)):
        l, r = rest[2 * i], rest[2 * i + 1]
        if l > r:
            l, r = r, l
        intervals.append((l, r))

    merged = merge_intervals(intervals)
    out = [
        str(max_non_overlap(intervals)),
        str(min_removed(intervals)),
        str(min_meeting_rooms_heap(intervals)),
        str(len(merged)),
        " ".join(f"{l} {r}" for l, r in merged),
        str(min_arrows(intervals)),
    ]
    sys.stdout.write("\n".join(out) + "\n")


# ---------------------------------------------------------------- 內置測試


def run_tests() -> None:
    # ---- 固定用例 ----
    a = [(1, 2), (2, 3), (3, 4), (1, 3)]
    assert max_non_overlap(a) == 3            # [1,2] [2,3] [3,4] 端點相接不算衝突
    assert min_removed(a) == 1
    assert min_meeting_rooms_heap(a) == 2
    assert min_meeting_rooms_sweep(a) == 2
    assert merge_intervals(a) == [(1, 4)]     # 閉區間合併：全部連成一段
    assert min_arrows(a) == 2                 # 箭射在 2（打掉前兩個與 [1,3]）、再射 4

    b = [(1, 4), (2, 3), (3, 5), (7, 9)]
    assert max_non_overlap(b) == 3            # [2,3] [3,5] [7,9]（端點相接不算衝突）
    assert min_removed(b) == 1
    assert min_meeting_rooms_heap(b) == 2
    assert merge_intervals(b) == [(1, 5), (7, 9)]
    assert min_arrows(b) == 2

    assert max_non_overlap([]) == 0
    assert min_removed([]) == 0
    assert min_meeting_rooms_heap([]) == 0
    assert merge_intervals([]) == []
    assert min_arrows([]) == 0
    assert min_arrows_bruteforce([]) == 0

    # 單點區間
    c = [(5, 5), (5, 5), (6, 6)]
    assert max_non_overlap(c) == 3            # [5,5] 與 [5,5] 半開語義下是空集，互不衝突
    assert min_meeting_rooms_heap(c) == 0     # 半開語義下 [l, l) 是空區間，不佔用會議室
    assert merge_intervals(c) == [(5, 5), (6, 6)]
    assert min_arrows(c) == 2                 # 一箭射 5、一箭射 6

    random.seed(20260930)
    for _ in range(600):
        n = random.randint(0, 8)
        iv = []
        for _ in range(n):
            l = random.randint(0, 8)
            r = l + random.randint(0, 4)
            iv.append((l, r))

        # 貪心 vs O(n^2) DP
        assert max_non_overlap(iv) == max_non_overlap_dp(iv)
        assert min_removed(iv) == n - max_non_overlap(iv)

        # 會議室：堆版 vs 掃描線版
        assert min_meeting_rooms_heap(iv) == min_meeting_rooms_sweep(iv)

        # 合併：結果必須有序、互不重疊、且覆蓋的點集與原區間完全一致
        mg = merge_intervals(iv)

        def covered(ivs: List[Interval], x: int) -> bool:
            return any(l <= x <= r for l, r in ivs)

        assert mg == sorted(mg)
        for i in range(1, len(mg)):
            assert mg[i - 1][1] < mg[i][0]     # 合併結果兩兩之間必須有空隙
        for x in range(-1, 15):
            assert covered(mg, x) == covered(iv, x)

        # 氣球：貪心 vs 子集枚舉（最優性驗證）
        ans = min_arrows(iv)
        assert min_arrows_bruteforce(iv) == ans
        # 貪心給出的箭數一定不多於「最多不重疊區間數」
        assert ans <= max_non_overlap(iv)
        # 用貪心構造出的箭位置確實能打掉全部氣球
        xs: List[int] = []
        last = None
        for l, r in sorted(iv, key=lambda x: (x[1], x[0])):
            if last is None or l > last:
                xs.append(r)
                last = r
        assert len(xs) == ans
        assert all(any(l <= x <= r for x in xs) for l, r in iv)

    print("all tests passed")


if __name__ == "__main__":
    raw = sys.stdin.read() if not sys.stdin.isatty() else ""
    if raw.strip():
        run_io(raw)
    else:
        run_tests()
