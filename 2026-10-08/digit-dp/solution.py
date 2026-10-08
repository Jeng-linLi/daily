"""數位 DP（Digit DP）：數位和計數 / 禁用數字 / 吉利數字 / 數位和總和 / Windy Number / 不要 62

題意：
    給定上界 n（0 ≤ n ≤ 10^15）、目標數位和 k、以及一個禁用數字 d，統計 [0, n] 內滿足
    各種「十進位表示上的性質」的整數個數：
      1. **數位和等於 k** 的數有幾個；
      2. **不含數字 d** 的數有幾個（經典題：不含 4 的門牌號）；
      3. **Lucky Number**：每一位都只由 4 或 7 組成的正整數有幾個；
      4. **數位和總和**：Σ_{i=0}^{n} digitsum(i)；
      5. **Windy Number**：相鄰兩位數字差的絕對值都 ≥ 2 的正整數有幾個；
      6. **不要 62**：不含數字 4、且不含連續子串 "62" 的數有幾個（經典題 HDU 2089）。

思路：
    ### 為什麼需要數位 DP
    直接枚舉 [0, n] 是 O(n)，n = 10^15 時不可能。數位 DP 的核心洞察是：
    **按位從高到低構造數字**，並且只關心「前綴是否已經小於 n 的對應前綴」這一個狀態。

    ### 貼緊（tight）與鬆散
    把 n 補成 L 位（前導零補齊），每個 [0, n] 的數都對應唯一一個 L 位的數字串。
    構造到某一位時有兩種情形：
      - **貼緊（tight = 1）**：此前的每一位都與 n 完全相同 → 本位上限是 n 的這一位 `ds[pos]`；
      - **鬆散（tight = 0）**：此前某一位已經比 n 小 → 本位可以隨便填 0..9。
    於是把「計數」問題轉成在 O(L × 10 × 狀態數) 的表格上做轉移，**複雜度只和位數有關**。

    ### 前導零
    數字 0 的合法表示是 "0"，更短的數要用前導零補齊（例如 n = 1000 時數字 7 要寫成 "0007"）。
    是否需要「是否已經開始（started）」這個狀態，取決於題目怎麼看待 0：
      - 數位和、數位和總和：前導零貢獻 0，無影響 → **不需要** started；
      - Lucky / Windy：題目只統計正整數 → 用「上一位是什麼」或 started 把全零串排除；
      - 不要 62：0 本身合法 → 全零串必須保留。

    ### 各題的狀態設計
    | 題目 | 狀態 | 轉移時的限制 |
    |---|---|---|
    | 數位和 = k | `(sum, tight)` | `sum + d ≤ k` 才轉移（剪枝） |
    | 不含數字 d | 用「集合型」DP：`(started, tight)` | `d` 不在允許集合就跳過 |
    | Lucky | `(started, tight)`，允許集合 = {4, 7} | 只填 4 / 7 |
    | 數位和總和 | `(tight)`，每格存 `(個數, 數位和總和)` | 加一位 d：總和 += d × 個數 |
    | Windy | `(上一位, tight)`，上一位 = 10 代表尚未開始 | `abs(d - 上一位) ≥ 2` |
    | 不要 62 | `(上一位是否為 6, tight)` | 跳過 d = 4；上一位是 6 時跳過 d = 2 |

    ### 數位和總和為什麼不用按 sum 分桶
    加一位數字 d 時，所有數的數位和都 + d，因此總和的更新是
        total' = total + d × count
    只要同時維護「個數」與「數位和總和」兩個量即可，不必枚舉 sum → 狀態數從 O(9L) 降到 1。

    ### 確定性（保證 Python 與 C++ 輸出逐字節一致）
    狀態轉移是純計數（加法），與遍歷順序無關，兩個語言結果必然相同。

應用場景：
    門牌 / 車牌 / 手機號的吉利號碼統計（不要 62、不要 4）、
    密碼學與編碼中的數字範圍計數、區間 [L, R] 統計（答案 = f(R) - f(L-1)，本題的 f 形式）、
    競賽中的「第 k 個滿足數位性質的數」（數位 DP + 二分）、大範圍的數字性質篩選。

複雜度：
    所有六問皆為 O(L × 10 × 狀態數) 時間、O(狀態數) 空間，L = 位數 ≤ 16（n ≤ 10^15）。
    暴力對拍 O(n × L)，僅測試用（n ≤ 3000）。

輸入格式（stdin，全部以空白分隔）：
    n   上界（負數視為 0；上限截到 10^15）
    k   目標數位和（負數視為 0）
    d   禁用數字（0..9；超出範圍表示不禁用任何數字）
    讀不到時缺的部分補 0；遇到非整數 token 視為輸入結束。
輸出格式（stdout）：
    第 1 行：[0, n] 中數位和等於 k 的數的個數
    第 2 行：[0, n] 中不含數字 d 的數的個數
    第 3 行：[1, n] 中的 Lucky Number 個數
    第 4 行：Σ_{i=0}^{n} digitsum(i)
    第 5 行：[1, n] 中的 Windy Number 個數
    第 6 行：[0, n] 中不含 4 且不含 "62" 的數的個數
無 stdin 輸入時運行內置斷言測試並輸出 `all tests passed`。
"""

import random
import re
import sys
from typing import Dict, List, Optional, Tuple

NMAX = 10 ** 15                              # n 的上限（與 C++ 一致，避免溢位）
_INT_RE = re.compile(r"^[+-]?[0-9]+$")      # 嚴格整數規則，與 C++ 的 tryLL 完全一致
NO_DIGIT = 10                                # 「尚未開始（全是前導零）」的哨兵值


def parse_int(tok: Optional[str]) -> Optional[int]:
    """把 token 轉成整數；非十進制整數則返回 None（視為輸入到此為止）。"""
    if tok is None or not _INT_RE.match(tok):
        return None
    return int(tok)


def digits_of(n: int) -> List[int]:
    """把非負整數拆成十進位數字串（高位在前）。"""
    return [int(c) for c in str(n)]


def digit_sum(x: int) -> int:
    """十進位數位和。"""
    return sum(int(c) for c in str(x))


# ---------------------------------------------------------------- 數位 DP


def count_digit_sum(n: int, k: int) -> int:
    """統計 [0, n] 中數位和等於 k 的數的個數。狀態 (sum, tight)。"""
    if n < 0 or k < 0:
        return 0
    ds = digits_of(n)
    L = len(ds)
    if k > 9 * L:
        return 0
    dp: Dict[Tuple[int, bool], int] = {(0, True): 1}
    for pos in range(L):
        nd: Dict[Tuple[int, bool], int] = {}
        for (s, tight), cnt in dp.items():
            lim = ds[pos] if tight else 9
            for d in range(lim + 1):
                ns = s + d
                if ns > k:
                    continue                 # 剪枝：數位和已超過目標
                key = (ns, tight and d == lim)
                nd[key] = nd.get(key, 0) + cnt
        dp = nd
    return sum(cnt for (s, _), cnt in dp.items() if s == k)


def count_digit_set(n: int, allowed: Tuple[int, ...]) -> int:
    """統計 [1, n] 中每一位都屬於 allowed 的正整數個數（不含 0）。狀態 (started, tight)。"""
    if n <= 0:
        return 0
    ds = digits_of(n)
    aset = set(allowed)
    dp: Dict[Tuple[bool, bool], int] = {(False, True): 1}
    for pos in range(len(ds)):
        nd: Dict[Tuple[bool, bool], int] = {}
        for (started, tight), cnt in dp.items():
            lim = ds[pos] if tight else 9
            for d in range(lim + 1):
                ntight = tight and d == lim
                if not started and d == 0:
                    key = (False, ntight)    # 仍是前導零，數字尚未開始
                    nd[key] = nd.get(key, 0) + cnt
                elif d in aset:
                    key = (True, ntight)
                    nd[key] = nd.get(key, 0) + cnt
        dp = nd
    return sum(cnt for (started, _), cnt in dp.items() if started)


def count_lucky(n: int) -> int:
    """[1, n] 中的 Lucky Number（每一位都是 4 或 7）個數。"""
    return count_digit_set(n, (4, 7))


def count_avoid_digit(n: int, forbid: int) -> int:
    """[0, n] 中不含數字 forbid 的數的個數。

    數字 0 的表示是 "0"：只有當禁用的是 0 時才把 0 排除，其餘情況 0 都要計入。
    """
    if n < 0:
        return 0
    allowed = tuple(x for x in range(10) if x != forbid)
    return count_digit_set(n, allowed) + (0 if forbid == 0 else 1)


def digit_sum_total(n: int) -> int:
    """Σ_{i=0}^{n} digitsum(i)。狀態只需 tight，每格存 (個數, 數位和總和)。"""
    if n < 0:
        return 0
    ds = digits_of(n)
    dp: Dict[bool, Tuple[int, int]] = {True: (1, 0)}
    for pos in range(len(ds)):
        nd: Dict[bool, Tuple[int, int]] = {}
        for tight, (cnt, tot) in dp.items():
            lim = ds[pos] if tight else 9
            for d in range(lim + 1):
                ntight = tight and d == lim
                ocnt, otot = nd.get(ntight, (0, 0))
                nd[ntight] = (ocnt + cnt, otot + tot + d * cnt)
        dp = nd
    return sum(tot for (_, tot) in dp.values())


def count_windy(n: int) -> int:
    """[1, n] 中的 Windy Number（相鄰位數字差 ≥ 2）個數。狀態 (上一位, tight)。"""
    if n < 0:
        return 0
    ds = digits_of(n)
    dp: Dict[Tuple[int, bool], int] = {(NO_DIGIT, True): 1}
    for pos in range(len(ds)):
        nd: Dict[Tuple[int, bool], int] = {}
        for (prev, tight), cnt in dp.items():
            lim = ds[pos] if tight else 9
            for d in range(lim + 1):
                ntight = tight and d == lim
                if prev == NO_DIGIT and d == 0:
                    key = (NO_DIGIT, ntight)     # 尚未開始，仍是前導零
                else:
                    if prev != NO_DIGIT and abs(d - prev) < 2:
                        continue
                    key = (d, ntight)
                nd[key] = nd.get(key, 0) + cnt
        dp = nd
    return sum(cnt for (prev, _), cnt in dp.items() if prev != NO_DIGIT)


def count_auspicious(n: int) -> int:
    """[0, n] 中不含數字 4、且不含連續 "62" 的數的個數（經典題「不要 62」）。"""
    if n < 0:
        return 0
    ds = digits_of(n)
    dp: Dict[Tuple[bool, bool], int] = {(False, True): 1}
    for pos in range(len(ds)):
        nd: Dict[Tuple[bool, bool], int] = {}
        for (prev6, tight), cnt in dp.items():
            lim = ds[pos] if tight else 9
            for d in range(lim + 1):
                if d == 4:
                    continue                     # 禁用數字 4
                if prev6 and d == 2:
                    continue                     # 禁用連續子串 "62"
                ntight = tight and d == lim
                key = (d == 6, ntight)
                nd[key] = nd.get(key, 0) + cnt
        dp = nd
    return sum(dp.values())


# ---------------------------------------------------------------- 暴力對拍


def bf_digit_sum(n: int, k: int) -> int:
    return sum(1 for i in range(n + 1) if digit_sum(i) == k)


def bf_avoid(n: int, forbid: int) -> int:
    if not (0 <= forbid <= 9):
        return n + 1                          # 禁用值不是一個數字 → 不禁用任何數字
    return sum(1 for i in range(n + 1) if str(forbid) not in str(i))


def bf_lucky(n: int) -> int:
    return sum(1 for i in range(1, n + 1) if all(c in "47" for c in str(i)))


def bf_sum_digits(n: int) -> int:
    return sum(digit_sum(i) for i in range(n + 1))


def bf_windy(n: int) -> int:
    cnt = 0
    for i in range(1, n + 1):
        st = str(i)
        if all(abs(int(st[j]) - int(st[j + 1])) >= 2 for j in range(len(st) - 1)):
            cnt += 1
    return cnt


def bf_auspicious(n: int) -> int:
    return sum(1 for i in range(n + 1) if "4" not in str(i) and "62" not in str(i))


# ---------------------------------------------------------------- IO 與測試


def run_io(raw: str) -> None:
    toks = raw.split()
    pos = 0

    def nxt() -> int:
        nonlocal pos
        v = parse_int(toks[pos]) if pos < len(toks) else None
        pos += 1
        if v is None:
            return 0
        if v > NMAX:
            return NMAX
        if v < -NMAX:
            return -NMAX
        return v

    n = nxt()
    k = nxt()
    d = nxt()
    if n < 0:
        n = 0
    if k < 0:
        k = 0

    out = [
        str(count_digit_sum(n, k)),
        str(count_avoid_digit(n, d)),
        str(count_lucky(n)),
        str(digit_sum_total(n)),
        str(count_windy(n)),
        str(count_auspicious(n)),
    ]
    sys.stdout.write("\n".join(out) + "\n")


def run_tests() -> None:
    # ---- 退化情形 ----
    assert count_digit_sum(0, 0) == 1                 # 只有數字 0
    assert count_digit_sum(0, 1) == 0
    assert count_digit_sum(-1, 0) == 0
    assert count_digit_sum(9, 5) == 1
    assert count_avoid_digit(0, 4) == 1               # [0,0] 不含 4 → 1 個
    assert count_avoid_digit(0, 0) == 0               # 0 含數字 0 → 0 個
    assert count_avoid_digit(-1, 4) == 0
    assert count_lucky(0) == 0
    assert count_lucky(3) == 0
    assert count_lucky(4) == 1
    assert count_lucky(7) == 2
    assert count_lucky(47) == 4                       # 4, 7, 44, 47
    assert digit_sum_total(0) == 0
    assert digit_sum_total(9) == 45
    assert count_windy(0) == 0
    assert count_windy(9) == 9                        # 一位數都算 windy
    assert count_auspicious(0) == 1                   # 0 不含 4 也不含 62

    # ---- 固定用例 ----
    assert count_digit_sum(100, 1) == 3               # 1, 10, 100
    assert count_digit_sum(100, 5) == 6
    assert count_avoid_digit(100, 4) == 82            # 0..100 去掉含 4 的（19 個）
    assert digit_sum_total(100) == 901
    assert count_windy(100) == 73
    assert count_auspicious(100) == 81                # 去掉含 4 與含 62 的

    # ---- 與暴力對拍（n ≤ 3000）----
    for n in list(range(0, 120)) + [199, 200, 999, 1000, 1999, 2000, 2999, 3000]:
        for k in range(0, 29):
            assert count_digit_sum(n, k) == bf_digit_sum(n, k), (n, k)
        for d in range(0, 11):                        # 含 d = 10（不禁用任何數字）
            assert count_avoid_digit(n, d) == bf_avoid(n, d), (n, d)
        assert count_lucky(n) == bf_lucky(n), n
        assert digit_sum_total(n) == bf_sum_digits(n), n
        assert count_windy(n) == bf_windy(n), n
        assert count_auspicious(n) == bf_auspicious(n), n

    # ---- 隨機對拍 ----
    rnd = random.Random(20261008)
    for _ in range(60):
        n = rnd.randint(0, 3000)
        k = rnd.randint(0, 30)
        d = rnd.randint(0, 10)
        assert count_digit_sum(n, k) == bf_digit_sum(n, k)
        assert count_avoid_digit(n, d) == bf_avoid(n, d)
        assert count_lucky(n) == bf_lucky(n)
        assert digit_sum_total(n) == bf_sum_digits(n)
        assert count_windy(n) == bf_windy(n)
        assert count_auspicious(n) == bf_auspicious(n)

    # ---- 大數自洽性（無法暴力，檢查單調性與區間可減性）----
    big = 10 ** 15
    assert count_digit_sum(big, 0) == 1               # 只有 0
    assert count_digit_sum(big, 9 * 15) == 1          # 只有 999999999999999
    assert count_lucky(big) == 2 ** 16 - 2            # 長度 1..15 的 4/7 串
    assert digit_sum_total(big - 1) == 15 * 45 * 10 ** 14   # 15 位補齊後每位平均 4.5
    assert digit_sum_total(big) == 15 * 45 * 10 ** 14 + 1   # 再加上 10^15 本身的數位和 1
    assert count_avoid_digit(big, 4) == count_digit_set(big, tuple(x for x in range(10) if x != 4)) + 1
    # 區間可減性：[L, R] = f(R) - f(L-1)
    L0, R0 = 12345, 98765
    assert count_windy(R0) - count_windy(L0 - 1) == sum(
        1 for i in range(L0, R0 + 1)
        if all(abs(int(str(i)[j]) - int(str(i)[j + 1])) >= 2 for j in range(len(str(i)) - 1))
    )
    assert count_auspicious(R0) - count_auspicious(L0 - 1) == sum(
        1 for i in range(L0, R0 + 1) if "4" not in str(i) and "62" not in str(i)
    )


if __name__ == "__main__":
    data = sys.stdin.read()
    if data.strip():
        run_io(data)
    else:
        run_tests()
        print("all tests passed")
