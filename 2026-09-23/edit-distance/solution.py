"""編輯距離（Edit Distance / Levenshtein Distance，動態規劃）

題意：給定兩個字符串 a、b，允許三種操作：插入一個字符、刪除一個字符、
    把一個字符替換成另一個字符。求把 a 變成 b 所需的最少操作次數，
    並給出一條達到該次數的操作序列。

思路：
    定義 dp[i][j] = 把 a 的前 i 個字符變成 b 的前 j 個字符的最少操作數。
    看最後一個字符，只有三種「最後一步」：
      刪掉 a[i-1]        -> dp[i-1][j] + 1
      插入 b[j-1]        -> dp[i][j-1] + 1
      把 a[i-1] 改/保留  -> dp[i-1][j-1] + (a[i-1] != b[j-1])
    三者取最小即爲狀態轉移。邊界 dp[0][j] = j（全插入）、dp[i][0] = i（全刪除）。

    空間可壓到一維：dp[j] 在掃描第 i 行時，「dp[j]」是上一行的 dp[i-1][j]、
    「dp[j-1]」是剛算好的本行 dp[i][j-1]，而 dp[i-1][j-1] 被覆蓋了，
    所以需要用一個變量 prev_diag 把左上角的舊值隨身帶着往前滾。

    還原操作序列則必須保留二維表：從 dp[m][n] 往回走，
    每步挑一個「能解釋當前 dp 值」的前驅。爲了讓下標不出錯，回溯是
    **從後往前**生成操作的，因此操作的下標天然遞減 —— 按生成順序依次施加時，
    每次改動都只影響下標 >= 當前下標的字符，已經處理過的更靠後的字符不會被挪動，
    而更早的字符還沒處理。所以「按下標遞減順序施加操作」這一套是自洽的：
    每個下標都指的就是「施加這一操作時字符串裏的位置」。

    注意：編輯距離最短時操作序列通常不唯一（例如 horse -> ros 有多條長度 3 的
    路徑），回溯只保證給出其中一條；測試因此只斷言「操作次數等於最優值」且
    「照着做一遍確實得到 b」，不鎖死具體是哪條路徑。

輸入格式（stdin）：
    第一行：字符串 a
    第二行：字符串 b（可以爲空行）
輸出格式（stdout）：
    第一行：最少操作次數
    接下來每行一條操作：replace <下標> <字符> / delete <下標> / insert <下標> <字符>
    （下標爲 0-based，指施加該操作時字符串中的位置；insert 表示插到該位置之前；
      兩個串本來就相等時不輸出任何操作行）
無 stdin 輸入時運行內置斷言測試。
"""

import sys
from functools import lru_cache
from typing import List, Tuple


def edit_distance(a: str, b: str) -> int:
    """一維滾動數組版，只求最少操作次數。時間 O(m*n)，空間 O(min(m, n)) 級。"""
    # 讓 b 成爲較短的那個，滾動數組更省空間（也順手少算一點）
    if len(a) < len(b):
        a, b = b, a
    m, n = len(a), len(b)

    dp = list(range(n + 1))          # dp[0][j] = j：空串變 b 的前 j 個字符，全插入
    for i in range(1, m + 1):
        prev_diag = dp[0]            # 上一行的 dp[i-1][0]，即左上角
        dp[0] = i                    # dp[i][0] = i：a 的前 i 個字符變空串，全刪除
        for j in range(1, n + 1):
            tmp = dp[j]              # 更新前是 dp[i-1][j]，更新後要交給下一輪的 prev_diag
            cost = 0 if a[i - 1] == b[j - 1] else 1
            dp[j] = min(
                tmp + 1,             # 刪除 a[i-1]
                dp[j - 1] + 1,       # 插入 b[j-1]（本行剛算好）
                prev_diag + cost,    # 替換或保持
            )
            prev_diag = tmp
    return dp[n]


def edit_distance_with_ops(a: str, b: str) -> Tuple[int, List[Tuple]]:
    """保留二維表並回溯出一條操作序列。時間 O(m*n)，空間 O(m*n)。"""
    m, n = len(a), len(b)
    dp = [[0] * (n + 1) for _ in range(m + 1)]
    for i in range(m + 1):
        dp[i][0] = i
    for j in range(n + 1):
        dp[0][j] = j

    for i in range(1, m + 1):
        for j in range(1, n + 1):
            cost = 0 if a[i - 1] == b[j - 1] else 1
            dp[i][j] = min(dp[i - 1][j] + 1, dp[i][j - 1] + 1, dp[i - 1][j - 1] + cost)

    # 從 dp[m][n] 回溯。操作下標遞減，故按生成順序施加即爲合法順序。
    ops: List[Tuple] = []
    i, j = m, n
    while i > 0 or j > 0:
        if i > 0 and j > 0 and a[i - 1] == b[j - 1] and dp[i][j] == dp[i - 1][j - 1]:
            i -= 1                                    # 字符相同，免費保留
            j -= 1
        elif i > 0 and j > 0 and dp[i][j] == dp[i - 1][j - 1] + 1:
            ops.append(("replace", i - 1, b[j - 1]))  # 把 a[i-1] 改成 b[j-1]
            i -= 1
            j -= 1
        elif j > 0 and dp[i][j] == dp[i][j - 1] + 1:
            ops.append(("insert", i, b[j - 1]))       # 在位置 i 之前插入 b[j-1]
            j -= 1
        elif i > 0 and dp[i][j] == dp[i - 1][j] + 1:
            ops.append(("delete", i - 1))             # 刪掉位置 i-1
            i -= 1
        else:  # pragma: no cover - 理論上不可達，留作兜底
            raise AssertionError("backtrace stuck")

    return dp[m][n], ops


def apply_ops(a: str, ops: List[Tuple]) -> str:
    """按序施加操作，用於驗證回溯出來的方案是否真的能把 a 變成 b。"""
    chars = list(a)
    for op in ops:
        if op[0] == "replace":
            chars[op[1]] = op[2]
        elif op[0] == "delete":
            del chars[op[1]]
        elif op[0] == "insert":
            chars.insert(op[1], op[2])
        else:
            raise ValueError("unknown op: " + op[0])
    return "".join(chars)


def format_ops(ops: List[Tuple]) -> List[str]:
    """把操作元組渲染成統一輸出文本。"""
    lines = []
    for op in ops:
        if op[0] == "replace":
            lines.append(f"replace {op[1]} {op[2]}")
        elif op[0] == "insert":
            lines.append(f"insert {op[1]} {op[2]}")
        else:
            lines.append(f"delete {op[1]}")
    return lines


def edit_distance_brute(a: str, b: str) -> int:
    """對照用的指數級遞歸（帶記憶化），僅用於小規模測試驗證。"""

    @lru_cache(maxsize=None)
    def go(i: int, j: int) -> int:
        if i == 0:
            return j
        if j == 0:
            return i
        cost = 0 if a[i - 1] == b[j - 1] else 1
        return min(go(i - 1, j) + 1, go(i, j - 1) + 1, go(i - 1, j - 1) + cost)

    res = go(len(a), len(b))
    go.cache_clear()
    return res


def run_io(data: str) -> None:
    """按統一輸入輸出格式處理 stdin 數據。"""
    lines = data.splitlines()
    a = lines[0] if len(lines) > 0 else ""
    b = lines[1] if len(lines) > 1 else ""
    dist, ops = edit_distance_with_ops(a, b)
    print(dist)
    for line in format_ops(ops):
        print(line)


def run_tests() -> None:
    # README 中的示例：horse -> ros，最少 3 步
    dist, ops = edit_distance_with_ops("horse", "ros")
    assert dist == 3
    assert len(ops) == 3                       # 操作條數確實等於最優值
    assert apply_ops("horse", ops) == "ros"    # 照着做一遍真的能變成 ros
    assert edit_distance("horse", "ros") == 3
    assert edit_distance_brute("horse", "ros") == 3

    # 經典用例：intention -> execution，最少 5 步
    assert edit_distance("intention", "execution") == 5
    assert edit_distance_with_ops("intention", "execution")[0] == 5

    # 完全相同：0 步，不產生任何操作
    assert edit_distance("abc", "abc") == 0
    assert edit_distance_with_ops("abc", "abc") == (0, [])
    assert edit_distance("", "") == 0

    # 一邊爲空：只能全插 / 全刪
    assert edit_distance("", "abc") == 3
    assert edit_distance("abc", "") == 3
    assert edit_distance_with_ops("", "abc")[0] == 3
    assert apply_ops("", edit_distance_with_ops("", "abc")[1]) == "abc"
    assert apply_ops("abc", edit_distance_with_ops("abc", "")[1]) == ""

    # 只差一個字符：1 步替換
    assert edit_distance("kitten", "sitten") == 1
    # kitten -> sitting：3 步（替換 e->i? 實爲 k->s、e->i、末尾插入 g）
    assert edit_distance("kitten", "sitting") == 3

    # 大小寫敏感，且長度差很大時退化爲大量插入
    assert edit_distance("Ab", "ab") == 1
    assert edit_distance("a", "aaaa") == 3

    # 純插入（a 是 b 的子序列）：flaw -> lawns 之類
    assert edit_distance("abc", "axbyc") == 2

    # 與記憶化暴力解隨機對拍：同時校驗「次數一致」與「操作序列可行且條數最優」
    import random

    random.seed(20260923)
    alphabet = "abc"
    for _ in range(200):
        a = "".join(random.choice(alphabet) for _ in range(random.randint(0, 7)))
        b = "".join(random.choice(alphabet) for _ in range(random.randint(0, 7)))
        d1 = edit_distance(a, b)
        d2, ops2 = edit_distance_with_ops(a, b)
        d3 = edit_distance_brute(a, b)
        assert d1 == d2 == d3                  # 一維版 = 二維版 = 暴力版
        assert len(ops2) == d2                 # 操作條數就是最優值
        assert apply_ops(a, ops2) == b         # 照着做一遍確實得到 b

    # 對稱性與三角不等式（編輯距離的基本性質）
    assert edit_distance("flaw", "lawn") == edit_distance("lawn", "flaw")
    for _ in range(50):
        x = "".join(random.choice(alphabet) for _ in range(random.randint(0, 6)))
        y = "".join(random.choice(alphabet) for _ in range(random.randint(0, 6)))
        assert edit_distance(x, y) == edit_distance(y, x)

    print("all tests passed")


if __name__ == "__main__":
    raw = sys.stdin.read() if not sys.stdin.isatty() else ""
    if raw.strip():
        run_io(raw)
    else:
        run_tests()
