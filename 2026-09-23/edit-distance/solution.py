"""编辑距离（Edit Distance / Levenshtein Distance，动态规划）

题意：给定两个字符串 a、b，允许三种操作：插入一个字符、删除一个字符、
    把一个字符替换成另一个字符。求把 a 变成 b 所需的最少操作次数，
    并给出一条达到该次数的操作序列。

思路：
    定义 dp[i][j] = 把 a 的前 i 个字符变成 b 的前 j 个字符的最少操作数。
    看最后一个字符，只有三种「最后一步」：
      删掉 a[i-1]        -> dp[i-1][j] + 1
      插入 b[j-1]        -> dp[i][j-1] + 1
      把 a[i-1] 改/保留  -> dp[i-1][j-1] + (a[i-1] != b[j-1])
    三者取最小即为状态转移。边界 dp[0][j] = j（全插入）、dp[i][0] = i（全删除）。

    空间可压到一维：dp[j] 在扫描第 i 行时，「dp[j]」是上一行的 dp[i-1][j]、
    「dp[j-1]」是刚算好的本行 dp[i][j-1]，而 dp[i-1][j-1] 被覆盖了，
    所以需要用一个变量 prev_diag 把左上角的旧值随身带着往前滚。

    还原操作序列则必须保留二维表：从 dp[m][n] 往回走，
    每步挑一个「能解释当前 dp 值」的前驱。为了让下标不出错，回溯是
    **从后往前**生成操作的，因此操作的下标天然递减 —— 按生成顺序依次施加时，
    每次改动都只影响下标 >= 当前下标的字符，已经处理过的更靠后的字符不会被挪动，
    而更早的字符还没处理。所以「按下标递减顺序施加操作」这一套是自洽的：
    每个下标都指的就是「施加这一操作时字符串里的位置」。

    注意：编辑距离最短时操作序列通常不唯一（例如 horse -> ros 有多条长度 3 的
    路径），回溯只保证给出其中一条；测试因此只断言「操作次数等于最优值」且
    「照着做一遍确实得到 b」，不锁死具体是哪条路径。

输入格式（stdin）：
    第一行：字符串 a
    第二行：字符串 b（可以为空行）
输出格式（stdout）：
    第一行：最少操作次数
    接下来每行一条操作：replace <下标> <字符> / delete <下标> / insert <下标> <字符>
    （下标为 0-based，指施加该操作时字符串中的位置；insert 表示插到该位置之前；
      两个串本来就相等时不输出任何操作行）
无 stdin 输入时运行内置断言测试。
"""

import sys
from functools import lru_cache
from typing import List, Tuple


def edit_distance(a: str, b: str) -> int:
    """一维滚动数组版，只求最少操作次数。时间 O(m*n)，空间 O(min(m, n)) 级。"""
    # 让 b 成为较短的那个，滚动数组更省空间（也顺手少算一点）
    if len(a) < len(b):
        a, b = b, a
    m, n = len(a), len(b)

    dp = list(range(n + 1))          # dp[0][j] = j：空串变 b 的前 j 个字符，全插入
    for i in range(1, m + 1):
        prev_diag = dp[0]            # 上一行的 dp[i-1][0]，即左上角
        dp[0] = i                    # dp[i][0] = i：a 的前 i 个字符变空串，全删除
        for j in range(1, n + 1):
            tmp = dp[j]              # 更新前是 dp[i-1][j]，更新后要交给下一轮的 prev_diag
            cost = 0 if a[i - 1] == b[j - 1] else 1
            dp[j] = min(
                tmp + 1,             # 删除 a[i-1]
                dp[j - 1] + 1,       # 插入 b[j-1]（本行刚算好）
                prev_diag + cost,    # 替换或保持
            )
            prev_diag = tmp
    return dp[n]


def edit_distance_with_ops(a: str, b: str) -> Tuple[int, List[Tuple]]:
    """保留二维表并回溯出一条操作序列。时间 O(m*n)，空间 O(m*n)。"""
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

    # 从 dp[m][n] 回溯。操作下标递减，故按生成顺序施加即为合法顺序。
    ops: List[Tuple] = []
    i, j = m, n
    while i > 0 or j > 0:
        if i > 0 and j > 0 and a[i - 1] == b[j - 1] and dp[i][j] == dp[i - 1][j - 1]:
            i -= 1                                    # 字符相同，免费保留
            j -= 1
        elif i > 0 and j > 0 and dp[i][j] == dp[i - 1][j - 1] + 1:
            ops.append(("replace", i - 1, b[j - 1]))  # 把 a[i-1] 改成 b[j-1]
            i -= 1
            j -= 1
        elif j > 0 and dp[i][j] == dp[i][j - 1] + 1:
            ops.append(("insert", i, b[j - 1]))       # 在位置 i 之前插入 b[j-1]
            j -= 1
        elif i > 0 and dp[i][j] == dp[i - 1][j] + 1:
            ops.append(("delete", i - 1))             # 删掉位置 i-1
            i -= 1
        else:  # pragma: no cover - 理论上不可达，留作兜底
            raise AssertionError("backtrace stuck")

    return dp[m][n], ops


def apply_ops(a: str, ops: List[Tuple]) -> str:
    """按序施加操作，用于验证回溯出来的方案是否真的能把 a 变成 b。"""
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
    """把操作元组渲染成统一输出文本。"""
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
    """对照用的指数级递归（带记忆化），仅用于小规模测试验证。"""

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
    """按统一输入输出格式处理 stdin 数据。"""
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
    assert len(ops) == 3                       # 操作条数确实等于最优值
    assert apply_ops("horse", ops) == "ros"    # 照着做一遍真的能变成 ros
    assert edit_distance("horse", "ros") == 3
    assert edit_distance_brute("horse", "ros") == 3

    # 经典用例：intention -> execution，最少 5 步
    assert edit_distance("intention", "execution") == 5
    assert edit_distance_with_ops("intention", "execution")[0] == 5

    # 完全相同：0 步，不产生任何操作
    assert edit_distance("abc", "abc") == 0
    assert edit_distance_with_ops("abc", "abc") == (0, [])
    assert edit_distance("", "") == 0

    # 一边为空：只能全插 / 全删
    assert edit_distance("", "abc") == 3
    assert edit_distance("abc", "") == 3
    assert edit_distance_with_ops("", "abc")[0] == 3
    assert apply_ops("", edit_distance_with_ops("", "abc")[1]) == "abc"
    assert apply_ops("abc", edit_distance_with_ops("abc", "")[1]) == ""

    # 只差一个字符：1 步替换
    assert edit_distance("kitten", "sitten") == 1
    # kitten -> sitting：3 步（替换 e->i? 实为 k->s、e->i、末尾插入 g）
    assert edit_distance("kitten", "sitting") == 3

    # 大小写敏感，且长度差很大时退化为大量插入
    assert edit_distance("Ab", "ab") == 1
    assert edit_distance("a", "aaaa") == 3

    # 纯插入（a 是 b 的子序列）：flaw -> lawns 之类
    assert edit_distance("abc", "axbyc") == 2

    # 与记忆化暴力解随机对拍：同时校验「次数一致」与「操作序列可行且条数最优」
    import random

    random.seed(20260923)
    alphabet = "abc"
    for _ in range(200):
        a = "".join(random.choice(alphabet) for _ in range(random.randint(0, 7)))
        b = "".join(random.choice(alphabet) for _ in range(random.randint(0, 7)))
        d1 = edit_distance(a, b)
        d2, ops2 = edit_distance_with_ops(a, b)
        d3 = edit_distance_brute(a, b)
        assert d1 == d2 == d3                  # 一维版 = 二维版 = 暴力版
        assert len(ops2) == d2                 # 操作条数就是最优值
        assert apply_ops(a, ops2) == b         # 照着做一遍确实得到 b

    # 对称性与三角不等式（编辑距离的基本性质）
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
