# 最长递增子序列（LIS）

| 语言 | 文件 |
|---|---|
| Python 3 | [`solution.py`](solution.py) |
| C++17 | [`solution.cpp`](solution.cpp) |

## 题意

给定整数数组 `nums`，求最长的**严格递增**子序列（Longest Increasing Subsequence）的长度，并还原出一条具体方案。

子序列不要求元素在原数组中连续，但相对顺序必须保持一致。

## 思路

### 方法一：动态规划 · O(n²)

定义 `dp[i]` 为「以 `nums[i]` 作为结尾」的最长递增子序列长度：

```
dp[i] = 1 + max{ dp[j] | j < i 且 nums[j] < nums[i] }
```

不存在这样的 `j` 时 `dp[i] = 1`。答案为 `max(dp)`。
另用 `pre[i]` 记录最优前驱下标，从 `dp` 最大的位置沿 `pre` 回跳再反转，即可还原一条具体序列。

### 方法二：贪心 + 二分 · O(n log n)

维护 `tails` 数组：`tails[k]` = **长度为 k+1 的递增子序列的结尾元素的最小可能值**。

`tails` 严格递增，因此对 `x = nums[i]`，二分找到第一个 `>= x` 的位置 `pos`：

- `pos == len(tails)`：`x` 比所有已有结尾都大，可以接长，长度 +1；
- 否则：用 `x` 覆盖 `tails[pos]`。结尾越小、后续接长的潜力越大——这是贪心正确性的关键（用 `bisect_left` 而非 `bisect_right`，保证严格递增）。

> 注意：`tails` 本身**不一定**是一条合法的子序列，只有它的长度 `len(tails)` 才是正确答案。
> 需要还原具体方案时请用方法一。

## 复杂度

| 方法 | 时间复杂度 | 空间复杂度 |
|---|---|---|
| 动态规划 | O(n²) | O(n)（`dp` + `pre`） |
| 贪心 + 二分 | O(n log n) | O(n)（`tails`） |

## 输入输出

```text
输入：                      输出：
8                           4
10 9 2 5 3 7 101 18         2 5 7 101
```

- 第一行 `n`，第二行 `n` 个整数。
- 输出第一行为 LIS 长度，第二行为一条 LIS（空格分隔）。
- Python 与 C++ 版本的输入输出格式完全一致。
- 无 stdin 输入时，两个版本都会运行内置断言测试并输出 `all tests passed`。

## 运行

```bash
python solution.py
g++ -std=c++17 -O2 solution.cpp -o solution && ./solution
```
