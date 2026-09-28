# 最長遞增子序列（LIS）

| 語言 | 文件 |
|---|---|
| Python 3 | [`solution.py`](solution.py) |
| C++17 | [`solution.cpp`](solution.cpp) |

## 題意

給定整數數組 `nums`，求最長的**嚴格遞增**子序列（Longest Increasing Subsequence）的長度，並還原出一條具體方案。

子序列不要求元素在原數組中連續，但相對順序必須保持一致。

## 思路

### 方法一：動態規劃 · O(n²)

定義 `dp[i]` 爲「以 `nums[i]` 作爲結尾」的最長遞增子序列長度：

```
dp[i] = 1 + max{ dp[j] | j < i 且 nums[j] < nums[i] }
```

不存在這樣的 `j` 時 `dp[i] = 1`。答案爲 `max(dp)`。
另用 `pre[i]` 記錄最優前驅下標，從 `dp` 最大的位置沿 `pre` 回跳再反轉，即可還原一條具體序列。

### 方法二：貪心 + 二分 · O(n log n)

維護 `tails` 數組：`tails[k]` = **長度爲 k+1 的遞增子序列的結尾元素的最小可能值**。

`tails` 嚴格遞增，因此對 `x = nums[i]`，二分找到第一個 `>= x` 的位置 `pos`：

- `pos == len(tails)`：`x` 比所有已有結尾都大，可以接長，長度 +1；
- 否則：用 `x` 覆蓋 `tails[pos]`。結尾越小、後續接長的潛力越大——這是貪心正確性的關鍵（用 `bisect_left` 而非 `bisect_right`，保證嚴格遞增）。

> 注意：`tails` 本身**不一定**是一條合法的子序列，只有它的長度 `len(tails)` 才是正確答案。
> 需要還原具體方案時請用方法一。

## 複雜度

| 方法 | 時間複雜度 | 空間複雜度 |
|---|---|---|
| 動態規劃 | O(n²) | O(n)（`dp` + `pre`） |
| 貪心 + 二分 | O(n log n) | O(n)（`tails`） |

## 輸入輸出

```text
輸入：                      輸出：
8                           4
10 9 2 5 3 7 101 18         2 5 7 101
```

- 第一行 `n`，第二行 `n` 個整數。
- 輸出第一行爲 LIS 長度，第二行爲一條 LIS（空格分隔）。
- Python 與 C++ 版本的輸入輸出格式完全一致。
- 無 stdin 輸入時，兩個版本都會運行內置斷言測試並輸出 `all tests passed`。

## 運行

```bash
python solution.py
g++ -std=c++17 -O2 solution.cpp -o solution && ./solution
```
