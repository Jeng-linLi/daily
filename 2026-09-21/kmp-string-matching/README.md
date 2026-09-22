# KMP 字符串匹配

在 `text` 中查找 `pattern` 出现的位置。

## 思路

朴素匹配在失配时会把 `text` 指针回退，导致最坏 `O(n·m)`。KMP 观察到：**失配时 `text` 中已经匹配的那段后缀，可以直接拿来当 `pattern` 的前缀复用**，因此 `text` 指针永不回退。

为此预处理出 `pattern` 的前缀函数 `pi`：

- `pi[i]` = `pattern[0..i]` 中「最长的、既是真前缀又是后缀」的子串长度
- 失配时执行 `j = pi[j-1]`，把 `pattern` 指针滑到下一个仍可能匹配的位置

匹配阶段与构造 `pi` 的过程几乎同构，都是「双指针 + 失配回退」。

## 复杂度

| 阶段 | 时间 | 空间 |
|---|---|---|
| 构造前缀函数 | `O(m)` | `O(m)` |
| 匹配 | `O(n)` | — |
| **合计** | **`O(n+m)`** | **`O(m)`** |

`n = len(text)`，`m = len(pattern)`。虽然回退是 while 循环，但 `j` 的总增加量不超过 `n`，故摊还为线性。

## 运行

```bash
python solution.py                                  # -> all tests passed
g++ -std=c++17 -O2 solution.cpp -o solution && ./solution
```
