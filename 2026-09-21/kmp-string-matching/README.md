# KMP 字符串匹配

在 `text` 中查找 `pattern` 出現的位置。

## 思路

樸素匹配在失配時會把 `text` 指針回退，導致最壞 `O(n·m)`。KMP 觀察到：**失配時 `text` 中已經匹配的那段後綴，可以直接拿來當 `pattern` 的前綴復用**，因此 `text` 指針永不回退。

爲此預處理出 `pattern` 的前綴函數 `pi`：

- `pi[i]` = `pattern[0..i]` 中「最長的、既是真前綴又是後綴」的子串長度
- 失配時執行 `j = pi[j-1]`，把 `pattern` 指針滑到下一個仍可能匹配的位置

匹配階段與構造 `pi` 的過程幾乎同構，都是「雙指針 + 失配回退」。

## 複雜度

| 階段 | 時間 | 空間 |
|---|---|---|
| 構造前綴函數 | `O(m)` | `O(m)` |
| 匹配 | `O(n)` | — |
| **合計** | **`O(n+m)`** | **`O(m)`** |

`n = len(text)`，`m = len(pattern)`。雖然回退是 while 循環，但 `j` 的總增加量不超過 `n`，故攤還爲線性。

## 運行

```bash
python solution.py                                  # -> all tests passed
g++ -std=c++17 -O2 solution.cpp -o solution && ./solution
```
