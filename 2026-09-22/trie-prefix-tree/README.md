# 前綴樹（Trie / Prefix Tree）

| 語言 | 文件 |
|---|---|
| Python 3 | [`solution.py`](solution.py) |
| C++17 | [`solution.cpp`](solution.cpp) |

## 題意

維護一個字符串集合，支持四種操作：

| 操作 | 含義 |
|---|---|
| `insert(word)` | 插入單詞 |
| `search(word)` | 查詢 `word` 是否被**完整**插入過 |
| `startsWith(p)` | 查詢是否存在以 `p` 爲前綴的單詞 |
| `countPrefix(p)` | 統計已插入單詞中以 `p` 爲前綴的個數（重複插入計多次） |

## 思路

Trie 是一棵「按字符分叉」的多叉樹：從根到任一節點的路徑拼起來就是一個前綴。
每個節點維護四個字段：

- `children`：字符 → 子節點（Python 用字典，C++ 用 `unordered_map`，均支持任意字符集）；
- `isEnd`：是否有單詞在此節點結束——用來區分 `"app"` 與 `"apple"`；
- `passCount`：有多少個已插入單詞**經過**該節點，即該前綴的出現次數；
- `endCount`：有多少個單詞恰好在此結束，支持重複插入計數。

**插入**：從根出發逐字符走，缺節點就新建，沿途 `passCount` 全部 +1，終點置 `isEnd = true` 且 `endCount += 1`。
**查詢**：同樣逐字符走，中途斷掉即不存在；`search` 還要求終點 `isEnd` 爲真。

所有操作的複雜度只與字符串長度有關，與集合中已有單詞數無關——這是 Trie 相對哈希表最大的優勢，也是它能高效做「前綴統計 / 自動補全」的原因。

> 代價是空間：若字符集很大且字符串稀疏，Trie 會相當喫內存（可改用壓縮 Trie / 雙數組 Trie 優化）。

## 複雜度

設字符串長度爲 `L`：

| 操作 | 時間複雜度 | 說明 |
|---|---|---|
| `insert` | O(L) | 逐字符建/走節點 |
| `search` | O(L) | 走到終點並檢查 `isEnd` |
| `startsWith` | O(L) | 能走通即可 |
| `countPrefix` | O(L) | 直接讀終點 `passCount`，無需遍歷子樹 |

空間複雜度 O(ΣL)，即所有已插入字符串的總字符數。

## 輸入輸出

```text
輸入：                       輸出：
6                            false
insert app                   true
insert apple                 2
insert application
search appl
startsWith appl
countPrefix appl
```

- 第一行爲操作條數 `n`，接下來 `n` 行每行 `操作 參數`。
- `search` / `startsWith` 各輸出一行 `true` / `false`；`countPrefix` 輸出一行整數；`insert` 無輸出。
- Python 與 C++ 版本的輸入輸出格式完全一致。
- 無 stdin 輸入時，兩個版本都會運行內置斷言測試並輸出 `all tests passed`。

## 運行

```bash
python solution.py
g++ -std=c++17 -O2 solution.cpp -o solution && ./solution
```
