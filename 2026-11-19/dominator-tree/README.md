# 支配樹與 Lengauer–Tarjan 演算法（Dominator Tree）

在有向圖中（從入口 `s` 出發），若「所有從 `s` 到 `v` 的路徑都必經 `u`」，則 **u 支配 v**。
`v` 的直接支配點 `idom[v]` 是離 `v` 最近的嚴格支配點；所有 `idom` 邊構成一棵以 `s` 為根的
**支配樹**。Lengauer–Tarjan 演算法 near-linear 時間求出所有點的 `idom`。

## 應用
- **編譯器控制流分析**：SSA（靜態單一賦值）形式、支配前向樹、自然循環識別。
- **程式最佳化**：不變數代碼外提（LICM）、可用表達式分析都建立在支配關係之上。

## 思路
1. **DFS 編號**：從 `s` 做 DFS，得到每個點的 dfn 與 DFS 樹父節點 `parent`。
2. **半支配點**：`semi[w]` = 所有指向 `w` 的點 `v` 中，`semi[eval(v)]` 的最小值（以 dfn 表示）。
   `eval` 用一個帶路徑壓縮的並查集維護：把 `w` 到其 DFS 祖先的鏈壓縮，鏈上 `semi` 最小者
   即為代表。這一步用 `link(parent[w], w)` 逐輪掛接。
3. **bucket 處理**：把 `w` 放進 `bucket[vertex[semi[w]]]`；當處理 DFS 父 `p` 時，對
   `bucket[p]` 中每個 `v`，令 `idom[v] = (semi[eval(v)] < semi[v]) ? eval(v) : p`。
4. **第二遍修正**：按 dfn 順序，若 `idom[w] != vertex[semi[w]]`，則
   `idom[w] = idom[idom[w]]`。
5. 根 `s` 與所有不可達點的 `idom` 設為 `-1`。

## 複雜度
- 時間：O((V + E) α(V))（近乎線性，α 為反 Ackermann 函數）。
- 空間：O(V + E)。

## 輸入輸出格式
```
n m s
u1 v1
...
um vm
```
- `n`：頂點數（0-based）；`m`：有向邊數；`s`：入口頂點。
- 之後 `m` 行 `u v`：有向邊 `u → v`。
- 輸出一行 `n` 個整數 `idom[0..n-1]`：第 `i` 個為 `i` 的直接支配點；
  `i == s` 或 `i` 從 `s` 不可達時為 `-1`。

### 範例輸入
```
6 6 0
0 1
0 2
1 3
2 3
3 4
4 5
```
### 範例輸出
```
-1 0 0 0 3 4
```

## 驗證
`python solution.py`（無 stdin）與 `solution.cpp` 無輸入執行時，皆與 O(n³) 迭代支配集合
參考演算法對拍 500 組隨機圖（含重邊），通過後輸出 `all tests passed`。
