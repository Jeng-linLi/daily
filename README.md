# daily-code

每天 1–3 個（加練日更多）**算法與數據結構**題目，每題同時提供 **Python 3** 與 **C++17** 兩個可直接運行/編譯的版本。

目錄按日期組織：`YYYY-MM-DD/<題目>/`，內含 `solution.py`、`solution.cpp`、`README.md`（題意、思路、複雜度）。

## 索引

| 日期 | 主題 | 語言 | 目錄 |
|---|---|---|---|
| 2026-09-28 | Bellman-Ford 最短路（負權邊 + 負環檢測） | Python / C++ | [`2026-09-28/bellman-ford`](2026-09-28/bellman-ford) |
| 2026-09-28 | 單調棧（下一個更大元素 / 每日溫度 / 最大矩形） | Python / C++ | [`2026-09-28/monotonic-stack`](2026-09-28/monotonic-stack) |
| 2026-09-28 | 快速排序與快速選擇（三路分區 + Quickselect） | Python / C++ | [`2026-09-28/quick-sort`](2026-09-28/quick-sort) |
| 2026-09-27 | 最長公共子序列（動態規劃 + 回溯） | Python / C++ | [`2026-09-27/longest-common-subsequence`](2026-09-27/longest-common-subsequence) |
| 2026-09-27 | 最小生成樹（Kruskal + 併查集） | Python / C++ | [`2026-09-27/kruskal-mst`](2026-09-27/kruskal-mst) |
| 2026-09-27 | 二分查找與二分答案（lower_bound / upper_bound） | Python / C++ | [`2026-09-27/binary-search`](2026-09-27/binary-search) |
| 2026-09-23 | 拓撲排序（Kahn + DFS 逆後序） | Python / C++ | [`2026-09-23/topological-sort`](2026-09-23/topological-sort) |
| 2026-09-23 | 編輯距離（動態規劃） | Python / C++ | [`2026-09-23/edit-distance`](2026-09-23/edit-distance) |
| 2026-09-23 | 歸併排序與逆序對計數（分治） | Python / C++ | [`2026-09-23/merge-sort-inversion-count`](2026-09-23/merge-sort-inversion-count) |
| 2026-09-22 | 0-1 背包（動態規劃） | Python / C++ | [`2026-09-22/knapsack-01`](2026-09-22/knapsack-01) |
| 2026-09-22 | 滑動窗口最大值（單調隊列） | Python / C++ | [`2026-09-22/sliding-window-maximum`](2026-09-22/sliding-window-maximum) |
| 2026-09-22 | 線段樹（區間和 + 懶標記） | Python / C++ | [`2026-09-22/segment-tree-range-sum`](2026-09-22/segment-tree-range-sum) |
| 2026-09-22 | 前綴樹（Trie） | Python / C++ | [`2026-09-22/trie-prefix-tree`](2026-09-22/trie-prefix-tree) |
| 2026-09-22 | 最長遞增子序列（LIS） | Python / C++ | [`2026-09-22/longest-increasing-subsequence`](2026-09-22/longest-increasing-subsequence) |
| 2026-09-21 | KMP 字符串匹配 | Python / C++ | [`2026-09-21/kmp-string-matching`](2026-09-21/kmp-string-matching) |
| 2026-09-21 | 併查集（Union-Find） | Python / C++ | [`2026-09-21/union-find`](2026-09-21/union-find) |
| 2026-09-21 | Dijkstra 單源最短路 | Python / C++ | [`2026-09-21/dijkstra-shortest-path`](2026-09-21/dijkstra-shortest-path) |

## 怎麼跑

```bash
# Python
python 2026-09-21/union-find/solution.py

# C++
g++ -std=c++17 -O2 2026-09-21/union-find/solution.cpp -o solution && ./solution
```

> Windows / MinGW 下不要把編譯產物輸出到 `/tmp`（鏈接器會報 No such file or directory），輸出到當前目錄即可。

每個文件都自帶斷言測試，通過會輸出 `all tests passed`。
