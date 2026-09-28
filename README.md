# daily-code

每天 1–3 个（加练日更多）**算法与数据结构**题目，每题同时提供 **Python 3** 与 **C++17** 两个可直接运行/编译的版本。

目录按日期组织：`YYYY-MM-DD/<题目>/`，内含 `solution.py`、`solution.cpp`、`README.md`（题意、思路、复杂度）。

## 索引

| 日期 | 主题 | 语言 | 目录 |
|---|---|---|---|
| 2026-09-28 | Bellman-Ford 最短路（负权边 + 负环检测） | Python / C++ | [`2026-09-28/bellman-ford`](2026-09-28/bellman-ford) |
| 2026-09-28 | 单调栈（下一个更大元素 / 每日温度 / 最大矩形） | Python / C++ | [`2026-09-28/monotonic-stack`](2026-09-28/monotonic-stack) |
| 2026-09-28 | 快速排序与快速选择（三路分区 + Quickselect） | Python / C++ | [`2026-09-28/quick-sort`](2026-09-28/quick-sort) |
| 2026-09-27 | 最长公共子序列（动态规划 + 回溯） | Python / C++ | [`2026-09-27/longest-common-subsequence`](2026-09-27/longest-common-subsequence) |
| 2026-09-27 | 最小生成树（Kruskal + 并查集） | Python / C++ | [`2026-09-27/kruskal-mst`](2026-09-27/kruskal-mst) |
| 2026-09-27 | 二分查找与二分答案（lower_bound / upper_bound） | Python / C++ | [`2026-09-27/binary-search`](2026-09-27/binary-search) |
| 2026-09-23 | 拓扑排序（Kahn + DFS 逆后序） | Python / C++ | [`2026-09-23/topological-sort`](2026-09-23/topological-sort) |
| 2026-09-23 | 编辑距离（动态规划） | Python / C++ | [`2026-09-23/edit-distance`](2026-09-23/edit-distance) |
| 2026-09-23 | 归并排序与逆序对计数（分治） | Python / C++ | [`2026-09-23/merge-sort-inversion-count`](2026-09-23/merge-sort-inversion-count) |
| 2026-09-22 | 0-1 背包（动态规划） | Python / C++ | [`2026-09-22/knapsack-01`](2026-09-22/knapsack-01) |
| 2026-09-22 | 滑动窗口最大值（单调队列） | Python / C++ | [`2026-09-22/sliding-window-maximum`](2026-09-22/sliding-window-maximum) |
| 2026-09-22 | 线段树（区间和 + 懒标记） | Python / C++ | [`2026-09-22/segment-tree-range-sum`](2026-09-22/segment-tree-range-sum) |
| 2026-09-22 | 前缀树（Trie） | Python / C++ | [`2026-09-22/trie-prefix-tree`](2026-09-22/trie-prefix-tree) |
| 2026-09-22 | 最长递增子序列（LIS） | Python / C++ | [`2026-09-22/longest-increasing-subsequence`](2026-09-22/longest-increasing-subsequence) |
| 2026-09-21 | KMP 字符串匹配 | Python / C++ | [`2026-09-21/kmp-string-matching`](2026-09-21/kmp-string-matching) |
| 2026-09-21 | 并查集（Union-Find） | Python / C++ | [`2026-09-21/union-find`](2026-09-21/union-find) |
| 2026-09-21 | Dijkstra 单源最短路 | Python / C++ | [`2026-09-21/dijkstra-shortest-path`](2026-09-21/dijkstra-shortest-path) |

## 怎么跑

```bash
# Python
python 2026-09-21/union-find/solution.py

# C++
g++ -std=c++17 -O2 2026-09-21/union-find/solution.cpp -o solution && ./solution
```

> Windows / MinGW 下不要把编译产物输出到 `/tmp`（链接器会报 No such file or directory），输出到当前目录即可。

每个文件都自带断言测试，通过会输出 `all tests passed`。
