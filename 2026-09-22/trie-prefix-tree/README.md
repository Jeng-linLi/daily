# 前缀树（Trie / Prefix Tree）

| 语言 | 文件 |
|---|---|
| Python 3 | [`solution.py`](solution.py) |
| C++17 | [`solution.cpp`](solution.cpp) |

## 题意

维护一个字符串集合，支持四种操作：

| 操作 | 含义 |
|---|---|
| `insert(word)` | 插入单词 |
| `search(word)` | 查询 `word` 是否被**完整**插入过 |
| `startsWith(p)` | 查询是否存在以 `p` 为前缀的单词 |
| `countPrefix(p)` | 统计已插入单词中以 `p` 为前缀的个数（重复插入计多次） |

## 思路

Trie 是一棵「按字符分叉」的多叉树：从根到任一节点的路径拼起来就是一个前缀。
每个节点维护四个字段：

- `children`：字符 → 子节点（Python 用字典，C++ 用 `unordered_map`，均支持任意字符集）；
- `isEnd`：是否有单词在此节点结束——用来区分 `"app"` 与 `"apple"`；
- `passCount`：有多少个已插入单词**经过**该节点，即该前缀的出现次数；
- `endCount`：有多少个单词恰好在此结束，支持重复插入计数。

**插入**：从根出发逐字符走，缺节点就新建，沿途 `passCount` 全部 +1，终点置 `isEnd = true` 且 `endCount += 1`。
**查询**：同样逐字符走，中途断掉即不存在；`search` 还要求终点 `isEnd` 为真。

所有操作的复杂度只与字符串长度有关，与集合中已有单词数无关——这是 Trie 相对哈希表最大的优势，也是它能高效做「前缀统计 / 自动补全」的原因。

> 代价是空间：若字符集很大且字符串稀疏，Trie 会相当吃内存（可改用压缩 Trie / 双数组 Trie 优化）。

## 复杂度

设字符串长度为 `L`：

| 操作 | 时间复杂度 | 说明 |
|---|---|---|
| `insert` | O(L) | 逐字符建/走节点 |
| `search` | O(L) | 走到终点并检查 `isEnd` |
| `startsWith` | O(L) | 能走通即可 |
| `countPrefix` | O(L) | 直接读终点 `passCount`，无需遍历子树 |

空间复杂度 O(ΣL)，即所有已插入字符串的总字符数。

## 输入输出

```text
输入：                       输出：
6                            false
insert app                   true
insert apple                 2
insert application
search appl
startsWith appl
countPrefix appl
```

- 第一行为操作条数 `n`，接下来 `n` 行每行 `操作 参数`。
- `search` / `startsWith` 各输出一行 `true` / `false`；`countPrefix` 输出一行整数；`insert` 无输出。
- Python 与 C++ 版本的输入输出格式完全一致。
- 无 stdin 输入时，两个版本都会运行内置断言测试并输出 `all tests passed`。

## 运行

```bash
python solution.py
g++ -std=c++17 -O2 solution.cpp -o solution && ./solution
```
