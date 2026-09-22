"""前缀树（Trie / Prefix Tree）

题意：维护一个字符串集合，支持三种操作：
    insert(word)     插入一个单词
    search(word)     查询单词是否被完整插入过
    starts_with(p)   查询是否存在以 p 为前缀的单词
    count_prefix(p)  统计已插入单词中以 p 为前缀的个数

思路：
    Trie 是一棵「按字符分叉」的多叉树，从根到某个节点的路径拼成一个前缀。
    每个节点维护：
      children：字符 -> 子节点（用字典，天然支持任意字符集）
      is_end：  是否有单词在此节点结束（区分 "app" 与 "apple"）
      pass：    有多少个已插入单词经过该节点（即该前缀出现次数）
      end：     有多少个单词恰好在此节点结束（支持重复插入计数）
    插入时沿字符逐层走，不存在就新建节点，沿途 pass 都 +1；结束时置 is_end。
    查找复杂度只与单词长度有关，与集合规模无关，这是它相对哈希表的最大优势。

输入格式（stdin）：
    第一行 n（操作条数）
    接下来 n 行，每行：insert <word> | search <word> | startsWith <prefix> | countPrefix <prefix>
输出格式（stdout）：
    search / startsWith 各输出一行 true / false；countPrefix 输出一行整数；insert 无输出
无 stdin 输入时运行内置断言测试。
"""

import sys
from typing import Dict, List


class TrieNode:
    __slots__ = ("children", "is_end", "pass_count", "end_count")

    def __init__(self) -> None:
        self.children: Dict[str, "TrieNode"] = {}  # 字符 -> 子节点
        self.is_end: bool = False                  # 是否有单词在此结束
        self.pass_count: int = 0                   # 经过该节点的单词数
        self.end_count: int = 0                    # 在此结束的单词数（含重复插入）


class Trie:
    def __init__(self) -> None:
        self.root = TrieNode()

    def insert(self, word: str) -> None:
        """插入单词，时间 O(|word|)。"""
        node = self.root
        node.pass_count += 1
        for ch in word:
            if ch not in node.children:
                node.children[ch] = TrieNode()
            node = node.children[ch]
            node.pass_count += 1
        node.is_end = True
        node.end_count += 1

    def _walk(self, s: str):
        """沿字符串走到对应节点；中途断掉返回 None。"""
        node = self.root
        for ch in s:
            if ch not in node.children:
                return None
            node = node.children[ch]
        return node

    def search(self, word: str) -> bool:
        """完整单词是否存在，时间 O(|word|)。"""
        node = self._walk(word)
        return node is not None and node.is_end

    def starts_with(self, prefix: str) -> bool:
        """是否存在该前缀，时间 O(|prefix|)。"""
        return self._walk(prefix) is not None

    def count_prefix(self, prefix: str) -> int:
        """以 prefix 为前缀的单词个数（重复插入计多次）。"""
        node = self._walk(prefix)
        return 0 if node is None else node.pass_count

    def count_word(self, word: str) -> int:
        """word 被插入了几次。"""
        node = self._walk(word)
        return 0 if node is None else node.end_count


def run_io(data: str) -> None:
    """按统一输入输出格式处理 stdin 数据。"""
    lines = [ln for ln in data.splitlines()]
    n = int(lines[0].split()[0])
    trie = Trie()
    out: List[str] = []
    for line in lines[1:1 + n]:
        parts = line.split()
        if len(parts) < 2:
            continue
        op, arg = parts[0], parts[1]
        if op == "insert":
            trie.insert(arg)
        elif op == "search":
            out.append("true" if trie.search(arg) else "false")
        elif op == "startsWith":
            out.append("true" if trie.starts_with(arg) else "false")
        elif op == "countPrefix":
            out.append(str(trie.count_prefix(arg)))
    print("\n".join(out))


def run_tests() -> None:
    trie = Trie()
    trie.insert("app")
    trie.insert("apple")
    trie.insert("application")
    trie.insert("banana")

    assert trie.search("app") is True                 # 完整单词，确实插入过
    assert trie.search("appl") is False               # 只是前缀，不是完整单词
    assert trie.search("apple") is True
    assert trie.search("applex") is False
    assert trie.starts_with("app") is True
    assert trie.starts_with("ban") is True
    assert trie.starts_with("cat") is False
    assert trie.count_prefix("app") == 3              # app / apple / application
    assert trie.count_prefix("appl") == 2             # apple / application
    assert trie.count_prefix("apple") == 1            # 只有 apple（application 是 appli...）
    assert trie.count_prefix("b") == 1
    assert trie.count_prefix("z") == 0

    trie.insert("app")                                 # 允许重复插入
    assert trie.count_word("app") == 2
    assert trie.count_prefix("app") == 4

    # 空串：根节点既是起点也是终点
    t2 = Trie()
    t2.insert("")
    assert t2.search("") is True
    assert t2.starts_with("") is True
    assert t2.count_prefix("") == 1

    print("all tests passed")


if __name__ == "__main__":
    raw = sys.stdin.read() if not sys.stdin.isatty() else ""
    if raw.strip():
        run_io(raw)
    else:
        run_tests()
