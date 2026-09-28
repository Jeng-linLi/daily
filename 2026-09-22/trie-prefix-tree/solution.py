"""前綴樹（Trie / Prefix Tree）

題意：維護一個字符串集合，支持三種操作：
    insert(word)     插入一個單詞
    search(word)     查詢單詞是否被完整插入過
    starts_with(p)   查詢是否存在以 p 爲前綴的單詞
    count_prefix(p)  統計已插入單詞中以 p 爲前綴的個數

思路：
    Trie 是一棵「按字符分叉」的多叉樹，從根到某個節點的路徑拼成一個前綴。
    每個節點維護：
      children：字符 -> 子節點（用字典，天然支持任意字符集）
      is_end：  是否有單詞在此節點結束（區分 "app" 與 "apple"）
      pass：    有多少個已插入單詞經過該節點（即該前綴出現次數）
      end：     有多少個單詞恰好在此節點結束（支持重複插入計數）
    插入時沿字符逐層走，不存在就新建節點，沿途 pass 都 +1；結束時置 is_end。
    查找複雜度只與單詞長度有關，與集合規模無關，這是它相對哈希表的最大優勢。

輸入格式（stdin）：
    第一行 n（操作條數）
    接下來 n 行，每行：insert <word> | search <word> | startsWith <prefix> | countPrefix <prefix>
輸出格式（stdout）：
    search / startsWith 各輸出一行 true / false；countPrefix 輸出一行整數；insert 無輸出
無 stdin 輸入時運行內置斷言測試。
"""

import sys
from typing import Dict, List


class TrieNode:
    __slots__ = ("children", "is_end", "pass_count", "end_count")

    def __init__(self) -> None:
        self.children: Dict[str, "TrieNode"] = {}  # 字符 -> 子節點
        self.is_end: bool = False                  # 是否有單詞在此結束
        self.pass_count: int = 0                   # 經過該節點的單詞數
        self.end_count: int = 0                    # 在此結束的單詞數（含重複插入）


class Trie:
    def __init__(self) -> None:
        self.root = TrieNode()

    def insert(self, word: str) -> None:
        """插入單詞，時間 O(|word|)。"""
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
        """沿字符串走到對應節點；中途斷掉返回 None。"""
        node = self.root
        for ch in s:
            if ch not in node.children:
                return None
            node = node.children[ch]
        return node

    def search(self, word: str) -> bool:
        """完整單詞是否存在，時間 O(|word|)。"""
        node = self._walk(word)
        return node is not None and node.is_end

    def starts_with(self, prefix: str) -> bool:
        """是否存在該前綴，時間 O(|prefix|)。"""
        return self._walk(prefix) is not None

    def count_prefix(self, prefix: str) -> int:
        """以 prefix 爲前綴的單詞個數（重複插入計多次）。"""
        node = self._walk(prefix)
        return 0 if node is None else node.pass_count

    def count_word(self, word: str) -> int:
        """word 被插入了幾次。"""
        node = self._walk(word)
        return 0 if node is None else node.end_count


def run_io(data: str) -> None:
    """按統一輸入輸出格式處理 stdin 數據。"""
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

    assert trie.search("app") is True                 # 完整單詞，確實插入過
    assert trie.search("appl") is False               # 只是前綴，不是完整單詞
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

    trie.insert("app")                                 # 允許重複插入
    assert trie.count_word("app") == 2
    assert trie.count_prefix("app") == 4

    # 空串：根節點既是起點也是終點
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
