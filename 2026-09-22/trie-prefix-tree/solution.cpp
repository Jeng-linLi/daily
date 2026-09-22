// 前缀树（Trie / Prefix Tree）
// 编译：g++ -std=c++17 -O2 solution.cpp -o solution && ./solution
//
// 思路：Trie 是「按字符分叉」的多叉树，根到某节点的路径拼成一个前缀。
//   每个节点维护 children（字符 -> 子节点）、isEnd（是否有单词在此结束，
//   用来区分 "app" 与 "apple"）、pass（经过该节点的单词数）、endCnt（在此结束的次数）。
//   插入时逐字符走，缺节点就新建，沿途 pass 全部 +1，结束时置 isEnd。
//   复杂度只与字符串长度有关，与集合规模无关——这是它相对哈希表最大的优势。
//
// 输入：第一行 n；接下来 n 行：insert w | search w | startsWith p | countPrefix p
// 输出：search / startsWith 输出 true / false；countPrefix 输出整数；insert 无输出
// 无 stdin 输入时运行内置断言测试。
#include <cassert>
#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>

using namespace std;

struct TrieNode {
    unordered_map<char, int> children;  // 字符 -> 子节点在 pool 中的下标
    bool isEnd = false;                 // 是否有单词在此结束
    int passCount = 0;                  // 经过该节点的单词数
    int endCount = 0;                   // 在此结束的单词数（含重复插入）
};

class Trie {
public:
    Trie() { pool.emplace_back(); }  // pool[0] 为根

    // 插入单词，时间 O(|word|)
    void insert(const string& word) {
        int u = 0;
        pool[u].passCount++;
        for (char ch : word) {
            auto it = pool[u].children.find(ch);
            if (it == pool[u].children.end()) {
                int v = static_cast<int>(pool.size());
                pool[u].children[ch] = v;
                pool.emplace_back();
                u = v;
            } else {
                u = it->second;
            }
            pool[u].passCount++;
        }
        pool[u].isEnd = true;
        pool[u].endCount++;
    }

    // 完整单词是否存在，时间 O(|word|)
    bool search(const string& word) const {
        int u = walk(word);
        return u != -1 && pool[u].isEnd;
    }

    // 是否存在该前缀，时间 O(|prefix|)
    bool startsWith(const string& prefix) const { return walk(prefix) != -1; }

    // 以 prefix 为前缀的单词个数（重复插入计多次）
    int countPrefix(const string& prefix) const {
        int u = walk(prefix);
        return u == -1 ? 0 : pool[u].passCount;
    }

    // word 被插入了几次
    int countWord(const string& word) const {
        int u = walk(word);
        return u == -1 ? 0 : pool[u].endCount;
    }

private:
    vector<TrieNode> pool;

    // 沿字符串走到对应节点；中途断掉返回 -1
    int walk(const string& s) const {
        int u = 0;
        for (char ch : s) {
            auto it = pool[u].children.find(ch);
            if (it == pool[u].children.end()) return -1;
            u = it->second;
        }
        return u;
    }
};

int main() {
    int n;
    if (cin >> n) {  // IO 模式
        Trie trie;
        string op, arg;
        for (int i = 0; i < n; ++i) {
            cin >> op >> arg;
            if (op == "insert") {
                trie.insert(arg);
            } else if (op == "search") {
                cout << (trie.search(arg) ? "true" : "false") << "\n";
            } else if (op == "startsWith") {
                cout << (trie.startsWith(arg) ? "true" : "false") << "\n";
            } else if (op == "countPrefix") {
                cout << trie.countPrefix(arg) << "\n";
            }
        }
        return 0;
    }

    Trie trie;
    trie.insert("app");
    trie.insert("apple");
    trie.insert("application");
    trie.insert("banana");

    assert(trie.search("app") == true);        // 完整单词，确实插入过
    assert(trie.search("appl") == false);      // 只是前缀，不是完整单词
    assert(trie.search("apple") == true);
    assert(trie.search("applex") == false);
    assert(trie.startsWith("app") == true);
    assert(trie.startsWith("ban") == true);
    assert(trie.startsWith("cat") == false);
    assert(trie.countPrefix("app") == 3);      // app / apple / application
    assert(trie.countPrefix("appl") == 2);     // apple / application
    assert(trie.countPrefix("apple") == 1);    // 只有 apple（application 是 appli...）
    assert(trie.countPrefix("b") == 1);
    assert(trie.countPrefix("z") == 0);

    trie.insert("app");                        // 允许重复插入
    assert(trie.countWord("app") == 2);
    assert(trie.countPrefix("app") == 4);

    // 空串：根节点既是起点也是终点
    Trie t2;
    t2.insert("");
    assert(t2.search("") == true);
    assert(t2.startsWith("") == true);
    assert(t2.countPrefix("") == 1);

    cout << "all tests passed" << endl;
    return 0;
}
