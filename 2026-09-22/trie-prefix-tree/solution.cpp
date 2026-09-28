// 前綴樹（Trie / Prefix Tree）
// 編譯：g++ -std=c++17 -O2 solution.cpp -o solution && ./solution
//
// 思路：Trie 是「按字符分叉」的多叉樹，根到某節點的路徑拼成一個前綴。
//   每個節點維護 children（字符 -> 子節點）、isEnd（是否有單詞在此結束，
//   用來區分 "app" 與 "apple"）、pass（經過該節點的單詞數）、endCnt（在此結束的次數）。
//   插入時逐字符走，缺節點就新建，沿途 pass 全部 +1，結束時置 isEnd。
//   複雜度只與字符串長度有關，與集合規模無關——這是它相對哈希表最大的優勢。
//
// 輸入：第一行 n；接下來 n 行：insert w | search w | startsWith p | countPrefix p
// 輸出：search / startsWith 輸出 true / false；countPrefix 輸出整數；insert 無輸出
// 無 stdin 輸入時運行內置斷言測試。
#include <cassert>
#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>

using namespace std;

struct TrieNode {
    unordered_map<char, int> children;  // 字符 -> 子節點在 pool 中的下標
    bool isEnd = false;                 // 是否有單詞在此結束
    int passCount = 0;                  // 經過該節點的單詞數
    int endCount = 0;                   // 在此結束的單詞數（含重複插入）
};

class Trie {
public:
    Trie() { pool.emplace_back(); }  // pool[0] 爲根

    // 插入單詞，時間 O(|word|)
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

    // 完整單詞是否存在，時間 O(|word|)
    bool search(const string& word) const {
        int u = walk(word);
        return u != -1 && pool[u].isEnd;
    }

    // 是否存在該前綴，時間 O(|prefix|)
    bool startsWith(const string& prefix) const { return walk(prefix) != -1; }

    // 以 prefix 爲前綴的單詞個數（重複插入計多次）
    int countPrefix(const string& prefix) const {
        int u = walk(prefix);
        return u == -1 ? 0 : pool[u].passCount;
    }

    // word 被插入了幾次
    int countWord(const string& word) const {
        int u = walk(word);
        return u == -1 ? 0 : pool[u].endCount;
    }

private:
    vector<TrieNode> pool;

    // 沿字符串走到對應節點；中途斷掉返回 -1
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

    assert(trie.search("app") == true);        // 完整單詞，確實插入過
    assert(trie.search("appl") == false);      // 只是前綴，不是完整單詞
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

    trie.insert("app");                        // 允許重複插入
    assert(trie.countWord("app") == 2);
    assert(trie.countPrefix("app") == 4);

    // 空串：根節點既是起點也是終點
    Trie t2;
    t2.insert("");
    assert(t2.search("") == true);
    assert(t2.startsWith("") == true);
    assert(t2.countPrefix("") == 1);

    cout << "all tests passed" << endl;
    return 0;
}
