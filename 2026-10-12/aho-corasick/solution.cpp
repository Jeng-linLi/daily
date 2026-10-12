// Aho-Corasick 多模式串匹配（Trie + fail 指針 + 完整轉移表）
//
// 與 solution.py 完全同構：同樣的演算法、同樣的確定性規則、同樣的輸入輸出格式。
// 編譯：g++ -std=c++17 -O2 -Wall solution.cpp -o solution
// 無 stdin 輸入時運行內置斷言測試並輸出 `all tests passed`。

#include <algorithm>
#include <cassert>
#include <cctype>
#include <deque>
#include <iostream>
#include <map>
#include <random>
#include <set>
#include <sstream>
#include <string>
#include <vector>

using namespace std;

// ---------------------------------------------------------------- 建自動機

struct AC {
    vector<map<char, int>> nxt;   // Trie 的轉移邊
    vector<int> fail;             // 失配指針
    vector<vector<int>> out;      // 合併後的輸出表（在此結尾的模式串編號）
    vector<int> order;            // BFS 序，建轉移表要用
};

// 建立 Trie + fail 指針 + 合併後的輸出表。
static AC buildTrie(const vector<string> &patterns) {
    AC ac;
    ac.nxt.push_back(map<char, int>());
    ac.fail.push_back(0);
    ac.out.push_back(vector<int>());

    for (int pid = 0; pid < (int)patterns.size(); pid++) {
        const string &p = patterns[pid];
        if (p.empty()) continue;              // 空模式串不進 Trie
        int v = 0;
        for (char ch : p) {
            auto it = ac.nxt[v].find(ch);
            if (it != ac.nxt[v].end()) {
                v = it->second;
            } else {
                int id = (int)ac.nxt.size();
                ac.nxt.push_back(map<char, int>());
                ac.fail.push_back(0);
                ac.out.push_back(vector<int>());
                ac.nxt[v][ch] = id;
                v = id;
            }
        }
        ac.out[v].push_back(pid);
    }

    ac.order.push_back(0);
    deque<int> q;
    for (const auto &kv : ac.nxt[0]) {        // map 已按鍵排序，與 Python 的 sorted() 一致
        int u = kv.second;
        ac.fail[u] = 0;
        q.push_back(u);
        ac.order.push_back(u);
    }

    while (!q.empty()) {
        int v = q.front();
        q.pop_front();
        // fail[v] 的深度更小，必定已先出隊並合併完畢 → 直接拼接
        vector<int> merged = ac.out[v];
        const vector<int> &fout = ac.out[ac.fail[v]];
        merged.insert(merged.end(), fout.begin(), fout.end());
        ac.out[v] = merged;

        for (const auto &kv : ac.nxt[v]) {
            char ch = kv.first;
            int u = kv.second;
            int f = ac.fail[v];
            while (f != 0 && ac.nxt[f].find(ch) == ac.nxt[f].end()) f = ac.fail[f];
            auto it = ac.nxt[f].find(ch);
            ac.fail[u] = (it != ac.nxt[f].end()) ? it->second : 0;
            q.push_back(u);
            ac.order.push_back(u);
        }
    }
    return ac;
}

// 建完整轉移表（Trie 圖）：go[v][c] = 從 v 吃掉字元 c 後到達的狀態。
static vector<vector<int>> buildGoto(const AC &ac, const vector<char> &alphabet,
                                     map<char, int> &idx) {
    for (int i = 0; i < (int)alphabet.size(); i++) idx[alphabet[i]] = i;
    int A = (int)alphabet.size();
    int S = (int)ac.nxt.size();
    vector<vector<int>> go(S, vector<int>(A, 0));
    for (int ci = 0; ci < A; ci++) {
        auto it = ac.nxt[0].find(alphabet[ci]);
        go[0][ci] = (it != ac.nxt[0].end()) ? it->second : 0;
    }
    for (int v : ac.order) {
        if (v == 0) continue;
        for (int ci = 0; ci < A; ci++) {
            char ch = alphabet[ci];
            auto it = ac.nxt[v].find(ch);
            if (it != ac.nxt[v].end())
                go[v][ci] = it->second;
            else
                go[v][ci] = go[ac.fail[v]][ci];
        }
    }
    return go;
}

// 沿 fail 鏈爬求轉移，O(深度)，僅供對拍驗證轉移表。
static int gotoNaive(const AC &ac, int v, char ch) {
    while (v != 0 && ac.nxt[v].find(ch) == ac.nxt[v].end()) v = ac.fail[v];
    auto it = ac.nxt[v].find(ch);
    return (it != ac.nxt[v].end()) ? it->second : 0;
}

// ---------------------------------------------------------------- 匹配

// 一次掃描找出所有模式串的所有出現位置（允許重疊）。
static vector<vector<int>> acMatch(const string &text, const vector<string> &patterns,
                                   const vector<vector<int>> &go, const map<char, int> &idx,
                                   const AC &ac) {
    vector<vector<int>> pos(patterns.size());
    int v = 0;
    for (int i = 0; i < (int)text.size(); i++) {
        auto it = idx.find(text[i]);
        if (it == idx.end()) {                // 字母表外的字元 → 回到根
            v = 0;
            continue;
        }
        v = go[v][it->second];
        for (int pid : ac.out[v])
            pos[pid].push_back(i - (int)patterns[pid].size() + 1);
    }
    return pos;
}

// 暴力匹配，O(n·|pat|)，僅供對拍。空模式串視為永不命中（與自動機一致）。
static vector<int> bfMatch(const string &text, const string &pat) {
    vector<int> res;
    if (pat.empty()) return res;
    int L = (int)pat.size();
    for (int i = 0; i + L <= (int)text.size(); i++)
        if (text.compare(i, L, pat) == 0) res.push_back(i);
    return res;
}

// ---------------------------------------------------------------- 輸出

static string run(const string &text, const vector<string> &patterns) {
    AC ac = buildTrie(patterns);
    set<char> cs;
    for (char c : text) cs.insert(c);
    for (const string &p : patterns)
        for (char c : p) cs.insert(c);
    vector<char> alphabet(cs.begin(), cs.end());

    map<char, int> idx;
    vector<vector<int>> go = buildGoto(ac, alphabet, idx);
    vector<vector<int>> pos = acMatch(text, patterns, go, idx, ac);

    int total = 0;
    for (const auto &p : pos) total += (int)p.size();

    // 與 Python 的 "\n".join(lines) 完全一致：行與行之間放 "\n"，最後一行不加
    ostringstream oss;
    oss << ac.nxt.size() << "\n" << total << "\n";
    for (int i = 0; i < (int)pos.size(); i++) {
        if (i) oss << " ";
        oss << pos[i].size();
    }
    for (int i = 0; i < (int)pos.size(); i++) {
        oss << "\n";
        for (int j = 0; j < (int)pos[i].size(); j++) {
            if (j) oss << " ";
            oss << pos[i][j];
        }
    }
    return oss.str();
}

// ---------------------------------------------------------------- 測試

static void check(const string &text, const vector<string> &patterns) {
    AC ac = buildTrie(patterns);
    set<char> cs;
    for (char c : text) cs.insert(c);
    for (const string &p : patterns)
        for (char c : p) cs.insert(c);
    vector<char> alphabet(cs.begin(), cs.end());

    map<char, int> idx;
    vector<vector<int>> go = buildGoto(ac, alphabet, idx);

    // 轉移表與 fail 鏈必須完全一致
    for (int v = 0; v < (int)ac.nxt.size(); v++)
        for (char c : alphabet) assert(go[v][idx[c]] == gotoNaive(ac, v, c));

    // fail 指標的基本性質
    assert(ac.fail[0] == 0);
    for (int v = 1; v < (int)ac.nxt.size(); v++) assert(ac.fail[v] != v);

    vector<vector<int>> pos = acMatch(text, patterns, go, idx, ac);
    assert(pos.size() == patterns.size());
    for (int pid = 0; pid < (int)patterns.size(); pid++)
        assert(pos[pid] == bfMatch(text, patterns[pid]));
}

static string randStr(mt19937 &rng, int n, const string &alpha) {
    string s;
    for (int i = 0; i < n; i++) s += alpha[rng() % alpha.size()];
    return s;
}

static void runTests() {
    mt19937 rng(20261012);

    // ---- 固定用例 ----
    check("ababa", {"aba", "ba", "a"});
    check("aaaa", {"aa", "aaa"});
    check("abc", {});
    check("", {"a", "ab"});
    check("xyz", {"q"});
    check("aaa", {""});
    check("she", {"he", "she", "hers", "his"});

    assert(run("ababa", {"aba", "ba", "a"}) == "6\n7\n2 2 3\n0 2\n1 3\n0 2 4");
    assert(run("aaaa", {"aa", "aaa"}) == "4\n5\n3 2\n0 1 2\n0 1");
    assert(run("abc", {}) == "1\n0\n");
    assert(run("", {"a", "ab"}) == "3\n0\n0 0\n\n");

    // ---- 隨機對拍：二字元字母表 ----
    for (int iter = 0; iter < 400; iter++) {
        int n = (int)(rng() % 31);
        string text = randStr(rng, n, "ab");
        int k = (int)(rng() % 6);
        vector<string> pats;
        for (int i = 0; i < k; i++) pats.push_back(randStr(rng, (int)(rng() % 5), "ab"));
        check(text, pats);
    }

    // ---- 隨機對拍：三字元字母表 + 較長模式串 ----
    for (int iter = 0; iter < 300; iter++) {
        int n = (int)(rng() % 41);
        string text = randStr(rng, n, "abc");
        int k = 1 + (int)(rng() % 6);
        vector<string> pats;
        for (int i = 0; i < k; i++) pats.push_back(randStr(rng, 1 + (int)(rng() % 6), "abc"));
        check(text, pats);
    }

    // ---- 隨機對拍：模式串取自文字的真實子串 ----
    for (int iter = 0; iter < 300; iter++) {
        int n = 1 + (int)(rng() % 25);
        string text = randStr(rng, n, "abcd");
        int k = 1 + (int)(rng() % 5);
        vector<string> pats;
        for (int i = 0; i < k; i++) {
            int a = (int)(rng() % n);
            int b = min(n, a + 1 + (int)(rng() % 5));
            pats.push_back(text.substr(a, b - a));
        }
        check(text, pats);
    }
}

// stdin 解析：嚴格只吃 `^[+-]?[0-9]+$` 的數字 token，其餘皆當作字串。
static bool isInt(const string &s) {
    if (s.empty()) return false;
    int i = 0;
    if (s[0] == '+' || s[0] == '-') i = 1;
    if (i >= (int)s.size()) return false;
    for (int j = i; j < (int)s.size(); j++)
        if (!isdigit((unsigned char)s[j])) return false;
    return true;
}

static void runIO(const vector<string> &tokens) {
    string text = tokens.size() >= 1 ? tokens[0] : "";
    int k = 0;
    if (tokens.size() >= 2 && isInt(tokens[1])) {
        k = stoi(tokens[1]);
        if (k < 0) k = 0;
    }
    vector<string> pats;
    for (int i = 0; i < k; i++) pats.push_back(tokens.size() >= 3 + (size_t)i ? tokens[2 + i] : "");
    cout << run(text, pats) << "\n";
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    vector<string> tokens;
    bool any = false;
    string cur;
    int c;
    while ((c = cin.get()) != EOF) {
        if (isspace(c)) {
            if (!cur.empty()) { tokens.push_back(cur); cur.clear(); }
        } else {
            cur.push_back((char)c);
            any = true;
        }
    }
    if (!cur.empty()) tokens.push_back(cur);

    if (any) {
        runIO(tokens);
    } else {
        runTests();
        cout << "all tests passed" << endl;
    }
    return 0;
}
