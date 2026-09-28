// KMP 字符串匹配（Knuth-Morris-Pratt）
// 編譯：g++ -std=c++17 -O2 solution.cpp -o solution && ./solution
#include <iostream>
#include <string>
#include <vector>
#include <cassert>

using namespace std;

// 構造前綴函數 pi：pi[i] = pattern[0..i] 的最長真前綴同時也是後綴的長度
vector<int> buildPrefix(const string& pattern) {
    vector<int> pi(pattern.size(), 0);
    int j = 0;  // 當前已匹配的前綴長度
    for (int i = 1; i < (int)pattern.size(); ++i) {
        // 失配時回退到更短的前綴，直到能接上或回到 0
        while (j > 0 && pattern[i] != pattern[j]) j = pi[j - 1];
        if (pattern[i] == pattern[j]) ++j;
        pi[i] = j;
    }
    return pi;
}

// 返回 pattern 在 text 中首次出現的下標，不存在返回 -1；空 pattern 返回 0
int kmpSearch(const string& text, const string& pattern) {
    if (pattern.empty()) return 0;
    vector<int> pi = buildPrefix(pattern);
    int j = 0;  // 當前在 pattern 上匹配到的長度
    for (int i = 0; i < (int)text.size(); ++i) {
        // 關鍵：text 指針永不回退，只回退 pattern 指針
        while (j > 0 && text[i] != pattern[j]) j = pi[j - 1];
        if (text[i] == pattern[j]) ++j;
        if (j == (int)pattern.size()) return i - (int)pattern.size() + 1;
    }
    return -1;
}

// 返回所有出現位置（允許重疊）
vector<int> kmpSearchAll(const string& text, const string& pattern) {
    vector<int> hits;
    if (pattern.empty()) {
        for (int i = 0; i <= (int)text.size(); ++i) hits.push_back(i);
        return hits;
    }
    vector<int> pi = buildPrefix(pattern);
    int j = 0;
    for (int i = 0; i < (int)text.size(); ++i) {
        while (j > 0 && text[i] != pattern[j]) j = pi[j - 1];
        if (text[i] == pattern[j]) ++j;
        if (j == (int)pattern.size()) {
            hits.push_back(i - (int)pattern.size() + 1);
            j = pi[j - 1];  // 繼續尋找下一次匹配
        }
    }
    return hits;
}

int main() {
    assert(kmpSearch("ababcabcabababd", "ababd") == 10);
    assert(kmpSearch("hello", "ll") == 2);
    assert(kmpSearch("aaaaa", "bba") == -1);
    assert(kmpSearch("abc", "") == 0);

    vector<int> all = kmpSearchAll("ababab", "aba");
    assert(all.size() == 2 && all[0] == 0 && all[1] == 2);

    vector<int> pi = buildPrefix("ababaca");
    vector<int> expect = {0, 0, 1, 2, 3, 0, 1};
    assert(pi == expect);

    cout << "all tests passed" << endl;
    return 0;
}
