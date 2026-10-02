// Manacher 算法：最長回文子串 / 回文子串計數（O(n)）
//
// 題意：
//     給定字串 s，求：1. 最長回文子串長度；2. 最長回文子串（同長度取最靠左）；
//     3. 回文子串總數；4. 奇回文半徑陣列 d1；5. 偶回文半徑陣列 d2。全部 O(n)。
//
// 思路：
//     中心擴展最直觀但最壞 O(n^2)（"aaaa…a"）。Manacher 的洞見是「已算過的回文可以複用」：
//     維護當前右端最遠的回文區間 [l, r]，處理中心 i 時，若 i 在區間內，
//     取其鏡像位置 j 的半徑作為初始值（並用 r − i + 1 截斷，超出右界的部分不保證），
//     然後才繼續往外擴。因為 [l, r] 的右端只會單調右移，while 的總執行次數是 O(n)。
//       - d1[i]：以 i 為中心的奇回文半徑（含中心），最長長度 = 2*d1[i] − 1；
//                同時 d1[i] 正好等於「以 i 為中心的奇回文個數」。
//       - d2[i]：以 i−1 與 i 之間為中心的偶回文半徑，最長長度 = 2*d2[i]；同理即個數。
//     因此回文子串總數 = sum(d1) + sum(d2)。
//     最長回文取最靠左，所以更新條件必須是嚴格 >（先掃到的先佔住）。
//
// 輸入格式（stdin）：
//     s          （第一行；空字串請給一個只有換行的輸入）
// 輸出格式（stdout）：
//     第 1 行：最長回文子串長度
//     第 2 行：最長回文子串（長度為 0 時輸出空行）
//     第 3 行：回文子串總數
//     第 4 行：d1 陣列，空格分隔（n = 0 時輸出空行）
//     第 5 行：d2 陣列，空格分隔（n = 0 時輸出空行）
// 無 stdin 輸入時運行內置斷言測試並輸出 `all tests passed`。

#include <algorithm>
#include <cassert>
#include <iostream>
#include <random>
#include <string>
#include <vector>

using namespace std;

// ---------------------------------------------------------------- Manacher
// d1[i] = 以 i 為中心的奇回文半徑（含中心）；最長奇回文長度 = 2*d1[i] − 1
static vector<int> manacherOdd(const string& s) {
    int n = (int)s.size();
    vector<int> d1(n, 0);
    int l = 0, r = -1;                       // 目前右端最遠的回文區間 [l, r]
    for (int i = 0; i < n; ++i) {
        int k = (i > r) ? 1 : min(d1[l + r - i], r - i + 1);
        while (i - k >= 0 && i + k < n && s[i - k] == s[i + k]) k++;
        d1[i] = k;
        if (i + k - 1 > r) {                 // 右端推進，更新當前最右回文
            l = i - k + 1;
            r = i + k - 1;
        }
    }
    return d1;
}

// d2[i] = 以 i−1 與 i 之間為中心的偶回文半徑；最長偶回文長度 = 2*d2[i]
static vector<int> manacherEven(const string& s) {
    int n = (int)s.size();
    vector<int> d2(n, 0);
    int l = 0, r = -1;
    for (int i = 0; i < n; ++i) {
        int k = (i > r) ? 0 : min(d2[l + r - i + 1], r - i + 1);
        while (i - k - 1 >= 0 && i + k < n && s[i - k - 1] == s[i + k]) k++;
        d2[i] = k;
        if (i + k - 1 > r) {
            l = i - k;
            r = i + k - 1;
        }
    }
    return d2;
}

// ---------------------------------------------------------------- 中心擴展（暴力基準）
static vector<int> centerExpansionOdd(const string& s) {
    int n = (int)s.size();
    vector<int> d1(n, 0);
    for (int i = 0; i < n; ++i) {
        int k = 1;
        while (i - k >= 0 && i + k < n && s[i - k] == s[i + k]) k++;
        d1[i] = k;
    }
    return d1;
}

static vector<int> centerExpansionEven(const string& s) {
    int n = (int)s.size();
    vector<int> d2(n, 0);
    for (int i = 0; i < n; ++i) {
        int k = 0;
        while (i - k - 1 >= 0 && i + k < n && s[i - k - 1] == s[i + k]) k++;
        d2[i] = k;
    }
    return d2;
}

// ---------------------------------------------------------------- 由半徑推出答案
// 最長回文子串的 (長度, 內容)；同長度取最靠左（嚴格 > 更新）
static pair<int, string> longestPalindrome(const string& s, const vector<int>& d1,
                                           const vector<int>& d2) {
    int n = (int)s.size();
    int best = 0, start = 0;
    for (int i = 0; i < n; ++i) {
        if (2 * d1[i] - 1 > best) {
            best = 2 * d1[i] - 1;
            start = i - d1[i] + 1;
        }
        if (2 * d2[i] > best) {
            best = 2 * d2[i];
            start = i - d2[i];
        }
    }
    return {best, s.substr((size_t)start, (size_t)best)};
}

// 回文子串總數 = sum(d1) + sum(d2)
static long long countPalindromes(const vector<int>& d1, const vector<int>& d2) {
    long long ans = 0;
    for (int x : d1) ans += x;
    for (int x : d2) ans += x;
    return ans;
}

// O(n^3) 枚舉所有子串求最長回文，只在測試裡當基準
static pair<int, string> longestPalindromeNaive(const string& s) {
    int n = (int)s.size();
    int best = 0;
    string ans;
    for (int i = 0; i < n; ++i)
        for (int j = i; j < n; ++j) {
            int len = j - i + 1;
            if (len <= best) continue;
            bool ok = true;
            for (int t = 0; t < len / 2; ++t)
                if (s[i + t] != s[j - t]) { ok = false; break; }
            if (ok) {
                best = len;
                ans = s.substr((size_t)i, (size_t)len);
            }
        }
    return {best, ans};
}

// O(n^3) 統計回文子串個數，只在測試裡當基準
static long long countPalindromesNaive(const string& s) {
    int n = (int)s.size();
    long long cnt = 0;
    for (int i = 0; i < n; ++i)
        for (int j = i; j < n; ++j) {
            bool ok = true;
            for (int t = 0; t < (j - i + 1) / 2; ++t)
                if (s[i + t] != s[j - t]) { ok = false; break; }
            if (ok) cnt++;
        }
    return cnt;
}

// ---------------------------------------------------------------- 小工具
static void printVec(const vector<int>& v) {
    for (size_t i = 0; i < v.size(); ++i) {
        if (i) cout << ' ';
        cout << v[i];
    }
    cout << '\n';
}

static mt19937 rngEngine;
static int rndInt(int lo, int hi) { return lo + (int)(rngEngine() % (unsigned)(hi - lo + 1)); }

static bool isPal(const string& s, int l, int r) {
    while (l < r)
        if (s[l++] != s[r--]) return false;
    return true;
}

// ---------------------------------------------------------------- IO 模式
static void runIo(const string& data) {
    string s;
    for (char c : data) {                    // 只取第一行
        if (c == '\n') break;
        s += c;
    }
    if (!s.empty() && s.back() == '\r') s.pop_back();

    vector<int> d1 = manacherOdd(s);
    vector<int> d2 = manacherEven(s);
    pair<int, string> lp = longestPalindrome(s, d1, d2);
    cout << lp.first << '\n';
    cout << lp.second << '\n';
    cout << countPalindromes(d1, d2) << '\n';
    printVec(d1);
    printVec(d2);
}

// ---------------------------------------------------------------- 測試
static void runTests() {
    // 空字串
    assert(manacherOdd("").empty());
    assert(manacherEven("").empty());
    assert(longestPalindrome("", {}, {}) == make_pair(0, string("")));
    assert(countPalindromes({}, {}) == 0);

    // 單字元
    assert((manacherOdd("a") == vector<int>{1}));
    assert((manacherEven("a") == vector<int>{0}));
    assert(longestPalindrome("a", {1}, {0}) == make_pair(1, string("a")));
    assert(countPalindromes({1}, {0}) == 1);

    // README 示例："cbbd" —— 最長回文是偶回文 "bb"
    {
        string s = "cbbd";
        vector<int> d1 = manacherOdd(s), d2 = manacherEven(s);
        assert((d1 == vector<int>{1, 1, 1, 1}));
        assert((d2 == vector<int>{0, 0, 1, 0}));
        assert(d1 == centerExpansionOdd(s));
        assert(d2 == centerExpansionEven(s));
        assert(longestPalindrome(s, d1, d2) == make_pair(2, string("bb")));
        assert(countPalindromes(d1, d2) == 5);
        assert(countPalindromes(d1, d2) == countPalindromesNaive(s));
    }

    // 同長度取最左："abacdc" 的 "aba" 與 "cdc" 都是 3，取 "aba"
    {
        string s = "abacdc";
        vector<int> d1 = manacherOdd(s), d2 = manacherEven(s);
        assert((d1 == vector<int>{1, 2, 1, 1, 2, 1}));
        assert((d2 == vector<int>{0, 0, 0, 0, 0, 0}));
        assert(d1 == centerExpansionOdd(s));
        assert(d2 == centerExpansionEven(s));
        assert(longestPalindrome(s, d1, d2) == make_pair(3, string("aba")));
        assert(countPalindromes(d1, d2) == 8);
        assert(countPalindromes(d1, d2) == countPalindromesNaive(s));
    }

    // 全同字元："aaaa" 的回文子串數 = 4*5/2 = 10
    {
        string s = "aaaa";
        vector<int> d1 = manacherOdd(s), d2 = manacherEven(s);
        assert((d1 == vector<int>{1, 2, 2, 1}));
        assert((d2 == vector<int>{0, 1, 2, 1}));       // d2[3] 只到 1（右邊沒有第 5 個 a）
        assert(d1 == centerExpansionOdd(s));
        assert(d2 == centerExpansionEven(s));
        assert(longestPalindrome(s, d1, d2) == make_pair(4, string("aaaa")));
        assert(countPalindromes(d1, d2) == 10);
        assert(countPalindromes(d1, d2) == countPalindromesNaive(s));
    }

    // 完全無回文（長度 > 1）："abcde"
    {
        string s = "abcde";
        vector<int> d1 = manacherOdd(s), d2 = manacherEven(s);
        assert(longestPalindrome(s, d1, d2) == make_pair(1, string("a")));
        assert(countPalindromes(d1, d2) == 5);
    }

    // 偶回文最長："abba"
    {
        string s = "abba";
        vector<int> d1 = manacherOdd(s), d2 = manacherEven(s);
        assert(longestPalindrome(s, d1, d2) == make_pair(4, string("abba")));
        assert(countPalindromes(d1, d2) == countPalindromesNaive(s));
        assert(countPalindromes(d1, d2) == 6);
    }

    // 隨機對拍：Manacher vs 中心擴展 vs 暴力枚舉
    rngEngine.seed(20261002);
    for (int round = 0; round < 2; ++round) {
        string alphabet = (round == 0) ? "ab" : "abc";
        int rounds = (round == 0) ? 400 : 200;
        int maxLen = (round == 0) ? 12 : 14;
        for (int t = 0; t < rounds; ++t) {
            int n = rndInt(0, maxLen);
            string s;
            for (int i = 0; i < n; ++i) s += alphabet[(size_t)rndInt(0, (int)alphabet.size() - 1)];

            vector<int> d1 = manacherOdd(s), d2 = manacherEven(s);
            assert(d1 == centerExpansionOdd(s));
            assert(d2 == centerExpansionEven(s));

            // 每個半徑都要真的對應一個回文，且再往外一格就不回文
            for (int i = 0; i < n; ++i) {
                assert(d1[i] >= 1);
                assert(isPal(s, i - d1[i] + 1, i + d1[i] - 1));
                if (i - d1[i] >= 0 && i + d1[i] < n) assert(s[i - d1[i]] != s[i + d1[i]]);
                if (d2[i] > 0) assert(isPal(s, i - d2[i], i + d2[i] - 1));
                if (i - d2[i] - 1 >= 0 && i + d2[i] < n) assert(s[i - d2[i] - 1] != s[i + d2[i]]);
            }

            assert(longestPalindrome(s, d1, d2) == longestPalindromeNaive(s));
            assert(countPalindromes(d1, d2) == countPalindromesNaive(s));
        }
    }

    cout << "all tests passed" << '\n';
}

int main() {
    string data, line;
    bool hasInput = false;
    while (getline(cin, line)) {
        data += line;
        data += '\n';
        (void)hasInput;
    }
    // 只要 stdin 有任何內容（哪怕只有一個換行，代表空字串）就走 IO 模式
    if (!data.empty()) runIo(data);
    else runTests();
    return 0;
}
