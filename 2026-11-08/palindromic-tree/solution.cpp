// 回文自動機 Palindromic Tree（Eertree）
//
// 題意：給定字串 s，統計其所有回文子字串：
//   1) 相異回文子字串個數 k（經典結論：<= |s|）；
//   2) 每個相異回文子字串的出現次數（按出現位置計，可重疊）；
//   3) 每個相異回文子字串首次出現的起始下標。
//
// 性質：每個節點代表一個相異回文；節點 0 為虛根（len = -1）、
//       節點 1 為空串根（len = 0）。每加一個字元至多新增一個相異回文，
//       其餘以該位置結尾的回文由 fail（最長真回文後綴）鏈給出；
//       cnt 沿 fail 樹自底向上累加即為實際出現次數。
//
// IO 模式（stdin 有內容時）：
//   一行字串 s（僅小寫字母）
//   輸出：第一行 k；第二行 total（回文子字串總出現次數）；
//         接下來 k 行按 (長度, 首次起始下標) 升冪輸出：len start count
// 無 stdin 時跑內建斷言（與 O(n^2) 暴力對拍）並輸出 all tests passed。
#include <cstdio>
#include <cctype>
#include <cstring>
#include <sstream>
#include <iostream>
#include <string>
#include <vector>
#include <array>
#include <algorithm>
#include <cassert>

using namespace std;
typedef long long ll;

static const int ALPHA = 26;

struct Eertree {
    vector<int> len, link, cnt, endpos;
    vector<array<int, ALPHA> > nxt;
    int last, size;
    const string* s;

    void build(const string& str) {
        s = &str;
        int n = (int)str.size();
        len.assign(n + 3, 0);
        link.assign(n + 3, 0);
        cnt.assign(n + 3, 0);
        endpos.assign(n + 3, 0);
        nxt.assign(n + 3, array<int, ALPHA>());
        for (int i = 0; i < n + 3; i++)
            for (int j = 0; j < ALPHA; j++) nxt[i][j] = -1;
        len[0] = -1; link[0] = 0;   // 虛根
        len[1] = 0;  link[1] = 0;   // 空串根
        last = 1; size = 2;

        for (int i = 0; i < n; i++) {
            int c = str[i] - 'a';
            int cur = getLink(last, i, c);
            if (nxt[cur][c] != -1) {
                last = nxt[cur][c];
                cnt[last]++;
            } else {
                int now = size++;
                len[now] = len[cur] + 2;
                endpos[now] = i;
                cnt[now] = 1;
                nxt[cur][c] = now;
                if (len[now] == 1) link[now] = 1;  // 特判：避免 fail 指回自己
                else link[now] = nxt[getLink(link[cur], i, c)][c];
                last = now;
            }
        }
        for (int v = size - 1; v >= 2; v--) cnt[link[v]] += cnt[v];
    }

    int getLink(int v, int i, int c) const {
        while (true) {
            int l = len[v];
            if (i - l - 1 >= 0 && (*s)[i - l - 1] - 'a' == c) return v;
            v = link[v];
        }
    }
};

static Eertree et;

static vector<array<ll, 3> > solve(const string& str) {
    et.build(str);
    vector<array<ll, 3> > nodes;  // (len, first_start, count)
    for (int v = 2; v < et.size; v++) {
        array<ll, 3> t;
        t[0] = et.len[v];
        t[1] = (ll)et.endpos[v] - et.len[v] + 1;
        t[2] = et.cnt[v];
        nodes.push_back(t);
    }
    sort(nodes.begin(), nodes.end());
    return nodes;
}

static bool isPal(const string& s, int i, int j) {  // [i, j)
    while (i < j - 1) {
        if (s[i] != s[j - 1]) return false;
        i++; j--;
    }
    return true;
}

static vector<array<ll, 3> > brute(const string& s) {
    int n = (int)s.size();
    vector<string> keys;
    vector<ll> cc, fs;
    for (int i = 0; i < n; i++)
        for (int j = i + 1; j <= n; j++) {
            if (!isPal(s, i, j)) continue;
            string sub = s.substr(i, j - i);
            size_t k = 0;
            while (k < keys.size() && keys[k] != sub) k++;
            if (k == keys.size()) {
                keys.push_back(sub); cc.push_back(0); fs.push_back(i);
            }
            cc[k]++;
        }
    vector<array<ll, 3> > nodes;
    for (size_t k = 0; k < keys.size(); k++) {
        array<ll, 3> t;
        t[0] = (ll)keys[k].size(); t[1] = fs[k]; t[2] = cc[k];
        nodes.push_back(t);
    }
    sort(nodes.begin(), nodes.end());
    return nodes;
}

static bool eq(const vector<array<ll, 3> >& a, const vector<array<ll, 3> >& b) {
    if (a.size() != b.size()) return false;
    for (size_t i = 0; i < a.size(); i++)
        for (int j = 0; j < 3; j++)
            if (a[i][j] != b[i][j]) return false;
    return true;
}

static void runAllTests() {
    const char* fixed[] = {"", "a", "aa", "ab", "aba", "aaaa", "abacaba", "abcba",
                           "aabaa", "banana", "abcabc", "abab", "racecar",
                           "aaaaaa", "abaaba"};
    for (const char* f : fixed) {
        string s(f);
        assert(eq(solve(s), brute(s)));
    }
    // 自製 LCG 隨機對拍（與 Python 版同一亂數流）
    unsigned long long seed = 20261108ULL;
    auto rnd = [&](int n) -> int {
        seed = seed * 1103515245ULL + 12345ULL;
        seed &= 0x7fffffffULL;
        return (int)(seed % (unsigned)n);
    };
    for (int tc = 0; tc < 300; tc++) {
        int n = rnd(24);
        string s;
        for (int i = 0; i < n; i++) s += (char)('a' + rnd(3));
        assert(eq(solve(s), brute(s)));
    }
}

int main() {
    stringstream buf;
    buf << cin.rdbuf();
    string input = buf.str();
    bool anyTok = false;
    for (size_t i = 0; i < input.size(); i++)
        if (!isspace((unsigned char)input[i])) { anyTok = true; break; }
    if (!anyTok) { runAllTests(); printf("all tests passed\n"); return 0; }
    istringstream it(input);
    string s;
    it >> s;
    vector<array<ll, 3> > nodes = solve(s);
    ll total = 0;
    for (size_t i = 0; i < nodes.size(); i++) total += nodes[i][2];
    ostringstream oss;
    oss << nodes.size() << '\n' << total << '\n';
    for (size_t i = 0; i < nodes.size(); i++)
        oss << nodes[i][0] << ' ' << nodes[i][1] << ' ' << nodes[i][2] << '\n';
    fputs(oss.str().c_str(), stdout);
    return 0;
}
