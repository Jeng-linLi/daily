// Lyndon 分解（Duval 演算法）與循環串最小／最大表示法（Booth 演算法）。
//
// 命令（IO 模式，從 stdin 讀入；S 為非空 token）：
//   lyndon    S   -> Lyndon 分解，各因子以單個空格分隔
//   minrot    S   -> 最小表示的起始下標（0-based）與該循環串
//   maxrot    S   -> 最大表示的起始下標（0-based）與該循環串
//
// stdin 無內容時跑內建斷言測試，輸出 all tests passed。
//
// 核心要點：
//   * Lyndon 詞：嚴格小於它所有真後綴的非空串。任何串都能唯一分解為
//     c1 c2 ... ck，每個 ci 都是 Lyndon 詞，且 c1 >= c2 >= ... >= ck（字典序）。
//   * Duval：從 i 開始延伸「週期性主鏈」——比較 s[k] 與 s[j]，相等則 k++
//     （週期不變），s[k] < s[j] 則 k 回到 i（更有希望的新週期），
//     s[k] > s[j] 則主鏈斷開。主鏈長度 j-i、週期 j-k，於是從 i 起連續切出
//     (j-k) 長度的因子直到超過 k。全程 i、j 單調前進，O(n)。
//   * Booth：把 s 複製成 s+s，維護兩個候選起點 i、j 與已匹配長度 k。
//     s[i+k] > s[j+k] 時 i 這一段全部不可能（i 到 i+k 都輸給對應的 j 段），
//     直接跳到 i+k+1；反之跳 j。k 從不倒退，O(n)。
//   * 最大表示只是把 Booth 的比較方向反過來（s[i+k] < s[j+k] 時淘汰 i）。

#include <cassert>
#include <cstdio>
#include <cstring>
#include <string>
#include <vector>
#include <iostream>
#include <sstream>

#ifdef _WIN32
#include <fcntl.h>
#include <io.h>
#endif

using namespace std;

// ------------------------------ 核心演算法 ------------------------------

// Lyndon 分解：回傳因子列表（非空串）。
static vector<string> duval(const string &s) {
    int n = (int)s.size();
    vector<string> res;
    int i = 0;
    while (i < n) {
        int j = i + 1, k = i;  // k 是主鏈的週期對應位置
        while (j < n && s[k] <= s[j]) {
            if (s[k] < s[j]) k = i;   // 出現更大的字元：週期重設
            else k++;                 // 相等：週期延續
            j++;
        }
        // 主鏈 break 掉，長度為 j-i 的週期節，可切出若干個長 (j-k) 的因子
        while (i <= k) {
            res.push_back(s.substr(i, j - k));
            i += j - k;
        }
    }
    return res;
}

// 回傳循環串最小（或最大）表示的起始下標。空串回傳 0。
static long long booth(const string &s, bool want_max = false) {
    int n = (int)s.size();
    if (n == 0) return 0;
    string ss = s + s;
    int i = 0, j = 1, k = 0;
    while (i < n && j < n && k < n) {
        char a = ss[i + k], b = ss[j + k];
        if (a == b) {
            k++;
            continue;
        }
        // want_max 時「a 較小」才淘汰 i，其餘與最小版對稱
        if (want_max ? (a < b) : (a > b)) {
            i += k + 1;
            if (i == j) i++;
        } else {
            j += k + 1;
            if (i == j) j++;
        }
        k = 0;
    }
    return (long long)min(i, j);
}

// 取循環串從 start 開始、長度 len(s) 的那一段。
static string rotation(const string &s, int start) {
    int n = (int)s.size();
    if (n == 0) return "";
    if (start < n) return s.substr(start) + s.substr(0, start);
    return s;
}

// ----------------------------- 暴力參考實作 -----------------------------

static bool is_lyndon(const string &t) {
    // t 嚴格小於所有真後綴。
    for (int p = 1; p < (int)t.size(); p++)
        if (!(t < t.substr(p))) return false;
    return true;
}

// 貪心：每次取最長的 Lyndon 前綴（O(n^3)，僅供對拍）。
static vector<string> brute_decompose(const string &s) {
    vector<string> res;
    int n = (int)s.size();
    int i = 0;
    while (i < n) {
        int best = -1;
        for (int j = n; j > i; j--) {
            if (is_lyndon(s.substr(i, j - i))) {
                best = j;
                break;
            }
        }
        if (best < 0) best = i + 1;  // 單個字元一定是 Lyndon 詞，正常不會走到
        res.push_back(s.substr(i, best - i));
        i = best;
    }
    return res;
}

// 枚舉所有循環移位取字典序最小（最大）者，平手取最小下標。
static pair<int, string> brute_rotation(const string &s, bool want_max = false) {
    int n = (int)s.size();
    if (n == 0) return {0, ""};
    int best_i = 0;
    string best;
    for (int p = 0; p < n; p++) {
        string cur = s.substr(p) + s.substr(0, p);
        if (p == 0 || (want_max ? cur > best : cur < best)) {
            best_i = p;
            best = cur;
        }
    }
    return {best_i, best};
}

// ------------------------------- IO 模式 -------------------------------

static string run_io(const vector<string> &toks) {
    size_t pos = 0;
    string out;
    auto nxt = [&]() -> const string * {
        if (pos >= toks.size()) return nullptr;
        return &toks[pos++];
    };
    while (true) {
        const string *cp = nxt();
        if (cp == nullptr) break;
        string cmd = *cp;
        if (cmd != "lyndon" && cmd != "minrot" && cmd != "maxrot") break;
        const string *sp = nxt();
        if (sp == nullptr) break;
        string s = *sp;
        if (cmd == "lyndon") {
            vector<string> f = duval(s);
            for (size_t t = 0; t < f.size(); t++) {
                if (t) out += ' ';
                out += f[t];
            }
            out += '\n';
        } else {
            bool want_max = (cmd == "maxrot");
            int p = (int)booth(s, want_max);
            out += to_string(p);
            out += ' ';
            out += rotation(s, p);
            out += '\n';
        }
    }
    fputs(out.c_str(), stdout);
    return out;
}

// ------------------------------ 內建測試 ------------------------------

static unsigned long long g_state = 20261114ULL;

static int rnd(int m) {
    g_state = (g_state * 1103515245ULL + 12345ULL) & 0x7FFFFFFFULL;
    return (int)(g_state % (unsigned)m);
}

static void run_tests() {
    // 1) 固定用例
    {
        vector<string> w1 = {"abcde"};
        assert(duval("abcde") == w1);
        vector<string> w2 = {"b", "an", "an", "a"};
        assert(duval("banana") == w2);
        vector<string> w3 = {"a", "a", "a", "a"};
        assert(duval("aaaa") == w3);
        vector<string> w4 = {"ababb"};
        assert(duval("ababb") == w4);
        assert(booth("banana") == 5);              // "abanan"
        assert(rotation("banana", 5) == "abanan");
        assert(booth("cba") == 2);                 // "acb"
        assert(booth("a") == 0);
        assert(booth("banana", true) == 2);        // "nanaba"
        assert(booth("abab") == 0);                // 平手取最小下標
        vector<string> w5 = {"c", "abc", "a"};
        assert(duval("cabca") == w5);
    }

    // 2) 隨機對拍：Duval vs 貪心暴力；Booth vs 全枚舉
    for (int tc = 0; tc < 400; tc++) {
        int n = rnd(24) + 1;
        int alpha = rnd(3) + 1;
        string s;
        for (int t = 0; t < n; t++) s += (char)('a' + rnd(alpha));
        vector<string> got = duval(s);
        vector<string> want = brute_decompose(s);
        assert(got == want);
        // Lyndon 分解的性質：每項是 Lyndon 詞且非遞增
        for (const string &t : got) assert(is_lyndon(t));
        for (size_t t = 0; t + 1 < got.size(); t++) assert(got[t] >= got[t + 1]);
        {
            string acc;
            for (const string &t : got) acc += t;
            assert(acc == s);
        }
        pair<int, string> pi = brute_rotation(s, false);
        assert(booth(s, false) == pi.first);
        assert(rotation(s, (int)booth(s, false)) == pi.second);
        pair<int, string> qi = brute_rotation(s, true);
        assert(booth(s, true) == qi.first);
        assert(rotation(s, (int)booth(s, true)) == qi.second);
    }

    // 3) 週期串與長串（考驗 k 不倒退的前提）
    for (int rep = 1; rep < 40; rep++) {
        string s;
        for (int t = 0; t < rep; t++) s += "ab";
        assert(booth(s) == 0);
        vector<string> want(rep, "ab");
        assert(duval(s) == want);
    }
    {
        string big(5000, 'z');
        big += 'a';
        big += string(4999, 'z');
        assert(booth(big) == 5000);        // 最小表示從那個 'a' 開始
        string allA(3000, 'a');
        vector<string> want(3000, "a");
        assert(duval(allA) == want);
    }
}

int main() {
#ifdef _WIN32
    _setmode(_fileno(stdout), _O_BINARY);
#endif
    // stdin 有內容走 IO 模式，否則跑內建測試
    ostringstream buf;
    buf << cin.rdbuf();
    string input = buf.str();
    if (!input.empty() && input.find_first_not_of(" \t\r\n") != string::npos) {
        istringstream iss(input);
        vector<string> toks;
        string t;
        while (iss >> t) toks.push_back(t);
        run_io(toks);
    } else {
        run_tests();
        fputs("all tests passed\n", stdout);
    }
    return 0;
}
