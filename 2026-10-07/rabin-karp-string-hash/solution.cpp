// Rabin-Karp 與滾動哈希（子串匹配 / 最長公共子串 / 最長回文子串 / 不同子串計數 / k-gram 相似度）
//
// 題意：
//     給定兩個字符串 text 與 pattern（ASCII），用「多項式滾動哈希（Rabin-Karp fingerprint）」解決五個問題：
//       1. pattern 在 text 中的所有出現位置（與 KMP、暴力法三方對拍）；
//       2. 兩串的最長公共子串（長度 + 最靠左的那一個）；
//       3. text 的最長回文子串（長度 + 最靠左的那一個）；
//       4. 兩串各自的不同子串個數；
//       5. 兩串的 3-gram Jaccard 相似度（**應用**：文件指紋 / 抄襲檢測），以最簡分數 `p/q` 輸出。
//
// 思路：
//     ### 多項式滾動哈希
//     把字符串看成 base 進制的多項式：
//         H(s) = ( ord(s[0])*base^(n-1) + ord(s[1])*base^(n-2) + ... + ord(s[n-1]) ) mod M
//     等價的遞推寫法是 `h[i+1] = h[i]*base + ord(s[i])`，預處理 O(n)。
//     則任意子串 `s[l:r]` 的指紋為
//         h[r] - h[l] * base^(r-l)      (mod M)
//     也就是「把前綴左移對齊後相減」，O(1) 得到任意子串指紋。
//     這就是 Rabin-Karp 的核心：窗口滑動時指紋可以增量更新，不必重新掃一遍子串。
//
//     ### 為什麼用雙哈希 + 真實比對
//     單模哈希有生日碰撞風險：n 個子串時碰撞概率約 n^2 / (2M)。
//     這裡同時取兩個質數模數 1e9+7 與 1e9+9（碰撞概率降到 ~1e-18），
//     並且在每個哈希命中的位置**再真實比對一次字符**（verify），
//     因此輸出結果與暴力法逐字節一致，不存在哈希碰撞導致的誤判。
//
//     ### 五個子問題
//       - **子串匹配**：滑動長度 |pat| 的窗口比指紋 → 平均 O(n+m)，最壞與暴力同階 O(n·m) 但常數極小。
//       - **最長公共子串**：「存在長度 L 的公共子串」對 L 單調 ⟹ 二分 L，
//         判定用哈希集合 O(n+m)，總 O((n+m) log min(n,m))，比 O(n·m) 的二維 DP 省空間。
//       - **最長回文子串**：枚舉 2n-1 個回文中心，**對每個中心二分回文半徑**
//         （半徑單調：半徑 r 是回文 ⟹ 半徑 r-1 也是），把「正串窗口」與「反串對應窗口」的指紋比對
//         → O(n log n)。注意不能直接對「長度」二分：存在長度 4 的回文不代表存在長度 3 的回文
//         （例如 `baab`），「存在長度恰為 L 的回文」對 L **不單調**。
//       - **不同子串個數**：枚舉長度 L，把 n-L+1 個窗口指紋塞進集合 → O(n^2) 時間。
//       - **k-gram Jaccard（應用）**：把文本切成長度 k 的滑動片段集合（shingle），
//         相似度 = |A∩B| / |A∪B|。這是抄襲檢測、網頁去重（Broder's shingling）的經典做法，
//         也是「把字符串問題轉成集合問題」的典型應用。
//
// 應用場景：
//     抄襲 / 重複內容檢測（shingling + Jaccard / MinHash）、生物序列比對（DNA 片段 fingerprint）、
//     大文件差分同步（rsync 的弱滾動校驗和）、編譯器與 IDE 的增量字符串搜索。
//
// 複雜度：
//     記 n = len(text)，m = len(pattern)。
//       預處理            O(n + m) 時間、O(n + m) 空間
//       子串匹配          平均 O(n + m)，最壞 O(n·m)
//       最長公共子串      O((n + m) log min(n, m))
//       最長回文子串      O(n log n)
//       不同子串個數      O(n^2)（小規模驗證用）
//       k-gram 相似度     O(n + m)
//
// 輸入格式（stdin）：
//     第 1 行：text
//     第 2 行：pattern
//     （兩行都可以是空行，代表空串；不足兩行時缺的部分視為空串；多餘的行忽略）
// 輸出格式（stdout）：
//     第 1 行：pattern 在 text 中的出現次數
//     第 2 行：所有起始位置（0-indexed，空格分隔；無則輸出空行）
//     第 3 行：最長公共子串長度
//     第 4 行：最長公共子串（最靠左的那個；長度為 0 時輸出空行）
//     第 5 行：text 的最長回文子串長度
//     第 6 行：最長回文子串（最靠左的那個；長度為 0 時輸出空行）
//     第 7 行：text 的不同子串個數
//     第 8 行：pattern 的不同子串個數
//     第 9 行：3-gram Jaccard 相似度，最簡分數 `p/q`
// 無 stdin 輸入時運行內置斷言測試並輸出 `all tests passed`。

#include <algorithm>
#include <cassert>
#include <iostream>
#include <numeric>
#include <random>
#include <set>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

using namespace std;

static const long long MOD1 = 1000000007LL;
static const long long MOD2 = 1000000009LL;
static const long long BASE = 911382323LL;   // 對兩個模數都是奇數且小於模數
static const int K_SHINGLE = 3;              // 應用：k-gram 的窗口長度

typedef pair<long long, long long> Hash2;

// ---------------------------------------------------------------- 滾動哈希

struct RollingHash {
    int n;
    vector<long long> h1, h2, p1, p2;

    explicit RollingHash(const string &s) {
        n = (int)s.size();
        h1.assign(n + 1, 0);
        h2.assign(n + 1, 0);
        p1.assign(n + 1, 1);
        p2.assign(n + 1, 1);
        for (int i = 0; i < n; ++i) {
            long long c = (long long)(unsigned char)s[i];
            h1[i + 1] = (h1[i] * BASE + c) % MOD1;
            h2[i + 1] = (h2[i] * BASE + c) % MOD2;
            p1[i + 1] = (p1[i] * BASE) % MOD1;
            p2[i + 1] = (p2[i] * BASE) % MOD2;
        }
    }

    Hash2 get(int l, int r) const {           // 子串 s[l:r] 的雙指紋
        long long x1 = (h1[r] - h1[l] * p1[r - l]) % MOD1;
        if (x1 < 0) x1 += MOD1;
        long long x2 = (h2[r] - h2[l] * p2[r - l]) % MOD2;
        if (x2 < 0) x2 += MOD2;
        return Hash2(x1, x2);
    }
};

static Hash2 polynomial_hash(const string &s) {
    long long a1 = 0, a2 = 0;
    for (char ch : s) {
        long long c = (long long)(unsigned char)ch;
        a1 = (a1 * BASE + c) % MOD1;
        a2 = (a2 * BASE + c) % MOD2;
    }
    return Hash2(a1, a2);
}

// ---------------------------------------------------------------- 子串匹配

static vector<int> rabin_karp(const string &text, const string &pat) {
    int n = (int)text.size(), m = (int)pat.size();
    vector<int> hits;
    if (m == 0 || m > n) return hits;
    RollingHash rt(text), rp(pat);
    Hash2 target = rp.get(0, m);
    for (int i = 0; i + m <= n; ++i) {
        if (rt.get(i, i + m) == target && text.compare(i, m, pat) == 0) {
            hits.push_back(i);
        }
    }
    return hits;
}

static vector<int> kmp_search(const string &text, const string &pat) {
    int n = (int)text.size(), m = (int)pat.size();
    vector<int> hits;
    if (m == 0 || m > n) return hits;
    vector<int> pi(m, 0);
    for (int i = 1; i < m; ++i) {
        int j = pi[i - 1];
        while (j > 0 && pat[i] != pat[j]) j = pi[j - 1];
        if (pat[i] == pat[j]) ++j;
        pi[i] = j;
    }
    int j = 0;
    for (int i = 0; i < n; ++i) {
        while (j > 0 && text[i] != pat[j]) j = pi[j - 1];
        if (text[i] == pat[j]) ++j;
        if (j == m) {
            hits.push_back(i - m + 1);
            j = pi[j - 1];
        }
    }
    return hits;
}

static vector<int> naive_search(const string &text, const string &pat) {
    int n = (int)text.size(), m = (int)pat.size();
    vector<int> hits;
    if (m == 0 || m > n) return hits;
    for (int i = 0; i + m <= n; ++i) {
        if (text.compare(i, m, pat) == 0) hits.push_back(i);
    }
    return hits;
}

// ---------------------------------------------------------------- 最長公共子串

static bool has_common_len(const string &a, const string &b,
                           const RollingHash &rha, const RollingHash &rhb, int length) {
    if (length <= 0) return true;
    if (length > (int)a.size() || length > (int)b.size()) return false;
    set<Hash2> seen;
    for (int j = 0; j + length <= (int)b.size(); ++j) seen.insert(rhb.get(j, j + length));
    for (int i = 0; i + length <= (int)a.size(); ++i) {
        if (seen.count(rha.get(i, i + length))) {
            if (b.find(a.substr(i, length)) != string::npos) return true;   // 真實比對
        }
    }
    return false;
}

static pair<int, string> longest_common_substring(const string &a, const string &b) {
    if (a.empty() || b.empty()) return make_pair(0, string(""));
    RollingHash rha(a), rhb(b);
    int lo = 0, hi = min((int)a.size(), (int)b.size());
    while (lo < hi) {
        int mid = (lo + hi + 1) / 2;
        if (has_common_len(a, b, rha, rhb, mid)) lo = mid;
        else hi = mid - 1;
    }
    if (lo == 0) return make_pair(0, string(""));
    for (int i = 0; i + lo <= (int)a.size(); ++i) {      // 取最靠左的 i
        if (b.find(a.substr(i, lo)) != string::npos) return make_pair(lo, a.substr(i, lo));
    }
    return make_pair(0, string(""));
}

static pair<int, string> lcs_bruteforce(const string &a, const string &b) {
    int best = 0, bi = 0;
    for (int i = 0; i < (int)a.size(); ++i) {
        for (int j = 0; j < (int)b.size(); ++j) {
            int k = 0;
            while (i + k < (int)a.size() && j + k < (int)b.size() && a[i + k] == b[j + k]) ++k;
            if (k > best) { best = k; bi = i; }
        }
    }
    return make_pair(best, best ? a.substr(bi, best) : string(""));
}

// ---------------------------------------------------------------- 最長回文子串

static bool is_pal_window(const string &s, const RollingHash &rhs, const RollingHash &rhr,
                          int l, int r) {
    int n = (int)s.size();
    if (l >= r) return true;                             // 空窗口（偶中心半徑 0）視為回文
    if (rhs.get(l, r) != rhr.get(n - r, n - l)) return false;
    string seg = s.substr(l, r - l);                     // 真實比對，杜絕哈希碰撞誤判
    string rev = seg;
    reverse(rev.begin(), rev.end());
    return seg == rev;
}

static pair<int, string> longest_palindrome_hash(const string &s) {
    int n = (int)s.size();
    if (n == 0) return make_pair(0, string(""));
    string rs = s;
    reverse(rs.begin(), rs.end());
    RollingHash rhs(s), rhr(rs);
    int best_len = 1, best_i = 0;

    for (int c = 0; c < n; ++c) {                        // 奇數長度：窗口 [c-r, c+r+1)
        int lo = 0, hi = min(c, n - 1 - c);
        while (lo < hi) {
            int mid = (lo + hi + 1) / 2;
            if (is_pal_window(s, rhs, rhr, c - mid, c + mid + 1)) lo = mid;
            else hi = mid - 1;
        }
        int length = 2 * lo + 1, start = c - lo;
        if (length > best_len || (length == best_len && start < best_i)) {
            best_len = length;
            best_i = start;
        }
    }
    for (int c = 1; c < n; ++c) {                        // 偶數長度：窗口 [c-r, c+r)
        int lo = 0, hi = min(c, n - c);
        while (lo < hi) {
            int mid = (lo + hi + 1) / 2;
            if (is_pal_window(s, rhs, rhr, c - mid, c + mid)) lo = mid;
            else hi = mid - 1;
        }
        if (lo >= 1) {
            int length = 2 * lo, start = c - lo;
            if (length > best_len || (length == best_len && start < best_i)) {
                best_len = length;
                best_i = start;
            }
        }
    }
    return make_pair(best_len, s.substr(best_i, best_len));
}

static pair<int, string> palindrome_bruteforce(const string &s) {
    int n = (int)s.size();
    int best = 0, bi = 0;
    for (int c = 0; c < n; ++c) {
        for (int t = 0; t < 2; ++t) {
            int i = c, j = (t == 0 ? c : c + 1);
            while (i >= 0 && j < n && s[i] == s[j]) { --i; ++j; }
            int len = j - i - 1;
            if (len > best || (len == best && i + 1 < bi)) { best = len; bi = i + 1; }
        }
    }
    return make_pair(best, best ? s.substr(bi, best) : string(""));
}

// ---------------------------------------------------------------- 不同子串 / 應用

static long long count_distinct_substrings(const string &s) {
    int n = (int)s.size();
    set<Hash2> seen;
    RollingHash rh(s);
    for (int length = 1; length <= n; ++length) {
        for (int i = 0; i + length <= n; ++i) seen.insert(rh.get(i, i + length));
    }
    return (long long)seen.size();
}

static long long distinct_substrings_bruteforce(const string &s) {
    set<string> seen;
    int n = (int)s.size();
    for (int i = 0; i < n; ++i)
        for (int j = i + 1; j <= n; ++j) seen.insert(s.substr(i, j - i));
    return (long long)seen.size();
}

static set<string> shingles(const string &s, int k = K_SHINGLE) {
    set<string> out;
    int n = (int)s.size();
    if (n == 0) return out;
    if (n <= k) { out.insert(s); return out; }
    for (int i = 0; i + k <= n; ++i) out.insert(s.substr(i, k));
    return out;
}

static string shingle_jaccard(const string &a, const string &b, int k = K_SHINGLE) {
    set<string> sa = shingles(a, k), sb = shingles(b, k);
    int inter = 0;
    for (const string &x : sa) if (sb.count(x)) ++inter;
    int union_sz = (int)(sa.size() + sb.size()) - inter;
    if (union_sz == 0) return "1/1";
    int g = (int)std::gcd(inter, union_sz);
    return to_string(inter / g) + "/" + to_string(union_sz / g);
}

// ---------------------------------------------------------------- IO 與測試

static string join_ints(const vector<int> &v) {
    string out;
    for (size_t i = 0; i < v.size(); ++i) {
        if (i) out += " ";
        out += to_string(v[i]);
    }
    return out;
}

static void run_io(const string &raw) {
    vector<string> lines;
    {
        string cur;
        for (char ch : raw) {
            if (ch == '\n') { lines.push_back(cur); cur.clear(); }
            else cur.push_back(ch);
        }
        if (!cur.empty()) lines.push_back(cur);   // 等價於 Python 的 split('\n') 後去掉尾端空串
    }
    for (string &ln : lines) {
        if (!ln.empty() && ln.back() == '\r') ln.pop_back();
    }
    string text = lines.size() > 0 ? lines[0] : string("");
    string pat = lines.size() > 1 ? lines[1] : string("");

    vector<int> pos = rabin_karp(text, pat);
    pair<int, string> lc = longest_common_substring(text, pat);
    pair<int, string> lp = longest_palindrome_hash(text);

    vector<string> out;
    out.push_back(to_string(pos.size()));
    out.push_back(join_ints(pos));
    out.push_back(to_string(lc.first));
    out.push_back(lc.second);
    out.push_back(to_string(lp.first));
    out.push_back(lp.second);
    out.push_back(to_string(count_distinct_substrings(text)));
    out.push_back(to_string(count_distinct_substrings(pat)));
    out.push_back(shingle_jaccard(text, pat));
    for (size_t i = 0; i < out.size(); ++i) cout << out[i] << (i + 1 == out.size() ? "" : "\n");
    cout << "\n";
}

static void run_tests() {
    // ---- 固定用例 ----
    assert(rabin_karp("abracadabra", "abra") == vector<int>({0, 7}));
    assert(rabin_karp("aaaaa", "aa") == vector<int>({0, 1, 2, 3}));
    assert(rabin_karp("abc", "").empty());
    assert(rabin_karp("", "a").empty());
    assert(rabin_karp("abc", "abcd").empty());
    assert(rabin_karp("a", "a") == vector<int>({0}));
    assert(kmp_search("abracadabra", "abra") == vector<int>({0, 7}));
    assert(naive_search("abracadabra", "abra") == vector<int>({0, 7}));

    assert(longest_common_substring("abcdef", "zcdemf") == make_pair(3, string("cde")));
    assert(longest_common_substring("", "abc") == make_pair(0, string("")));
    assert(longest_common_substring("abc", "def") == make_pair(0, string("")));

    assert(longest_palindrome_hash("babad") == make_pair(3, string("bab")));
    assert(longest_palindrome_hash("cbbd") == make_pair(2, string("bb")));
    assert(longest_palindrome_hash("") == make_pair(0, string("")));
    assert(longest_palindrome_hash("a") == make_pair(1, string("a")));

    {   // 整串指紋與前綴哈希口徑一致
        RollingHash rh("abracadabra");
        assert(polynomial_hash("abracadabra") == rh.get(0, 11));
        assert(polynomial_hash("") == rh.get(0, 0));
    }

    assert(count_distinct_substrings("aaa") == 3);
    assert(count_distinct_substrings("") == 0);
    assert(shingle_jaccard("abcde", "abcde") == "1/1");
    assert(shingle_jaccard("abcde", "fghij") == "0/1");

    mt19937 rng(20261007);
    const string alphabet = "ab";

    // ---- 隨機對拍：三種匹配實現必須一致 ----
    for (int t = 0; t < 400; ++t) {
        int n = (int)(rng() % 31), m = (int)(rng() % 6);
        string text, pat;
        for (int i = 0; i < n; ++i) text.push_back(alphabet[rng() % 2]);
        for (int i = 0; i < m; ++i) pat.push_back(alphabet[rng() % 2]);
        assert(rabin_karp(text, pat) == kmp_search(text, pat));
        assert(kmp_search(text, pat) == naive_search(text, pat));
    }

    // ---- 隨機對拍：最長公共子串 / 最長回文子串 / 不同子串 ----
    for (int t = 0; t < 300; ++t) {
        int n = (int)(rng() % 15), m = (int)(rng() % 15);
        string a, b;
        for (int i = 0; i < n; ++i) a.push_back(alphabet[rng() % 2]);
        for (int i = 0; i < m; ++i) b.push_back(alphabet[rng() % 2]);

        pair<int, string> got = longest_common_substring(a, b);
        pair<int, string> want = lcs_bruteforce(a, b);
        assert(got.first == want.first);
        if (want.first > 0) {
            string rev = got.second;
            reverse(rev.begin(), rev.end());
            assert(b.find(got.second) != string::npos);
        }

        pair<int, string> gp = longest_palindrome_hash(a);
        pair<int, string> wp = palindrome_bruteforce(a);
        assert(gp.first == wp.first);
        {
            string rev = gp.second;
            reverse(rev.begin(), rev.end());
            assert(gp.second == rev);
            assert(a.find(gp.second) != string::npos);
            assert((int)gp.second.size() == gp.first);
        }

        assert(count_distinct_substrings(a) == distinct_substrings_bruteforce(a));
        assert(count_distinct_substrings(b) == distinct_substrings_bruteforce(b));
    }

    // ---- 應用：k-gram Jaccard 與直接集合計算對拍 ----
    const string ab3 = "abc";
    for (int t = 0; t < 200; ++t) {
        int n = (int)(rng() % 21), m = (int)(rng() % 21);
        string a, b;
        for (int i = 0; i < n; ++i) a.push_back(ab3[rng() % 3]);
        for (int i = 0; i < m; ++i) b.push_back(ab3[rng() % 3]);
        set<string> sa = shingles(a), sb = shingles(b);
        int inter = 0;
        for (const string &x : sa) if (sb.count(x)) ++inter;
        int union_sz = (int)(sa.size() + sb.size()) - inter;
        if (union_sz == 0) {
            assert(shingle_jaccard(a, b) == "1/1");
        } else {
            int g = (int)std::gcd(inter, union_sz);
            assert(shingle_jaccard(a, b) == to_string(inter / g) + "/" + to_string(union_sz / g));
        }
    }

    cout << "all tests passed" << endl;
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    ostringstream oss;
    oss << cin.rdbuf();
    string raw = oss.str();

    bool has_input = raw.find_first_not_of(" \t\r\n") != string::npos;
    if (has_input) run_io(raw);
    else run_tests();
    return 0;
}
