// 後綴陣列（Suffix Array）：倍增構造 + Kasai LCP + 最長重複子串 / 不同子串數 / BWT / 子串搜尋
//
// 與 solution.py 完全同構：同樣的演算法、同樣的確定性規則、同樣的輸入輸出格式。
// 編譯：g++ -std=c++17 -O2 -Wall solution.cpp -o solution
// 無 stdin 輸入時運行內置斷言測試並輸出 `all tests passed`。

#include <algorithm>
#include <cassert>
#include <cctype>
#include <iostream>
#include <numeric>
#include <random>
#include <sstream>
#include <string>
#include <vector>

using namespace std;

static const int NEG = -1;   // 越界的排名鍵（-1 比任何合法排名都小，對應「空串最小」）

// ---------------------------------------------------------------- 構造與 LCP

// 倍增法構造後綴陣列，返回 sa（把 n 個後綴按字典序排序後的起始下標）。
static vector<int> buildSA(const string &s) {
    int n = (int)s.size();
    if (n == 0) return {};
    vector<int> sa(n);
    iota(sa.begin(), sa.end(), 0);
    vector<int> rank(n);
    for (int i = 0; i < n; i++) rank[i] = (unsigned char)s[i];
    for (int k = 1; k < n; k <<= 1) {
        // 排序鍵 (rank[i], rank[i+k] 或 -1)，與 Python 的 key 完全一致
        auto second = [&](int i) -> int { return (i + k < n) ? rank[i + k] : NEG; };
        sort(sa.begin(), sa.end(), [&](int a, int b) {
            if (rank[a] != rank[b]) return rank[a] < rank[b];
            int ra = second(a), rb = second(b);
            return ra < rb;
        });
        vector<int> nrank(n);
        nrank[sa[0]] = 0;
        for (int i = 1; i < n; i++) {
            bool same = (rank[sa[i]] == rank[sa[i - 1]]) && (second(sa[i]) == second(sa[i - 1]));
            nrank[sa[i]] = nrank[sa[i - 1]] + (same ? 0 : 1);
        }
        rank = nrank;
        if (rank[sa[n - 1]] == n - 1) break;     // 排名已兩兩不同，排序完成
    }
    return sa;
}

// 暴力構造：直接對所有後綴做字串排序，O(n² log n)，僅供對拍。
static vector<int> saBruteforce(const string &s) {
    int n = (int)s.size();
    vector<int> v(n);
    iota(v.begin(), v.end(), 0);
    sort(v.begin(), v.end(), [&](int a, int b) { return s.compare(a, string::npos, s, b, string::npos) < 0; });
    return v;
}

// 由 sa 得到名次陣列 rank（sa 的逆置換）。
static vector<int> buildRank(const vector<int> &sa) {
    vector<int> rank(sa.size(), 0);
    for (int i = 0; i < (int)sa.size(); i++) rank[sa[i]] = i;
    return rank;
}

// 直接求 s[i:] 與 s[j:] 的最長公共前綴長度，O(n)，僅供對拍。
static int lcpNaive(const string &s, int i, int j) {
    int n = (int)s.size(), h = 0;
    while (i + h < n && j + h < n && s[i + h] == s[j + h]) h++;
    return h;
}

// Kasai 演算法，O(n)：lcp[i] = LCP(sa[i], sa[i-1])，lcp[0] = 0。
static vector<int> kasaiLCP(const string &s, const vector<int> &sa) {
    int n = (int)s.size();
    vector<int> lcp(n, 0);
    if (n == 0) return lcp;
    vector<int> rank = buildRank(sa);
    int h = 0;
    for (int i = 0; i < n; i++) {
        int r = rank[i];
        if (r == 0) {
            h = 0;
            continue;
        }
        int j = sa[r - 1];
        while (i + h < n && j + h < n && s[i + h] == s[j + h]) h++;
        lcp[r] = h;
        if (h > 0) h--;                          // 下一個後綴最多只會少 1，這是線性的關鍵
    }
    return lcp;
}

// ---------------------------------------------------------------- 應用

struct LRS {
    int len;
    int start;
    string sub;
};

// 最長重複子串，返回 (長度, 起始下標, 子串)；不存在則 (0, -1, "")。
// 平手時取 sa 序最前者（即第一次出現位置在 sa 中較靠前）。
static LRS longestRepeatedSubstring(const string &s, const vector<int> &sa, const vector<int> &lcp) {
    int best = 0, start = -1;
    for (int i = 0; i < (int)lcp.size(); i++) {
        if (lcp[i] > best) {                     // 嚴格更優才更新 → 平手取最前，保證確定性
            best = lcp[i];
            start = sa[i];
        }
    }
    LRS res{best, start, best > 0 ? s.substr(start, best) : ""};
    return res;
}

// 暴力求最長重複子串（n ≤ 12）：枚舉所有 (i, j) 求 LCP，平手取起始下標較小者。
static pair<int, int> longestRepeatedBruteforce(const string &s) {
    int n = (int)s.size(), best = 0, start = -1;
    for (int i = 0; i < n; i++)
        for (int j = i + 1; j < n; j++) {
            int h = lcpNaive(s, i, j);
            if (h == 0) continue;
            int si = min(i, j);
            if (h > best || (h == best && (start == -1 || si < start))) {
                best = h;
                start = si;
            }
        }
    return {best, start};
}

// 本質不同的子串數 = n(n+1)/2 - Σ lcp。
static long long countDistinctSubstrings(const string &s, const vector<int> &lcp) {
    long long n = (long long)s.size();
    long long total = n * (n + 1) / 2;
    long long sum = 0;
    for (int x : lcp) sum += x;
    return total - sum;
}

// 暴力統計不同子串數：直接建集合，僅供對拍。
static long long countDistinctBruteforce(const string &s) {
    vector<string> seen;
    for (int i = 0; i < (int)s.size(); i++)
        for (int j = i + 1; j <= (int)s.size(); j++) seen.push_back(s.substr(i, j - i));
    sort(seen.begin(), seen.end());
    seen.erase(unique(seen.begin(), seen.end()), seen.end());
    return (long long)seen.size();
}

// BWT：bwt[i] = s[(sa[i] - 1 + n) % n]，把重複子串聚成連續相同字元。
static string burrowsWheeler(const string &s, const vector<int> &sa) {
    int n = (int)s.size();
    if (n == 0) return "";
    string out;
    out.reserve(n);
    for (int p : sa) out.push_back(s[(p - 1 + n) % n]);
    return out;
}

// 在 sa 上二分搜尋 pat，返回 (出現次數, 出現位置升序列表)。
// 所有以 pat 開頭的後綴在 sa 中是一段連續區間，兩次二分即可定位。
static pair<int, vector<int>> saSearch(const string &s, const vector<int> &sa, const string &pat) {
    int n = (int)s.size();
    if (pat.empty()) return {0, {}};
    int L = (int)pat.size();
    auto prefix = [&](int p) -> string { return s.substr(p, min(L, n - p < 0 ? 0 : n - p)); };
    // 下界：第一個「前綴 ≥ pat」的後綴
    int lo = 0, hi = n;
    while (lo < hi) {
        int mid = (lo + hi) / 2;
        if (prefix(sa[mid]) < pat) lo = mid + 1;
        else hi = mid;
    }
    int left = lo;
    // 上界：第一個「前綴 > pat」的後綴（s 的尾端截斷視為較小，故 == pat 時往右走）
    lo = left;
    hi = n;
    while (lo < hi) {
        int mid = (lo + hi) / 2;
        if (prefix(sa[mid]) > pat) hi = mid;
        else lo = mid + 1;
    }
    int right = lo;
    vector<int> occ(sa.begin() + left, sa.begin() + right);
    sort(occ.begin(), occ.end());
    return {right - left, occ};
}

// 暴力搜尋所有出現位置，僅供對拍。
static pair<int, vector<int>> searchBruteforce(const string &s, const string &pat) {
    if (pat.empty()) return {0, {}};
    vector<int> pos;
    for (int i = 0; i + (int)pat.size() <= (int)s.size(); i++)
        if (s.compare(i, pat.size(), pat) == 0) pos.push_back(i);
    return {(int)pos.size(), pos};
}

// ---------------------------------------------------------------- IO 與測試

static string joinInts(const vector<int> &v) {
    ostringstream oss;
    for (size_t i = 0; i < v.size(); i++) {
        if (i) oss << " ";
        oss << v[i];
    }
    return oss.str();
}

static void runIO(const string &raw) {
    istringstream iss(raw);
    vector<string> toks;
    string tk;
    while (iss >> tk) toks.push_back(tk);
    string s = toks.size() > 0 ? toks[0] : "";
    string pat = toks.size() > 1 ? toks[1] : "";

    vector<int> sa = buildSA(s);
    vector<int> rank = buildRank(sa);
    vector<int> lcp = kasaiLCP(s, sa);
    LRS lrs = longestRepeatedSubstring(s, sa, lcp);
    long long distinct = countDistinctSubstrings(s, lcp);
    string bwt = burrowsWheeler(s, sa);
    auto sr = saSearch(s, sa, pat);

    ostringstream oss;
    oss << joinInts(sa) << "\n" << joinInts(rank) << "\n" << joinInts(lcp) << "\n";
    oss << lrs.len << "\n" << lrs.sub << "\n" << lrs.start << "\n";
    oss << distinct << "\n" << bwt << "\n";
    oss << sr.first << "\n" << joinInts(sr.second) << "\n";
    cout << oss.str();
}

static void runTests() {
    // ---- 空串 ----
    assert(buildSA("").empty());
    assert(kasaiLCP("", {}).empty());
    {
        LRS r = longestRepeatedSubstring("", {}, {});
        assert(r.len == 0 && r.start == -1 && r.sub == "");
    }
    assert(countDistinctSubstrings("", {}) == 0);
    assert(burrowsWheeler("", {}) == "");
    assert(saSearch("", {}, "a").first == 0);

    // ---- 單字元 ----
    assert((buildSA("a") == vector<int>{0}));
    assert((kasaiLCP("a", {0}) == vector<int>{0}));
    {
        LRS r = longestRepeatedSubstring("a", {0}, {0});
        assert(r.len == 0 && r.start == -1 && r.sub == "");
    }
    assert(countDistinctSubstrings("a", {0}) == 1);
    assert(burrowsWheeler("a", {0}) == "a");
    assert(saSearch("a", {0}, "a").first == 1);
    assert((saSearch("a", {0}, "a").second == vector<int>{0}));
    assert(saSearch("a", {0}, "b").first == 0);

    // ---- 經典：banana ----
    {
        string b = "banana";
        vector<int> sa = buildSA(b);
        assert((sa == vector<int>{5, 3, 1, 0, 4, 2}));
        vector<int> lcp = kasaiLCP(b, sa);
        assert((lcp == vector<int>{0, 1, 3, 0, 0, 2}));
        assert((buildRank(sa) == vector<int>{3, 2, 5, 1, 4, 0}));
        LRS r = longestRepeatedSubstring(b, sa, lcp);
        assert(r.len == 3 && r.start == 1 && r.sub == "ana");
        assert(countDistinctSubstrings(b, lcp) == 15);
        assert(countDistinctSubstrings(b, lcp) == countDistinctBruteforce(b));
        assert(saSearch(b, sa, "ana").first == 2);
        assert((saSearch(b, sa, "ana").second == vector<int>{1, 3}));
        assert((saSearch(b, sa, "na").second == vector<int>{2, 4}));
        assert(saSearch(b, sa, "banana").first == 1);
        assert(saSearch(b, sa, "bananas").first == 0);
    }

    // ---- 經典：mississippi ----
    {
        string m = "mississippi";
        vector<int> sa = buildSA(m);
        vector<int> lcp = kasaiLCP(m, sa);
        LRS r = longestRepeatedSubstring(m, sa, lcp);
        assert(r.len == 4 && r.sub == "issi");
        assert(countDistinctSubstrings(m, lcp) == 53);
        assert(countDistinctSubstrings(m, lcp) == countDistinctBruteforce(m));
        assert((saSearch(m, sa, "ssi").second == vector<int>{2, 5}));
        assert((saSearch(m, sa, "issi").second == vector<int>{1, 4}));
    }

    // ---- 全同字元：aaaa ----
    {
        string a4 = "aaaa";
        vector<int> sa = buildSA(a4);
        vector<int> lcp = kasaiLCP(a4, sa);
        assert((lcp == vector<int>{0, 1, 2, 3}));
        LRS r = longestRepeatedSubstring(a4, sa, lcp);
        assert(r.len == 3 && r.start == 0 && r.sub == "aaa");
        assert(countDistinctSubstrings(a4, lcp) == 4);
        assert(burrowsWheeler(a4, sa) == "aaaa");
        assert((saSearch(a4, sa, "aa").second == vector<int>{0, 1, 2}));
    }

    // ---- BWT 經典：abracadabra ----
    {
        string ab = "abracadabra";
        vector<int> sa = buildSA(ab);
        assert((sa == vector<int>{10, 7, 0, 3, 5, 8, 1, 4, 6, 9, 2}));
        vector<int> lcp = kasaiLCP(ab, sa);
        assert((lcp == vector<int>{0, 1, 4, 1, 1, 0, 3, 0, 0, 0, 2}));
        assert(burrowsWheeler(ab, sa) == "rdarcaaaabb");
        LRS r = longestRepeatedSubstring(ab, sa, lcp);
        assert(r.len == 4 && r.start == 0 && r.sub == "abra");
        assert(countDistinctSubstrings(ab, lcp) == 54);
    }

    // ---- 隨機對拍 ----
    mt19937 rng(20261008);
    auto randStr = [&](int n, const string &alpha) {
        string s;
        for (int i = 0; i < n; i++) s.push_back(alpha[rng() % alpha.size()]);
        return s;
    };
    for (int iter = 0; iter < 200; iter++) {
        int n = (int)(rng() % 13);
        string s = randStr(n, "ab");
        vector<int> sa = buildSA(s);
        assert(sa == saBruteforce(s));
        assert((int)sa.size() == n);
        vector<int> sortedSa = sa;
        sort(sortedSa.begin(), sortedSa.end());
        for (int i = 0; i < n; i++) assert(sortedSa[i] == i);
        vector<int> rank = buildRank(sa);
        for (int i = 0; i < n; i++) assert(rank[sa[i]] == i);
        vector<int> lcp = kasaiLCP(s, sa);
        for (int i = 1; i < n; i++) assert(lcp[i] == lcpNaive(s, sa[i], sa[i - 1]));
        LRS r = longestRepeatedSubstring(s, sa, lcp);
        pair<int, int> r2 = longestRepeatedBruteforce(s);
        assert(r.len == r2.first);
        if (r.len > 0) {
            assert(s.substr(r.start, r.len) == r.sub);
            int occ = 0;
            for (int j = 0; j + r.len <= n; j++)
                if (s.substr(j, r.len) == r.sub) occ++;
            assert(occ >= 2);
        }
        assert(countDistinctSubstrings(s, lcp) == countDistinctBruteforce(s));
        string bwt = burrowsWheeler(s, sa);
        assert((int)bwt.size() == n);
        string a1 = bwt, a2 = s;
        sort(a1.begin(), a1.end());
        sort(a2.begin(), a2.end());
        assert(a1 == a2);                        // BWT 是 s 的一個排列
        for (int plen = 1; plen <= 3; plen++) {
            string pat = randStr(plen, "ab");
            auto x1 = saSearch(s, sa, pat);
            auto x2 = searchBruteforce(s, pat);
            assert(x1.first == x2.first);
            assert(x1.second == x2.second);
        }
    }

    // ---- 三元字母表，提高多樣性 ----
    for (int iter = 0; iter < 100; iter++) {
        int n = (int)(rng() % 15);
        string s = randStr(n, "abc");
        vector<int> sa = buildSA(s);
        assert(sa == saBruteforce(s));
        vector<int> lcp = kasaiLCP(s, sa);
        for (int i = 1; i < n; i++) assert(lcp[i] == lcpNaive(s, sa[i], sa[i - 1]));
        assert(countDistinctSubstrings(s, lcp) == countDistinctBruteforce(s));
        LRS r = longestRepeatedSubstring(s, sa, lcp);
        assert(r.len == longestRepeatedBruteforce(s).first);
    }
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    string raw, line;
    bool any = false;
    while (getline(cin, line)) {
        raw += line;
        raw += "\n";
        for (char ch : line)
            if (!isspace((unsigned char)ch)) any = true;
    }
    if (any) {
        runIO(raw);
    } else {
        runTests();
        cout << "all tests passed" << endl;
    }
    return 0;
}
