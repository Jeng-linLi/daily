// 後綴自動機（Suffix Automaton, SAM）
//
// 題意：
//     給定字串 s，建立後綴自動機並回答五種查詢：
//       1. s 有幾個不同的子串        2. t 是否為 s 的子串
//       3. t 在 s 中出現幾次          4. s 與 t 的最長公共子串
//       5. 字典序第 k 小的子串
//
// 思路：
//     SAM 是接受 s 所有子串的最小 DFA，用 O(n) 個狀態裝下 O(n^2) 個子串
//     （狀態數 ≤ 2n-1、轉移數 ≤ 3n-4）。
//     每個狀態 v 代表一個 endpos 等價類，記 length[v]（最長串長度）與
//     link[v]（後綴鏈），則 v 代表的串長區間為 (length[link[v]], length[v]]，
//     故不同子串數 = Σ (length[v] - length[link[v]])。
//
//     出現次數：構造時每個新建狀態 cur 記 cnt=1（克隆記 0），
//     再按 length 由大到小把 cnt[v] 累加到 cnt[link[v]]，即完成拓撲累積。
//     第 k 小子串：由根出發的路徑與子串一一對應，令
//     dp[v] = Σ_c (1 + dp[nxt[v][c]])，按 length 由大到小算出後從根往下走即可。
//     按字元排序走訪 = 前綴樹前序走訪 = 字典序。
//
// 複雜度：構造 O(n log α)；判定 / 出現次數 / LCS 皆 O(|t|)；第 k 小 O(長度 · log α)；空間 O(n)。
//
// 輸入（stdin，空白分隔）：s  t  k       （k 非數字或缺失視為 0；字串截斷到 3000 位元組）
// 輸出（stdout，6 行）：
//     1: 不同子串數      2: t 是否為子串(1/0)   3: t 出現次數
//     4: LCS 長度        5: LCS 本身             6: 字典序第 k 小子串（越界為空行）
// 無 stdin 時跑內置斷言測試並輸出 all tests passed。

#include <cassert>
#include <cerrno>
#include <cstdlib>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <map>
#include <random>
#include <set>
#include <string>
#include <vector>

#ifdef _WIN32
#include <fcntl.h>
#include <io.h>
#endif

namespace {

constexpr std::size_t MAX_LEN = 3000;   // 輸入字串截斷長度（與 Python 版一致）

// ---------------------------------------------------------------- 後綴自動機
class SuffixAutomaton {
public:
    explicit SuffixAutomaton(const std::string& s = "") {
        nxt_.emplace_back();
        link_.push_back(-1);
        length_.push_back(0);
        last_ = 0;
        for (unsigned char ch : s) extend(ch);
        buildCnt();
        buildDp();
    }

    long long distinctCount() const {
        long long byLink = 0;
        for (std::size_t v = 1; v < nxt_.size(); ++v)
            byLink += static_cast<long long>(length_[v]) - length_[link_[v]];
        assert(byLink == dp_[0]);        // 兩種算法互相對拍
        return dp_[0];
    }

    int walk(const std::string& t) const {
        int v = 0;
        for (unsigned char ch : t) {
            auto it = nxt_[v].find(ch);
            if (it == nxt_[v].end()) return -1;
            v = it->second;
        }
        return v;
    }

    bool isSubstring(const std::string& t) const { return walk(t) != -1; }

    long long countOccurrences(const std::string& t) const {
        int v = walk(t);
        return v < 0 ? 0 : cnt_[v];
    }

    // 回傳 (最長公共子串長度, 字串本身)；並列時取 t 中最早出現者
    std::pair<int, std::string> longestCommonSubstring(const std::string& t) const {
        int v = 0, matched = 0, best = 0, bestPos = -1;
        for (int i = 0; i < static_cast<int>(t.size()); ++i) {
            unsigned char ch = static_cast<unsigned char>(t[i]);
            while (v != 0 && nxt_[v].find(ch) == nxt_[v].end()) {
                v = link_[v];            // 沿後綴鏈退，縮短目前匹配
                matched = length_[v];
            }
            auto it = nxt_[v].find(ch);
            if (it == nxt_[v].end()) {
                v = 0;
                matched = 0;
            } else {
                v = it->second;
                matched += 1;
            }
            if (matched > best) {
                best = matched;
                bestPos = i;
            }
        }
        std::string text;
        if (best > 0) text = t.substr(bestPos - best + 1, static_cast<std::size_t>(best));
        return {best, text};
    }

    std::string kthSubstring(long long k) const {
        if (k < 1 || k > dp_[0]) return "";
        std::string out;
        int v = 0;
        while (k > 0) {
            bool moved = false;
            for (const auto& kv : nxt_[v]) {          // std::map 已按 unsigned char 排序
                long long block = 1 + dp_[kv.second];
                if (k <= block) {
                    out.push_back(static_cast<char>(kv.first));
                    k -= 1;
                    v = kv.second;
                    moved = true;
                    break;
                }
                k -= block;
            }
            if (!moved) break;                        // 理論上不會發生
        }
        return out;
    }

    std::size_t stateCount() const { return nxt_.size(); }
    const std::vector<int>& link() const { return link_; }
    const std::vector<int>& length() const { return length_; }
    long long cnt(int v) const { return cnt_[v]; }

private:
    std::vector<std::map<unsigned char, int>> nxt_;
    std::vector<int> link_, length_, created_;
    std::vector<long long> cnt_, dp_;
    int last_;

    void extend(unsigned char ch) {
        int cur = static_cast<int>(nxt_.size());
        nxt_.emplace_back();
        link_.push_back(0);
        length_.push_back(length_[last_] + 1);
        int p = last_;
        while (p != -1 && nxt_[p].find(ch) == nxt_[p].end()) {
            nxt_[p][ch] = cur;
            p = link_[p];
        }
        if (p == -1) {
            link_[cur] = 0;
        } else {
            int q = nxt_[p][ch];
            if (length_[p] + 1 == length_[q]) {
                link_[cur] = q;
            } else {
                // q 代表的串長短不一 → 克隆出只代表「短的那些」的狀態
                int clone = static_cast<int>(nxt_.size());
                nxt_.push_back(nxt_[q]);
                length_.push_back(length_[p] + 1);
                link_.push_back(link_[q]);
                while (p != -1) {
                    auto it = nxt_[p].find(ch);
                    if (it == nxt_[p].end() || it->second != q) break;
                    it->second = clone;
                    p = link_[p];
                }
                link_[q] = clone;
                link_[cur] = clone;
            }
        }
        last_ = cur;
        created_.push_back(cur);
    }

    std::vector<int> order() const {
        std::vector<int> ord(nxt_.size());
        for (std::size_t i = 0; i < nxt_.size(); ++i) ord[i] = static_cast<int>(i);
        std::sort(ord.begin(), ord.end(), [this](int a, int b) {
            if (length_[a] != length_[b]) return length_[a] > length_[b];
            return a < b;
        });
        return ord;
    }

    void buildCnt() {
        cnt_.assign(nxt_.size(), 0);
        for (int v : created_) cnt_[v] = 1;
        for (int v : order()) {
            int p = link_[v];
            if (p >= 0) cnt_[p] += cnt_[v];
        }
    }

    void buildDp() {
        dp_.assign(nxt_.size(), 0);
        for (int v : order()) {
            long long total = 0;
            for (const auto& kv : nxt_[v]) total += 1 + dp_[kv.second];
            dp_[v] = total;
        }
    }
};

// ---------------------------------------------------------------- 暴力對拍
std::vector<std::string> allSubstrings(const std::string& s) {
    std::vector<std::string> v;
    for (std::size_t i = 0; i < s.size(); ++i)
        for (std::size_t j = i + 1; j <= s.size(); ++j) v.push_back(s.substr(i, j - i));
    return v;
}

long long bfDistinct(const std::string& s) {
    std::vector<std::string> v = allSubstrings(s);
    std::set<std::string> st(v.begin(), v.end());
    return static_cast<long long>(st.size());
}

long long bfOccurrences(const std::string& s, const std::string& t) {
    if (t.empty()) return static_cast<long long>(s.size());
    long long ans = 0;
    std::size_t pos = s.find(t);
    while (pos != std::string::npos) {
        ++ans;
        pos = s.find(t, pos + 1);
    }
    return ans;
}

int bfLcsLen(const std::string& s, const std::string& t) {
    std::vector<std::string> v = allSubstrings(s);
    std::set<std::string> subs(v.begin(), v.end());
    int best = 0;
    for (std::size_t i = 0; i < t.size(); ++i)
        for (std::size_t j = i + 1; j <= t.size(); ++j)
            if (subs.count(t.substr(i, j - i)) && static_cast<int>(j - i) > best)
                best = static_cast<int>(j - i);
    return best;
}

std::string bfKth(const std::string& s, long long k) {
    std::vector<std::string> v = allSubstrings(s);
    std::sort(v.begin(), v.end());
    v.erase(std::unique(v.begin(), v.end()), v.end());
    if (k < 1 || k > static_cast<long long>(v.size())) return "";
    return v[static_cast<std::size_t>(k - 1)];
}

// ---------------------------------------------------------------- 輸出
std::string run(const std::string& s, const std::string& t, long long k) {
    SuffixAutomaton sam(s);
    std::pair<int, std::string> lcs = sam.longestCommonSubstring(t);
    std::string out;
    out += std::to_string(sam.distinctCount());
    out += "\n";
    out += sam.isSubstring(t) ? "1" : "0";
    out += "\n";
    out += std::to_string(sam.countOccurrences(t));
    out += "\n";
    out += std::to_string(lcs.first);
    out += "\n";
    out += lcs.second;
    out += "\n";
    out += sam.kthSubstring(k);
    return out;
}

void check(const std::string& s, const std::string& t, long long k) {
    SuffixAutomaton sam(s);
    assert(sam.distinctCount() == bfDistinct(s));
    assert(sam.isSubstring(t) == (s.find(t) != std::string::npos));
    assert(sam.countOccurrences(t) == bfOccurrences(s, t));
    std::pair<int, std::string> lcs = sam.longestCommonSubstring(t);
    assert(lcs.first == bfLcsLen(s, t));
    if (lcs.first > 0) {
        assert(s.find(lcs.second) != std::string::npos);
        assert(t.find(lcs.second) != std::string::npos);
    }
    assert(sam.kthSubstring(k) == bfKth(s, k));
}

// ---------------------------------------------------------------- 解析
bool isSpaceCh(unsigned char c) {
    return c == ' ' || c == '\t' || c == '\n' || c == '\r' || c == '\v' || c == '\f';
}

std::vector<std::string> splitWs(const std::string& raw) {
    std::vector<std::string> out;
    std::string cur;
    for (unsigned char c : raw) {
        if (isSpaceCh(c)) {
            if (!cur.empty()) { out.push_back(cur); cur.clear(); }
        } else {
            cur.push_back(static_cast<char>(c));
        }
    }
    if (!cur.empty()) out.push_back(cur);
    return out;
}

bool isNumToken(const std::string& tok) {
    if (tok.empty()) return false;
    std::size_t i = (tok[0] == '+' || tok[0] == '-') ? 1 : 0;
    if (i >= tok.size()) return false;
    for (; i < tok.size(); ++i)
        if (tok[i] < '0' || tok[i] > '9') return false;
    return true;
}

long long parseLLClamp(const std::string& tok) {
    errno = 0;
    char* end = nullptr;
    long long v = std::strtoll(tok.c_str(), &end, 10);
    if (end == tok.c_str()) return 0;
    return v;                      // strtoll 溢位時已飽和到 LLONG_MAX / LLONG_MIN
}

// ---------------------------------------------------------------- 測試
void runTests() {
    std::mt19937 rng(20261013);

    // ---- 固定用例（與 Python 版逐字節一致）----
    assert(run("aba", "babab", 4) == "5\n0\n0\n3\naba\nb");
    assert(run("aaaa", "aa", 3) == "4\n1\n3\n2\naa\naaa");
    assert(run("", "", 1) == "0\n1\n0\n0\n\n");
    assert(run("abc", "", 1) == "6\n1\n3\n0\n\na");
    assert(run("a", "a", 1) == "1\n1\n1\n1\na\na");
    assert(run("abab", "baba", 5) == "7\n0\n0\n3\nbab\nb");

    check("aba", "babab", 4);
    check("aaaa", "aa", 3);
    check("", "", 1);
    check("abc", "", 1);
    check("abab", "baba", 5);
    check("mississippi", "issipp", 20);

    // ---- 隨機對拍（小字串，暴力列舉）----
    for (int iter = 0; iter < 500; ++iter) {
        std::string alpha = (rng() % 3 == 0) ? "ab" : ((rng() % 2 == 0) ? "abc" : "abcd");
        int ns = static_cast<int>(rng() % 13);
        int nt = static_cast<int>(rng() % 9);
        std::string s, t;
        for (int i = 0; i < ns; ++i) s.push_back(alpha[rng() % alpha.size()]);
        for (int i = 0; i < nt; ++i) t.push_back(alpha[rng() % alpha.size()]);
        long long k = static_cast<long long>(rng() % 33) - 2;
        check(s, t, k);
    }

    // ---- 隨機對拍（大字串，只比「不同子串數 / 第 k 小 / 出現次數」）----
    for (int iter = 0; iter < 60; ++iter) {
        int n = 60 + static_cast<int>(rng() % 141);
        std::string s;
        for (int i = 0; i < n; ++i) s.push_back(static_cast<char>('a' + rng() % 5));
        SuffixAutomaton sam(s);
        assert(sam.distinctCount() == bfDistinct(s));

        std::vector<std::string> v = allSubstrings(s);
        std::sort(v.begin(), v.end());
        v.erase(std::unique(v.begin(), v.end()), v.end());
        for (int q = 0; q < 15; ++q) {
            long long k = 1 + static_cast<long long>(rng() % v.size());
            assert(sam.kthSubstring(k) == v[static_cast<std::size_t>(k - 1)]);
        }
        for (int q = 0; q < 15; ++q) {
            int i = static_cast<int>(rng() % s.size());
            int j = i + 1 + static_cast<int>(rng() % static_cast<unsigned>(s.size() - i));
            std::string t = s.substr(i, static_cast<std::size_t>(j - i));
            assert(sam.countOccurrences(t) == bfOccurrences(s, t));
            assert(sam.isSubstring(t));
        }
        std::pair<int, std::string> lcs = sam.longestCommonSubstring(s);
        assert(lcs.first == static_cast<int>(s.size()) && lcs.second == s);
    }

    // ---- 邊界：單字元重複、全相異、長度上限 ----
    assert(SuffixAutomaton(std::string(300, 'a')).distinctCount() == 300);
    std::string allDiff;
    for (int i = 1; i <= 200; ++i) allDiff.push_back(static_cast<char>(i));
    assert(SuffixAutomaton(allDiff).distinctCount() == 200LL * 201 / 2);

    // 長度上限：3000 只跑 SAM 本身（暴力要列舉 O(n^2) 個子串，3000 會吃掉數 GB），
    // 暴力對拍改用 400：子串總量約 8×10^4，兩個語言都在一秒內結束。
    std::string big;
    for (std::size_t i = 0; i < MAX_LEN; ++i) big.push_back((rng() % 2) ? 'a' : 'b');
    long long bigCnt = SuffixAutomaton(big).distinctCount();
    long long cap = static_cast<long long>(MAX_LEN) * static_cast<long long>(MAX_LEN + 1) / 2;
    assert(bigCnt > 1 && bigCnt <= cap);
    std::string mid;
    for (int i = 0; i < 400; ++i) mid.push_back(static_cast<char>('a' + rng() % 3));
    assert(SuffixAutomaton(mid).distinctCount() == bfDistinct(mid));

    // ---- 結構性質 ----
    for (int iter = 0; iter < 200; ++iter) {
        int n = 1 + static_cast<int>(rng() % 40);
        std::string s;
        for (int i = 0; i < n; ++i) s.push_back((rng() % 2) ? 'a' : 'b');
        SuffixAutomaton sam(s);
        assert(sam.stateCount() <= 2 * static_cast<std::size_t>(s.size()));
        long long rootCnt = 0;
        for (std::size_t v = 0; v < sam.stateCount(); ++v) {
            if (sam.link()[v] == -1) {
                rootCnt += sam.cnt(static_cast<int>(v));
                continue;
            }
            assert(sam.link()[v] >= 0 && sam.link()[v] < static_cast<int>(sam.stateCount()));
            assert(sam.length()[sam.link()[v]] < sam.length()[v]);
        }
        assert(rootCnt == static_cast<long long>(s.size()));
    }
}

void runIo(const std::string& raw) {
    std::vector<std::string> tok = splitWs(raw);
    std::string s = tok.size() >= 1 ? tok[0] : std::string();
    std::string t = tok.size() >= 2 ? tok[1] : std::string();
    if (s.size() > MAX_LEN) s = s.substr(0, MAX_LEN);
    if (t.size() > MAX_LEN) t = t.substr(0, MAX_LEN);
    long long k = 0;
    if (tok.size() >= 3 && isNumToken(tok[2])) k = parseLLClamp(tok[2]);

    std::string out = run(s, t, k);
    out.push_back('\n');
    std::cout << out;
}

}  // namespace

int main() {
#ifdef _WIN32
    // Windows 文字模式會把 \n 轉成 \r\n，與 Python 的二進位輸出對不起來
    _setmode(_fileno(stdout), _O_BINARY);
#endif
    std::string raw((std::istreambuf_iterator<char>(std::cin)),
                    std::istreambuf_iterator<char>());
    bool blank = true;
    for (unsigned char c : raw)
        if (!isSpaceCh(c)) { blank = false; break; }
    if (!blank) {
        runIo(raw);
    } else {
        runTests();
        std::cout << "all tests passed" << std::endl;
    }
    return 0;
}
