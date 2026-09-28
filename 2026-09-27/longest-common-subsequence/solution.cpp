// 最長公共子序列（LCS，動態規劃 + 回溯還原）
// 編譯：g++ -std=c++17 -O2 solution.cpp -o solution && ./solution
//
// 思路：dp[i][j] = a 的前 i 個字符與 b 的前 j 個字符的 LCS 長度。
//   看最後一對字符 a[i-1] 與 b[j-1]，只有兩種情況：
//     - 相等：這個字符一定可以接在 a[:i-1] 與 b[:j-1] 的 LCS 後面，
//             故 dp[i][j] = dp[i-1][j-1] + 1；
//     - 不等：它倆不可能同時出現在同一個匹配裏，dp[i][j] = max(dp[i-1][j], dp[i][j-1])。
//   邊界 dp[0][j] = dp[i][0] = 0。
//
//   只求長度時空間可壓到兩行（dp[i][*] 只依賴 dp[i-1][*]）；但要還原方案必須保留
//   整張表，再從 dp[m][n] 往回走：字符相等就收下並同時後退一步，否則往 dp 值大的
//   方向走（相等時優先走 i，即丟棄 a[i-1]）。回溯倒着走，收集到的字符最後要反轉。
//
//   注意 LCS 通常不唯一，回溯只保證給出其中一條。
//
// 輸入：第一行字符串 a；第二行字符串 b（可爲空行）
// 輸出：第一行 LCS 長度；第二行一條達到該長度的公共子序列（長度爲 0 時輸出空行）
// 無 stdin 輸入時運行內置斷言測試。
#include <algorithm>
#include <cassert>
#include <cctype>
#include <iostream>
#include <iterator>
#include <string>
#include <utility>
#include <vector>

using namespace std;

// 兩行滾動數組，只求長度。時間 O(m*n)，空間 O(min(m, n))
int lcsLength(const string& aIn, const string& bIn) {
    const string *pa = &aIn, *pb = &bIn;
    if (pa->size() < pb->size()) swap(pa, pb);   // 讓 b 成爲較短的那個，滾動數組更省空間
    const string& a = *pa;
    const string& b = *pb;
    int m = static_cast<int>(a.size()), n = static_cast<int>(b.size());

    vector<int> prev(n + 1, 0), cur(n + 1, 0);
    for (int i = 1; i <= m; ++i) {
        for (int j = 1; j <= n; ++j) {
            if (a[i - 1] == b[j - 1]) cur[j] = prev[j - 1] + 1;
            else cur[j] = max(prev[j], cur[j - 1]);
        }
        prev.swap(cur);                          // 交換，下一行復用上一行的空間
        cur[0] = 0;
    }
    return prev[n];
}

// 保留完整 dp 表並回溯出一條 LCS。時間 O(m*n)，空間 O(m*n)
pair<int, string> lcsWithString(const string& a, const string& b) {
    int m = static_cast<int>(a.size()), n = static_cast<int>(b.size());
    vector<vector<int>> dp(m + 1, vector<int>(n + 1, 0));
    for (int i = 1; i <= m; ++i) {
        for (int j = 1; j <= n; ++j) {
            if (a[i - 1] == b[j - 1]) dp[i][j] = dp[i - 1][j - 1] + 1;
            else dp[i][j] = max(dp[i - 1][j], dp[i][j - 1]);
        }
    }

    string sub;
    int i = m, j = n;
    while (i > 0 && j > 0) {
        if (a[i - 1] == b[j - 1]) {
            sub.push_back(a[i - 1]);             // 這個字符屬於 LCS
            --i; --j;
        } else if (dp[i - 1][j] >= dp[i][j - 1]) {
            --i;                                 // 丟棄 a[i-1]（相等時優先走 i）
        } else {
            --j;                                 // 丟棄 b[j-1]
        }
    }
    reverse(sub.begin(), sub.end());
    return {dp[m][n], sub};
}

// 判斷 sub 是否爲 s 的子序列
bool isSubsequence(const string& sub, const string& s) {
    size_t pos = 0;
    for (char ch : sub) {
        pos = s.find(ch, pos);
        if (pos == string::npos) return false;
        ++pos;
    }
    return true;
}

// 對照用的帶記憶化遞歸（指數級搜索 + 剪枝），僅用於小規模測試驗證
int lcsBrute(const string& a, const string& b) {
    int m = static_cast<int>(a.size()), n = static_cast<int>(b.size());
    vector<vector<int>> memo(m + 1, vector<int>(n + 1, -1));

    // 遞歸 lambda：go(i, j) = a 的前 i 個與 b 的前 j 個的 LCS 長度
    auto go = [&](auto&& self, int i, int j) -> int {
        if (i == 0 || j == 0) return 0;
        int& res = memo[i][j];
        if (res != -1) return res;
        if (a[i - 1] == b[j - 1]) return res = self(self, i - 1, j - 1) + 1;
        return res = max(self(self, i - 1, j), self(self, i, j - 1));
    };

    return go(go, m, n);
}

// 與 Python 版同規模的固定隨機序列（LCG），兩版各自獨立與暴力解對拍
struct LCG {
    unsigned long long s;
    LCG(unsigned long long seed) : s(seed) {}
    int next(int lo, int hi) {
        s = s * 6364136223846793005ULL + 1442695040888963407ULL;
        return lo + static_cast<int>((s >> 33) % static_cast<unsigned long long>(hi - lo + 1));
    }
};

int main() {
    // 一次性讀完整份 stdin，與 Python 版 `raw.strip()` 的判定保持完全一致：
    // 只有存在非空白字符時才進入 IO 模式
    string data((istreambuf_iterator<char>(cin)), istreambuf_iterator<char>());
    bool hasContent = false;
    for (char c : data) {
        if (!isspace(static_cast<unsigned char>(c))) { hasContent = true; break; }
    }
    if (hasContent) {
        vector<string> lines;
        string line;
        for (size_t p = 0; p <= data.size();) {          // 手工按 \n 切分，與 splitlines 對齊
            size_t q = data.find('\n', p);
            if (q == string::npos) {
                line = data.substr(p);
                if (!line.empty()) lines.push_back(line);
                break;
            }
            line = data.substr(p, q - p);
            if (!line.empty() && line.back() == '\r') line.pop_back();  // 兼容 CRLF
            lines.push_back(line);                        // 空行也要保留，與 Python splitlines 一致
            p = q + 1;
        }
        string a = lines.size() > 0 ? lines[0] : "";
        string b = lines.size() > 1 ? lines[1] : "";
        auto res = lcsWithString(a, b);
        cout << res.first << "\n" << res.second << "\n";
        return 0;
    }

    // README 示例：abcde 與 ace 的 LCS 是 ace，長度 3
    {
        auto res = lcsWithString("abcde", "ace");
        assert(res.first == 3);
        assert(res.second == "ace");
        assert(lcsLength("abcde", "ace") == 3);
        assert(lcsBrute("abcde", "ace") == 3);
    }

    // 完全相同：LCS 就是自身
    assert(lcsLength("abc", "abc") == 3);
    assert(lcsWithString("abc", "abc") == make_pair(3, string("abc")));

    // 沒有公共字符：空串
    assert(lcsLength("abc", "def") == 0);
    assert(lcsWithString("abc", "def") == make_pair(0, string("")));

    // 一邊爲空
    assert(lcsLength("", "abc") == 0);
    assert(lcsWithString("", "abc") == make_pair(0, string("")));
    assert(lcsWithString("abc", "") == make_pair(0, string("")));
    assert(lcsWithString("", "") == make_pair(0, string("")));

    // 經典用例：長度爲 4（"bdab" / "bcba" 等都算對，只斷言長度與合法性）
    {
        auto res = lcsWithString("abcbdab", "bdcaba");
        assert(res.first == 4);
        assert(static_cast<int>(res.second.size()) == 4);
        assert(isSubsequence(res.second, "abcbdab"));
        assert(isSubsequence(res.second, "bdcaba"));
        assert(lcsBrute("abcbdab", "bdcaba") == 4);
    }

    // a 是 b 的子序列：LCS 就是 a
    assert(lcsLength("ace", "abcde") == 3);
    assert(lcsWithString("ace", "abcde").second == "ace");

    // 大小寫敏感
    assert(lcsLength("Abc", "abc") == 2);

    // 重複字符
    assert(lcsLength("aaaa", "aa") == 2);
    assert(lcsLength("aab", "aba") == 2);

    LCG rng(20260927ULL);
    const string alphabet = "abc";

    // 隨機對拍：兩行滾動版 = 完整表版 = 記憶化暴力版，且還原出的串確實合法
    for (int t = 0; t < 300; ++t) {
        string a, b;
        int la = rng.next(0, 8), lb = rng.next(0, 8);
        for (int i = 0; i < la; ++i) a.push_back(alphabet[rng.next(0, 2)]);
        for (int i = 0; i < lb; ++i) b.push_back(alphabet[rng.next(0, 2)]);

        int d1 = lcsLength(a, b);
        auto res = lcsWithString(a, b);
        int d3 = lcsBrute(a, b);
        assert(d1 == res.first);                          // 滾動版 = 完整表版
        assert(res.first == d3);                          // 完整表版 = 暴力版
        assert(static_cast<int>(res.second.size()) == res.first);  // 串長就是最優值
        assert(isSubsequence(res.second, a));             // 確實是 a 的子序列
        assert(isSubsequence(res.second, b));             // 確實是 b 的子序列
        assert(0 <= res.first && res.first <= min(la, lb));  // 長度落在合理區間內
    }

    // 對稱性與上界性質
    for (int t = 0; t < 100; ++t) {
        string a, b;
        int la = rng.next(0, 6), lb = rng.next(0, 6);
        for (int i = 0; i < la; ++i) a.push_back(alphabet[rng.next(0, 2)]);
        for (int i = 0; i < lb; ++i) b.push_back(alphabet[rng.next(0, 2)]);
        assert(lcsLength(a, b) == lcsLength(b, a));
        assert(lcsWithString(a, b).first == lcsWithString(b, a).first);
        assert(0 <= lcsLength(a, b) && lcsLength(a, b) <= min(la, lb));
    }

    cout << "all tests passed" << endl;
    return 0;
}
