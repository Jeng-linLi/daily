// 最长公共子序列（LCS，动态规划 + 回溯还原）
// 编译：g++ -std=c++17 -O2 solution.cpp -o solution && ./solution
//
// 思路：dp[i][j] = a 的前 i 个字符与 b 的前 j 个字符的 LCS 长度。
//   看最后一对字符 a[i-1] 与 b[j-1]，只有两种情况：
//     - 相等：这个字符一定可以接在 a[:i-1] 与 b[:j-1] 的 LCS 后面，
//             故 dp[i][j] = dp[i-1][j-1] + 1；
//     - 不等：它俩不可能同时出现在同一个匹配里，dp[i][j] = max(dp[i-1][j], dp[i][j-1])。
//   边界 dp[0][j] = dp[i][0] = 0。
//
//   只求长度时空间可压到两行（dp[i][*] 只依赖 dp[i-1][*]）；但要还原方案必须保留
//   整张表，再从 dp[m][n] 往回走：字符相等就收下并同时后退一步，否则往 dp 值大的
//   方向走（相等时优先走 i，即丢弃 a[i-1]）。回溯倒着走，收集到的字符最后要反转。
//
//   注意 LCS 通常不唯一，回溯只保证给出其中一条。
//
// 输入：第一行字符串 a；第二行字符串 b（可为空行）
// 输出：第一行 LCS 长度；第二行一条达到该长度的公共子序列（长度为 0 时输出空行）
// 无 stdin 输入时运行内置断言测试。
#include <algorithm>
#include <cassert>
#include <cctype>
#include <iostream>
#include <iterator>
#include <string>
#include <utility>
#include <vector>

using namespace std;

// 两行滚动数组，只求长度。时间 O(m*n)，空间 O(min(m, n))
int lcsLength(const string& aIn, const string& bIn) {
    const string *pa = &aIn, *pb = &bIn;
    if (pa->size() < pb->size()) swap(pa, pb);   // 让 b 成为较短的那个，滚动数组更省空间
    const string& a = *pa;
    const string& b = *pb;
    int m = static_cast<int>(a.size()), n = static_cast<int>(b.size());

    vector<int> prev(n + 1, 0), cur(n + 1, 0);
    for (int i = 1; i <= m; ++i) {
        for (int j = 1; j <= n; ++j) {
            if (a[i - 1] == b[j - 1]) cur[j] = prev[j - 1] + 1;
            else cur[j] = max(prev[j], cur[j - 1]);
        }
        prev.swap(cur);                          // 交换，下一行复用上一行的空间
        cur[0] = 0;
    }
    return prev[n];
}

// 保留完整 dp 表并回溯出一条 LCS。时间 O(m*n)，空间 O(m*n)
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
            sub.push_back(a[i - 1]);             // 这个字符属于 LCS
            --i; --j;
        } else if (dp[i - 1][j] >= dp[i][j - 1]) {
            --i;                                 // 丢弃 a[i-1]（相等时优先走 i）
        } else {
            --j;                                 // 丢弃 b[j-1]
        }
    }
    reverse(sub.begin(), sub.end());
    return {dp[m][n], sub};
}

// 判断 sub 是否为 s 的子序列
bool isSubsequence(const string& sub, const string& s) {
    size_t pos = 0;
    for (char ch : sub) {
        pos = s.find(ch, pos);
        if (pos == string::npos) return false;
        ++pos;
    }
    return true;
}

// 对照用的带记忆化递归（指数级搜索 + 剪枝），仅用于小规模测试验证
int lcsBrute(const string& a, const string& b) {
    int m = static_cast<int>(a.size()), n = static_cast<int>(b.size());
    vector<vector<int>> memo(m + 1, vector<int>(n + 1, -1));

    // 递归 lambda：go(i, j) = a 的前 i 个与 b 的前 j 个的 LCS 长度
    auto go = [&](auto&& self, int i, int j) -> int {
        if (i == 0 || j == 0) return 0;
        int& res = memo[i][j];
        if (res != -1) return res;
        if (a[i - 1] == b[j - 1]) return res = self(self, i - 1, j - 1) + 1;
        return res = max(self(self, i - 1, j), self(self, i, j - 1));
    };

    return go(go, m, n);
}

// 与 Python 版同规模的固定随机序列（LCG），两版各自独立与暴力解对拍
struct LCG {
    unsigned long long s;
    LCG(unsigned long long seed) : s(seed) {}
    int next(int lo, int hi) {
        s = s * 6364136223846793005ULL + 1442695040888963407ULL;
        return lo + static_cast<int>((s >> 33) % static_cast<unsigned long long>(hi - lo + 1));
    }
};

int main() {
    // 一次性读完整份 stdin，与 Python 版 `raw.strip()` 的判定保持完全一致：
    // 只有存在非空白字符时才进入 IO 模式
    string data((istreambuf_iterator<char>(cin)), istreambuf_iterator<char>());
    bool hasContent = false;
    for (char c : data) {
        if (!isspace(static_cast<unsigned char>(c))) { hasContent = true; break; }
    }
    if (hasContent) {
        vector<string> lines;
        string line;
        for (size_t p = 0; p <= data.size();) {          // 手工按 \n 切分，与 splitlines 对齐
            size_t q = data.find('\n', p);
            if (q == string::npos) {
                line = data.substr(p);
                if (!line.empty()) lines.push_back(line);
                break;
            }
            line = data.substr(p, q - p);
            if (!line.empty() && line.back() == '\r') line.pop_back();  // 兼容 CRLF
            lines.push_back(line);                        // 空行也要保留，与 Python splitlines 一致
            p = q + 1;
        }
        string a = lines.size() > 0 ? lines[0] : "";
        string b = lines.size() > 1 ? lines[1] : "";
        auto res = lcsWithString(a, b);
        cout << res.first << "\n" << res.second << "\n";
        return 0;
    }

    // README 示例：abcde 与 ace 的 LCS 是 ace，长度 3
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

    // 没有公共字符：空串
    assert(lcsLength("abc", "def") == 0);
    assert(lcsWithString("abc", "def") == make_pair(0, string("")));

    // 一边为空
    assert(lcsLength("", "abc") == 0);
    assert(lcsWithString("", "abc") == make_pair(0, string("")));
    assert(lcsWithString("abc", "") == make_pair(0, string("")));
    assert(lcsWithString("", "") == make_pair(0, string("")));

    // 经典用例：长度为 4（"bdab" / "bcba" 等都算对，只断言长度与合法性）
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

    // 大小写敏感
    assert(lcsLength("Abc", "abc") == 2);

    // 重复字符
    assert(lcsLength("aaaa", "aa") == 2);
    assert(lcsLength("aab", "aba") == 2);

    LCG rng(20260927ULL);
    const string alphabet = "abc";

    // 随机对拍：两行滚动版 = 完整表版 = 记忆化暴力版，且还原出的串确实合法
    for (int t = 0; t < 300; ++t) {
        string a, b;
        int la = rng.next(0, 8), lb = rng.next(0, 8);
        for (int i = 0; i < la; ++i) a.push_back(alphabet[rng.next(0, 2)]);
        for (int i = 0; i < lb; ++i) b.push_back(alphabet[rng.next(0, 2)]);

        int d1 = lcsLength(a, b);
        auto res = lcsWithString(a, b);
        int d3 = lcsBrute(a, b);
        assert(d1 == res.first);                          // 滚动版 = 完整表版
        assert(res.first == d3);                          // 完整表版 = 暴力版
        assert(static_cast<int>(res.second.size()) == res.first);  // 串长就是最优值
        assert(isSubsequence(res.second, a));             // 确实是 a 的子序列
        assert(isSubsequence(res.second, b));             // 确实是 b 的子序列
        assert(0 <= res.first && res.first <= min(la, lb));  // 长度落在合理区间内
    }

    // 对称性与上界性质
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
