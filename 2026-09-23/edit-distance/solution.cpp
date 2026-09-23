// 编辑距离（Edit Distance / Levenshtein Distance，动态规划）
// 编译：g++ -std=c++17 -O2 solution.cpp -o solution && ./solution
//
// 思路：dp[i][j] = 把 a 的前 i 个字符变成 b 的前 j 个字符的最少操作数。
//   看最后一个字符，只有三种「最后一步」：
//     删掉 a[i-1]        -> dp[i-1][j] + 1
//     插入 b[j-1]        -> dp[i][j-1] + 1
//     把 a[i-1] 改/保留  -> dp[i-1][j-1] + (a[i-1] != b[j-1])
//   三者取最小。边界：dp[0][j] = j（全插入）、dp[i][0] = i（全删除）。
//
//   空间压一维：扫描第 i 行时 dp[j] 是上一行的 dp[i-1][j]、dp[j-1] 是本行刚算好的
//   dp[i][j-1]，而 dp[i-1][j-1] 已被覆盖，所以用 prev_diag 随身带着左上角往前滚。
//
//   还原操作序列要保留二维表，从 dp[m][n] 往回走，每步挑一个能解释当前 dp 值的前驱。
//   回溯是**从后往前**生成操作的，下标天然递减：按生成顺序依次施加时，每次改动只影响
//   下标 >= 当前下标的字符，更靠后的已处理完、更靠前的还没动，因此下标始终指的是
//   「施加这一操作时字符串里的位置」，这套顺序是自洽的。
//
// 输入：第一行字符串 a；第二行字符串 b（可为空行）
// 输出：第一行最少操作次数；接下来每行一条操作
//       replace <下标> <字符> / delete <下标> / insert <下标> <字符>
// 无 stdin 输入时运行内置断言测试。
#include <algorithm>
#include <cassert>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

using namespace std;

// 一维滚动数组版，只求最少操作次数。时间 O(m*n)，空间 O(n)
int editDistance(const string& aIn, const string& bIn) {
    // 让 b 成为较短的那个，滚动数组更省空间
    const string *pa = &aIn, *pb = &bIn;
    if (pa->size() < pb->size()) swap(pa, pb);
    const string& a = *pa;
    const string& b = *pb;
    int m = static_cast<int>(a.size()), n = static_cast<int>(b.size());

    vector<int> dp(n + 1);
    for (int j = 0; j <= n; ++j) dp[j] = j;   // dp[0][j] = j：空串变 b 的前 j 个字符
    for (int i = 1; i <= m; ++i) {
        int prevDiag = dp[0];                 // 上一行的 dp[i-1][0]，即左上角
        dp[0] = i;                            // dp[i][0] = i：全删除
        for (int j = 1; j <= n; ++j) {
            int tmp = dp[j];                  // 更新前是 dp[i-1][j]，之后交给 prevDiag
            int cost = (a[i - 1] == b[j - 1]) ? 0 : 1;
            int val = tmp + 1;                          // 删除 a[i-1]
            if (dp[j - 1] + 1 < val) val = dp[j - 1] + 1;   // 插入 b[j-1]
            if (prevDiag + cost < val) val = prevDiag + cost;  // 替换或保持
            dp[j] = val;
            prevDiag = tmp;
        }
    }
    return dp[n];
}

// 操作：type 为 "replace" / "delete" / "insert"；idx 为施加时的下标；ch 仅插入/替换用
struct Op {
    string type;
    int idx;
    char ch;
};

// 保留二维表并回溯出一条操作序列。时间 O(m*n)，空间 O(m*n)
pair<int, vector<Op>> editDistanceWithOps(const string& a, const string& b) {
    int m = static_cast<int>(a.size()), n = static_cast<int>(b.size());
    vector<vector<int>> dp(m + 1, vector<int>(n + 1, 0));
    for (int i = 0; i <= m; ++i) dp[i][0] = i;
    for (int j = 0; j <= n; ++j) dp[0][j] = j;

    for (int i = 1; i <= m; ++i) {
        for (int j = 1; j <= n; ++j) {
            int cost = (a[i - 1] == b[j - 1]) ? 0 : 1;
            int val = dp[i - 1][j] + 1;                        // 删除
            if (dp[i][j - 1] + 1 < val) val = dp[i][j - 1] + 1;   // 插入
            if (dp[i - 1][j - 1] + cost < val) val = dp[i - 1][j - 1] + cost;  // 替换/保持
            dp[i][j] = val;
        }
    }

    vector<Op> ops;
    int i = m, j = n;
    while (i > 0 || j > 0) {
        if (i > 0 && j > 0 && a[i - 1] == b[j - 1] && dp[i][j] == dp[i - 1][j - 1]) {
            --i; --j;                                       // 字符相同，免费保留
        } else if (i > 0 && j > 0 && dp[i][j] == dp[i - 1][j - 1] + 1) {
            ops.push_back({"replace", i - 1, b[j - 1]});    // 把 a[i-1] 改成 b[j-1]
            --i; --j;
        } else if (j > 0 && dp[i][j] == dp[i][j - 1] + 1) {
            ops.push_back({"insert", i, b[j - 1]});         // 在位置 i 之前插入 b[j-1]
            --j;
        } else if (i > 0 && dp[i][j] == dp[i - 1][j] + 1) {
            ops.push_back({"delete", i - 1, '\0'});         // 删掉位置 i-1
            --i;
        } else {
            assert(false && "backtrace stuck");             // 理论上不可达
        }
    }

    return {dp[m][n], ops};
}

// 按序施加操作，用于验证回溯出来的方案是否真的能把 a 变成 b
string applyOps(const string& a, const vector<Op>& ops) {
    string cur = a;
    for (const Op& op : ops) {
        if (op.type == "replace") {
            cur[op.idx] = op.ch;
        } else if (op.type == "delete") {
            cur.erase(op.idx, 1);
        } else if (op.type == "insert") {
            cur.insert(op.idx, 1, op.ch);   // 插到位置 idx 之前
        } else {
            assert(false && "unknown op");
        }
    }
    return cur;
}

// 与 Python 版一致的统一输出文本
vector<string> formatOps(const vector<Op>& ops) {
    vector<string> lines;
    for (const Op& op : ops) {
        if (op.type == "replace" || op.type == "insert") {
            lines.push_back(op.type + " " + to_string(op.idx) + " " + string(1, op.ch));
        } else {
            lines.push_back(op.type + " " + to_string(op.idx));
        }
    }
    return lines;
}

// 对照用的带记忆化递归（指数级搜索 + 剪枝），仅用于小规模测试验证
int editDistanceBrute(const string& a, const string& b) {
    int m = static_cast<int>(a.size()), n = static_cast<int>(b.size());
    vector<vector<int>> memo(m + 1, vector<int>(n + 1, -1));

    // 递归 lambda：go(i, j) = a 的前 i 个变 b 的前 j 个的最少操作数
    auto go = [&](auto&& self, int i, int j) -> int {
        if (i == 0) return j;
        if (j == 0) return i;
        int& res = memo[i][j];
        if (res != -1) return res;
        int cost = (a[i - 1] == b[j - 1]) ? 0 : 1;
        int val = self(self, i - 1, j) + 1;                              // 删除
        int ins = self(self, i, j - 1) + 1;                              // 插入
        int rep = self(self, i - 1, j - 1) + cost;                       // 替换/保持
        if (ins < val) val = ins;
        if (rep < val) val = rep;
        return res = val;
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
    string a, b;
    if (getline(cin, a)) {  // IO 模式
        if (!a.empty() && a.back() == '\r') a.pop_back();
        if (getline(cin, b)) {
            if (!b.empty() && b.back() == '\r') b.pop_back();
        } else {
            b.clear();
        }
        auto res = editDistanceWithOps(a, b);
        cout << res.first << "\n";
        for (const string& line : formatOps(res.second)) cout << line << "\n";
        return 0;
    }

    // README 中的示例：horse -> ros，最少 3 步
    {
        auto res = editDistanceWithOps("horse", "ros");
        assert(res.first == 3);
        assert(static_cast<int>(res.second.size()) == 3);   // 操作条数确实等于最优值
        assert(applyOps("horse", res.second) == "ros");     // 照着做一遍真的得到 ros
        assert(editDistance("horse", "ros") == 3);
        assert(editDistanceBrute("horse", "ros") == 3);
    }

    // 经典用例：intention -> execution，最少 5 步
    assert(editDistance("intention", "execution") == 5);
    assert(editDistanceWithOps("intention", "execution").first == 5);

    // 完全相同：0 步，不产生任何操作
    assert(editDistance("abc", "abc") == 0);
    assert(editDistanceWithOps("abc", "abc").first == 0);
    assert(editDistanceWithOps("abc", "abc").second.empty());
    assert(editDistance("", "") == 0);

    // 一边为空：只能全插 / 全删
    assert(editDistance("", "abc") == 3);
    assert(editDistance("abc", "") == 3);
    {
        auto r1 = editDistanceWithOps("", "abc");
        auto r2 = editDistanceWithOps("abc", "");
        assert(r1.first == 3);
        assert(applyOps("", r1.second) == "abc");
        assert(r2.first == 3);
        assert(applyOps("abc", r2.second) == "");
    }

    // 只差一个字符：1 步替换
    assert(editDistance("kitten", "sitten") == 1);
    // kitten -> sitting：3 步（k->s、e->i、末尾插入 g）
    assert(editDistance("kitten", "sitting") == 3);

    // 大小写敏感；长度差很大时退化为大量插入
    assert(editDistance("Ab", "ab") == 1);
    assert(editDistance("a", "aaaa") == 3);

    // 纯插入（a 是 b 的子序列）
    assert(editDistance("abc", "axbyc") == 2);

    // 与记忆化暴力解随机对拍：同时校验次数一致、操作条数最优、照着做一遍确实得到 b
    LCG rng(20260923ULL);
    const string alphabet = "abc";
    for (int t = 0; t < 200; ++t) {
        string x, y;
        int lx = rng.next(0, 7), ly = rng.next(0, 7);
        for (int i = 0; i < lx; ++i) x.push_back(alphabet[rng.next(0, 2)]);
        for (int i = 0; i < ly; ++i) y.push_back(alphabet[rng.next(0, 2)]);

        int d1 = editDistance(x, y);
        auto res = editDistanceWithOps(x, y);
        int d3 = editDistanceBrute(x, y);
        assert(d1 == res.first);            // 一维版 = 二维版
        assert(res.first == d3);            // 二维版 = 暴力版
        assert(static_cast<int>(res.second.size()) == res.first);  // 操作条数就是最优值
        assert(applyOps(x, res.second) == y);                      // 照着做一遍确实得到 y
    }

    // 对称性（编辑距离的基本性质）
    assert(editDistance("flaw", "lawn") == editDistance("lawn", "flaw"));
    for (int t = 0; t < 50; ++t) {
        string x, y;
        int lx = rng.next(0, 6), ly = rng.next(0, 6);
        for (int i = 0; i < lx; ++i) x.push_back(alphabet[rng.next(0, 2)]);
        for (int i = 0; i < ly; ++i) y.push_back(alphabet[rng.next(0, 2)]);
        assert(editDistance(x, y) == editDistance(y, x));
    }

    cout << "all tests passed" << endl;
    return 0;
}
