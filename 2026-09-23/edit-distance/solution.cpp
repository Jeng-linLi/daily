// 編輯距離（Edit Distance / Levenshtein Distance，動態規劃）
// 編譯：g++ -std=c++17 -O2 solution.cpp -o solution && ./solution
//
// 思路：dp[i][j] = 把 a 的前 i 個字符變成 b 的前 j 個字符的最少操作數。
//   看最後一個字符，只有三種「最後一步」：
//     刪掉 a[i-1]        -> dp[i-1][j] + 1
//     插入 b[j-1]        -> dp[i][j-1] + 1
//     把 a[i-1] 改/保留  -> dp[i-1][j-1] + (a[i-1] != b[j-1])
//   三者取最小。邊界：dp[0][j] = j（全插入）、dp[i][0] = i（全刪除）。
//
//   空間壓一維：掃描第 i 行時 dp[j] 是上一行的 dp[i-1][j]、dp[j-1] 是本行剛算好的
//   dp[i][j-1]，而 dp[i-1][j-1] 已被覆蓋，所以用 prev_diag 隨身帶着左上角往前滾。
//
//   還原操作序列要保留二維表，從 dp[m][n] 往回走，每步挑一個能解釋當前 dp 值的前驅。
//   回溯是**從後往前**生成操作的，下標天然遞減：按生成順序依次施加時，每次改動只影響
//   下標 >= 當前下標的字符，更靠後的已處理完、更靠前的還沒動，因此下標始終指的是
//   「施加這一操作時字符串裏的位置」，這套順序是自洽的。
//
// 輸入：第一行字符串 a；第二行字符串 b（可爲空行）
// 輸出：第一行最少操作次數；接下來每行一條操作
//       replace <下標> <字符> / delete <下標> / insert <下標> <字符>
// 無 stdin 輸入時運行內置斷言測試。
#include <algorithm>
#include <cassert>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

using namespace std;

// 一維滾動數組版，只求最少操作次數。時間 O(m*n)，空間 O(n)
int editDistance(const string& aIn, const string& bIn) {
    // 讓 b 成爲較短的那個，滾動數組更省空間
    const string *pa = &aIn, *pb = &bIn;
    if (pa->size() < pb->size()) swap(pa, pb);
    const string& a = *pa;
    const string& b = *pb;
    int m = static_cast<int>(a.size()), n = static_cast<int>(b.size());

    vector<int> dp(n + 1);
    for (int j = 0; j <= n; ++j) dp[j] = j;   // dp[0][j] = j：空串變 b 的前 j 個字符
    for (int i = 1; i <= m; ++i) {
        int prevDiag = dp[0];                 // 上一行的 dp[i-1][0]，即左上角
        dp[0] = i;                            // dp[i][0] = i：全刪除
        for (int j = 1; j <= n; ++j) {
            int tmp = dp[j];                  // 更新前是 dp[i-1][j]，之後交給 prevDiag
            int cost = (a[i - 1] == b[j - 1]) ? 0 : 1;
            int val = tmp + 1;                          // 刪除 a[i-1]
            if (dp[j - 1] + 1 < val) val = dp[j - 1] + 1;   // 插入 b[j-1]
            if (prevDiag + cost < val) val = prevDiag + cost;  // 替換或保持
            dp[j] = val;
            prevDiag = tmp;
        }
    }
    return dp[n];
}

// 操作：type 爲 "replace" / "delete" / "insert"；idx 爲施加時的下標；ch 僅插入/替換用
struct Op {
    string type;
    int idx;
    char ch;
};

// 保留二維表並回溯出一條操作序列。時間 O(m*n)，空間 O(m*n)
pair<int, vector<Op>> editDistanceWithOps(const string& a, const string& b) {
    int m = static_cast<int>(a.size()), n = static_cast<int>(b.size());
    vector<vector<int>> dp(m + 1, vector<int>(n + 1, 0));
    for (int i = 0; i <= m; ++i) dp[i][0] = i;
    for (int j = 0; j <= n; ++j) dp[0][j] = j;

    for (int i = 1; i <= m; ++i) {
        for (int j = 1; j <= n; ++j) {
            int cost = (a[i - 1] == b[j - 1]) ? 0 : 1;
            int val = dp[i - 1][j] + 1;                        // 刪除
            if (dp[i][j - 1] + 1 < val) val = dp[i][j - 1] + 1;   // 插入
            if (dp[i - 1][j - 1] + cost < val) val = dp[i - 1][j - 1] + cost;  // 替換/保持
            dp[i][j] = val;
        }
    }

    vector<Op> ops;
    int i = m, j = n;
    while (i > 0 || j > 0) {
        if (i > 0 && j > 0 && a[i - 1] == b[j - 1] && dp[i][j] == dp[i - 1][j - 1]) {
            --i; --j;                                       // 字符相同，免費保留
        } else if (i > 0 && j > 0 && dp[i][j] == dp[i - 1][j - 1] + 1) {
            ops.push_back({"replace", i - 1, b[j - 1]});    // 把 a[i-1] 改成 b[j-1]
            --i; --j;
        } else if (j > 0 && dp[i][j] == dp[i][j - 1] + 1) {
            ops.push_back({"insert", i, b[j - 1]});         // 在位置 i 之前插入 b[j-1]
            --j;
        } else if (i > 0 && dp[i][j] == dp[i - 1][j] + 1) {
            ops.push_back({"delete", i - 1, '\0'});         // 刪掉位置 i-1
            --i;
        } else {
            assert(false && "backtrace stuck");             // 理論上不可達
        }
    }

    return {dp[m][n], ops};
}

// 按序施加操作，用於驗證回溯出來的方案是否真的能把 a 變成 b
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

// 與 Python 版一致的統一輸出文本
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

// 對照用的帶記憶化遞歸（指數級搜索 + 剪枝），僅用於小規模測試驗證
int editDistanceBrute(const string& a, const string& b) {
    int m = static_cast<int>(a.size()), n = static_cast<int>(b.size());
    vector<vector<int>> memo(m + 1, vector<int>(n + 1, -1));

    // 遞歸 lambda：go(i, j) = a 的前 i 個變 b 的前 j 個的最少操作數
    auto go = [&](auto&& self, int i, int j) -> int {
        if (i == 0) return j;
        if (j == 0) return i;
        int& res = memo[i][j];
        if (res != -1) return res;
        int cost = (a[i - 1] == b[j - 1]) ? 0 : 1;
        int val = self(self, i - 1, j) + 1;                              // 刪除
        int ins = self(self, i, j - 1) + 1;                              // 插入
        int rep = self(self, i - 1, j - 1) + cost;                       // 替換/保持
        if (ins < val) val = ins;
        if (rep < val) val = rep;
        return res = val;
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
        assert(static_cast<int>(res.second.size()) == 3);   // 操作條數確實等於最優值
        assert(applyOps("horse", res.second) == "ros");     // 照着做一遍真的得到 ros
        assert(editDistance("horse", "ros") == 3);
        assert(editDistanceBrute("horse", "ros") == 3);
    }

    // 經典用例：intention -> execution，最少 5 步
    assert(editDistance("intention", "execution") == 5);
    assert(editDistanceWithOps("intention", "execution").first == 5);

    // 完全相同：0 步，不產生任何操作
    assert(editDistance("abc", "abc") == 0);
    assert(editDistanceWithOps("abc", "abc").first == 0);
    assert(editDistanceWithOps("abc", "abc").second.empty());
    assert(editDistance("", "") == 0);

    // 一邊爲空：只能全插 / 全刪
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

    // 只差一個字符：1 步替換
    assert(editDistance("kitten", "sitten") == 1);
    // kitten -> sitting：3 步（k->s、e->i、末尾插入 g）
    assert(editDistance("kitten", "sitting") == 3);

    // 大小寫敏感；長度差很大時退化爲大量插入
    assert(editDistance("Ab", "ab") == 1);
    assert(editDistance("a", "aaaa") == 3);

    // 純插入（a 是 b 的子序列）
    assert(editDistance("abc", "axbyc") == 2);

    // 與記憶化暴力解隨機對拍：同時校驗次數一致、操作條數最優、照着做一遍確實得到 b
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
        assert(d1 == res.first);            // 一維版 = 二維版
        assert(res.first == d3);            // 二維版 = 暴力版
        assert(static_cast<int>(res.second.size()) == res.first);  // 操作條數就是最優值
        assert(applyOps(x, res.second) == y);                      // 照着做一遍確實得到 y
    }

    // 對稱性（編輯距離的基本性質）
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
