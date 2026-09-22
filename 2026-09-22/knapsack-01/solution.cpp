// 0-1 背包（Knapsack 0-1，动态规划）
// 编译：g++ -std=c++17 -O2 solution.cpp -o solution && ./solution
//
// 思路：定义 dp[i][c] = 只考虑前 i 个物品、容量为 c 时的最大价值。对第 i 个物品：
//     不选 dp[i][c] = dp[i-1][c]；选 dp[i][c] = dp[i-1][c-w[i]] + v[i]（需 c >= w[i]），取最大。
//   物品只能选一次，转移只依赖 i-1 层，不会同层自我叠加（区别于完全背包）。
//
//   空间可压到一维 dp[c]，但**容量 c 必须倒序遍历**（W -> w[i]）：
//   倒序保证 dp[c-w[i]] 读到的仍是 i-1 层旧值；正序则会把同一物品重复放入，退化成完全背包。
//   要还原方案必须保留二维表：dp[i][c] != dp[i-1][c] 说明第 i 个物品被选了，反推即可。
//
// 输入：第一行 n W；接下来 n 行，每行 w_i v_i
// 输出：第一行最大总价值；第二行被选中物品下标（0-based，升序，空格分隔）
// 无 stdin 输入时运行内置断言测试。
#include <algorithm>
#include <cassert>
#include <iostream>
#include <utility>
#include <vector>

using namespace std;

// 只求最大价值：一维滚动数组，时间 O(n*W)，空间 O(W)
int knapsackMax(const vector<int>& w, const vector<int>& v, int capacity) {
    vector<int> dp(capacity + 1, 0);
    for (size_t i = 0; i < w.size(); ++i) {
        // 容量倒序；正序会让同一物品被重复选取，退化为完全背包
        for (int c = capacity; c >= w[i]; --c) {
            int cand = dp[c - w[i]] + v[i];
            if (cand > dp[c]) dp[c] = cand;
        }
    }
    return dp[capacity];
}

// 求最大价值并还原一种选取方案：保留二维表，空间 O(n*W)
pair<int, vector<int>> knapsackWithItems(const vector<int>& w, const vector<int>& v, int capacity) {
    int n = static_cast<int>(w.size());
    vector<vector<int>> dp(n + 1, vector<int>(capacity + 1, 0));

    for (int i = 1; i <= n; ++i) {
        for (int c = 0; c <= capacity; ++c) {
            int best = dp[i - 1][c];                       // 不选第 i-1 个物品
            if (c >= w[i - 1]) {
                int take = dp[i - 1][c - w[i - 1]] + v[i - 1];  // 选它
                if (take > best) best = take;
            }
            dp[i][c] = best;
        }
    }

    // 反向还原：dp[i][c] 比 dp[i-1][c] 大，说明第 i-1 个物品被选中了
    vector<int> chosen;
    int c = capacity;
    for (int i = n; i >= 1; --i) {
        if (dp[i][c] != dp[i - 1][c]) {
            chosen.push_back(i - 1);
            c -= w[i - 1];
        }
    }
    reverse(chosen.begin(), chosen.end());
    return {dp[n][capacity], chosen};
}

// 对照用的指数级枚举，仅用于小规模测试验证
int bruteForce(const vector<int>& w, const vector<int>& v, int capacity) {
    int n = static_cast<int>(w.size());
    int best = 0;
    for (int mask = 0; mask < (1 << n); ++mask) {
        int tw = 0, tv = 0;
        for (int i = 0; i < n; ++i) {
            if (mask >> i & 1) { tw += w[i]; tv += v[i]; }
        }
        if (tw <= capacity && tv > best) best = tv;
    }
    return best;
}

// 与 Python 版完全一致的固定随机序列，保证两版跑同一批用例
struct LCG {
    unsigned long long s;
    LCG(unsigned long long seed) : s(seed) {}
    int next(int lo, int hi) {
        s = s * 6364136223846793005ULL + 1442695040888963407ULL;
        return lo + static_cast<int>((s >> 33) % static_cast<unsigned long long>(hi - lo + 1));
    }
};

int main() {
    int n, W;
    if (cin >> n >> W) {  // IO 模式
        vector<int> w(n), v(n);
        for (int i = 0; i < n; ++i) cin >> w[i] >> v[i];
        auto res = knapsackWithItems(w, v, W);
        vector<int>& chosen = res.second;
        int total = 0;
        for (int i : chosen) total += v[i];
        cout << total << "\n";
        for (size_t i = 0; i < chosen.size(); ++i) {
            if (i) cout << " ";
            cout << chosen[i];
        }
        cout << "\n";
        return 0;
    }

    // 经典用例：容量 10 的最优值为 12，但存在多解
    //   (a) 2 号 + 3 号：w=4+6=10, v=5+7=12
    //   (b) 0 号 + 1 号 + 2 号：w=2+3+4=9, v=3+4+5=12
    // 反推时「dp[i][c] == dp[i-1][c] 视为未选」，会优先得到 (b)，故只断言值最优、方案合法
    {
        vector<int> w = {2, 3, 4, 6};
        vector<int> v = {3, 4, 5, 7};
        auto res = knapsackWithItems(w, v, 10);
        assert(res.first == 12);
        assert(knapsackMax(w, v, 10) == 12);
        int tw = 0, tv = 0;
        for (int i : res.second) { tw += w[i]; tv += v[i]; }
        assert(tw <= 10);
        assert(tv == 12);
    }

    // 容量为 0 / 物品为空
    assert(knapsackMax({}, {}, 10) == 0);
    assert(knapsackMax({5}, {9}, 0) == 0);
    assert(knapsackWithItems({5}, {9}, 0).first == 0);
    assert(knapsackWithItems({5}, {9}, 0).second.empty());

    // 单件装不下
    assert(knapsackMax({5}, {9}, 4) == 0);

    // 全部都能装下
    assert(knapsackMax({1, 2, 3}, {1, 2, 3}, 10) == 6);

    // 按价值密度贪心会选错：密度最高的是 0 号(2.0)，但最优解是 1+2 号 w=10 v=18
    assert(knapsackMax({5, 4, 6}, {10, 7, 11}, 10) == 18);
    assert(knapsackMax({5, 4, 6}, {10, 7, 11}, 10) == bruteForce({5, 4, 6}, {10, 7, 11}, 10));

    // 零重量物品：价值白拿，且不会造成死循环
    assert(knapsackMax({0, 3}, {5, 4}, 3) == 9);

    // 与暴力枚举随机对照：同时校验最优值一致、方案可行且达到最优值
    LCG rng(20260922ULL);
    for (int t = 0; t < 300; ++t) {
        int len = rng.next(1, 10);
        int cap = rng.next(0, 20);
        vector<int> ws(len), vs(len);
        for (int i = 0; i < len; ++i) {
            ws[i] = rng.next(0, 8);
            vs[i] = rng.next(0, 20);
        }
        auto res = knapsackWithItems(ws, vs, cap);
        assert(res.first == bruteForce(ws, vs, cap));   // 最优值与暴力一致
        assert(knapsackMax(ws, vs, cap) == res.first);  // 一维版与二维版一致
        int tw = 0, tv = 0;
        for (int i : res.second) { tw += ws[i]; tv += vs[i]; }
        assert(tw <= cap);                              // 方案不超容量
        assert(tv == res.first);                        // 方案确实达到最优值
    }

    cout << "all tests passed" << endl;
    return 0;
}
