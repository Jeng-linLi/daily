// 最长递增子序列（Longest Increasing Subsequence, LIS）
// 编译：g++ -std=c++17 -O2 solution.cpp -o solution && ./solution
//
// 思路：
//   方法一 · 动态规划 O(n^2)
//       dp[i] 表示「以 nums[i] 作为结尾」的最长递增子序列长度。
//       dp[i] = 1 + max{ dp[j] | j < i 且 nums[j] < nums[i] }，取不到则 dp[i] = 1。
//       pre[i] 记前驱下标，从 dp 最大处往回跳即可还原序列。
//   方法二 · 贪心 + 二分 O(n log n)
//       tails[k] = 长度为 k+1 的递增子序列的最小结尾值，tails 本身严格递增。
//       对每个 x 二分找到第一个 >= x 的位置：能接长就追加，否则用 x 覆盖（结尾越小潜力越大）。
//       tails 未必是合法子序列，但 len(tails) 一定是答案。
//
// 输入：第一行 n；第二行 n 个整数
// 输出：第一行 LIS 长度；第二行一条 LIS（空格分隔）
// 无 stdin 输入时运行内置断言测试。
#include <algorithm>
#include <cassert>
#include <iostream>
#include <sstream>
#include <vector>

using namespace std;

// 贪心 + 二分，O(n log n)，只求长度
int lengthOfLIS(const vector<int>& nums) {
    vector<int> tails;  // tails[k] = 长度 k+1 的递增子序列的最小结尾
    for (int x : nums) {
        auto it = lower_bound(tails.begin(), tails.end(), x);  // 第一个 >= x
        if (it == tails.end()) {
            tails.push_back(x);  // 能接在所有已有序列后面
        } else {
            *it = x;             // 用更小的结尾替换，留出增长空间
        }
    }
    return static_cast<int>(tails.size());
}

// 动态规划，O(n^2)，返回一条具体的最长递增子序列
vector<int> lisDP(const vector<int>& nums) {
    int n = static_cast<int>(nums.size());
    if (n == 0) return {};
    vector<int> dp(n, 1), pre(n, -1);
    int best = 0;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < i; ++j) {
            if (nums[j] < nums[i] && dp[j] + 1 > dp[i]) {
                dp[i] = dp[j] + 1;
                pre[i] = j;
            }
        }
        if (dp[i] > dp[best]) best = i;
    }
    vector<int> seq;
    for (int k = best; k != -1; k = pre[k]) seq.push_back(nums[k]);
    reverse(seq.begin(), seq.end());
    return seq;
}

static string joinInts(const vector<int>& v) {
    ostringstream oss;
    for (size_t i = 0; i < v.size(); ++i) {
        if (i) oss << ' ';
        oss << v[i];
    }
    return oss.str();
}

int main() {
    int n;
    if (cin >> n) {  // IO 模式
        vector<int> nums(n);
        for (int i = 0; i < n; ++i) cin >> nums[i];
        vector<int> seq = lisDP(nums);
        cout << seq.size() << "\n" << joinInts(seq) << "\n";
        return 0;
    }

    // 经典用例
    assert(lengthOfLIS({10, 9, 2, 5, 3, 7, 101, 18}) == 4);
    assert(lisDP({10, 9, 2, 5, 3, 7, 101, 18}) == (vector<int>{2, 5, 7, 101}));
    // 相等元素不算递增
    assert(lengthOfLIS({7, 7, 7, 7}) == 1);
    assert(lisDP({7, 7, 7, 7}) == (vector<int>{7}));
    // 含重复但仍能取更长
    assert(lengthOfLIS({0, 1, 0, 3, 2, 3}) == 4);
    assert(lisDP({0, 1, 0, 3, 2, 3}) == (vector<int>{0, 1, 2, 3}));
    // 边界
    assert(lengthOfLIS({}) == 0);
    assert(lisDP({}).empty());
    assert(lengthOfLIS({1}) == 1);
    assert(lisDP({1}) == (vector<int>{1}));
    // 完全递减 / 完全递增
    assert(lengthOfLIS({5, 4, 3, 2, 1}) == 1);
    assert(lengthOfLIS({1, 2, 3, 4, 5}) == 5);
    assert(lisDP({1, 2, 3, 4, 5}) == (vector<int>{1, 2, 3, 4, 5}));

    cout << "all tests passed" << endl;
    return 0;
}
