// 最長遞增子序列（Longest Increasing Subsequence, LIS）
// 編譯：g++ -std=c++17 -O2 solution.cpp -o solution && ./solution
//
// 思路：
//   方法一 · 動態規劃 O(n^2)
//       dp[i] 表示「以 nums[i] 作爲結尾」的最長遞增子序列長度。
//       dp[i] = 1 + max{ dp[j] | j < i 且 nums[j] < nums[i] }，取不到則 dp[i] = 1。
//       pre[i] 記前驅下標，從 dp 最大處往回跳即可還原序列。
//   方法二 · 貪心 + 二分 O(n log n)
//       tails[k] = 長度爲 k+1 的遞增子序列的最小結尾值，tails 本身嚴格遞增。
//       對每個 x 二分找到第一個 >= x 的位置：能接長就追加，否則用 x 覆蓋（結尾越小潛力越大）。
//       tails 未必是合法子序列，但 len(tails) 一定是答案。
//
// 輸入：第一行 n；第二行 n 個整數
// 輸出：第一行 LIS 長度；第二行一條 LIS（空格分隔）
// 無 stdin 輸入時運行內置斷言測試。
#include <algorithm>
#include <cassert>
#include <iostream>
#include <sstream>
#include <vector>

using namespace std;

// 貪心 + 二分，O(n log n)，只求長度
int lengthOfLIS(const vector<int>& nums) {
    vector<int> tails;  // tails[k] = 長度 k+1 的遞增子序列的最小結尾
    for (int x : nums) {
        auto it = lower_bound(tails.begin(), tails.end(), x);  // 第一個 >= x
        if (it == tails.end()) {
            tails.push_back(x);  // 能接在所有已有序列後面
        } else {
            *it = x;             // 用更小的結尾替換，留出增長空間
        }
    }
    return static_cast<int>(tails.size());
}

// 動態規劃，O(n^2)，返回一條具體的最長遞增子序列
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

    // 經典用例
    assert(lengthOfLIS({10, 9, 2, 5, 3, 7, 101, 18}) == 4);
    assert(lisDP({10, 9, 2, 5, 3, 7, 101, 18}) == (vector<int>{2, 5, 7, 101}));
    // 相等元素不算遞增
    assert(lengthOfLIS({7, 7, 7, 7}) == 1);
    assert(lisDP({7, 7, 7, 7}) == (vector<int>{7}));
    // 含重複但仍能取更長
    assert(lengthOfLIS({0, 1, 0, 3, 2, 3}) == 4);
    assert(lisDP({0, 1, 0, 3, 2, 3}) == (vector<int>{0, 1, 2, 3}));
    // 邊界
    assert(lengthOfLIS({}) == 0);
    assert(lisDP({}).empty());
    assert(lengthOfLIS({1}) == 1);
    assert(lisDP({1}) == (vector<int>{1}));
    // 完全遞減 / 完全遞增
    assert(lengthOfLIS({5, 4, 3, 2, 1}) == 1);
    assert(lengthOfLIS({1, 2, 3, 4, 5}) == 5);
    assert(lisDP({1, 2, 3, 4, 5}) == (vector<int>{1, 2, 3, 4, 5}));

    cout << "all tests passed" << endl;
    return 0;
}
