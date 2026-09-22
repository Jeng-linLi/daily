// 滑动窗口最大值（单调队列 / Monotonic Queue）
// 编译：g++ -std=c++17 -O2 solution.cpp -o solution && ./solution
//
// 思路：朴素做法每个窗口扫一遍取最大是 O(n*k)，k 大时超时。
//   单调队列优化到 O(n)：队列存"下标"，对应值严格递减，队首永远是当前窗口最大值。
//   1) 入队前从队尾弹出所有 <= 当前值的下标——它们既更小、又更早出窗口，永远不可能成为答案；
//   2) 当前下标入队；
//   3) 从队首弹出所有已滑出窗口的下标（下标 <= i-k）；
//   4) i >= k-1 时，队首下标对应的值即为当前窗口最大值。
//   每个元素恰好入队一次、出队一次，故总时间线性。
//
// 输入：第一行 n k；第二行 n 个整数
// 输出：一行 n-k+1 个整数，空格分隔，为各窗口最大值
// 无 stdin 输入时运行内置断言测试。
#include <cassert>
#include <deque>
#include <iostream>
#include <vector>

using namespace std;

// 返回长度为 k 的滑动窗口在每个位置的最大值，时间 O(n)，空间 O(k)
vector<int> maxSlidingWindow(const vector<int>& nums, int k) {
    vector<int> ans;
    int n = static_cast<int>(nums.size());
    if (n == 0 || k <= 0) return ans;
    if (k == 1) return vector<int>(nums.begin(), nums.end());
    if (k >= n) {
        int mx = nums[0];
        for (int v : nums) mx = max(mx, v);
        ans.push_back(mx);
        return ans;
    }

    deque<int> q;  // 存下标，保证 nums[q.front()] > nums[q.back()]
    for (int i = 0; i < n; ++i) {
        // 1) 队尾所有不大于当前值的下标都不可能再成为答案
        while (!q.empty() && nums[q.back()] <= nums[i]) q.pop_back();
        // 2) 当前下标入队
        q.push_back(i);
        // 3) 队首已滑出窗口的下标出队
        while (!q.empty() && q.front() <= i - k) q.pop_front();
        // 4) 窗口成型后，队首即最大值
        if (i >= k - 1) ans.push_back(nums[q.front()]);
    }
    return ans;
}

// 对照用的朴素实现，O(n*k)，仅用于测试验证
vector<int> bruteForce(const vector<int>& nums, int k) {
    vector<int> ans;
    int n = static_cast<int>(nums.size());
    for (int i = 0; i + k <= n; ++i) {
        int mx = nums[i];
        for (int j = i; j < i + k; ++j) mx = max(mx, nums[j]);
        ans.push_back(mx);
    }
    return ans;
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
    int n, k;
    if (cin >> n >> k) {  // IO 模式
        vector<int> nums(n);
        for (int i = 0; i < n; ++i) cin >> nums[i];
        vector<int> ans = maxSlidingWindow(nums, k);
        for (size_t i = 0; i < ans.size(); ++i) {
            if (i) cout << " ";
            cout << ans[i];
        }
        cout << "\n";
        return 0;
    }

    assert((maxSlidingWindow({1, 3, -1, -3, 5, 3, 6, 7}, 3) == vector<int>{3, 3, 5, 5, 6, 7}));
    assert((maxSlidingWindow({1}, 1) == vector<int>{1}));
    assert((maxSlidingWindow({1, -1}, 1) == vector<int>{1, -1}));
    assert((maxSlidingWindow({9, 8, 7, 6, 5}, 3) == vector<int>{9, 8, 7}));   // 递减：队首不断被挤出
    assert((maxSlidingWindow({1, 2, 3, 4, 5}, 3) == vector<int>{3, 4, 5}));   // 递增：队尾不断被弹出
    assert((maxSlidingWindow({5, 5, 5, 5}, 2) == vector<int>{5, 5, 5}));      // 全相等
    assert((maxSlidingWindow({-7, -8, -7, -6, -5}, 3) == vector<int>{-7, -6, -5}));
    assert((maxSlidingWindow({1, 3, 1, 2, 0, 5}, 3) == vector<int>{3, 3, 2, 5}));
    assert((maxSlidingWindow({4, 3, 2, 1}, 4) == vector<int>{4}));            // k == n
    assert((maxSlidingWindow({4, 3, 2, 1}, 5) == vector<int>{4}));            // k > n
    assert((maxSlidingWindow({}, 3) == vector<int>{}));                       // 空数组

    // 与朴素实现随机对照：确保单调队列没有边界错误
    LCG rng(20260922ULL);
    for (int t = 0; t < 200; ++t) {
        int len = rng.next(1, 40);
        int kk = rng.next(1, len);
        vector<int> arr(len);
        for (int i = 0; i < len; ++i) arr[i] = rng.next(-50, 50);
        assert(maxSlidingWindow(arr, kk) == bruteForce(arr, kk));
    }

    cout << "all tests passed" << endl;
    return 0;
}
