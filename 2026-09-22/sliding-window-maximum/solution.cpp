// 滑動窗口最大值（單調隊列 / Monotonic Queue）
// 編譯：g++ -std=c++17 -O2 solution.cpp -o solution && ./solution
//
// 思路：樸素做法每個窗口掃一遍取最大是 O(n*k)，k 大時超時。
//   單調隊列優化到 O(n)：隊列存"下標"，對應值嚴格遞減，隊首永遠是當前窗口最大值。
//   1) 入隊前從隊尾彈出所有 <= 當前值的下標——它們既更小、又更早出窗口，永遠不可能成爲答案；
//   2) 當前下標入隊；
//   3) 從隊首彈出所有已滑出窗口的下標（下標 <= i-k）；
//   4) i >= k-1 時，隊首下標對應的值即爲當前窗口最大值。
//   每個元素恰好入隊一次、出隊一次，故總時間線性。
//
// 輸入：第一行 n k；第二行 n 個整數
// 輸出：一行 n-k+1 個整數，空格分隔，爲各窗口最大值
// 無 stdin 輸入時運行內置斷言測試。
#include <cassert>
#include <deque>
#include <iostream>
#include <vector>

using namespace std;

// 返回長度爲 k 的滑動窗口在每個位置的最大值，時間 O(n)，空間 O(k)
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

    deque<int> q;  // 存下標，保證 nums[q.front()] > nums[q.back()]
    for (int i = 0; i < n; ++i) {
        // 1) 隊尾所有不大於當前值的下標都不可能再成爲答案
        while (!q.empty() && nums[q.back()] <= nums[i]) q.pop_back();
        // 2) 當前下標入隊
        q.push_back(i);
        // 3) 隊首已滑出窗口的下標出隊
        while (!q.empty() && q.front() <= i - k) q.pop_front();
        // 4) 窗口成型後，隊首即最大值
        if (i >= k - 1) ans.push_back(nums[q.front()]);
    }
    return ans;
}

// 對照用的樸素實現，O(n*k)，僅用於測試驗證
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

// 與 Python 版完全一致的固定隨機序列，保證兩版跑同一批用例
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
    assert((maxSlidingWindow({9, 8, 7, 6, 5}, 3) == vector<int>{9, 8, 7}));   // 遞減：隊首不斷被擠出
    assert((maxSlidingWindow({1, 2, 3, 4, 5}, 3) == vector<int>{3, 4, 5}));   // 遞增：隊尾不斷被彈出
    assert((maxSlidingWindow({5, 5, 5, 5}, 2) == vector<int>{5, 5, 5}));      // 全相等
    assert((maxSlidingWindow({-7, -8, -7, -6, -5}, 3) == vector<int>{-7, -6, -5}));
    assert((maxSlidingWindow({1, 3, 1, 2, 0, 5}, 3) == vector<int>{3, 3, 2, 5}));
    assert((maxSlidingWindow({4, 3, 2, 1}, 4) == vector<int>{4}));            // k == n
    assert((maxSlidingWindow({4, 3, 2, 1}, 5) == vector<int>{4}));            // k > n
    assert((maxSlidingWindow({}, 3) == vector<int>{}));                       // 空數組

    // 與樸素實現隨機對照：確保單調隊列沒有邊界錯誤
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
