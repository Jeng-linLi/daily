// 歸併排序與逆序對計數（Merge Sort & Inversion Count）
// 編譯：g++ -std=c++17 -O2 solution.cpp -o solution && ./solution
//
// 思路：對區間 [lo, hi) 以 mid 切分後，逆序對 (i, j) 恰好分三類且不重不漏：
//     都在左半邊 / 都在右半邊（各自遞歸統計）/ i 在左、j 在右（跨中線）。
//   跨中線的那部分在「合併兩個已有序子數組」時批量結算：若 arr[i] > arr[j]，
//   左半邊 arr[i..mid) 全都 > arr[j]，一次比較就貢獻 (mid - i) 個逆序對。
//
//   兩個易錯點：
//     - 比較寫成 arr[i] <= arr[j] 才走左半邊（取等號），否則相等元素被誤判爲逆序對；
//     - 逆序對上界是 n*(n-1)/2（完全逆序），計數必須用 long long，int 會溢出。
//
// 輸入：第一行 n；第二行 n 個整數（可跨行）
// 輸出：第一行逆序對個數；第二行升序排序後的序列（空格分隔；n = 0 時輸出空行）
// 無 stdin 輸入時運行內置斷言測試。
#include <algorithm>
#include <cassert>
#include <iostream>
#include <sstream>
#include <vector>

using namespace std;

// 對 a[lo, hi) 歸併排序，返回其中的逆序對個數（結果寫回 a）
long long mergeSortCount(vector<int>& a, vector<int>& buf, int lo, int hi) {
    if (hi - lo <= 1) return 0;

    int mid = lo + (hi - lo) / 2;
    // 左右兩半內部的逆序對各自遞歸統計
    long long inv = mergeSortCount(a, buf, lo, mid);
    inv += mergeSortCount(a, buf, mid, hi);

    int i = lo, j = mid, k = lo;
    while (i < mid && j < hi) {
        if (a[i] <= a[j]) {
            // 取等號：相等元素不構成逆序對，同時保證排序穩定
            buf[k++] = a[i++];
        } else {
            buf[k++] = a[j++];
            // 左半邊 a[i..mid) 全部 > a[j]，一次性結算 mid - i 個逆序對
            inv += mid - i;
        }
    }
    while (i < mid) buf[k++] = a[i++];
    while (j < hi) buf[k++] = a[j++];

    for (int t = lo; t < hi; ++t) a[t] = buf[t];
    return inv;
}

// 返回 {逆序對個數, 升序排序後的新數組}。時間 O(n log n)，空間 O(n)
pair<long long, vector<int>> sortAndCount(const vector<int>& nums) {
    vector<int> a = nums;
    vector<int> buf(a.size());
    long long inv = mergeSortCount(a, buf, 0, static_cast<int>(a.size()));
    return {inv, a};
}

// 對照用的 O(n^2) 暴力枚舉，僅用於小規模測試驗證
long long countInversionsBrute(const vector<int>& nums) {
    long long cnt = 0;
    int n = static_cast<int>(nums.size());
    for (int i = 0; i < n; ++i)
        for (int j = i + 1; j < n; ++j)
            if (nums[i] > nums[j]) ++cnt;
    return cnt;
}

// 用標準庫排序的結果做參照
vector<int> stdSorted(vector<int> v) {
    sort(v.begin(), v.end());
    return v;
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
    int n;
    if (cin >> n) {  // IO 模式
        vector<int> nums(n);
        for (int i = 0; i < n; ++i) cin >> nums[i];
        auto res = sortAndCount(nums);
        const vector<int>& sortedNums = res.second;
        cout << res.first << "\n";
        for (size_t i = 0; i < sortedNums.size(); ++i) {
            if (i) cout << " ";
            cout << sortedNums[i];
        }
        cout << "\n";
        return 0;
    }

    // README 中的示例：2 3 8 6 1 -> 逆序對 5 個
    //   (2,1) (3,1) (8,6) (8,1) (6,1)
    {
        auto res = sortAndCount({2, 3, 8, 6, 1});
        assert(res.first == 5);
        assert(res.second == (vector<int>{1, 2, 3, 6, 8}));
        assert(countInversionsBrute({2, 3, 8, 6, 1}) == 5);
    }

    // 空序列與單元素
    assert(sortAndCount({}).first == 0);
    assert(sortAndCount({}).second.empty());
    assert(sortAndCount({42}).first == 0);
    assert(sortAndCount({42}).second == (vector<int>{42}));

    // 已升序：0 個逆序對
    {
        auto res = sortAndCount({1, 2, 3, 4, 5});
        assert(res.first == 0);
        assert(res.second == (vector<int>{1, 2, 3, 4, 5}));
    }

    // 完全逆序：n*(n-1)/2 個逆序對，驗證 64 位計數不溢出
    {
        auto res = sortAndCount({5, 4, 3, 2, 1});
        assert(res.first == 10);
        assert(res.second == (vector<int>{1, 2, 3, 4, 5}));

        vector<int> big(2000);
        for (int i = 0; i < 2000; ++i) big[i] = 2000 - i;
        assert(sortAndCount(big).first == 2000LL * 1999 / 2);
    }

    // 相等元素不算逆序對（這裡最容易把 <= 寫成 < 而數多）
    {
        auto res = sortAndCount({2, 2, 1});
        assert(res.first == 2);
        assert(res.second == (vector<int>{1, 2, 2}));
        assert(sortAndCount({1, 1, 1}).first == 0);
        auto res2 = sortAndCount({3, 1, 3, 1});
        assert(res2.first == 3);
        assert(res2.second == (vector<int>{1, 1, 3, 3}));
    }

    // 負數與零：逆序對爲 (-1,-3) (-1,-2) (0,-2) (2,-2)，共 4 個
    {
        auto res = sortAndCount({-1, -3, 0, 2, -2});
        assert(res.first == 4);
        assert(res.second == (vector<int>{-3, -2, -1, 0, 2}));
    }

    // 與暴力解隨機對拍：校驗逆序對數一致、排序結果正確、長度不變
    LCG rng(20260923ULL);
    for (int t = 0; t < 300; ++t) {
        int len = rng.next(0, 40);
        vector<int> nums(len);
        for (int i = 0; i < len; ++i) nums[i] = rng.next(-20, 20);
        auto res = sortAndCount(nums);
        assert(res.first == countInversionsBrute(nums));   // 與暴力枚舉一致
        assert(res.second == stdSorted(nums));              // 排序結果正確
        assert(static_cast<int>(res.second.size()) == len);  // 元素一個不多一個不少
    }

    // 排序不應改動調用方傳入的原數組
    {
        vector<int> original = {3, 1, 2};
        sortAndCount(original);
        assert(original == (vector<int>{3, 1, 2}));
    }

    cout << "all tests passed" << endl;
    return 0;
}
