// 归并排序与逆序对计数（Merge Sort & Inversion Count）
// 编译：g++ -std=c++17 -O2 solution.cpp -o solution && ./solution
//
// 思路：对区间 [lo, hi) 以 mid 切分后，逆序对 (i, j) 恰好分三类且不重不漏：
//     都在左半边 / 都在右半边（各自递归统计）/ i 在左、j 在右（跨中线）。
//   跨中线的那部分在「合并两个已有序子数组」时批量结算：若 arr[i] > arr[j]，
//   左半边 arr[i..mid) 全都 > arr[j]，一次比较就贡献 (mid - i) 个逆序对。
//
//   两个易错点：
//     - 比较写成 arr[i] <= arr[j] 才走左半边（取等号），否则相等元素被误判为逆序对；
//     - 逆序对上界是 n*(n-1)/2（完全逆序），计数必须用 long long，int 会溢出。
//
// 输入：第一行 n；第二行 n 个整数（可跨行）
// 输出：第一行逆序对个数；第二行升序排序后的序列（空格分隔；n = 0 时输出空行）
// 无 stdin 输入时运行内置断言测试。
#include <algorithm>
#include <cassert>
#include <iostream>
#include <sstream>
#include <vector>

using namespace std;

// 对 a[lo, hi) 归并排序，返回其中的逆序对个数（结果写回 a）
long long mergeSortCount(vector<int>& a, vector<int>& buf, int lo, int hi) {
    if (hi - lo <= 1) return 0;

    int mid = lo + (hi - lo) / 2;
    // 左右两半内部的逆序对各自递归统计
    long long inv = mergeSortCount(a, buf, lo, mid);
    inv += mergeSortCount(a, buf, mid, hi);

    int i = lo, j = mid, k = lo;
    while (i < mid && j < hi) {
        if (a[i] <= a[j]) {
            // 取等号：相等元素不构成逆序对，同时保证排序稳定
            buf[k++] = a[i++];
        } else {
            buf[k++] = a[j++];
            // 左半边 a[i..mid) 全部 > a[j]，一次性结算 mid - i 个逆序对
            inv += mid - i;
        }
    }
    while (i < mid) buf[k++] = a[i++];
    while (j < hi) buf[k++] = a[j++];

    for (int t = lo; t < hi; ++t) a[t] = buf[t];
    return inv;
}

// 返回 {逆序对个数, 升序排序后的新数组}。时间 O(n log n)，空间 O(n)
pair<long long, vector<int>> sortAndCount(const vector<int>& nums) {
    vector<int> a = nums;
    vector<int> buf(a.size());
    long long inv = mergeSortCount(a, buf, 0, static_cast<int>(a.size()));
    return {inv, a};
}

// 对照用的 O(n^2) 暴力枚举，仅用于小规模测试验证
long long countInversionsBrute(const vector<int>& nums) {
    long long cnt = 0;
    int n = static_cast<int>(nums.size());
    for (int i = 0; i < n; ++i)
        for (int j = i + 1; j < n; ++j)
            if (nums[i] > nums[j]) ++cnt;
    return cnt;
}

// 用标准库排序的结果做参照
vector<int> stdSorted(vector<int> v) {
    sort(v.begin(), v.end());
    return v;
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

    // README 中的示例：2 3 8 6 1 -> 逆序对 5 个
    //   (2,1) (3,1) (8,6) (8,1) (6,1)
    {
        auto res = sortAndCount({2, 3, 8, 6, 1});
        assert(res.first == 5);
        assert(res.second == (vector<int>{1, 2, 3, 6, 8}));
        assert(countInversionsBrute({2, 3, 8, 6, 1}) == 5);
    }

    // 空序列与单元素
    assert(sortAndCount({}).first == 0);
    assert(sortAndCount({}).second.empty());
    assert(sortAndCount({42}).first == 0);
    assert(sortAndCount({42}).second == (vector<int>{42}));

    // 已升序：0 个逆序对
    {
        auto res = sortAndCount({1, 2, 3, 4, 5});
        assert(res.first == 0);
        assert(res.second == (vector<int>{1, 2, 3, 4, 5}));
    }

    // 完全逆序：n*(n-1)/2 个逆序对，验证 64 位计数不溢出
    {
        auto res = sortAndCount({5, 4, 3, 2, 1});
        assert(res.first == 10);
        assert(res.second == (vector<int>{1, 2, 3, 4, 5}));

        vector<int> big(2000);
        for (int i = 0; i < 2000; ++i) big[i] = 2000 - i;
        assert(sortAndCount(big).first == 2000LL * 1999 / 2);
    }

    // 相等元素不算逆序对（这里最容易把 <= 写成 < 而数多）
    {
        auto res = sortAndCount({2, 2, 1});
        assert(res.first == 2);
        assert(res.second == (vector<int>{1, 2, 2}));
        assert(sortAndCount({1, 1, 1}).first == 0);
        auto res2 = sortAndCount({3, 1, 3, 1});
        assert(res2.first == 3);
        assert(res2.second == (vector<int>{1, 1, 3, 3}));
    }

    // 负数与零：逆序对为 (-1,-3) (-1,-2) (0,-2) (2,-2)，共 4 个
    {
        auto res = sortAndCount({-1, -3, 0, 2, -2});
        assert(res.first == 4);
        assert(res.second == (vector<int>{-3, -2, -1, 0, 2}));
    }

    // 与暴力解随机对拍：校验逆序对数一致、排序结果正确、长度不变
    LCG rng(20260923ULL);
    for (int t = 0; t < 300; ++t) {
        int len = rng.next(0, 40);
        vector<int> nums(len);
        for (int i = 0; i < len; ++i) nums[i] = rng.next(-20, 20);
        auto res = sortAndCount(nums);
        assert(res.first == countInversionsBrute(nums));   // 与暴力枚举一致
        assert(res.second == stdSorted(nums));              // 排序结果正确
        assert(static_cast<int>(res.second.size()) == len);  // 元素一个不多一个不少
    }

    // 排序不应改动调用方传入的原数组
    {
        vector<int> original = {3, 1, 2};
        sortAndCount(original);
        assert(original == (vector<int>{3, 1, 2}));
    }

    cout << "all tests passed" << endl;
    return 0;
}
