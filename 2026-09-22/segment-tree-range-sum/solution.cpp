// 线段树（区间和 + 懒标记区间加）
// 编译：g++ -std=c++17 -O2 solution.cpp -o solution && ./solution
//
// 题意：给定数组 a，支持（0-based，闭区间 [l, r]）：
//     query(l, r) 求 a[l..r] 的和
//     add(l, r, v) 把 a[l..r] 每个元素加 v
//   暴力单次 O(n)，m 次操作 O(nm)；线段树把两种操作都降到 O(log n)。
//
// 思路：线段树是描述「区间」的二叉树——叶子存单个元素，内部节点存左右儿子区间的合并和。
//   区间查询：把目标区间拆成 O(log n) 个节点，完全覆盖的节点直接返回。
//   区间修改：整段被覆盖时只更新该节点 sum 并打懒标记 lazy，不再下递归；
//             下次必须深入时再 push_down 把标记分给两个儿子。
//   懒标记的本质是「延迟执行」：没人问细节就不必真的改到叶子。
//   堆式存储：节点 p 的左儿子 2p、右儿子 2p+1，空间开 4n。
//
// 输入：第一行 n q；第二行 n 个整数；接下来 q 行：1 l r（查询）/ 2 l r v（区间加）
// 输出：每个查询操作输出一行区间和
// 无 stdin 输入时运行内置断言测试。
#include <cassert>
#include <iostream>
#include <vector>

using namespace std;
using ll = long long;

class SegmentTree {
public:
    explicit SegmentTree(const vector<ll>& arr) : n(static_cast<int>(arr.size())) {
        sum.assign(4 * n + 5, 0);   // sum[p]：节点 p 对应区间的元素和
        lazy.assign(4 * n + 5, 0);  // lazy[p]：待下推给子孙的「整体加多少」
        if (n > 0) build(1, 0, n - 1, arr);
    }

    int size() const { return n; }

    // 区间加：a[ql..qr] += v，O(log n)
    void add(int ql, int qr, ll v) {
        if (n == 0) return;
        add(1, 0, n - 1, ql, qr, v);
    }

    // 区间求和：返回 a[ql..qr] 的和，O(log n)
    ll query(int ql, int qr) {
        if (n == 0) return 0;
        return query(1, 0, n - 1, ql, qr);
    }

private:
    int n;
    vector<ll> sum, lazy;

    // 自底向上建树，O(n)
    void build(int p, int l, int r, const vector<ll>& arr) {
        if (l == r) {
            sum[p] = arr[l];
            return;
        }
        int mid = (l + r) / 2;
        build(2 * p, l, mid, arr);
        build(2 * p + 1, mid + 1, r, arr);
        sum[p] = sum[2 * p] + sum[2 * p + 1];
    }

    // 给节点 p 整段加 v：更新 sum，并累积懒标记
    void apply(int p, int l, int r, ll v) {
        sum[p] += v * (r - l + 1);
        lazy[p] += v;
    }

    // 把 p 的懒标记下推给两个儿子，只在需要深入时使用
    void pushDown(int p, int l, int r) {
        if (lazy[p] == 0 || l == r) return;
        int mid = (l + r) / 2;
        ll v = lazy[p];
        apply(2 * p, l, mid, v);
        apply(2 * p + 1, mid + 1, r, v);
        lazy[p] = 0;
    }

    void add(int p, int l, int r, int ql, int qr, ll v) {
        if (ql <= l && r <= qr) {  // 完全覆盖，打标记后直接返回
            apply(p, l, r, v);
            return;
        }
        pushDown(p, l, r);
        int mid = (l + r) / 2;
        if (ql <= mid) add(2 * p, l, mid, ql, qr, v);
        if (qr > mid) add(2 * p + 1, mid + 1, r, ql, qr, v);
        sum[p] = sum[2 * p] + sum[2 * p + 1];
    }

    ll query(int p, int l, int r, int ql, int qr) {
        if (ql <= l && r <= qr) return sum[p];  // 完全覆盖，直接返回整段和
        pushDown(p, l, r);
        int mid = (l + r) / 2;
        ll ans = 0;
        if (ql <= mid) ans += query(2 * p, l, mid, ql, qr);
        if (qr > mid) ans += query(2 * p + 1, mid + 1, r, ql, qr);
        return ans;
    }
};

int main() {
    int n, q;
    if (cin >> n >> q) {  // IO 模式
        vector<ll> arr(n);
        for (int i = 0; i < n; ++i) cin >> arr[i];
        SegmentTree st(arr);
        for (int i = 0; i < q; ++i) {
            int kind, l, r;
            cin >> kind >> l >> r;
            if (kind == 1) {
                cout << st.query(l, r) << "\n";
            } else {
                ll v;
                cin >> v;
                st.add(l, r, v);
            }
        }
        return 0;
    }

    {
        vector<ll> arr = {1, 3, 5, 7, 9, 11};
        SegmentTree st(arr);
        assert(st.query(0, 5) == 36);
        assert(st.query(1, 3) == 15);  // 3 + 5 + 7
        assert(st.query(2, 2) == 5);   // 单点查询
        st.add(1, 3, 10);              // [1, 13, 15, 17, 9, 11]
        assert(st.query(1, 3) == 45);
        assert(st.query(0, 5) == 66);  // 36 + 3 * 10
        assert(st.query(0, 0) == 1);   // 区间外的点不受影响
        assert(st.query(4, 5) == 20);
    }

    // 懒标记叠加与部分覆盖交叉验证
    {
        SegmentTree st(vector<ll>(5, 0));
        // 结果数组为 [2, 5, 5, 5, 2]
        st.add(0, 4, 2);
        st.add(1, 3, 3);
        assert(st.query(0, 4) == 19);
        assert(st.query(0, 0) == 2);
        assert(st.query(1, 3) == 15);  // 5 + 5 + 5
        assert(st.query(4, 4) == 2);
    }

    // 单元素与空数组边界
    {
        SegmentTree st(vector<ll>{42});
        assert(st.query(0, 0) == 42);
        st.add(0, 0, -40);
        assert(st.query(0, 0) == 2);

        SegmentTree st0(vector<ll>{});
        assert(st0.query(0, 0) == 0);
    }

    cout << "all tests passed" << endl;
    return 0;
}
