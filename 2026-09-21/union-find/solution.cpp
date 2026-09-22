// 并查集（Disjoint Set Union / Union-Find）
// 编译：g++ -std=c++17 -O2 solution.cpp -o solution && ./solution
#include <iostream>
#include <vector>
#include <cassert>

using namespace std;

class UnionFind {
public:
    explicit UnionFind(int n) : parent(n), sz(n, 1), components(n) {
        for (int i = 0; i < n; ++i) parent[i] = i;
    }

    // 带路径压缩的查找：递归回溯时把沿途节点直接挂到根上
    int find(int x) {
        if (parent[x] != x) parent[x] = find(parent[x]);
        return parent[x];
    }

    // 合并 a、b 所在集合；返回是否真的发生了合并
    bool unite(int a, int b) {
        int ra = find(a), rb = find(b);
        if (ra == rb) return false;
        // 按大小合并：小树挂到大树下面，避免树退化成链
        if (sz[ra] < sz[rb]) swap(ra, rb);
        parent[rb] = ra;
        sz[ra] += sz[rb];
        --components;
        return true;
    }

    bool connected(int a, int b) { return find(a) == find(b); }
    int componentSize(int x) { return sz[find(x)]; }
    int count() const { return components; }

private:
    vector<int> parent;  // parent[i] = i 的父节点
    vector<int> sz;      // 以 i 为根的树的节点数
    int components;      // 连通分量个数
};

// 给定 n 个节点与若干无向边，返回连通分量个数
int countComponents(int n, const vector<pair<int, int>>& edges) {
    UnionFind uf(n);
    for (auto& e : edges) uf.unite(e.first, e.second);
    return uf.count();
}

int main() {
    UnionFind uf(7);
    assert(uf.count() == 7);
    assert(uf.unite(0, 1) == true);
    assert(uf.unite(1, 2) == true);
    assert(uf.unite(1, 2) == false);  // 重复合并不再生效
    assert(uf.connected(0, 2) == true);
    assert(uf.connected(0, 3) == false);
    assert(uf.componentSize(0) == 3);
    uf.unite(3, 4);
    uf.unite(5, 6);
    assert(uf.count() == 3);  // {0,1,2} {3,4} {5,6}

    assert(countComponents(5, {{0, 1}, {1, 2}, {3, 4}}) == 2);
    assert(countComponents(4, {}) == 4);

    cout << "all tests passed" << endl;
    return 0;
}
