"""并查集（Disjoint Set Union / Union-Find）

支持近乎 O(1) 的「合并」与「查询是否同组」。
两个优化缺一不可：路径压缩（find 时把节点直接挂到根上）+ 按秩/按大小合并（小树挂大树）。
"""

from typing import List


class UnionFind:
    def __init__(self, n: int) -> None:
        self.parent = list(range(n))  # parent[i] = i 的父节点
        self.size = [1] * n          # 以 i 为根的树的节点数
        self.count = n               # 连通分量个数

    def find(self, x: int) -> int:
        """带路径压缩的查找：递归回溯时把沿途节点直接挂到根上。"""
        if self.parent[x] != x:
            self.parent[x] = self.find(self.parent[x])
        return self.parent[x]

    def union(self, a: int, b: int) -> bool:
        """合并 a、b 所在集合；返回是否真的发生了合并（本来就在同组则返回 False）。"""
        ra, rb = self.find(a), self.find(b)
        if ra == rb:
            return False
        # 按大小合并：小树挂到大树下面，避免树退化成链
        if self.size[ra] < self.size[rb]:
            ra, rb = rb, ra
        self.parent[rb] = ra
        self.size[ra] += self.size[rb]
        self.count -= 1
        return True

    def connected(self, a: int, b: int) -> bool:
        return self.find(a) == self.find(b)

    def component_size(self, x: int) -> int:
        return self.size[self.find(x)]


def count_components(n: int, edges: List[List[int]]) -> int:
    """给定 n 个节点与若干无向边，返回连通分量个数。"""
    uf = UnionFind(n)
    for a, b in edges:
        uf.union(a, b)
    return uf.count


if __name__ == "__main__":
    uf = UnionFind(7)
    assert uf.count == 7
    assert uf.union(0, 1) is True
    assert uf.union(1, 2) is True
    assert uf.union(1, 2) is False          # 重复合并不再生效
    assert uf.connected(0, 2) is True
    assert uf.connected(0, 3) is False
    assert uf.component_size(0) == 3
    uf.union(3, 4)
    uf.union(5, 6)
    assert uf.count == 3                    # {0,1,2} {3,4} {5,6}

    assert count_components(5, [[0, 1], [1, 2], [3, 4]]) == 2
    assert count_components(4, []) == 4
    print("all tests passed")
