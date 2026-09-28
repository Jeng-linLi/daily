"""併查集（Disjoint Set Union / Union-Find）

支持近乎 O(1) 的「合併」與「查詢是否同組」。
兩個優化缺一不可：路徑壓縮（find 時把節點直接掛到根上）+ 按秩/按大小合併（小樹掛大樹）。
"""

from typing import List


class UnionFind:
    def __init__(self, n: int) -> None:
        self.parent = list(range(n))  # parent[i] = i 的父節點
        self.size = [1] * n          # 以 i 爲根的樹的節點數
        self.count = n               # 連通分量個數

    def find(self, x: int) -> int:
        """帶路徑壓縮的查找：遞歸回溯時把沿途節點直接掛到根上。"""
        if self.parent[x] != x:
            self.parent[x] = self.find(self.parent[x])
        return self.parent[x]

    def union(self, a: int, b: int) -> bool:
        """合併 a、b 所在集合；返回是否真的發生了合併（本來就在同組則返回 False）。"""
        ra, rb = self.find(a), self.find(b)
        if ra == rb:
            return False
        # 按大小合併：小樹掛到大樹下面，避免樹退化成鏈
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
    """給定 n 個節點與若干無向邊，返回連通分量個數。"""
    uf = UnionFind(n)
    for a, b in edges:
        uf.union(a, b)
    return uf.count


if __name__ == "__main__":
    uf = UnionFind(7)
    assert uf.count == 7
    assert uf.union(0, 1) is True
    assert uf.union(1, 2) is True
    assert uf.union(1, 2) is False          # 重複合併不再生效
    assert uf.connected(0, 2) is True
    assert uf.connected(0, 3) is False
    assert uf.component_size(0) == 3
    uf.union(3, 4)
    uf.union(5, 6)
    assert uf.count == 3                    # {0,1,2} {3,4} {5,6}

    assert count_components(5, [[0, 1], [1, 2], [3, 4]]) == 2
    assert count_components(4, []) == 4
    print("all tests passed")
