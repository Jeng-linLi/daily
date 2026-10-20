import sys

sys.path.insert(0, ".")
from solution import DLX, brute, solve_exact_cover  # noqa

mat1 = [
    [0, 0, 1, 0, 1, 1, 0],
    [1, 0, 0, 1, 0, 0, 1],
    [0, 1, 1, 0, 0, 1, 0],
    [1, 0, 0, 1, 0, 0, 0],
    [0, 1, 0, 0, 0, 0, 1],
    [0, 0, 0, 1, 1, 0, 1],
]

print("brute:", brute(6, 7, mat1))
print("dlx  :", solve_exact_cover(6, 7, mat1, 10))

rows = []
for i in range(6):
    cols = [j for j in range(7) if mat1[i][j]]
    rows.append((i, cols))
d = DLX(7, rows)
print("size[1..7]:", d.size[1:8])
print("nodes:", len(d.L))
# 欄 1 (A) 的 D 串列應該是 row1(1), row3(3)
i = d.D[1]
seq = []
while i != 1:
    seq.append(d.ROW[i])
    i = d.D[i]
print("col A rows:", seq)
# row 1 的水平串列
print("row1 cols:", sorted(d.C[j] for j in range(len(d.L)) if d.ROW[j] == 1))
