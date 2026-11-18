# Lyndon 分解（Duval）與循環串最小／最大表示（Booth）

## 題意

- **Lyndon 分解**：把任意字串唯一分解為 `c1 c2 ... ck`，其中每個 `ci` 都是 Lyndon 詞
  （嚴格小於它所有真後綴），且字典序 `c1 ≥ c2 ≥ ... ≥ ck`。
- **循環串最小 / 最大表示**：對 `s` 的所有 `n` 個循環移位，找出字典序最小（最大）者，
  輸出它的起始下標（0-based）與該移位字串；平手取最小下標。

## 思路

**Duval 演算法**：從位置 `i` 出發，用兩個指標 `j, k` 延伸「週期性主鏈」。
比較 `s[k]` 與 `s[j]`：相等則 `k++`（週期延續）；`s[k] < s[j]` 則 `k = i`（出現更大的字元，
週期重設）；`s[k] > s[j]` 則主鏈斷開。此時主鏈長度 `j - i`、週期 `j - k`，
從 `i` 起可連續切出若干個長度為 `j - k` 的 Lyndon 因子。`i`、`j` 單調前進，全程 O(n)。

**Booth 演算法**：把 `s` 複製成 `s + s`，維護兩個候選起點 `i, j` 與已匹配長度 `k`。
若 `s[i+k] > s[j+k]`，則 `i ... i+k` 這一段都不可能是最小表示起點，直接跳到 `i+k+1`；
反之跳 `j`。`k` 從不倒退，因此 O(n)。最大表示只需把比較方向反過來。

## 複雜度

| 項目 | 時間 | 空間 |
|---|---|---|
| Duval 分解 | O(n) | O(n)（輸出） |
| Booth 最小/最大表示 | O(n) | O(n)（`s + s`） |

## 輸入輸出

stdin 每筆命令兩個 token：`<cmd> <S>`，命令為 `lyndon` / `minrot` / `maxrot`。
讀到非法命令或 token 不足即結束。無 stdin 時跑內建隨機對拍（Duval vs 貪心暴力、
Booth vs 全枚舉），通過輸出 `all tests passed`。

### 範例

輸入：
```
lyndon banana
minrot banana
maxrot banana
lyndon aaaa
minrot cba
```

輸出：
```
b an an a
5 abanan
2 nanaba
a a a a
2 acb
```

## 應用

Lyndon 分解是 Burrows–Wheeler 變換與字串週期分析的基礎；Booth 演算法常用於
「最小表示法」問題（例如 Ring Rotation、等價字串正規化）。
