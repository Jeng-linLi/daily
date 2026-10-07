# Rabin-Karp 與滾動哈希

| 語言 | 檔案 |
|---|---|
| Python 3 | [`solution.py`](solution.py) |
| C++17 | [`solution.cpp`](solution.cpp) |

## 題意

給定兩個字符串 `text` 與 `pattern`（ASCII），用**多項式滾動哈希（Rabin-Karp fingerprint）**解決五個問題：

1. `pattern` 在 `text` 中的所有出現位置（與 KMP、暴力法三方對拍）；
2. 兩串的**最長公共子串**（長度 + 最靠左的那一個）；
3. `text` 的**最長回文子串**（長度 + 最靠左的那一個）；
4. 兩串各自的**不同子串個數**；
5. 兩串的 **3-gram Jaccard 相似度**（應用：文件指紋 / 抄襲檢測），以最簡分數 `p/q` 輸出。

## 思路

### 多項式滾動哈希

把字符串看成 `base` 進制的多項式：

```
H(s) = ( ord(s[0])·base^(n-1) + ord(s[1])·base^(n-2) + … + ord(s[n-1]) ) mod M
```

等價的遞推寫法是 `h[i+1] = h[i]·base + ord(s[i])`，預處理 O(n)。則任意子串 `s[l:r]` 的指紋為

```
h[r] − h[l]·base^(r−l)   (mod M)
```

也就是「把前綴左移對齊後相減」，O(1) 得到任意子串指紋。這就是 Rabin-Karp 的核心：窗口滑動時指紋可以增量更新，不必重新掃一遍子串。

### 為什麼用雙哈希 + 真實比對

單模哈希有生日碰撞風險：n 個子串時碰撞概率約 `n² / 2M`。這裡同時取兩個質數模數 `1e9+7` 與 `1e9+9`（碰撞概率降到 ~1e-18），並且在每個哈希命中的位置**再真實比對一次字符**（verify），因此輸出結果與暴力法逐字節一致，不存在哈希碰撞導致的誤判。

### 五個子問題

| 子問題 | 做法 | 複雜度 |
|---|---|---|
| 子串匹配 | 滑動長度 \|pat\| 的窗口比指紋 | 平均 O(n+m)，最壞 O(n·m) |
| 最長公共子串 | 「存在長度 L 的公共子串」對 L **單調** ⟹ 二分 L，判定用哈希集合 | O((n+m)·log min(n,m)) |
| 最長回文子串 | 枚舉 2n−1 個中心，**對半徑二分**（半徑單調） | O(n log n) |
| 不同子串個數 | 枚舉長度 + 指紋集合 | O(n²) |
| k-gram Jaccard | 滑動片段（shingle）集合的 Jaccard | O(n+m) |

> ⚠️ **回文不能直接對「長度」二分**：存在長度 4 的回文不代表存在長度 3 的回文（例如 `baab`），
> 所以「存在長度**恰為** L 的回文」對 L **不單調**。
> 但「以某個固定中心、半徑 r 的子串是回文」對 r **單調**（半徑 r 成立 ⟹ r−1 成立），
> 所以正確做法是枚舉中心、對半徑二分。這是本次實作踩到的第一個坑。

### 應用：k-gram Jaccard（抄襲檢測）

把文本切成長度 k 的滑動片段集合（shingle），相似度 `= |A∩B| / |A∪B|`。
這是抄襲檢測、網頁去重（Broder's shingling）、文件指紋的經典做法，
也是「把字符串問題轉成集合問題」的典型應用。程式中以**最簡分數**輸出（用 gcd 約分），
避免浮點誤差讓 Python 與 C++ 的輸出對不上。

其它應用場景：生物序列比對（DNA 片段 fingerprint）、大文件差分同步（rsync 的弱滾動校驗和）、
編譯器與 IDE 的增量字符串搜索。

## 複雜度

記 `n = len(text)`，`m = len(pattern)`。

| 項目 | 時間 | 空間 |
|---|---|---|
| 預處理 | O(n + m) | O(n + m) |
| 子串匹配 | 平均 O(n + m)，最壞 O(n·m) | O(n + m) |
| 最長公共子串 | O((n + m)·log min(n, m)) | O(n + m) |
| 最長回文子串 | O(n·log n) | O(n) |
| 不同子串個數 | O(n²) | O(n²)（指紋集合） |
| k-gram 相似度 | O(n + m) | O(n + m) |

## 輸入輸出

輸入（stdin）：

```
第 1 行：text
第 2 行：pattern
```

兩行都可以是空行（代表空串）；不足兩行時缺的部分視為空串；多餘的行忽略。

輸出（stdout）：

```
第 1 行：pattern 在 text 中的出現次數
第 2 行：所有起始位置（0-indexed，空格分隔；無則輸出空行）
第 3 行：最長公共子串長度
第 4 行：最長公共子串（最靠左的那個）
第 5 行：text 的最長回文子串長度
第 6 行：最長回文子串（最靠左的那個）
第 7 行：text 的不同子串個數
第 8 行：pattern 的不同子串個數
第 9 行：3-gram Jaccard 相似度，最簡分數 p/q
```

無 stdin 輸入時運行內置斷言測試並輸出 `all tests passed`。

### 範例 1

輸入：

```
abracadabra
abra
```

輸出：

```
2
0 7
4
abra
3
aca
54
9
2/7
```

（`abra` 出現在 0 與 7；最長公共子串就是 `abra` 本身，長 4；`abracadabra` 的最長回文是 `aca`；
`abracadabra` 有 54 個不同子串，`abra` 有 9 個；3-gram 集合 `{abr,bra,rac,aca,cad,ada,dab}` 與 `{abr,bra}` 的 Jaccard = 2/7。）

### 範例 2

輸入：

```
babad
abba
```

輸出：

```
0

2
ba
3
bab
12
8
0/1
```

（`abba` 不在 `babad` 中；最長公共子串 `ba`（長 2，最靠左）；`babad` 的最長回文 `bab`（`aba` 同長但位置較右）；
兩串的 3-gram 集合完全不相交，相似度 0/1。）

## 怎麼跑

```bash
python solution.py          # 無輸入 → 跑內置斷言
printf 'abracadabra\nabra\n' | python solution.py

g++ -std=c++17 -O2 solution.cpp -o solution && ./solution
```

> Windows / MinGW 下不要把編譯產物輸出到 `/tmp`（鏈接器會報 No such file or directory），輸出到當前目錄即可。
