// 插頭 DP（輪廓線 DP）：n×m 棋盤 1×2 骨牌完全覆蓋計數（C++17）
//
// 命令（IO 模式，從 stdin 讀入）：
//   count n m   -> 輸出 n×m 棋盤被 1×2 骨牌完全覆蓋的方案數（n, m <= 12）
//
// stdin 無內容時跑內建斷言測試，輸出 all tests passed。
//
// 核心要點（輪廓線 / 插頭 DP 的入門形態）：
//   * 逐格掃描，狀態是長度 m 的 bitmask：進入格子 (i, j) 時
//       - bit k (k < j)  = 第 i 行第 k 列向下伸的插頭（伸給 i+1 行）
//       - bit k (k >= j) = 上方伸下來的插頭（(i-1,k) 豎放下伸 = (i,k) 的上插頭；
//         (i,k-1) 橫放右伸 = (i,k) 的左插頭，編碼在同一個 bit 上）
//     骨牌每格度數為 2，同一輪廓位至多一個插頭，單 bit 即可。
//   * 轉移：輪廓位 j 已有插頭 -> 本格已被覆蓋，消費該插頭；
//     無插頭 -> 本格必須被覆蓋：豎放（bit j 置 1，伸向下一行）
//     或橫放（要求 bit j+1 為 0，置 bit j+1 = 1 給右鄰）。
//   * 最後一格無法豎放/橫放時該狀態自然死亡；終態要求 mask == 0
//     （殘留插頭意味著骨牌伸出棋盤）。
//   * 複雜度 O(n·m·2^m)，12×12 約 5.3e16 種，unsigned long long 足夠。

#include <cassert>
#include <cstdio>
#include <iostream>
#include <map>
#include <string>
#include <vector>

using namespace std;
typedef unsigned long long ull;

ull countDomino(int n, int m) {
    if (n <= 0 || m <= 0) return 0;
    map<int, ull> dp;
    dp[0] = 1;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < m; ++j) {
            map<int, ull> ndp;
            for (const auto& kv : dp) {
                int mask = kv.first;
                ull cnt = kv.second;
                if ((mask >> j) & 1) {
                    // 本格已被左鄰/上鄰的骨牌佔住：消費插頭
                    int nmask = mask & ~(1 << j);
                    ndp[nmask] += cnt;
                } else {
                    // 豎放：(i,j)-(i+1,j)，向 bit j 伸出
                    if (i + 1 < n) {
                        int nmask = mask | (1 << j);
                        ndp[nmask] += cnt;
                    }
                    // 橫放：(i,j)-(i,j+1)，bit j+1 必須空閒
                    if (j + 1 < m && !((mask >> (j + 1)) & 1)) {
                        int nmask = mask | (1 << (j + 1));
                        ndp[nmask] += cnt;
                    }
                }
            }
            dp.swap(ndp);
        }
    }
    auto it = dp.find(0);
    return it == dp.end() ? 0 : it->second;
}

ull bruteCount(int n, int m) {
    // 暴力 DFS 填格計數（僅用於小棋盤對拍）
    vector<bool> used((size_t)n * m, false);
    // 迭代式 DFS 用顯式棧太繞，遞歸即可（小棋盤深度 <= 14）
    struct Dfs {
        int n, m;
        vector<bool>* used;
        ull go(int pos) {
            while (pos < n * m && (*used)[pos]) ++pos;
            if (pos == n * m) return 1;
            ull total = 0;
            int r = pos / m, c = pos % m;
            (*used)[pos] = true;
            if (r + 1 < n && !(*used)[pos + m]) {  // 豎放
                (*used)[pos + m] = true;
                total += go(pos + 1);
                (*used)[pos + m] = false;
            }
            if (c + 1 < m && !(*used)[pos + 1]) {  // 橫放
                (*used)[pos + 1] = true;
                total += go(pos + 1);
                (*used)[pos + 1] = false;
            }
            (*used)[pos] = false;
            return total;
        }
    } dfs{n, m, &used};
    return dfs.go(0);
}

// 自製 LCG（與倉庫其他題目一致）
static ull g_state = 20261110;
ull rnd(ull mm) {
    g_state = (g_state * 1103515245 + 12345) & 0x7FFFFFFFULL;
    return g_state % mm;
}

void runTests() {
    // 1) 已知數列：1×m、2×n = Fib(n+1)、3×2k 數列
    assert(countDomino(1, 1) == 0);
    assert(countDomino(1, 2) == 1);
    assert(countDomino(2, 1) == 1);
    assert(countDomino(1, 7) == 0);
    assert(countDomino(1, 8) == 1);
    ull a = 1, b = 1;  // Fib
    for (int k = 1; k <= 10; ++k) {
        ull want = b;  // Fib(k+1)，初始 Fib(1)=Fib(2)=1
        ull next = a + b;
        a = b;
        b = next;
        assert(countDomino(2, k) == want);
    }
    ull seq3[5] = {3, 11, 41, 153, 571};  // 3×2, 3×4, 3×6, 3×8, 3×10
    for (int idx = 0; idx < 5; ++idx)
        assert(countDomino(3, 2 * (idx + 1)) == seq3[idx]);

    // 2) 對稱性 count(n,m) == count(m,n)
    for (int n = 1; n <= 7; ++n)
        for (int m = 1; m <= 7; ++m)
            assert(countDomino(n, m) == countDomino(m, n));

    // 3) 格數為奇數必為 0
    for (int n = 1; n <= 6; ++n)
        for (int m = 1; m <= 6; ++m)
            if ((n * m) % 2 == 1) assert(countDomino(n, m) == 0);

    // 4) 與暴力 DFS 對拍（小棋盤隨機）
    for (int t = 0; t < 60; ++t) {
        int n = (int)rnd(4) + 1;
        int m = (int)rnd(4) + 1;
        if (n * m > 14) continue;
        assert(countDomino(n, m) == bruteCount(n, m));
    }

    // 5) 大數值固定回歸（對照 OEIS A004003：n×n 骨牌鋪法數）
    assert(countDomino(8, 8) == 12988816ULL);
    assert(countDomino(10, 10) == 258584046368ULL);
    assert(countDomino(12, 12) == 53060477521960000ULL);
    assert(countDomino(7, 5) == 0);  // 奇數格數（35 格）
}

int main() {
    string cmd;
    if (cin >> cmd) {
        bool any = true;
        while (any) {
            if (cmd == "count") {
                long long n, m;
                if (!(cin >> n >> m)) break;
                cout << countDomino((int)n, (int)m) << '\n';
            } else {
                break;  // 未知命令
            }
            if (!(cin >> cmd)) break;
        }
    } else {
        runTests();
        printf("all tests passed\n");
    }
    return 0;
}
