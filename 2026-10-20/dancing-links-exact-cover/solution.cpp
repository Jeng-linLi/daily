// Dancing Links（DLX）與精確覆蓋：Knuth Algorithm X 的雙向環狀串列實作
//
// 題意：
//     給一個 0/1 矩陣，**精確覆蓋（Exact Cover）**要選出若干列（row），
//     使得每一欄（column）恰好被一個被選到的列覆蓋一次。
//
//     本程式用 Knuth 的 **Algorithm X + Dancing Links（DLX）**求解，並提供兩種模式：
//
//         mode 1  一般精確覆蓋：給定 n × m 的 0/1 矩陣，輸出找到的解數（上限 k）
//                 與字典序最小的解，並與 2^n 暴力枚舉對拍（n ≤ 14 時）
//         mode 2  數獨求解：把 B²×B² 的數獨轉成精確覆蓋（4·N² 個約束欄），
//                 N = B²，輸出填完的盤面
//
// 思路：
//     ### 1. Algorithm X：遞迴 + 回溯
//     - 選一個「還沒被覆蓋」的欄 c（啟發式：選**剩餘候選列最少**的那一欄，
//       這是所謂的 S-heuristic，能把分支因子壓到最小）。
//     - 枚舉所有在 c 欄為 1 的列 r：選 r，把 r 覆蓋到的所有欄都從矩陣中刪掉，
//       遞迴；回溯時再把它們放回去。
//
//     ### 2. Dancing Links：O(1) 的刪除與還原
//     樸素做法是每次複製一份矩陣，代價 O(nm)。Knuth 的巧思是把每個 1 做成一個節點，
//     用**雙向環狀串列**橫向（同一列）與縱向（同一欄）串起來，加上一個「欄表頭」，
//     再串成一條含 root 的環。刪除一個欄只要：
//
//         R[L[c]] = R[c]; L[R[c]] = L[c]            // 從欄表頭環上摘掉
//
//     而**還原**只需要反向做一次（因為被摘掉的節點指標根本沒被清掉，
//     它還記得左右的鄰居）——這就是 "dancing links" 名字的由來：
//
//         R[L[c]] = c; L[R[c]] = c
//
//     這個技巧讓「覆蓋 / 取消覆蓋」是 O(該欄的 1 的個數)，而不是 O(nm)，
//     實務上數獨這種規模（324 欄、729 列）能在毫秒級解完。
//
//     ### 3. 確定性：兩語言必須長出一模一樣的結構
//     - 節點插入順序固定：逐列（row id 遞增）逐欄加入，同一欄內的 D 串列
//       因此就是 row id 遞增的順序。
//     - 選欄用 S-heuristic，平手時取**欄編號最小**的那一個（從左往右掃、
//       只有嚴格小於才更新）。
//     - 因為 uncover 是 cover 的精確逆運算，欄表頭環的左右順序在回溯後完全一致，
//       所以「從左往右掃」的結果可重現，Python 與 C++ 能逐位元組對拍。
//
//     ### 4. 數獨怎麼轉成精確覆蓋
//     N = B²，盤面 N×N。對每個 (r, c, d) 建一列（共 N³ 列），約束欄分四組各 N² 個：
//         1. 格子 (r,c) 必須填一個數      2. 第 r 列必須有數字 d
//         3. 第 c 行必須有數字 d          4. 第 b 個宮必須有數字 d
//     已填的格子只保留對應的 d 那一列，其餘從候選中刪掉。
//     精確覆蓋的一組解 ⟺ 一組合法的數獨填法。
//
//     ### 5. 全 0 的列要事先剔除
//     一個全 0 的列覆蓋不到任何欄，選它或不選它都不影響「每欄恰好一次」，
//     會讓解的數量白白乘以 2。兩個版本都在建結構時**跳過全 0 的列**，
//     暴力對拍也一樣略過，兩邊語義才一致。
//
// 應用場景：
//     - 數獨 / KenKen / Star Battle / Hidato 等謎題求解器（DLX 是數獨的標準解法）
//     - N 皇后、八皇后變體（轉成精確覆蓋後可瞬間數出全部解）
//     - 多聯骨牌（Pentomino）鋪滿矩形棋盤
//     - 排班問題：護士輪班、機組指派、會議室分配（本質是「集合劃分」）
//     - VLSI 測試向量最小化、覆蓋陣列（covering array）生成
//     - 把 NP 問題建模成精確覆蓋後，用 DLX 當作小型的組合搜尋引擎
//
// 複雜度：
//     建結構                            O(n · m)（等於矩陣中 1 的個數）
//     cover / uncover 一個欄            O(該欄目前的列數)
//     Algorithm X 搜尋                  指數級（精確覆蓋是 NP-complete），
//                                       但 S-heuristic 在實務上極快
//     空間                              O(1 的個數 + m)

#include <algorithm>
#include <cerrno>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <string>
#include <array>
#include <vector>

#ifdef _WIN32
#include <fcntl.h>
#include <io.h>
#endif

namespace {

const int NMAX = 40;        // 矩陣列數上限
const int MMAX = 40;        // 矩陣欄數上限
const int KMAX = 10;        // 解數上限
const int BF_NMAX = 14;     // 超過這個列數就不做 2^n 暴力對拍
const int BMAX = 4;         // 數獨宮的邊長上限（B=3 → 9×9，B=4 → 16×16）

// ---------------------------------------------------------------- DLX

class DLX {
public:
    // ncol 個欄（編號 1..ncol），rows 是 (row id, 欄編號串列)
    DLX(int ncol, const std::vector<std::pair<int, std::vector<int>>>& rows) {
        for (int j = 0; j <= ncol; ++j) {
            L_.push_back(j - 1);
            R_.push_back(j + 1);
            U_.push_back(j);
            D_.push_back(j);
            C_.push_back(j);
            ROW_.push_back(-1);
            size_.push_back(0);
        }
        L_[0] = ncol;
        R_[ncol] = 0;

        for (const auto& pr : rows) {
            int rid = pr.first;
            const std::vector<int>& cols = pr.second;
            if (cols.empty()) continue;             // 全 0 的列直接略過
            int first = -1;
            for (int c : cols) {
                int idx = static_cast<int>(L_.size());
                L_.push_back(0);
                R_.push_back(0);
                U_.push_back(0);
                D_.push_back(0);
                C_.push_back(c);
                ROW_.push_back(rid);
                size_.push_back(0);
                // 垂直：插到欄 c 的環狀串列尾部（表頭上方）
                U_[idx] = U_[c];
                D_[idx] = c;
                D_[U_[c]] = idx;
                U_[c] = idx;
                size_[c] += 1;
                // 水平：串成這一列的環
                if (first < 0) {
                    first = idx;
                    L_[idx] = idx;
                    R_[idx] = idx;
                } else {
                    L_[idx] = L_[first];
                    R_[idx] = first;
                    R_[L_[first]] = idx;
                    L_[first] = idx;
                }
            }
        }
    }

    void cover(int c) {
        R_[L_[c]] = R_[c];
        L_[R_[c]] = L_[c];
        int i = D_[c];
        while (i != c) {
            int j = R_[i];
            while (j != i) {
                U_[D_[j]] = U_[j];
                D_[U_[j]] = D_[j];
                size_[C_[j]] -= 1;
                j = R_[j];
            }
            i = D_[i];
        }
    }

    void uncover(int c) {
        int i = U_[c];
        while (i != c) {
            int j = L_[i];
            while (j != i) {
                size_[C_[j]] += 1;
                U_[D_[j]] = j;
                D_[U_[j]] = j;
                j = L_[j];
            }
            i = U_[i];
        }
        R_[L_[c]] = c;
        L_[R_[c]] = c;
    }

    // Algorithm X：回傳 (找到的解數, 字典序最小的解)
    std::pair<int, std::vector<int>> solve(int limit) {
        limit_ = limit;
        cnt_ = 0;
        have_ = false;
        best_.clear();
        std::vector<int> cur;
        search(cur);
        return {cnt_, best_};
    }

    bool hasBest() const { return have_; }

private:
    void search(std::vector<int>& cur) {
        if (cnt_ >= limit_) return;
        if (R_[0] == 0) {
            std::vector<int> s = cur;
            std::sort(s.begin(), s.end());
            ++cnt_;
            if (!have_ || s < best_) {
                best_ = s;
                have_ = true;
            }
            return;
        }
        // S-heuristic：候選列最少的欄；平手取欄編號最小（從左往右掃，嚴格小於才更新）
        int c = -1;
        int bestSz = -1;
        int j = R_[0];
        while (j != 0) {
            if (bestSz < 0 || size_[j] < bestSz) {
                bestSz = size_[j];
                c = j;
            }
            j = R_[j];
        }
        if (bestSz == 0) return;                 // 有欄湊不出來，此路不通
        cover(c);
        int i = D_[c];
        while (i != c) {
            cur.push_back(ROW_[i]);
            j = R_[i];
            while (j != i) {
                cover(C_[j]);
                j = R_[j];
            }
            search(cur);
            cur.pop_back();
            j = L_[i];
            while (j != i) {
                uncover(C_[j]);
                j = L_[j];
            }
            i = D_[i];
            if (cnt_ >= limit_) break;
        }
        uncover(c);
    }

    std::vector<int> L_, R_, U_, D_, C_, ROW_, size_;
    int limit_ = 1;
    int cnt_ = 0;
    bool have_ = false;
    std::vector<int> best_;
};

// ---------------------------------------------------------------- 暴力對拍

std::pair<int, std::vector<int>> brute(int n, int m, const std::vector<std::vector<int>>& mat) {
    // 2^n 枚舉所有列的子集（n ≤ 14）。全 0 的列與 DLX 一樣略過。
    std::vector<unsigned long long> masks;
    std::vector<int> ids;
    for (int i = 0; i < n; ++i) {
        unsigned long long mk = 0;
        for (int j = 0; j < m; ++j) {
            if (mat[i][j]) mk |= (1ULL << j);
        }
        if (mk) {
            masks.push_back(mk);
            ids.push_back(i);
        }
    }
    unsigned long long full = (m >= 64) ? ~0ULL : ((1ULL << m) - 1ULL);
    int cnt = 0;
    bool have = false;
    std::vector<int> best;
    int k = static_cast<int>(masks.size());
    for (int mask = 0; mask < (1 << k); ++mask) {
        unsigned long long acc = 0;
        bool ok = true;
        std::vector<int> sel;
        for (int t = 0; t < k; ++t) {
            if ((mask >> t) & 1) {
                if (acc & masks[t]) {
                    ok = false;
                    break;
                }
                acc |= masks[t];
                sel.push_back(ids[t]);
            }
        }
        if (ok && acc == full) {
            ++cnt;
            if (!have || sel < best) {
                best = sel;
                have = true;
            }
        }
    }
    if (cnt == 0) return {0, std::vector<int>()};
    return {cnt, best};
}

std::pair<int, std::vector<int>> solveExactCover(int n, int m,
                                                 const std::vector<std::vector<int>>& mat,
                                                 int k) {
    std::vector<std::pair<int, std::vector<int>>> rows;
    for (int i = 0; i < n; ++i) {
        std::vector<int> cols;
        for (int j = 0; j < m; ++j) {
            // DLX 的欄編號是 1..m，節點 0 被 root 佔用
            if (mat[i][j]) cols.push_back(j + 1);
        }
        if (!cols.empty()) rows.push_back({i, cols});
    }
    DLX dlx(m, rows);
    return dlx.solve(k);
}

// ---------------------------------------------------------------- 數獨

struct SudokuModel {
    int ncol;
    std::vector<std::pair<int, std::vector<int>>> rows;
    std::vector<std::array<int, 3>> meta;
};

SudokuModel sudokuRows(int B, const std::vector<std::vector<int>>& grid) {
    int N = B * B;
    int cells = N * N;
    auto cellCol = [&](int r, int c) { return r * N + c + 1; };
    auto rowCol = [&](int r, int d) { return cells + r * N + d + 1; };
    auto colCol = [&](int c, int d) { return 2 * cells + c * N + d + 1; };
    auto boxCol = [&](int r, int c, int d) {
        int b = (r / B) * B + (c / B);
        return 3 * cells + b * N + d + 1;
    };

    SudokuModel model;
    model.ncol = 4 * cells;
    int rid = 0;
    for (int r = 0; r < N; ++r) {
        for (int c = 0; c < N; ++c) {
            int given = grid[r][c];
            for (int d = 0; d < N; ++d) {
                if (given != 0 && given != d + 1) continue;
                std::vector<int> cols;
                cols.push_back(cellCol(r, c));
                cols.push_back(rowCol(r, d));
                cols.push_back(colCol(c, d));
                cols.push_back(boxCol(r, c, d));
                model.rows.push_back({rid, cols});
                model.meta.push_back({r, c, d});
                ++rid;
            }
        }
    }
    return model;
}

std::vector<std::vector<int>> solveSudoku(int B, const std::vector<std::vector<int>>& grid) {
    int N = B * B;
    SudokuModel model = sudokuRows(B, grid);
    DLX dlx(model.ncol, model.rows);
    std::pair<int, std::vector<int>> res = dlx.solve(1);
    if (res.first == 0) return std::vector<std::vector<int>>();
    std::vector<std::vector<int>> out = grid;
    for (int rid : res.second) {
        const std::array<int, 3>& t = model.meta[rid];
        out[t[0]][t[1]] = t[2] + 1;
    }
    // 合法性檢查（自洽檢查，兩語言都做）
    for (int i = 0; i < N; ++i) {
        std::vector<int> seen(N + 1, 0);
        for (int j = 0; j < N; ++j) seen[out[i][j]] += 1;
        for (int d = 1; d <= N; ++d) {
            if (seen[d] != 1) {
                std::cout << "FATAL sudoku row " << i << "\n";
                std::abort();
            }
        }
        seen.assign(N + 1, 0);
        for (int j = 0; j < N; ++j) seen[out[j][i]] += 1;
        for (int d = 1; d <= N; ++d) {
            if (seen[d] != 1) {
                std::cout << "FATAL sudoku col " << i << "\n";
                std::abort();
            }
        }
        seen.assign(N + 1, 0);
        int br = (i / B) * B;
        int bc = (i % B) * B;
        for (int dr = 0; dr < B; ++dr) {
            for (int dc = 0; dc < B; ++dc) seen[out[br + dr][bc + dc]] += 1;
        }
        for (int d = 1; d <= N; ++d) {
            if (seen[d] != 1) {
                std::cout << "FATAL sudoku box " << i << "\n";
                std::abort();
            }
        }
    }
    return out;
}

// ---------------------------------------------------------------- 輸出組裝

std::string runMode1(int n, int m, const std::vector<std::vector<int>>& mat, int k) {
    std::pair<int, std::vector<int>> res = solveExactCover(n, m, mat, k);
    std::ostringstream oss;
    oss << res.first;
    for (int v : res.second) oss << " " << v;
    oss << "\n";
    if (n > BF_NMAX) {
        oss << "-1\n-1\n";
        return oss.str();
    }
    std::pair<int, std::vector<int>> bf = brute(n, m, mat);
    int ok = (res.first == std::min(bf.first, k)) ? 1 : 0;
    if (bf.first <= k && ok) ok = (res.second == bf.second) ? 1 : 0;
    oss << bf.first << "\n" << ok << "\n";
    return oss.str();
}

std::string runMode2(int B, const std::vector<std::vector<int>>& grid) {
    std::vector<std::vector<int>> sol = solveSudoku(B, grid);
    int N = B * B;
    if (sol.empty()) return "no solution\n";
    std::ostringstream oss;
    for (int r = 0; r < N; ++r) {
        for (int c = 0; c < N; ++c) {
            if (c) oss << " ";
            oss << sol[r][c];
        }
        oss << "\n";
    }
    return oss.str();
}

// ---------------------------------------------------------------- 自製 LCG

unsigned long long lcgState = 0x9E3779B97F4A7C15ULL;

unsigned long long lcgNext() {
    lcgState = lcgState * 6364136223846793005ULL + 1442695040888963407ULL;
    return lcgState;
}

long long lcgInt(long long lo, long long hi) {
    unsigned long long r = lcgNext() >> 11;
    long long span = hi - lo + 1;
    return lo + static_cast<long long>(r % static_cast<unsigned long long>(span));
}

// ---------------------------------------------------------------- 內置測試

void runTests() {
    // ---- 固定用例 1：經典的小精確覆蓋（Knuth 論文裡的例子）----
    {
        std::vector<std::vector<int>> mat = {
            {0, 0, 1, 0, 1, 1, 0},
            {1, 0, 0, 1, 0, 0, 1},
            {0, 1, 1, 0, 0, 1, 0},
            {1, 0, 0, 1, 0, 0, 0},
            {0, 1, 0, 0, 0, 0, 1},
            {0, 0, 0, 1, 1, 0, 1},
        };
        std::string out = runMode1(6, 7, mat, 10);
        if (out != "1 0 3 4\n1\n1\n") {
            std::cout << "FATAL fixed1 [" << out << "]\n";
            std::abort();
        }
    }

    // ---- 固定用例 2：無解（某一欄全 0）----
    {
        std::vector<std::vector<int>> mat = {{1, 0}, {0, 0}, {1, 0}};
        std::string out = runMode1(3, 2, mat, 10);
        if (out != "0\n0\n1\n") {
            std::cout << "FATAL fixed2 [" << out << "]\n";
            std::abort();
        }
    }

    // ---- 固定用例 3：欄數為 0（空集合也是一個解）----
    {
        std::vector<std::vector<int>> mat(3, std::vector<int>(0));
        std::string out = runMode1(3, 0, mat, 10);
        if (out != "1\n1\n1\n") {
            std::cout << "FATAL fixed3 [" << out << "]\n";
            std::abort();
        }
    }

    // ---- 固定用例 4：多解，驗證解數上限 k ----
    {
        std::vector<std::vector<int>> mat = {{1, 0}, {0, 1}, {1, 0}, {0, 1}};
        std::string out = runMode1(4, 2, mat, 10);
        if (out != "4 0 1\n4\n1\n") {
            std::cout << "FATAL fixed4a [" << out << "]\n";
            std::abort();
        }
        std::string out2 = runMode1(4, 2, mat, 2);
        if (out2 != "2 0 1\n4\n1\n") {
            std::cout << "FATAL fixed4b [" << out2 << "]\n";
            std::abort();
        }
    }

    // ---- 數獨：經典 9×9 ----
    {
        std::vector<std::vector<int>> puzzle = {
            {5, 3, 0, 0, 7, 0, 0, 0, 0},
            {6, 0, 0, 1, 9, 5, 0, 0, 0},
            {0, 9, 8, 0, 0, 0, 0, 6, 0},
            {8, 0, 0, 0, 6, 0, 0, 0, 3},
            {4, 0, 0, 8, 0, 3, 0, 0, 1},
            {7, 0, 0, 0, 2, 0, 0, 0, 6},
            {0, 6, 0, 0, 0, 0, 2, 8, 0},
            {0, 0, 0, 4, 1, 9, 0, 0, 5},
            {0, 0, 0, 0, 8, 0, 0, 7, 9},
        };
        std::string txt = runMode2(3, puzzle);
        std::vector<std::string> lines;
        std::istringstream iss(txt);
        std::string ln;
        while (std::getline(iss, ln)) lines.push_back(ln);
        if (lines.size() != 9) {
            std::cout << "FATAL sudoku lines " << lines.size() << "\n";
            std::abort();
        }
        if (lines[0] != "5 3 4 6 7 8 9 1 2") {
            std::cout << "FATAL sudoku l0 [" << lines[0] << "]\n";
            std::abort();
        }
        if (lines[8] != "3 4 5 2 8 6 1 7 9") {
            std::cout << "FATAL sudoku l8 [" << lines[8] << "]\n";
            std::abort();
        }
    }

    // ---- 數獨：4×4 ----
    {
        std::vector<std::vector<int>> p4(4, std::vector<int>(4, 0));
        std::string txt = runMode2(2, p4);
        if (txt.substr(0, txt.find('\n')) != "1 2 3 4") {
            std::cout << "FATAL sudoku4 [" << txt << "]\n";
            std::abort();
        }
    }

    // ---- 數獨：無解 ----
    {
        std::vector<std::vector<int>> bad(4, std::vector<int>(4, 0));
        bad[0][0] = 1;
        bad[0][1] = 2;
        bad[0][2] = 3;
        bad[0][3] = 4;
        bad[1][0] = 3;
        bad[1][1] = 1;                 // 與 (0,0)=1 同列重複 1 → 無解
        std::string txt = runMode2(2, bad);
        if (txt != "no solution\n") {
            std::cout << "FATAL sudoku-bad [" << txt << "]\n";
            std::abort();
        }
    }

    // ---- 隨機對拍：DLX vs 2^n 暴力 ----
    {
        lcgState = 20261020ULL;
        for (int t = 0; t < 600; ++t) {
            int n = static_cast<int>(lcgInt(0, 10));
            int m = static_cast<int>(lcgInt(0, 8));
            int pp = static_cast<int>(lcgInt(15, 70));       // 百分比
            std::vector<std::vector<int>> mat(n, std::vector<int>(m, 0));
            for (int i = 0; i < n; ++i) {
                for (int j = 0; j < m; ++j) {
                    if (lcgInt(0, 99) < pp) mat[i][j] = 1;
                }
            }
            std::string out = runMode1(n, m, mat, 10);
            size_t p = out.rfind('\n', out.size() - 2);
            std::string flag = out.substr(p + 1);
            if (flag != "1\n") {
                std::cout << "FATAL fuzz " << t << " [" << out << "]\n";
                std::abort();
            }
        }
    }

    // ---- 隨機對拍：稍微大一點的矩陣（n ≤ 14）----
    {
        for (int t = 0; t < 200; ++t) {
            int n = static_cast<int>(lcgInt(8, 14));
            int m = static_cast<int>(lcgInt(4, 10));
            int pp = static_cast<int>(lcgInt(20, 50));
            std::vector<std::vector<int>> mat(n, std::vector<int>(m, 0));
            for (int i = 0; i < n; ++i) {
                for (int j = 0; j < m; ++j) {
                    if (lcgInt(0, 99) < pp) mat[i][j] = 1;
                }
            }
            std::string out = runMode1(n, m, mat, 5);
            size_t p = out.rfind('\n', out.size() - 2);
            std::string flag = out.substr(p + 1);
            if (flag != "1\n") {
                std::cout << "FATAL fuzz2 " << t << " [" << out << "]\n";
                std::abort();
            }
        }
    }

    // ---- 隨機數獨（挖空 9×9）----
    {
        std::vector<std::vector<int>> empty(9, std::vector<int>(9, 0));
        std::vector<std::vector<int>> base = solveSudoku(3, empty);
        for (int t = 0; t < 6; ++t) {
            lcgState = 3000ULL + static_cast<unsigned long long>(t);
            std::vector<std::vector<int>> g = base;
            for (int h = 0; h < 55; ++h) {
                int pos = static_cast<int>(lcgInt(0, 80));
                g[pos / 9][pos % 9] = 0;
            }
            std::vector<std::vector<int>> sol = solveSudoku(3, g);
            if (sol.empty()) {
                std::cout << "FATAL rsudoku " << t << "\n";
                std::abort();
            }
            for (int r = 0; r < 9; ++r) {
                for (int c = 0; c < 9; ++c) {
                    if (g[r][c] != 0 && sol[r][c] != g[r][c]) {
                        std::cout << "FATAL rsudoku mismatch " << t << "\n";
                        std::abort();
                    }
                }
            }
        }
    }
}

// ---------------------------------------------------------------- IO

bool isInt(const std::string& s) {
    if (s.empty()) return false;
    size_t i = (s[0] == '+' || s[0] == '-') ? 1 : 0;
    if (i >= s.size()) return false;
    for (size_t j = i; j < s.size(); ++j) {
        if (s[j] < '0' || s[j] > '9') return false;
    }
    return true;
}

long long clampLL(const std::string& t, long long lo, long long hi, long long def) {
    if (!isInt(t)) return def;
    errno = 0;
    char* end = nullptr;
    long long v = std::strtoll(t.c_str(), &end, 10);
    if (errno == ERANGE || end != t.c_str() + t.size()) return def;
    if (v < lo) v = lo;
    if (v > hi) v = hi;
    return v;
}

void runIo(const std::string& raw) {
    std::istringstream iss(raw);
    std::vector<std::string> toks;
    std::string t;
    while (iss >> t) toks.push_back(t);
    size_t idx = 0;
    auto nextInt = [&](long long def) -> long long {
        if (idx >= toks.size()) return def;
        std::string s = toks[idx++];
        return clampLL(s, -1000000000000000000LL, 1000000000000000000LL, def);
    };

    long long mode = nextInt(1);
    if (mode == 2) {
        long long B = nextInt(3);
        if (B < 1) B = 1;
        if (B > BMAX) B = BMAX;
        int N = static_cast<int>(B) * static_cast<int>(B);
        std::vector<std::vector<int>> grid(N, std::vector<int>(N, 0));
        for (int r = 0; r < N; ++r) {
            for (int c = 0; c < N; ++c) {
                long long v = nextInt(0);
                if (v < 0) v = 0;
                if (v > N) v = N;
                grid[r][c] = static_cast<int>(v);
            }
        }
        std::cout << runMode2(static_cast<int>(B), grid);
        return;
    }

    long long n = nextInt(0);
    long long m = nextInt(0);
    if (n < 0) n = 0;
    if (m < 0) m = 0;
    if (n > NMAX) n = NMAX;
    if (m > MMAX) m = MMAX;
    std::vector<std::vector<int>> mat(static_cast<size_t>(n),
                                      std::vector<int>(static_cast<size_t>(m), 0));
    for (long long i = 0; i < n; ++i) {
        for (long long j = 0; j < m; ++j) {
            long long v = nextInt(0);
            mat[static_cast<size_t>(i)][static_cast<size_t>(j)] = (v != 0) ? 1 : 0;
        }
    }
    long long k = nextInt(1);
    if (k < 0) k = 0;
    if (k > KMAX) k = KMAX;
    std::cout << runMode1(static_cast<int>(n), static_cast<int>(m), mat, static_cast<int>(k));
}

}  // namespace

int main() {
#ifdef _WIN32
    _setmode(_fileno(stdout), _O_BINARY);
#endif
    std::string raw((std::istreambuf_iterator<char>(std::cin)), std::istreambuf_iterator<char>());
    bool blank = true;
    for (char c : raw) {
        unsigned char u = static_cast<unsigned char>(c);
        if (u > ' ') {
            blank = false;
            break;
        }
    }
    if (!blank) {
        runIo(raw);
    } else {
        runTests();
        std::cout << "all tests passed" << std::endl;
    }
    return 0;
}
