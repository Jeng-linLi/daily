// 樹上莫隊（Mo's algorithm on tree）與迴滾莫隊（rollback Mo's algorithm）
//
// 題意：
//     給一棵 n 個點的有根樹（根為 1 號點），每個點帶一個整數權值。
//     q 個詢問，每次給 (u, v)，問 u 到 v 的**路徑上**有幾個不同的權值
//     （路徑含端點；u = v 時路徑就是單一個點）。
//
//     IO 格式（有 stdin 時）：
//         第 1 行    n q
//         第 2 行    a_1 a_2 ... a_n         權值
//         第 3 行    p_2 p_3 ... p_n         點 i（2..n）的父節點
//         接下來 q 行  u v
//     輸出 q 行，每行一個整數。
//
// 思路：
//     ### 1. 把「路徑查詢」拍平成一維「區間查詢」
//     對樹做一次 DFS，記錄每個點的**進入時間** tin 與**離開時間** tout，
//     得到一條長度 2n 的歐拉序（每個點恰好出現兩次：進入、離開）。
//     對詢問 (u, v)（不妨設 tin[u] < tin[v]），路徑 u→v 在歐拉序上對應：
//
//         - 若 u 是 v 的祖先：一段連續區間 [tin[u], tin[v]]，
//           路徑上每個點在區間裡恰好出現一次，其餘點出現 0 或 2 次；
//         - 否則：區間 [tin[u], tin[v]] 裡，u 出現兩次（進+出）、
//           LCA(u, v) 出現零次，其餘路徑點各一次。
//           於是「區間 + 額外翻轉 u 與 LCA」就正好得到整條路徑。
//
//     翻轉用 vis[] 標記實現：區間裡出現兩次的點來回翻轉兩次等於沒翻，
//     出現一次的點留在「啟用」狀態——這正是莫隊維護「區間內點集」的標準技巧。
//
//     ### 2. 迴滾莫隊（rollback Mo）
//     普通莫隊在區間上移動左右端點，需要「加」與「減」兩種操作；
//     迴滾莫隊只保留「加」，左端點向左移產生的修改全部記在一棵**棧**上，
//     回答完當前詢問後統一彈棧撤銷。適用於「刪不掉」的資料結構
//     （例如只支援 insert 的平衡樹、二元組統計等）。
//
//     具體結構（塊大小 B ≈ √m，m 為歐拉序長度）：
//         - 詢問按 (l 所在塊, r) 排序，同一塊內 r 單調遞增；
//         - 進入新塊時清空全部狀態，左右端點都停在**塊尾**；
//         - r 只增不減（只加）；l 從塊尾往左加，記棧；
//         - 回答完把棧彈回檢查點，l 復原到塊尾，r 保留；
//         - 若某詢問的 r 落在塊內（r < 塊尾，整窗在塊中），
//           先整體回滾到空窗，再一次性加入 [l, r]，答完再回滾——
//           這類「短詢問」代價 O(B)，不破壞「只加」的不變式。
//
//     本題的「加」= 翻轉一個點（維護每種權值的出現次數 cnt 與不同權值數 d），
//     翻轉是自逆操作，棧裡記下被翻轉的點，反序再翻一次即精確撤銷。
//
//     ### 3. LCA 用倍增法
//     O(n log n) 預處理、O(log n) 查詢；也可以用歐拉序 + Sparse Table 做到 O(1)，
//     但莫隊本身是 O(√n) 級別，倍增常數小、好寫，兩語言行為也完全一致。
//
//     ### 4. 兩語言一致性
//     - 答案與 DFS 順序無關（路徑點集是樹的固有屬性），歐拉序不同不影響輸出；
//     - 莫隊的排序鍵 (l / B, r) 兩邊一致；平手不影響答案；
//     - 詢問按輸入順序輸出，與處理順序解耦。
//
// 複雜度：
//     預處理（DFS + 倍增）    O(n log n)
//     單次詢問                O(√n)（均攤，迴滾莫隊標準複雜度）
//     q 次詢問合計            O((n + q) √n)
//     空間                    O(n log n)
//
// 應用場景：
//     - 遊戲 / 社交網路：「共同好友」「幫派勢力範圍」等路徑聚合統計
//     - 版本控制系統：兩次 commit 之間修改過的檔案種類數（diff 統計）
//     - 資料庫查詢最佳化：樹形目錄下「某兩節點路徑上的標籤數」
//     - 離線演算法競賽經典題型：靜態樹路徑第 k 小 / 種類數 / 權值和

#include <algorithm>
#include <cctype>
#include <cerrno>
#include <cstdlib>
#include <iostream>
#include <iterator>
#include <random>
#include <set>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#ifdef _WIN32
#include <fcntl.h>
#include <io.h>
#endif

namespace {

const int NMAX = 20000;                   // IO 模式的節點數上限
const int QMAX = 20000;                   // IO 模式的詢問數上限
const long long VMAX = 1000000000LL;      // 權值絕對值上限
const long long CLAMP_ABS = 1000000000000000000LL;  // 通用整數夾取範圍

// ---------------------------------------------------------------- 樹與 LCA

struct TreeMoSolver {
    int n;
    std::vector<int> parent;
    int root;
    std::vector<int> first, tout, tour, depth;
    std::vector<std::vector<int>> up;
    int log_;

    TreeMoSolver(int n_, const std::vector<int>& values,
                 const std::vector<int>& parent_)
        : n(n_), parent(parent_) {
        (void)values;
        if (n == 0) {
            root = 0;
            log_ = 1;
            up.assign(1, std::vector<int>());
            return;
        }
        std::vector<std::vector<int>> children(n);
        root = 0;
        for (int v = 0; v < n; ++v) {
            int p = parent[v];
            if (p < 0) {
                root = v;
            } else if (p < n && p != v) {
                children[p].push_back(v);
            }
        }
        // 自環一律接回根，保證任何輸入都能終止
        for (int v = 0; v < n; ++v) {
            if (parent[v] == v) children[root].push_back(v);
        }

        // ---- 迭代 DFS 求歐拉序（進入 / 離開各記一次）----
        first.assign(n, 0);
        tout.assign(n, 0);
        tour.clear();
        tour.reserve(2 * n);
        std::vector<std::pair<int, bool>> st;
        st.reserve(2 * n);
        st.push_back(std::make_pair(root, true));
        while (!st.empty()) {
            int v = st.back().first;
            bool enter = st.back().second;
            st.pop_back();
            if (enter) {
                first[v] = static_cast<int>(tour.size());
                tour.push_back(v);
                st.push_back(std::make_pair(v, false));
                // 反序入棧，讓編號小的子樹先遍歷（決定性，方便對拍）
                for (int i = static_cast<int>(children[v].size()) - 1; i >= 0; --i) {
                    st.push_back(std::make_pair(children[v][i], true));
                }
            } else {
                tout[v] = static_cast<int>(tour.size());
                tour.push_back(v);
            }
        }

        // ---- BFS 求深度 ----
        depth.assign(n, 0);
        std::vector<int> order;
        order.push_back(root);
        std::vector<char> seen(n, 0);
        seen[root] = 1;
        for (size_t head = 0; head < order.size(); ++head) {
            int v = order[head];
            for (size_t i = 0; i < children[v].size(); ++i) {
                int c = children[v][i];
                if (!seen[c]) {
                    seen[c] = 1;
                    depth[c] = depth[v] + 1;
                    order.push_back(c);
                }
            }
        }

        // ---- 倍增 ----
        log_ = 1;
        while ((1LL << log_) < static_cast<long long>(n)) ++log_;
        up.assign(log_, std::vector<int>(n, 0));
        for (int v = 0; v < n; ++v) {
            int p = parent[v];
            up[0][v] = (p >= 0 && p < n && p != v) ? p : root;
        }
        up[0][root] = root;
        for (int k = 1; k < log_; ++k) {
            for (int v = 0; v < n; ++v) {
                up[k][v] = up[k - 1][up[k - 1][v]];
            }
        }
    }

    int lca(int u, int v) const {
        if (depth[u] < depth[v]) std::swap(u, v);
        int diff = depth[u] - depth[v];
        int k = 0;
        while (diff) {
            if (diff & 1) u = up[k][u];
            diff >>= 1;
            ++k;
        }
        if (u == v) return u;
        for (int k2 = log_ - 1; k2 >= 0; --k2) {
            if (up[k2][u] != up[k2][v]) {
                u = up[k2][u];
                v = up[k2][v];
            }
        }
        return parent[u] >= 0 ? parent[u] : u;
    }

    // 詢問 → (歐拉序區間 [l, r], 需要額外翻轉的點)
    void toWindow(int u, int v, int& l, int& r,
                  std::vector<int>& extras) const {
        extras.clear();
        int i = first[u], j = first[v];
        if (i > j) {
            std::swap(u, v);
            std::swap(i, j);
        }
        int a = lca(u, v);
        if (a == u) {
            l = i;
            r = j;
        } else {
            l = i;
            r = j;
            extras.push_back(u);
            extras.push_back(a);
        }
    }
};

// ---------------------------------------------------------------- 迴滾莫隊

struct RollbackMo {
    const std::vector<int>& tour;
    const std::vector<int>& values;  // 離散化後的權值
    int m;
    int B;
    int n;
    std::vector<int> cnt;
    std::vector<char> vis;
    int d;
    std::vector<int> stack;
    int curL, curR;

    RollbackMo(const std::vector<int>& tour_, const std::vector<int>& values_,
               int nvals, int n_)
        : tour(tour_), values(values_), n(n_) {
        m = static_cast<int>(tour_.size());
        long long b = 1;
        while (b * b < static_cast<long long>(m)) ++b;
        B = static_cast<int>(b);
        cnt.assign(std::max(nvals, 1), 0);
        vis.assign(std::max(n, 1), 0);
        d = 0;
        curL = 0;
        curR = -1;
    }

    void toggle(int x) {
        // 翻轉點 x 的啟用狀態（自逆操作，不記棧）
        int c = values[x];
        if (vis[x]) {
            if (--cnt[c] == 0) --d;
            vis[x] = 0;
        } else {
            if (cnt[c] == 0) ++d;
            ++cnt[c];
            vis[x] = 1;
        }
    }

    void addNode(int x) {
        toggle(x);
        stack.push_back(x);
    }

    void rollback(size_t snapSize, int curLv, int curRv) {
        while (stack.size() > snapSize) {
            toggle(stack.back());
            stack.pop_back();
        }
        curL = curLv;
        curR = curRv;
    }

    // windows[i] = {l, r, extra0, extra1, ...}；回傳與輸入同序的答案
    std::vector<int> solve(const std::vector<std::vector<int>>& windows) {
        int q = static_cast<int>(windows.size());
        std::vector<int> ans(q, 0);
        std::vector<int> order(q);
        for (int i = 0; i < q; ++i) order[i] = i;
        std::sort(order.begin(), order.end(), [&](int a, int b) {
            int ba = windows[a][0] / B, bb = windows[b][0] / B;
            if (ba != bb) return ba < bb;
            return windows[a][1] < windows[b][1];
        });
        int curBlk = -1;
        int blockR = 0;
        for (int oi = 0; oi < q; ++oi) {
            int idx = order[oi];
            int l = windows[idx][0];
            int r = windows[idx][1];
            const std::vector<int>& extras = windows[idx];
            int blk = l / B;
            if (blk != curBlk) {
                // 新塊：清空全部狀態，端點停在塊尾
                stack.clear();
                std::fill(cnt.begin(), cnt.end(), 0);
                std::fill(vis.begin(), vis.end(), 0);
                d = 0;
                long long br = (static_cast<long long>(blk) + 1) * B;
                if (br > m) br = m;
                blockR = static_cast<int>(br) - 1;
                curL = blockR;
                curR = blockR - 1;
                curBlk = blk;
            }
            if (r < blockR) {
                // 短詢問：整個窗口在塊內 → 回滾到空窗再一次加滿 [l, r]
                rollback(0, blockR, blockR - 1);
                curR = r;
                curL = r + 1;
                while (curL > l) {
                    --curL;
                    addNode(tour[curL]);
                }
                for (size_t t = 2; t < extras.size(); ++t) addNode(extras[t]);
                ans[idx] = d;
                rollback(0, blockR, blockR - 1);
            } else {
                // 長詢問：r 只增（只加），l 從塊尾往左加（記棧，答完撤銷）
                while (curR < r) {
                    ++curR;
                    addNode(tour[curR]);
                }
                size_t save = stack.size();
                while (curL > l) {
                    --curL;
                    addNode(tour[curL]);
                }
                for (size_t t = 2; t < extras.size(); ++t) addNode(extras[t]);
                ans[idx] = d;
                // 只撤銷 l 的部分；r 的部分保留給同一塊的後續詢問
                rollback(save, blockR, curR);
            }
        }
        return ans;
    }
};

// ---------------------------------------------------------------- 暴力對拍

int bfDistinct(int n, const std::vector<int>& values,
               const std::vector<int>& parent, int u, int v) {
    // 沿父指標爬行收集路徑點集，統計不同權值（完全獨立的參考實作）
    std::vector<char> anc(std::max(n, 1), 0);
    int x = u, steps = 0;
    while (x >= 0 && steps <= n) {
        anc[x] = 1;
        x = parent[x];
        ++steps;
    }
    int y = v;
    steps = 0;
    while (y >= 0 && steps <= n && !anc[y]) {
        y = parent[y];
        ++steps;
    }
    int l = (y >= 0) ? y : u;
    std::set<int> seen;
    x = u;
    steps = 0;
    while (x != l && steps <= n) {
        seen.insert(values[x]);
        x = (parent[x] >= 0) ? parent[x] : l;
        ++steps;
    }
    seen.insert(values[l]);
    x = v;
    steps = 0;
    while (x != l && steps <= n) {
        seen.insert(values[x]);
        x = (parent[x] >= 0) ? parent[x] : l;
        ++steps;
    }
    return static_cast<int>(seen.size());
}

// ---------------------------------------------------------------- 主流程

std::string run(int n, const std::vector<int>& values,
                const std::vector<int>& parent,
                const std::vector<std::pair<int, int>>& queries) {
    TreeMoSolver solver(n, values, parent);
    // 權值離散化
    std::vector<int> uniq(values);
    std::sort(uniq.begin(), uniq.end());
    uniq.erase(std::unique(uniq.begin(), uniq.end()), uniq.end());
    std::vector<int> cvals(n);
    for (int i = 0; i < n; ++i) {
        cvals[i] = static_cast<int>(
            std::lower_bound(uniq.begin(), uniq.end(), values[i]) - uniq.begin());
    }
    std::vector<std::vector<int>> windows(queries.size());
    for (size_t i = 0; i < queries.size(); ++i) {
        int l = 0, r = 0;
        std::vector<int> extras;
        solver.toWindow(queries[i].first, queries[i].second, l, r, extras);
        windows[i].push_back(l);
        windows[i].push_back(r);
        for (size_t t = 0; t < extras.size(); ++t) windows[i].push_back(extras[t]);
    }
    RollbackMo mo(solver.tour, cvals, static_cast<int>(uniq.size()), n);
    std::vector<int> ans = mo.solve(windows);
    std::string res;
    for (size_t i = 0; i < ans.size(); ++i) {
        if (i) res += "\n";
        res += std::to_string(ans[i]);
    }
    if (!ans.empty()) res += "\n";
    return res;
}

// ---------------------------------------------------------------- 內置測試

void fail(const std::string& msg) {
    std::cerr << "FAIL: " << msg << std::endl;
    std::exit(1);
}

std::vector<int> parseAnswers(const std::string& text) {
    std::vector<int> out;
    std::istringstream iss(text);
    int v;
    while (iss >> v) out.push_back(v);
    return out;
}

void checkAgainstBrute(int n, const std::vector<int>& values,
                       const std::vector<int>& parent,
                       const std::vector<std::pair<int, int>>& queries,
                       const std::vector<int>& got, const std::string& tag) {
    if (static_cast<int>(got.size()) != static_cast<int>(queries.size())) {
        fail(tag + ": answer count " + std::to_string(got.size()));
    }
    for (size_t i = 0; i < queries.size(); ++i) {
        int w = bfDistinct(n, values, parent, queries[i].first, queries[i].second);
        if (got[i] != w) {
            fail(tag + ": query " + std::to_string(i) + " got " +
                 std::to_string(got[i]) + " want " + std::to_string(w));
        }
    }
}

void runTests() {
    // ---- 固定用例：樹 1-2, 1-3, 2-4, 2-5（根 1）----
    // 權值（0 起算）：[10, 20, 30, 20, 10]
    {
        int n = 5;
        std::vector<int> values = {10, 20, 30, 20, 10};
        std::vector<int> parent = {-1, 0, 0, 1, 1};
        std::vector<std::pair<int, int>> qs = {
            {3, 4}, {2, 3}, {3, 3}, {0, 4}, {1, 2}, {4, 2}};
        // (3,4): 路徑 3-1-4      → {20,10}       → 2
        // (2,3): 路徑 2-0-1-3    → {30,10,20}    → 3
        // (3,3): 單點            → {20}          → 1
        // (0,4): 路徑 0-1-4      → {10,20}       → 2
        // (1,2): 路徑 1-0-2      → {20,10,30}    → 3
        // (4,2): 路徑 4-1-0-2    → {10,20,30}    → 3
        std::vector<int> want = {2, 3, 1, 2, 3, 3};
        std::string got = run(n, values, parent, qs);
        std::string exp;
        for (size_t i = 0; i < want.size(); ++i) {
            if (i) exp += "\n";
            exp += std::to_string(want[i]);
        }
        exp += "\n";
        if (got != exp) fail("fixed case: got [" + got + "]");
        checkAgainstBrute(n, values, parent, qs, want, "fixed case brute");
    }

    // ---- 單點樹 ----
    {
        std::vector<int> values = {7};
        std::vector<int> parent = {-1};
        std::vector<std::pair<int, int>> qs = {{0, 0}};
        if (run(1, values, parent, qs) != "1\n") fail("single node");
    }

    // ---- 沒有詢問 → 空輸出 ----
    {
        std::vector<int> values = {1, 2, 3};
        std::vector<int> parent = {-1, 0, 0};
        std::vector<std::pair<int, int>> qs;
        if (run(3, values, parent, qs) != "") fail("no query");
    }

    // ---- 全部權值相同 → 每個答案都是 1 ----
    {
        std::vector<int> values = {5, 5, 5, 5};
        std::vector<int> parent = {-1, 0, 0, 0};
        std::vector<std::pair<int, int>> qs = {{0, 1}, {1, 3}, {2, 3}};
        if (run(4, values, parent, qs) != "1\n1\n1\n") fail("all same values");
    }

    // ---- 鏈狀樹（最深路徑，壓力測試塊邏輯）----
    {
        int n = 1500;
        std::vector<int> values(n);
        for (int i = 0; i < n; ++i) values[i] = (i * 37) % 9;
        std::vector<int> parent(n, -1);
        for (int v = 1; v < n; ++v) parent[v] = v - 1;
        std::mt19937 rng(777);
        std::vector<std::pair<int, int>> qs(150);
        for (int i = 0; i < 150; ++i) {
            qs[i] = std::make_pair(static_cast<int>(rng() % n),
                                   static_cast<int>(rng() % n));
        }
        checkAgainstBrute(n, values, parent, qs,
                          parseAnswers(run(n, values, parent, qs)), "chain");
    }

    // ---- 隨機對拍：四種樹形 × 兩種值域 ----
    {
        std::mt19937 rng(20261023);
        for (int t = 0; t < 400; ++t) {
            int n = 1 + static_cast<int>(rng() % 60);
            int shape = t % 4;
            std::vector<int> parent(n, -1);
            if (shape == 0) {                    // 隨機樹
                for (int v = 1; v < n; ++v) parent[v] = static_cast<int>(rng() % v);
            } else if (shape == 1) {             // 鏈
                for (int v = 1; v < n; ++v) parent[v] = v - 1;
            } else if (shape == 2) {             // 星
                for (int v = 1; v < n; ++v) parent[v] = 0;
            } else {                             // 雙層
                for (int v = 1; v < n; ++v) parent[v] = 0;
                for (int v = 2; v < n; ++v) {
                    if (rng() % 2 == 0) parent[v] = 1;
                }
            }
            bool small = (t % 2 == 0);
            std::vector<int> values(n);
            for (int i = 0; i < n; ++i) {
                if (small) {
                    values[i] = static_cast<int>(rng() % 4);  // 強制大量重複
                } else {
                    values[i] = static_cast<int>(
                        static_cast<long long>(rng() % 2000000000LL) - 1000000000LL);
                }
            }
            int q = static_cast<int>(rng() % 21);
            std::vector<std::pair<int, int>> queries(q);
            for (int i = 0; i < q; ++i) {
                queries[i] = std::make_pair(static_cast<int>(rng() % n),
                                            static_cast<int>(rng() % n));
            }
            checkAgainstBrute(n, values, parent, queries,
                              parseAnswers(run(n, values, parent, queries)),
                              "random t=" + std::to_string(t));
        }
    }

    // ---- 隨機大樹：2000 點、400 詢問 ----
    {
        std::mt19937 rng(777);
        int n = 2000;
        std::vector<int> values(n);
        for (int i = 0; i < n; ++i) values[i] = static_cast<int>(rng() % 51);
        std::vector<int> parent(n, -1);
        for (int v = 1; v < n; ++v) parent[v] = static_cast<int>(rng() % v);
        std::vector<std::pair<int, int>> qs(400);
        for (int i = 0; i < 400; ++i) {
            qs[i] = std::make_pair(static_cast<int>(rng() % n),
                                   static_cast<int>(rng() % n));
        }
        checkAgainstBrute(n, values, parent, qs,
                          parseAnswers(run(n, values, parent, qs)), "big tree");
    }
}

// ---------------------------------------------------------------- IO

bool parseIntToken(const std::string& t, long long& out) {
    // 嚴格匹配 ^[+-]?[0-9]+$
    size_t i = 0;
    if (i < t.size() && (t[i] == '+' || t[i] == '-')) ++i;
    if (i >= t.size()) return false;
    for (; i < t.size(); ++i) {
        if (t[i] < '0' || t[i] > '9') return false;
    }
    errno = 0;
    out = strtoll(t.c_str(), nullptr, 10);
    if (errno == ERANGE) out = (t[0] == '-') ? -CLAMP_ABS : CLAMP_ABS;
    return true;
}

std::vector<std::string> tokenize(const std::string& raw) {
    std::vector<std::string> toks;
    std::string cur;
    for (size_t i = 0; i < raw.size(); ++i) {
        char c = raw[i];
        if (c == ' ' || c == '\n' || c == '\r' || c == '\t' ||
            c == '\f' || c == '\v') {
            if (!cur.empty()) {
                toks.push_back(cur);
                cur.clear();
            }
        } else {
            cur += c;
        }
    }
    if (!cur.empty()) toks.push_back(cur);
    return toks;
}

void runIo(const std::string& raw) {
    std::vector<std::string> toks = tokenize(raw);
    size_t idx = 0;
    // 與 Python 版一致：非數字 token 消費掉並回傳預設值；數值先夾到 ±1e18
    auto nextInt = [&](long long def) -> long long {
        if (idx >= toks.size()) return def;
        std::string t = toks[idx];
        ++idx;
        long long v = 0;
        if (!parseIntToken(t, v)) return def;
        if (v > CLAMP_ABS) v = CLAMP_ABS;
        if (v < -CLAMP_ABS) v = -CLAMP_ABS;
        return v;
    };
    auto cl = [](long long v, long long lo, long long hi) -> long long {
        if (v < lo) return lo;
        if (v > hi) return hi;
        return v;
    };

    long long n = cl(nextInt(0), 0, NMAX);
    long long q = cl(nextInt(0), 0, QMAX);
    std::vector<int> values(static_cast<size_t>(std::max(n, 0LL)), 0);
    for (long long i = 0; i < n; ++i) {
        values[static_cast<size_t>(i)] =
            static_cast<int>(cl(nextInt(0), -VMAX, VMAX));
    }
    std::vector<int> parent(static_cast<size_t>(std::max(n, 0LL)), -1);
    for (long long v = 1; v < n; ++v) {
        long long p = cl(nextInt(1), 1, n);
        // 自環一律當成指向根，保證後續爬行終止
        parent[static_cast<size_t>(v)] =
            (p - 1 != v) ? static_cast<int>(p - 1) : 0;
    }
    std::vector<std::pair<int, int>> queries;
    for (long long i = 0; i < q; ++i) {
        if (n == 0) break;
        long long u = cl(nextInt(1), 1, n) - 1;
        long long v = cl(nextInt(1), 1, n) - 1;
        queries.push_back(std::make_pair(static_cast<int>(u),
                                         static_cast<int>(v)));
    }
    std::cout << run(static_cast<int>(n), values, parent, queries);
}

}  // namespace

int main() {
#ifdef _WIN32
    _setmode(_fileno(stdout), _O_BINARY);
#endif
    std::string raw((std::istreambuf_iterator<char>(std::cin)),
                    std::istreambuf_iterator<char>());
    bool blank = true;
    for (size_t i = 0; i < raw.size(); ++i) {
        unsigned char u = static_cast<unsigned char>(raw[i]);
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




