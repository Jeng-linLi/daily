// 基環樹與函數圖（Functional Graph / Pseudo-Forest）
//
// 題意：
//     本題分兩部分，兩部分都同時給出「高效演算法」與「獨立實作的暴力版本」並互相對拍。
//
//     A. 有向函數圖：n 個點，每個點 i 恰有一條出邊 nxt[i]（0 ≤ nxt[i] < n）。
//        1. 找出所有環，每個環正規化成「從環上最小點號開始」的順序
//        2. belong[i]：i 沿出邊一直走，最終進入的那個環的代表（環上最小點號）
//        3. dist_to_cycle[i]：i 走到環上所需步數（環上點為 0）
//        4. cycle_len[i]：i 最終進入的那個環的長度
//        5. jump(i, K)：從 i 出發走 K 步後落在哪個點（K 可達 10^18，用二進制跳躍）
//
//     B. 無向基環樹森林：n2 個點、m2 條邊，每個連通分量恰有一個環（退化成樹時環長為 0）。
//        對每個分量輸出：最小點號、環長、直徑（任意兩點間最短路徑的最大值）。
//
// 思路：
//     ### A1. 找環：沿著唯一出邊走，一定會「走進一條尾巴，再進入一個環」
//     函數圖的每個連通分量恰有一個環，形狀是「很多棵內向樹掛在環上」。
//     用 0/1/2 三色標記（0 未訪問、1 在當前這條鏈上、2 已處理完）：
//       - 從未訪問點 s 出發一直走，把沿途節點記進 path 與 pos（節點 → 在 path 中的下標）。
//       - 若撞到「顏色 1」的點 v：path[pos[v]:] 就是那個環（v 在當前鏈上出現了兩次）。
//       - 若撞到「顏色 2」的點 v：這條尾巴接到了別人已經處理好的環上。
//     環上的點距離為 0；尾巴上的點用**倒序**傳遞 dist[u] = dist[nxt[u]] + 1，總共 O(n)。
//
//     ### A2. 走 K 步：二進制跳躍（倍增）
//     up[j][i] 表示從 i 出發走 2^j 步到達的點；up[j][i] = up[j-1][ up[j-1][i] ]。
//     把 K 拆成二進位，從低位往高位跳即可。O(n log K) 預處理、O(log K) 單次查詢。
//
//     ### B1. 找環：拓撲剝葉（砍掉所有度數 ≤ 1 的點）
//     反覆取出度數 ≤ 1 的點、從圖上刪掉、把鄰居度數減 1；最後**沒被刪掉的點**就是 2-core，
//     在基環樹裡恰好就是那個唯一的環（樹會被全部刪光，環長 0）。
//     剝葉同時維護 dep[u]（往已刪除那一側的最大深度）與 top1/top2（各分支深度前兩大）：
//         「經過 u 的最長路徑」= top1[u] + top2[u]
//
//     ### B2. 基環樹直徑 = max(每棵掛在環上的樹的直徑, 穿過環的最長路徑)
//     環上點 c_i 掛著深度 val[i] 的樹，環上兩點距離為 min(|i-j|, L-|i-j|)，
//     候選答案 val[i] + val[j] + min(|i-j|, L-|i-j|)。把 val 複製一份接在後面（破環成鏈），
//     對 (val[j]+j)（d ≤ L/2 時）與 (val[j]-j)（d > L/2 時）分別做滑動窗口最大值，O(L) 算完。
//
// 應用場景：
//     - 有向函數圖：置換的循環分解、映射迭代的收斂行為、雜湊碰撞（rho 形路徑）、
//       Pollard-Rho 因數分解的形狀本質上就是函數圖
//     - 二進制跳躍：置換的 k 次冪、自動機跑 k 步的落點、指標跳躍、LCA 倍增同理
//     - 無向基環樹：n 點 n 邊的圖（環狀道路 + 支線）、「基環樹森林直徑」經典題
//
// 複雜度：
//     A 找環                O(n) 時間 / O(n) 空間
//     A 二進制跳躍          O(n log K) 預處理，O(log K) 單次查詢
//     B 剝葉 + 環上 DP      O(n + m + L) 時間 / O(n + m) 空間（L 為環長）

#include <algorithm>
#include <array>
#include <cassert>
#include <cerrno>
#include <cstdio>
#include <cstdlib>
#include <deque>
#include <functional>
#include <iostream>
#include <iterator>
#include <map>
#include <numeric>
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

const int MAXN = 60;                    // 有向函數圖在 IO 模式下的最大點數
const int MAXN2 = 30;                   // 無向圖最大點數
const int MAXM2 = 200;                  // 無向圖最大邊數
const long long KMAX = 1000000000000000000LL;  // 跳躍步數上限（10^18）
const long long K_BRUTE = 3000;         // 超過就不跑暴力模擬
const int LOG = 62;                     // 2^62 > KMAX

// ------------------------------------------------------- A. 有向函數圖

struct FuncInfo {
    std::vector<int> belong;
    std::vector<int> dtc;
    std::vector<int> clen;
    std::vector<std::vector<int>> cycles;
};

FuncInfo analyzeFunctional(const std::vector<int>& nxt) {
    int n = static_cast<int>(nxt.size());
    FuncInfo info;
    info.belong.assign(n, -1);
    info.dtc.assign(n, 0);
    info.clen.assign(n, 0);
    std::vector<int> state(n, 0);       // 0 未訪問 / 1 在當前鏈上 / 2 已處理

    for (int s = 0; s < n; ++s) {
        if (state[s] != 0) continue;
        std::vector<int> path;
        std::vector<int> pos(n, -1);
        int v = s;
        while (state[v] == 0) {
            state[v] = 1;
            pos[v] = static_cast<int>(path.size());
            path.push_back(v);
            v = nxt[v];
        }
        int start;
        if (state[v] == 1) {            // 撞回當前鏈 → 找到新環
            std::vector<int> cyc(path.begin() + pos[v], path.end());
            int L = static_cast<int>(cyc.size());
            int mn = *std::min_element(cyc.begin(), cyc.end());
            int k = static_cast<int>(std::find(cyc.begin(), cyc.end(), mn) - cyc.begin());
            std::vector<int> norm;
            for (int i = 0; i < L; ++i) norm.push_back(cyc[(k + i) % L]);
            info.cycles.push_back(norm);
            for (int u : cyc) {
                info.belong[u] = mn;
                info.dtc[u] = 0;
                info.clen[u] = L;
                state[u] = 2;
            }
            start = pos[v] - 1;
        } else {                        // 接到已處理好的環上
            start = static_cast<int>(path.size()) - 1;
        }
        for (int j = start; j >= 0; --j) {   // 尾巴倒序傳遞
            int w = nxt[path[j]];
            info.belong[path[j]] = info.belong[w];
            info.dtc[path[j]] = info.dtc[w] + 1;
            info.clen[path[j]] = info.clen[w];
            state[path[j]] = 2;
        }
    }
    std::sort(info.cycles.begin(), info.cycles.end(),
              [](const std::vector<int>& a, const std::vector<int>& b) { return a[0] < b[0]; });
    return info;
}

std::vector<int> jumpAll(const std::vector<int>& nxt, long long k) {
    int n = static_cast<int>(nxt.size());
    if (n == 0) return {};
    if (k < 0) k = 0;
    if (k > KMAX) k = KMAX;
    std::vector<std::vector<int>> up(LOG, std::vector<int>(n, 0));
    up[0] = nxt;
    for (int j = 1; j < LOG; ++j) {
        for (int i = 0; i < n; ++i) up[j][i] = up[j - 1][up[j - 1][i]];
    }
    std::vector<int> res(n, 0);
    for (int i = 0; i < n; ++i) {
        int v = i;
        long long x = k;
        int b = 0;
        while (x) {
            if (x & 1LL) v = up[b][v];
            x >>= 1;
            ++b;
        }
        res[i] = v;
    }
    return res;
}

std::vector<int> jumpBrute(const std::vector<int>& nxt, long long k) {
    std::vector<int> res;
    for (int i = 0; i < static_cast<int>(nxt.size()); ++i) {
        int v = i;
        for (long long s = 0; s < k; ++s) v = nxt[v];
        res.push_back(v);
    }
    return res;
}

struct FuncBrute {
    std::vector<int> belong;
    std::vector<int> dtc;
    std::vector<int> clen;
    std::vector<std::pair<int, int>> pairs;   // (環上最小點號, 環長)，已排序
};

FuncBrute functionalBrute(const std::vector<int>& nxt) {
    int n = static_cast<int>(nxt.size());
    FuncBrute out;
    std::set<std::pair<int, int>> ps;
    for (int i = 0; i < n; ++i) {
        std::vector<int> seen(n, -1);
        int v = i, step = 0;
        while (seen[v] == -1) {
            seen[v] = step;
            v = nxt[v];
            ++step;
        }
        int first = seen[v];
        int mn = n, cnt = 0;
        for (int u = 0; u < n; ++u) {
            if (seen[u] >= first) {
                ++cnt;
                if (u < mn) mn = u;
            }
        }
        out.belong.push_back(mn);
        out.dtc.push_back(first);
        out.clen.push_back(step - first);
        ps.insert({mn, step - first});
    }
    out.pairs.assign(ps.begin(), ps.end());
    return out;
}

// ------------------------------------------------------- B. 無向基環樹

std::vector<int> slidingMax(const std::vector<int>& arr, int size) {
    int n = static_cast<int>(arr.size());
    std::vector<int> out;
    if (size <= 0 || size > n) return out;
    std::deque<int> dq;                 // 存下標，對應值遞減
    for (int i = 0; i < n; ++i) {
        while (!dq.empty() && arr[dq.back()] <= arr[i]) dq.pop_back();
        dq.push_back(i);
        while (dq.front() <= i - size) dq.pop_front();
        if (i >= size - 1) out.push_back(arr[dq.front()]);
    }
    return out;
}

int ringDiameterFast(const std::vector<int>& vals) {
    int L = static_cast<int>(vals.size());
    if (L <= 1) return 0;
    std::vector<int> v2(2 * L);
    for (int i = 0; i < 2 * L; ++i) v2[i] = vals[i % L];
    int half = L / 2;
    std::vector<int> keyA(2 * L), keyB(2 * L);
    for (int j = 0; j < 2 * L; ++j) {
        keyA[j] = v2[j] + j;
        keyB[j] = v2[j] - j;
    }
    std::vector<int> winA = slidingMax(keyA, half);          // 窗口 [i+1, i+half]
    int sizeB = L - half - 1;
    std::vector<int> winB = slidingMax(keyB, sizeB);         // 窗口 [i+half+1, i+L-1]
    int best = 0;
    for (int i = 0; i < L; ++i) {
        if (half > 0) {
            int k = i + 1;
            if (k < static_cast<int>(winA.size())) best = std::max(best, vals[i] - i + winA[k]);
        }
        if (sizeB > 0) {
            int k = i + half + 1;
            if (k < static_cast<int>(winB.size())) best = std::max(best, vals[i] + i + L + winB[k]);
        }
    }
    return best;
}

int ringDiameterBrute(const std::vector<int>& vals) {
    int L = static_cast<int>(vals.size());
    int best = 0;
    for (int i = 0; i < L; ++i) {
        for (int j = i + 1; j < L; ++j) {
            int d = j - i;
            best = std::max(best, vals[i] + vals[j] + std::min(d, L - d));
        }
    }
    return best;
}

std::vector<int> ringOrder(const std::vector<std::vector<int>>& adj,
                           const std::vector<char>& on,
                           const std::vector<int>& nodes) {
    std::vector<int> order;
    if (nodes.empty()) return order;
    int start = *std::min_element(nodes.begin(), nodes.end());
    order.push_back(start);
    int prev = -1, cur = start;
    while (static_cast<int>(order.size()) < static_cast<int>(nodes.size())) {
        int nxt = -1;
        for (int w : adj[cur]) {
            if (on[w] && w != prev) {
                if (nxt == -1 || w < nxt) nxt = w;
            }
        }
        if (nxt == -1 || nxt == start) break;
        order.push_back(nxt);
        prev = cur;
        cur = nxt;
    }
    return order;
}

// (最小點號, 環長, 直徑)
std::vector<std::array<int, 3>> analyzePseudoForest(int n, const std::vector<std::pair<int, int>>& edges) {
    std::vector<std::vector<int>> adj(n);
    for (auto& e : edges) {
        adj[e.first].push_back(e.second);
        adj[e.second].push_back(e.first);
    }
    std::vector<int> deg(n);
    for (int i = 0; i < n; ++i) deg[i] = static_cast<int>(adj[i].size());

    std::vector<char> on(n, 1);
    std::vector<int> dep(n, 0), top1(n, 0), top2(n, 0);
    std::deque<int> q;
    for (int i = 0; i < n; ++i) {
        if (deg[i] <= 1) q.push_back(i);
    }
    while (!q.empty()) {
        int u = q.front();
        q.pop_front();
        if (!on[u]) continue;
        on[u] = 0;
        for (int v : adj[u]) {
            if (!on[v]) continue;
            int branch = dep[u] + 1;
            if (branch > top1[v]) {
                top2[v] = top1[v];
                top1[v] = branch;
            } else if (branch > top2[v]) {
                top2[v] = branch;
            }
            if (branch > dep[v]) dep[v] = branch;
            --deg[v];
            if (deg[v] == 1) q.push_back(v);
        }
    }

    std::vector<int> comp(n, -1);
    std::vector<std::vector<int>> comps;
    for (int s = 0; s < n; ++s) {
        if (comp[s] != -1) continue;
        int cid = static_cast<int>(comps.size());
        comp[s] = cid;
        std::vector<int> nodes;
        std::vector<int> st(1, s);
        while (!st.empty()) {
            int x = st.back();
            st.pop_back();
            nodes.push_back(x);
            for (int w : adj[x]) {
                if (comp[w] == -1) {
                    comp[w] = cid;
                    st.push_back(w);
                }
            }
        }
        comps.push_back(nodes);
    }

    std::vector<std::array<int, 3>> out;
    for (auto& nodes : comps) {
        std::vector<int> cyc;
        int inner = 0;
        for (int x : nodes) {
            if (on[x]) cyc.push_back(x);
            inner = std::max(inner, top1[x] + top2[x]);
        }
        int ring = 0;
        if (static_cast<int>(cyc.size()) > 1) {
            std::vector<int> order = ringOrder(adj, on, cyc);
            std::vector<int> vals;
            for (int x : order) vals.push_back(dep[x]);
            ring = ringDiameterFast(vals);
        }
        int mn = *std::min_element(nodes.begin(), nodes.end());
        out.push_back({mn, static_cast<int>(cyc.size()), std::max(inner, ring)});
    }
    std::sort(out.begin(), out.end(),
              [](const std::array<int, 3>& a, const std::array<int, 3>& b) {
                  if (a[0] != b[0]) return a[0] < b[0];
                  if (a[1] != b[1]) return a[1] < b[1];
                  return a[2] < b[2];
              });
    return out;
}

std::vector<std::array<int, 3>> pseudoForestBrute(int n, const std::vector<std::pair<int, int>>& edges) {
    std::vector<std::vector<int>> adj(n);
    for (auto& e : edges) {
        adj[e.first].push_back(e.second);
        adj[e.second].push_back(e.first);
    }
    std::vector<int> comp(n, -1);
    std::vector<std::vector<int>> comps;
    for (int s = 0; s < n; ++s) {
        if (comp[s] != -1) continue;
        int cid = static_cast<int>(comps.size());
        comp[s] = cid;
        std::vector<int> nodes(1, s);
        for (size_t i = 0; i < nodes.size(); ++i) {
            int x = nodes[i];
            for (int w : adj[x]) {
                if (comp[w] == -1) {
                    comp[w] = cid;
                    nodes.push_back(w);
                }
            }
        }
        comps.push_back(nodes);
    }

    // 禁止走 (banU, banV) 這條邊時，x 能否到達 y
    auto reachable = [&](int x, int y, int banU, int banV) {
        std::vector<char> seen(n, 0);
        seen[x] = 1;
        std::vector<int> st(1, x);
        while (!st.empty()) {
            int a = st.back();
            st.pop_back();
            if (a == y) return true;
            for (int b : adj[a]) {
                if ((a == banU && b == banV) || (a == banV && b == banU)) continue;
                if (!seen[b]) {
                    seen[b] = 1;
                    st.push_back(b);
                }
            }
        }
        return false;
    };

    std::vector<std::array<int, 3>> out;
    for (auto& nodes : comps) {
        int cycCnt = 0;
        for (int v : nodes) {
            bool onCycle = false;
            std::set<int> nbrs(adj[v].begin(), adj[v].end());
            for (int u : nbrs) {
                if (u == v) { onCycle = true; break; }        // 自環
                if (reachable(u, v, u, v)) { onCycle = true; break; }
            }
            if (onCycle) ++cycCnt;
        }
        int dia = 0;
        for (int s : nodes) {
            std::vector<int> dist(n, -1);
            dist[s] = 0;
            std::vector<int> stack(1, s);
            for (size_t i = 0; i < stack.size(); ++i) {
                int x = stack[i];
                for (int w : adj[x]) {
                    if (dist[w] == -1) {
                        dist[w] = dist[x] + 1;
                        stack.push_back(w);
                    }
                }
            }
            for (int t : nodes) dia = std::max(dia, dist[t]);
        }
        int mn = *std::min_element(nodes.begin(), nodes.end());
        out.push_back({mn, cycCnt, dia});
    }
    std::sort(out.begin(), out.end(),
              [](const std::array<int, 3>& a, const std::array<int, 3>& b) {
                  if (a[0] != b[0]) return a[0] < b[0];
                  if (a[1] != b[1]) return a[1] < b[1];
                  return a[2] < b[2];
              });
    return out;
}

// ------------------------------------------------------- 輸出

std::string join(const std::vector<int>& v) {
    std::ostringstream oss;
    for (size_t i = 0; i < v.size(); ++i) {
        if (i) oss << " ";
        oss << v[i];
    }
    return oss.str();
}

std::string run(const std::vector<int>& nxt, long long k, int n2,
                const std::vector<std::pair<int, int>>& edges) {
    FuncInfo info = analyzeFunctional(nxt);
    FuncBrute bf = functionalBrute(nxt);
    std::vector<int> jumps = jumpAll(nxt, k);

    int ok = 1;
    if (info.belong != bf.belong || info.dtc != bf.dtc || info.clen != bf.clen) ok = 0;
    std::vector<std::pair<int, int>> algoPairs;
    for (auto& c : info.cycles) algoPairs.push_back({c[0], static_cast<int>(c.size())});
    std::sort(algoPairs.begin(), algoPairs.end());
    if (algoPairs != bf.pairs) ok = 0;
    if (k <= K_BRUTE && jumps != jumpBrute(nxt, k)) ok = 0;
    for (auto& c : info.cycles) {                        // 環的形狀自洽性
        std::set<int> ss(c.begin(), c.end());
        if (ss.size() != c.size() || *std::min_element(c.begin(), c.end()) != c[0]) ok = 0;
        for (size_t idx = 0; idx < c.size(); ++idx) {
            if (nxt[c[idx]] != c[(idx + 1) % c.size()]) ok = 0;
        }
    }

    std::vector<std::array<int, 3>> pf = analyzePseudoForest(n2, edges);
    std::vector<std::array<int, 3>> pbf = pseudoForestBrute(n2, edges);
    if (pf != pbf) ok = 0;

    std::ostringstream oss;
    oss << info.cycles.size() << "\n";
    for (size_t i = 0; i < info.cycles.size(); ++i) {
        if (i) oss << " ";
        oss << info.cycles[i][0] << ":" << info.cycles[i].size();
    }
    oss << "\n";
    oss << join(info.belong) << "\n";
    oss << join(info.dtc) << "\n";
    oss << join(info.clen) << "\n";
    oss << join(jumps) << "\n";
    oss << pf.size() << "\n";
    for (size_t i = 0; i < pf.size(); ++i) {
        if (i) oss << " ";
        oss << pf[i][0] << ":" << pf[i][1] << ":" << pf[i][2];
    }
    oss << "\n";
    oss << ok << "\n";
    return oss.str();
}

// ------------------------------------------------------- 內置測試

std::vector<std::pair<int, int>> genPseudoForest(std::mt19937& rng, int n) {
    std::vector<std::pair<int, int>> edges;
    std::set<std::pair<int, int>> eset;
    std::vector<int> uf(n);
    for (int i = 0; i < n; ++i) uf[i] = i;
    std::function<int(int)> find = [&](int x) {
        while (uf[x] != x) {
            uf[x] = uf[uf[x]];
            x = uf[x];
        }
        return x;
    };
    for (int i = 1; i < n; ++i) {
        if (rng() % 100 < 85) {
            int p = static_cast<int>(rng() % i);
            edges.push_back({p, i});
            eset.insert({p, i});
            eset.insert({i, p});
        }
    }
    for (auto& e : edges) {
        int a = find(e.first), b = find(e.second);
        if (a != b) uf[a] = b;
    }
    std::map<int, std::vector<int>> groups;
    for (int i = 0; i < n; ++i) groups[find(i)].push_back(i);
    for (auto& kv : groups) {
        const std::vector<int>& members = kv.second;
        if (members.size() < 2 || rng() % 100 < 25) continue;   // 保持是樹
        for (int t = 0; t < 20; ++t) {
            int u = members[rng() % members.size()];
            int v = members[rng() % members.size()];
            if (u != v && !eset.count({u, v})) {
                edges.push_back({u, v});
                eset.insert({u, v});
                eset.insert({v, u});
                break;
            }
        }
    }
    return edges;
}

void runTests() {
    // ---- 固定用例：有向函數圖 ----
    {
        std::vector<int> nxt = {1, 2, 0, 1};
        FuncInfo r = analyzeFunctional(nxt);
        assert(r.cycles.size() == 1 && r.cycles[0] == std::vector<int>({0, 1, 2}));
        assert(r.belong == std::vector<int>({0, 0, 0, 0}));
        assert(r.dtc == std::vector<int>({0, 0, 0, 1}));
        assert(r.clen == std::vector<int>({3, 3, 3, 3}));
        assert(jumpAll(nxt, 0) == std::vector<int>({0, 1, 2, 3}));
        assert(jumpAll(nxt, 1) == std::vector<int>({1, 2, 0, 1}));
        assert(jumpAll(nxt, 100) == jumpBrute(nxt, 100));
    }
    {
        std::vector<int> nxt = {1, 2, 3, 4, 2};    // 環 {2,3,4}
        FuncInfo r = analyzeFunctional(nxt);
        assert(r.cycles.size() == 1 && r.cycles[0] == std::vector<int>({2, 3, 4}));
        assert(r.dtc == std::vector<int>({2, 1, 0, 0, 0}));
        assert(r.clen == std::vector<int>({3, 3, 3, 3, 3}));
        assert(r.belong == std::vector<int>({2, 2, 2, 2, 2}));
    }
    {
        FuncInfo r = analyzeFunctional({});
        assert(r.cycles.empty() && r.belong.empty());
        assert(jumpAll({}, 5).empty());
        assert(analyzeFunctional({0}).cycles[0] == std::vector<int>({0}));
        assert(jumpAll({0}, KMAX) == std::vector<int>({0}));
    }
    {   // 純大環：jump 等價於模運算
        std::vector<int> big(50);
        std::vector<int> want(50);
        for (int i = 0; i < 50; ++i) {
            big[i] = (i + 1) % 50;
            want[i] = static_cast<int>((i + KMAX) % 50);
        }
        assert(jumpAll(big, KMAX) == want);
    }

    // ---- 固定用例：無向基環樹 ----
    {
        std::vector<std::pair<int, int>> e = {{0, 1}, {1, 2}, {2, 0}, {0, 3}, {3, 4}};
        std::vector<std::array<int, 3>> got = analyzePseudoForest(5, e);
        std::vector<std::array<int, 3>> want = {{0, 3, 3}};
        assert(got == want);
        assert(got == pseudoForestBrute(5, e));
    }
    {
        std::vector<std::pair<int, int>> e = {{0, 1}, {1, 2}, {2, 3}};
        std::vector<std::array<int, 3>> want = {{0, 0, 3}};
        assert(analyzePseudoForest(4, e) == want);
        assert(analyzePseudoForest(4, e) == pseudoForestBrute(4, e));
    }
    {
        std::vector<std::array<int, 3>> want = {{0, 0, 0}};
        assert(analyzePseudoForest(1, {}) == want);
        assert(analyzePseudoForest(0, {}).empty());
    }
    {
        std::vector<std::pair<int, int>> e = {{0, 1}, {1, 2}, {2, 3}, {3, 0}};
        std::vector<std::array<int, 3>> want = {{0, 4, 2}};
        assert(analyzePseudoForest(4, e) == want);
        assert(analyzePseudoForest(4, e) == pseudoForestBrute(4, e));
    }

    // ---- 環上 DP：滑動窗口版 vs O(L²) 暴力 ----
    {
        std::mt19937 rng0(1);
        for (int it = 0; it < 300; ++it) {
            int L = 1 + static_cast<int>(rng0() % 12);
            std::vector<int> vals(L);
            for (int i = 0; i < L; ++i) vals[i] = static_cast<int>(rng0() % 10);
            assert(ringDiameterFast(vals) == ringDiameterBrute(vals));
        }
    }

    // ---- 隨機對拍 ----
    {
        std::mt19937 rng(20261017);
        for (int it = 0; it < 400; ++it) {
            int n = static_cast<int>(rng() % 21);
            std::vector<int> nxt(n);
            for (int i = 0; i < n; ++i) nxt[i] = static_cast<int>(rng() % (n ? n : 1));
            long long k;
            switch (rng() % 6) {
                case 0: k = 0; break;
                case 1: k = 1; break;
                case 2: k = 2; break;
                case 3: k = 3; break;
                case 4: k = static_cast<int>(rng() % 61); break;
                default: k = static_cast<int>(rng() % 2001); break;
            }
            int n2 = static_cast<int>(rng() % 11);
            std::vector<std::pair<int, int>> edges = genPseudoForest(rng, n2);
            std::string out = run(nxt, k, n2, edges);
            assert(out.size() >= 2 && out[out.size() - 2] == '1');
        }
    }

    // ---- 長鏈只跑線性演算法 ----
    {
        std::vector<int> chain(2000);
        for (int i = 0; i < 1999; ++i) chain[i] = i + 1;
        chain[1999] = 1999;
        FuncInfo r = analyzeFunctional(chain);
        assert(r.cycles.size() == 1 && r.cycles[0] == std::vector<int>({1999}));
        assert(r.dtc[0] == 1999);
        assert(r.clen[0] == 1);
        assert(jumpAll(chain, KMAX)[0] == 1999);
    }
}

// ------------------------------------------------------- IO

bool isInt(const std::string& t) {
    if (t.empty()) return false;
    size_t i = 0;
    if (t[0] == '+' || t[0] == '-') i = 1;
    if (i >= t.size()) return false;
    for (; i < t.size(); ++i) {
        if (t[i] < '0' || t[i] > '9') return false;
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

    int n = static_cast<int>(nextInt(0));
    if (n < 0) n = 0;
    if (n > MAXN) n = MAXN;
    std::vector<int> nxt(n, 0);
    for (int i = 0; i < n; ++i) {
        long long v = nextInt(0);
        if (v < 0 || v >= n) v = 0;
        nxt[i] = static_cast<int>(v);
    }
    long long k = nextInt(0);
    if (k < 0) k = 0;
    int n2 = static_cast<int>(nextInt(0));
    if (n2 < 0) n2 = 0;
    if (n2 > MAXN2) n2 = MAXN2;
    int m2 = static_cast<int>(nextInt(0));
    if (m2 < 0) m2 = 0;
    if (m2 > MAXM2) m2 = MAXM2;
    std::vector<std::pair<int, int>> edges;
    for (int i = 0; i < m2; ++i) {
        long long u = nextInt(0), v = nextInt(0);
        if (u >= 0 && u < n2 && v >= 0 && v < n2) edges.push_back({static_cast<int>(u), static_cast<int>(v)});
    }
    std::cout << run(nxt, k, n2, edges);
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
        if (u > ' ') { blank = false; break; }
    }
    if (!blank) {
        runIo(raw);
    } else {
        runTests();
        std::cout << "all tests passed" << std::endl;
    }
    return 0;
}
