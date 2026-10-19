// 最大團與 Bron–Kerbosch 演算法（含貪心著色上界與補圖獨立集）
//
// 題意：
//     給一張 n 點的無向簡單圖（鄰接矩陣），求：
//
//     1. **最大團大小**（maximum clique）：最大的「兩兩相鄰」的點集有多大
//     2. **極大團的數量**（maximal cliques）：不能再加點的團有幾個
//     3. **字典序最小的那個最大團**的點集合
//     4. **最大獨立集**與**最小點覆蓋**
//     5. 全部結果都與「枚舉 2^n 個子集」的暴力版對拍
//
// 思路：
//     ### 1. Bron–Kerbosch：列舉所有極大團
//     維護三個集合：
//         R = 目前這條搜尋路徑上已經選進團裡的點
//         P = 還能加進來的候選點（與 R 中所有點都相鄰）
//         X = 已處理過、用來避免重複列舉同一個極大團的點
//     當 P 與 X 同時為空，R 就是一個**極大團**（加不進任何點了）。
//
//     ### 2. 樞軸（pivot）剪枝
//     沒有樞軸的版本會對同一個極大團重複搜尋很多次。
//     取 u ∈ P ∪ X 中「|P ∩ N(u)| 最大」的那個當樞軸（平手取編號最小的）：
//     任何極大團要嘛包含 u、要嘛包含某個「與 u 不相鄰」的點，
//     所以只需遞迴 P \ N(u) 裡的點，候選數目明顯變小，且不漏解。
//
//     ### 3. 貪心著色上界（branch and bound）
//     只求「大小」不求「列舉」時可以剪得更兇。對候選集合 P 做一次貪心著色：
//     依序給每個點塗上「鄰居沒用過的最小顏色」，得到 color[v]。
//     因為同一個團裡任兩點相鄰，顏色必互不相同，所以
//         |R| + color[v] ≤ 目前最佳解  →  從 v 出發不可能更好，直接剪掉
//     這就是 Tomita 的 MCQ 演算法，實務上比純樞軸版快非常多。
//     注意：貪心著色**不保證** color 隨下標遞增，上界必須取「前綴最大值」，
//     否則會把還有機會的分支誤剪掉（本題開發時就踩過這個坑）。
//
//     ### 4. 最大獨立集 = 補圖的最大團
//     原圖的獨立集（兩兩**不**相鄰）正好是補圖的團（兩兩相鄰）。
//     而對任意圖都有「α + τ = n」（最大獨立集 + 最小點覆蓋 = 點數），
//     因為點覆蓋的補集一定是獨立集，反之亦然。
//
//     ### 5. 跨語言一致的排序
//     著色的順序、樞軸的平手規則、遞迴枚舉順序全部寫死成「編號由小到大」，
//     否則 Python 與 C++ 的集合迭代順序不同，剪枝量就不同，輸出對不起來。
//
// 應用場景：
//     - **社群網路**：找出「一群互相認識的人」（clique 就是最緊密的社群結構）
//     - **生物資訊**：蛋白質交互作用網路裡找功能模組；基因共表現聚類
//     - **編譯器**：暫存器分配時的干涉圖、指令排程的衝突圖
//     - **化學資訊學**：用最大團做分子結構比對（兩個分子圖的最大公共子結構）
//     - 排班 / 會議安排：獨立集 = 可以同時進行的任務集合
//
// 複雜度：
//     樞軸版列舉所有極大團   O(3^{n/3}) 時間（最壞情況，且是緊的）
//     著色上界版求最大團     最壞同階，實務上遠快
//     空間                   O(n²) 鄰接矩陣 + O(n) 遞迴深度
//     暴力版（n ≤ 14）       O(2^n · n²) 時間

#include <algorithm>
#include <cerrno>
#include <cstdlib>
#include <functional>
#include <iostream>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#ifdef _WIN32
#include <fcntl.h>
#include <io.h>
#endif

namespace {

const int MAXN = 20;         // IO 模式下的最大點數
const int MAXM = 190;        // 最大邊數（簡單圖上限 n(n-1)/2）
const int BRUTE_MAX = 14;    // 頂點數不超過這個值才跑暴力對拍

// ---------------------------------------------------------------- 工具

std::vector<std::vector<int>> buildAdj(int n, const std::vector<std::pair<int, int>>& edges) {
    std::vector<std::vector<int>> adj(n, std::vector<int>(n, 0));
    for (const auto& e : edges) {
        int u = e.first, v = e.second;
        if (0 <= u && u < n && 0 <= v && v < n && u != v) {
            adj[u][v] = 1;
            adj[v][u] = 1;
        }
    }
    return adj;
}

std::vector<std::vector<int>> complement(const std::vector<std::vector<int>>& adj) {
    int n = static_cast<int>(adj.size());
    std::vector<std::vector<int>> c(n, std::vector<int>(n, 0));
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            c[i][j] = (i != j && !adj[i][j]) ? 1 : 0;
        }
    }
    return c;
}

// ------------------------------------------------- 演算法一：樞軸版（列舉極大團）

struct BkRes {
    int best;
    long long count;
    std::vector<int> bestSet;
};

BkRes bronKerbosch(const std::vector<std::vector<int>>& adj) {
    int n = static_cast<int>(adj.size());
    BkRes res;
    res.best = 0;
    res.count = 0;

    std::vector<int> r;
    std::vector<int> p;
    std::vector<int> x;
    for (int i = 0; i < n; ++i) p.push_back(i);

    // 用顯式遞迴 lambda；每層的 p / x 都是複本，語意與 Python 版一致
    std::function<void(std::vector<int>, std::vector<int>, std::vector<int>)> rec;
    rec = [&](std::vector<int> rr, std::vector<int> pp, std::vector<int> xx) {
        if (pp.empty() && xx.empty()) {
            ++res.count;
            std::vector<int> cand = rr;
            std::sort(cand.begin(), cand.end());
            if (static_cast<int>(cand.size()) > res.best) {
                res.best = static_cast<int>(cand.size());
                res.bestSet = cand;
            } else if (static_cast<int>(cand.size()) == res.best && res.best > 0) {
                if (res.bestSet.empty() || cand < res.bestSet) res.bestSet = cand;
            }
            return;
        }
        int pu = -1;
        int bestCover = -1;
        for (int u : pp) {
            int c = 0;
            for (int w : pp) {
                if (adj[u][w]) ++c;
            }
            if (c > bestCover || (c == bestCover && (pu < 0 || u < pu))) {
                bestCover = c;
                pu = u;
            }
        }
        for (int u : xx) {
            int c = 0;
            for (int w : pp) {
                if (adj[u][w]) ++c;
            }
            if (c > bestCover || (c == bestCover && (pu < 0 || u < pu))) {
                bestCover = c;
                pu = u;
            }
        }
        if (pu < 0) return;

        std::vector<int> cand;
        for (int v : pp) {
            if (!adj[pu][v]) cand.push_back(v);
        }
        for (int v : cand) {
            std::vector<int> np;
            for (int w : pp) {
                if (adj[v][w]) np.push_back(w);
            }
            std::vector<int> nx;
            for (int w : xx) {
                if (adj[v][w]) nx.push_back(w);
            }
            std::vector<int> nr = rr;
            nr.push_back(v);
            rec(nr, np, nx);
            // 把 v 從 P 移到 X
            std::vector<int> pp2;
            for (int w : pp) {
                if (w != v) pp2.push_back(w);
            }
            pp = pp2;
            xx.push_back(v);
        }
    };

    rec(r, p, x);
    return res;
}

// ------------------------------------------------- 演算法二：著色上界版

std::vector<int> greedyColor(const std::vector<std::vector<int>>& adj,
                             const std::vector<int>& p) {
    int n = static_cast<int>(adj.size());
    std::vector<int> colorOf(n, -1);
    std::vector<int> out(p.size(), 0);
    int m = static_cast<int>(p.size());
    for (int idx = 0; idx < m; ++idx) {
        int v = p[idx];
        std::vector<int> used(m + 2, 0);
        for (int u : p) {
            if (colorOf[u] > 0 && adj[v][u]) used[colorOf[u]] = 1;
        }
        int c = 1;
        while (c <= m + 1 && used[c]) ++c;
        colorOf[v] = c;
        out[idx] = c;
    }
    return out;
}

int maxCliqueColor(const std::vector<std::vector<int>>& adj) {
    int n = static_cast<int>(adj.size());
    int best = 0;

    std::function<void(int, std::vector<int>, std::vector<int>)> expand;
    expand = [&](int cur, std::vector<int> p, std::vector<int> colors) {
        if (p.empty()) {
            if (cur > best) best = cur;
            return;
        }
        // 貪心著色不保證 colors 隨下標遞增 → 上界取前綴最大值
        std::vector<int> pmax(p.size(), 0);
        int acc = 0;
        for (size_t i = 0; i < colors.size(); ++i) {
            if (colors[i] > acc) acc = colors[i];
            pmax[i] = acc;
        }
        for (int i = static_cast<int>(p.size()) - 1; i >= 0; --i) {
            if (cur + pmax[i] <= best) return;
            if (cur + static_cast<int>(p.size()) <= best) return;
            int v = p[i];
            std::vector<int> np;
            for (int j = 0; j < i; ++j) {
                if (adj[v][p[j]]) np.push_back(p[j]);
            }
            std::vector<int> ncolors = greedyColor(adj, np);
            expand(cur + 1, np, ncolors);
            p.erase(p.begin() + i);
        }
    };

    std::vector<int> all;
    for (int i = 0; i < n; ++i) all.push_back(i);
    expand(0, all, greedyColor(adj, all));
    return best;
}

// ------------------------------------------------- 演算法三：暴力枚舉子集

struct BruteRes {
    int best;
    long long cntMaximal;
    std::vector<int> bestSet;
    int bestInd;
};

BruteRes bruteAll(const std::vector<std::vector<int>>& adj) {
    int n = static_cast<int>(adj.size());
    BruteRes res;
    res.best = -1;
    res.cntMaximal = -1;
    res.bestInd = -1;
    if (n > BRUTE_MAX) return res;
    res.best = 0;
    res.cntMaximal = 0;
    res.bestInd = 0;

    long long total = 1LL << n;
    for (long long mask = 0; mask < total; ++mask) {
        std::vector<int> verts;
        for (int i = 0; i < n; ++i) {
            if ((mask >> i) & 1) verts.push_back(i);
        }
        bool isClique = true;
        for (size_t a = 0; a < verts.size() && isClique; ++a) {
            for (size_t b = a + 1; b < verts.size(); ++b) {
                if (!adj[verts[a]][verts[b]]) {
                    isClique = false;
                    break;
                }
            }
        }
        if (isClique) {
            if (static_cast<int>(verts.size()) > res.best) {
                res.best = static_cast<int>(verts.size());
                res.bestSet = verts;
            } else if (static_cast<int>(verts.size()) == res.best && verts < res.bestSet) {
                res.bestSet = verts;
            }
            bool maximal = true;
            for (int u = 0; u < n && maximal; ++u) {
                if ((mask >> u) & 1) continue;
                bool okAll = true;
                for (int w : verts) {
                    if (!adj[u][w]) {
                        okAll = false;
                        break;
                    }
                }
                if (okAll) maximal = false;
            }
            if (maximal) ++res.cntMaximal;
        }
        bool isInd = true;
        for (size_t a = 0; a < verts.size() && isInd; ++a) {
            for (size_t b = a + 1; b < verts.size(); ++b) {
                if (adj[verts[a]][verts[b]]) {
                    isInd = false;
                    break;
                }
            }
        }
        if (isInd && static_cast<int>(verts.size()) > res.bestInd) {
            res.bestInd = static_cast<int>(verts.size());
        }
    }
    return res;
}

// ---------------------------------------------------------------- 主流程

std::string run(int n, const std::vector<std::pair<int, int>>& edges) {
    std::vector<std::vector<int>> adj = buildAdj(n, edges);
    bool ok = true;

    BkRes bk = bronKerbosch(adj);
    int mc2 = maxCliqueColor(adj);
    std::vector<std::vector<int>> comp = complement(adj);
    int mis = maxCliqueColor(comp);
    int mvc = n - mis;

    BruteRes br = bruteAll(adj);
    if (br.best != -1) {
        if (bk.best != br.best || mc2 != br.best || bk.count != br.cntMaximal) ok = false;
        if (bk.bestSet != br.bestSet) ok = false;
        if (mis != br.bestInd) ok = false;
    }

    if (bk.best != mc2) ok = false;
    if (!(0 <= mis && mis <= n) || !(0 <= mvc && mvc <= n)) ok = false;
    if (mis + mvc != n) ok = false;
    if (bronKerbosch(comp).best != mis) ok = false;

    std::ostringstream oss;
    oss << bk.best << ' ' << mc2 << ' ' << br.best << '\n';
    oss << bk.count << ' ' << br.cntMaximal << '\n';
    oss << bk.bestSet.size();
    for (int v : bk.bestSet) oss << ' ' << v;
    oss << '\n';
    oss << mis << ' ' << mvc << ' ' << br.bestInd << '\n';
    oss << (ok ? '1' : '0') << '\n';
    return oss.str();
}

// ---------------------------------------------------------------- 內置測試

struct Lcg {
    unsigned long long s;
    explicit Lcg(unsigned long long seed) : s(seed) {}
    unsigned long long next() {
        s = s * 6364136223846793005ULL + 1442695040888963407ULL;
        return (s >> 32) & 0xFFFFFFFFULL;
    }
    long long below(long long n) {
        return static_cast<long long>(next() % static_cast<unsigned long long>(n));
    }
    double unit() { return static_cast<double>(next() & 0xFFFFFFULL) / 16777216.0; }
};

void expectEq(int n, const std::vector<std::pair<int, int>>& edges, const std::string& want,
              int line) {
    std::string got = run(n, edges);
    if (got != want) {
        std::cout << "FAIL line " << line << "\n--- got ---\n"
                  << got << "--- want ---\n" << want << std::endl;
        std::exit(1);
    }
}

void runTests() {
    expectEq(0, {}, "0 0 0\n1 1\n0\n0 0 0\n1\n", 1);
    expectEq(1, {}, "1 1 1\n1 1\n1 0\n1 0 1\n1\n", 2);
    expectEq(3, {{0, 1}, {1, 2}, {0, 2}}, "3 3 3\n1 1\n3 0 1 2\n1 2 1\n1\n", 3);
    expectEq(3, {{0, 1}, {1, 2}}, "2 2 2\n2 2\n2 0 1\n2 1 2\n1\n", 4);
    expectEq(3, {}, "1 1 1\n3 3\n1 0\n3 0 3\n1\n", 5);
    expectEq(4, {{0, 1}, {1, 2}, {0, 2}, {2, 3}}, "3 3 3\n2 2\n3 0 1 2\n2 2 2\n1\n", 6);

    // ---- 隨機對拍 ----
    {
        Lcg rng(20261019ULL);
        for (int t = 0; t < 400; ++t) {
            int n = static_cast<int>(rng.below(15));
            double dens = (rng.next() % 3 == 0) ? 0.2 : ((rng.next() % 2 == 0) ? 0.5 : 0.8);
            std::vector<std::pair<int, int>> edges;
            for (int i = 0; i < n; ++i) {
                for (int j = i + 1; j < n; ++j) {
                    if (rng.unit() < dens) edges.push_back({i, j});
                }
            }
            std::string out = run(n, edges);
            size_t pos = out.rfind('\n', out.size() - 2);
            std::string last = (pos == std::string::npos) ? out
                                                          : out.substr(pos + 1, out.size() - 2 - pos);
            if (last != "1") {
                std::cout << "FAIL random t=" << t << " n=" << n << "\n" << out << std::endl;
                std::exit(1);
            }
        }
    }

    // ---- 較大的圖（不跑暴力）----
    {
        Lcg rng(555ULL);
        for (int t = 0; t < 120; ++t) {
            int n = static_cast<int>(15 + rng.below(MAXN - 14));
            std::vector<std::pair<int, int>> edges;
            for (int i = 0; i < n; ++i) {
                for (int j = i + 1; j < n; ++j) {
                    if (rng.unit() < 0.3) edges.push_back({i, j});
                }
            }
            std::string out = run(n, edges);
            size_t pos = out.rfind('\n', out.size() - 2);
            std::string last = (pos == std::string::npos) ? out
                                                          : out.substr(pos + 1, out.size() - 2 - pos);
            if (last != "1") {
                std::cout << "FAIL big t=" << t << " n=" << n << "\n" << out << std::endl;
                std::exit(1);
            }
        }
    }

    // ---- 完全圖 Kn ----
    for (int n = 0; n <= 8; ++n) {
        std::vector<std::pair<int, int>> edges;
        for (int i = 0; i < n; ++i) {
            for (int j = i + 1; j < n; ++j) edges.push_back({i, j});
        }
        std::string out = run(n, edges);
        std::ostringstream want1;
        int mc = (n == 0) ? 0 : n;
        want1 << mc << ' ' << mc << ' ' << mc << '\n';
        if (out.compare(0, want1.str().size(), want1.str()) != 0) {
            std::cout << "FAIL Kn n=" << n << "\n" << out << std::endl;
            std::exit(1);
        }
        size_t pos = out.find('\n');
        std::string second = out.substr(pos + 1, out.find('\n', pos + 1) - pos - 1);
        if (second != "1 1") {
            std::cout << "FAIL Kn count n=" << n << " -> " << second << std::endl;
            std::exit(1);
        }
    }
}

// ---------------------------------------------------------------- IO

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

    long long n = nextInt(0);
    if (n < 0) n = 0;
    if (n > MAXN) n = MAXN;
    long long m = nextInt(0);
    if (m < 0) m = 0;
    if (m > MAXM) m = MAXM;
    std::vector<std::pair<int, int>> edges;
    for (long long i = 0; i < m; ++i) {
        long long u = nextInt(0);
        long long v = nextInt(0);
        if (0 <= u && u < n && 0 <= v && v < n && u != v) {
            edges.push_back({static_cast<int>(u), static_cast<int>(v)});
        }
    }
    std::cout << run(static_cast<int>(n), edges);
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
