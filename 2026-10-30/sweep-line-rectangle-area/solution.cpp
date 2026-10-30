// 掃描線 + 線段樹：求 n 個軸對齊矩形的聯集面積（Union Area）
//
// 題意：
//     給定 n 個**軸對齊**（axis-aligned）的矩形，第 i 個矩形的左下角與右上角
//     分別是 (x1[i], y1[i]) 與 (x2[i], y2[i])，其中 x1[i] < x2[i]、y1[i] < y2[i]。
//     求這 n 個矩形**聯集**（union）的面積，也就是把重疊部分只算一次的總面積。
//     限制：1 <= n <= 100000，|座標| <= 10^9。
//
// 思路：
//     ### 1. 為什麼不能直接加總
//     矩形彼此可能重疊，直接把面積相加會重複計算。
//     樸素做法是「格點化」：把所有 x、y 座標收集起來排序，
//     把平面切成 O(n^2) 個小格子再逐格判斷 —— 格子數太多，O(n^2) 不可行。
//
//     ### 2. 掃描線（sweep line）
//     讓一條垂直線從左往右掃過 x 軸。
//     對每個矩形 (x1, y1, x2, y2) 產生兩個事件：
//         (x1, +1, [y1, y2))  —— 進入矩形，把 y 區間 [y1, y2) 的覆蓋次數 +1
//         (x2, -1, [y1, y2))  —— 離開矩形，把 y 區間 [y1, y2) 的覆蓋次數 -1
//     相鄰兩個事件的 x 之間，覆蓋情況完全不變 —— 這就是一個「條帶」（slab）。
//     若條帶寬為 dx、此時 y 方向的聯集長度為 L，則該條帶的面積 = dx * L。
//     把所有條帶的面積加起來就是答案。
//     所以問題變成：**動態維護一組 y 區間的聯集長度**。
//
//     ### 3. 座標壓縮 + 線段樹
//     把所有出現過的 y 座標排序去重，得到 ys[0..m-1]。
//     線段樹的葉子不是「點」而是「基本區間」[ys[i], ys[i+1])，共 m-1 個。
//     節點維護兩個量：
//         cnt   —— 完整覆蓋此區間的「矩形數」（不下傳的覆蓋次數）
//         length—— 此區間內目前被覆蓋的總長度
//     更新規則（經典的「覆蓋次數」寫法）：
//         若 cnt > 0：length = 區間的實際長度（整段被蓋住）
//         否則若為葉子：length = 0
//         否則：length = 左子 + 右子
//     根節點的 length 就是當前整個 y 軸的聯集長度。
//
//     ### 4. 為什麼不需要 push-down
//     這裡刻意不做 lazy 下傳：cnt 記錄的是「有多少個矩形完整覆蓋了這個節點」，
//     而 deeper 的 cnt 記錄的是「在父節點未被完整覆蓋的前提下，
//     額外覆蓋此區間的矩形數」。當 cnt 由正變 0 時，
//     子節點的 length 早已包含了所有下層更新，直接相加即可。
//     這是矩形聯集面積的標準寫法。
//
//     ### 5. 精度與整數性
//     所有座標都是整數，條帶寬度與聯集長度也都是整數，
//     所以面積一定是整數 —— 直接用整數運算即可，不需要浮點數。
//
// 複雜度：
//     座標壓縮 + 排序    O(n log n)
//     2n 次線段樹更新    O(n log n)
//     總計              O(n log n)，空間 O(n)
//
// 應用場景（Applications）：
//     - **地圖 / GIS 覆蓋分析**：多個行政區、地塊、保護區的聯集面積，
//       是選址分析（site selection）的基本運算
//     - **資料中心機櫃佔地**：把機櫃視為矩形，求機房內被機櫃佔用的總面積
//     - **建築物輪廓**：多棟建物平面圖的聯集，估算基地使用率（容積率、覆蓋率）
//     - **影像處理 / UI 渲染**：多個有透明重疊的圖層合成後的實際覆蓋像素數
//     - **碰撞與覆蓋檢測**：機械臂多臂在平面上的作業區域聯集
//     - **廣告位出租**：一塊版面被多個廣告區塊佔用時的可售面積

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <map>
#include <optional>
#include <string>
#include <tuple>
#include <vector>

#ifdef _WIN32
#include <fcntl.h>
#include <io.h>
#endif

namespace {

constexpr int NMAX = 100000;
constexpr long long COORD_MAX = 1000000000LL;
constexpr long long PARSE_LIMIT = 1000000000000000000LL;

struct Rect {
    long long x1 = 0, y1 = 0, x2 = 0, y2 = 0;
};

// y 軸上的「覆蓋次數 + 聯集長度」線段樹。
// 葉子 i 代表基本區間 [ys[i], ys[i+1])，共 n = ys.size() - 1 個。
class SegmentTree {
public:
    SegmentTree(const std::vector<long long>& ys) : ys_(ys), n_(static_cast<int>(ys.size()) - 1) {
        int size = 4 * (n_ + 2);
        cnt_.assign(size, 0);
        length_.assign(size, 0);
    }

    void update(int ql, int qr, int delta) {
        if (n_ <= 0 || ql >= qr) return;
        rec(1, 0, n_, ql, qr, delta);
    }

    long long coveredLength() const { return length_[1]; }

private:
    const std::vector<long long>& ys_;
    int n_ = 0;
    std::vector<int> cnt_;
    std::vector<long long> length_;

    void rec(int node, int l, int r, int ql, int qr, int delta) {
        if (ql <= l && r <= qr) {
            cnt_[node] += delta;
        } else {
            int mid = (l + r) >> 1;
            if (ql < mid) rec(node << 1, l, mid, ql, qr, delta);
            if (qr > mid) rec((node << 1) | 1, mid, r, ql, qr, delta);
        }
        if (cnt_[node] > 0) {
            length_[node] = ys_[r] - ys_[l];
        } else if (r - l == 1) {
            length_[node] = 0;
        } else {
            length_[node] = length_[node << 1] + length_[(node << 1) | 1];
        }
    }
};

// 掃描線 + 線段樹求聯集面積
long long unionArea(const std::vector<Rect>& rects) {
    if (rects.empty()) return 0;
    std::vector<long long> ys;
    ys.reserve(rects.size() * 2);
    // 事件：(x, delta, y1, y2)
    std::vector<std::tuple<long long, int, long long, long long>> events;
    events.reserve(rects.size() * 2);
    for (const Rect& r : rects) {
        ys.push_back(r.y1);
        ys.push_back(r.y2);
        events.emplace_back(r.x1, 1, r.y1, r.y2);
        events.emplace_back(r.x2, -1, r.y1, r.y2);
    }
    std::sort(ys.begin(), ys.end());
    ys.erase(std::unique(ys.begin(), ys.end()), ys.end());
    std::map<long long, int> yIndex;
    for (std::size_t i = 0; i < ys.size(); ++i) {
        yIndex[ys[i]] = static_cast<int>(i);
    }
    std::sort(events.begin(), events.end());

    SegmentTree tree(ys);
    long long area = 0;
    long long prevX = std::get<0>(events[0]);
    std::size_t i = 0;
    while (i < events.size()) {
        long long x = std::get<0>(events[i]);
        // 條帶 [prevX, x) 的面積
        area += (x - prevX) * tree.coveredLength();
        // 同一個 x 的所有事件一起處理
        while (i < events.size() && std::get<0>(events[i]) == x) {
            int delta = std::get<1>(events[i]);
            long long y1 = std::get<2>(events[i]);
            long long y2 = std::get<3>(events[i]);
            tree.update(yIndex[y1], yIndex[y2], delta);
            ++i;
        }
        prevX = x;
    }
    return area;
}

// 參考實作：座標壓縮後逐格累加（O(n^2)，只給小測試用）
long long bruteArea(const std::vector<Rect>& rects) {
    if (rects.empty()) return 0;
    std::vector<long long> xs, ys;
    for (const Rect& r : rects) {
        xs.push_back(r.x1);
        xs.push_back(r.x2);
        ys.push_back(r.y1);
        ys.push_back(r.y2);
    }
    std::sort(xs.begin(), xs.end());
    xs.erase(std::unique(xs.begin(), xs.end()), xs.end());
    std::sort(ys.begin(), ys.end());
    ys.erase(std::unique(ys.begin(), ys.end()), ys.end());
    long long area = 0;
    for (std::size_t i = 0; i + 1 < xs.size(); ++i) {
        // 取區間內的一個整數代表點（用下取整避免負數時偏到右端）
        long long xm = xs[i] + (xs[i + 1] - xs[i]) / 2;
        for (std::size_t j = 0; j + 1 < ys.size(); ++j) {
            long long ym = ys[j] + (ys[j + 1] - ys[j]) / 2;
            for (const Rect& r : rects) {
                if (r.x1 <= xm && xm < r.x2 && r.y1 <= ym && ym < r.y2) {
                    area += (xs[i + 1] - xs[i]) * (ys[j + 1] - ys[j]);
                    break;
                }
            }
        }
    }
    return area;
}

// 簡易 LCG，兩種語言共用同一組常數以確保亂數序列一致
struct Lcg {
    unsigned long long s;
    explicit Lcg(unsigned long long seed) : s(seed) {}
    unsigned long long next() {
        s = (s * 6364136223846793005ULL + 1442695040888963407ULL) &
            0xFFFFFFFFFFFFFFFFULL;
        return s;
    }
    int range(int lo, int hi) {  // [lo, hi]
        return lo + static_cast<int>(
                          next() % static_cast<unsigned long long>(hi - lo + 1));
    }
};

void fail(const char* msg) {
    std::cerr << msg << "\n";
    std::exit(1);
}

void runTests() {
    // ---- 固定用例 ----
    if (unionArea({{0, 0, 4, 3}}) != 12) fail("case 1 failed");
    if (unionArea({{0, 0, 2, 2}, {5, 0, 7, 2}}) != 8) fail("case 2 failed");
    if (unionArea({{0, 0, 2, 2}, {1, 0, 3, 2}}) != 6) fail("case 3 failed");
    if (unionArea({{0, 0, 4, 4}, {1, 1, 3, 3}}) != 16) fail("case 4 failed");
    if (unionArea({{0, 0, 2, 2}, {0, 0, 2, 2}}) != 4) fail("case 5 failed");
    if (unionArea({{0, 1, 4, 3}, {1, 0, 3, 4}}) != 12) fail("case 6 failed");
    if (unionArea({{-2, -2, 2, 2}}) != 16) fail("case 7 failed");
    if (unionArea({{-2, -2, 0, 0}, {0, 0, 2, 2}}) != 8) fail("case 8 failed");
    if (unionArea({{0, 0, 2, 2}, {2, 0, 4, 2}}) != 8) fail("case 9 failed");
    if (unionArea({{0, 0, 2, 2}, {1, 1, 3, 3}, {2, 2, 4, 4}}) != 10) {
        fail("case 10 failed");
    }

    // ---- 隨機對拍 ----
    {
        Lcg rng(20261030ULL);
        for (int iter = 0; iter < 200; ++iter) {
            int n = rng.range(1, 6);
            std::vector<Rect> rects;
            for (int i = 0; i < n; ++i) {
                long long x1 = rng.range(-8, 8);
                long long y1 = rng.range(-8, 8);
                rects.push_back({x1, y1, x1 + rng.range(1, 6), y1 + rng.range(1, 6)});
            }
            if (unionArea(rects) != bruteArea(rects)) fail("random mismatch");
        }
    }
    {
        Lcg rng(99887766ULL);
        for (int iter = 0; iter < 20; ++iter) {
            int n = rng.range(5, 12);
            std::vector<Rect> rects;
            for (int i = 0; i < n; ++i) {
                long long x1 = rng.range(-20, 20);
                long long y1 = rng.range(-20, 20);
                rects.push_back({x1, y1, x1 + rng.range(1, 10), y1 + rng.range(1, 10)});
            }
            if (unionArea(rects) != bruteArea(rects)) fail("mid random mismatch");
        }
    }
    // 純網格
    {
        std::vector<Rect> rects;
        for (int i = 0; i < 10; ++i) {
            for (int j = 0; j < 10; ++j) {
                rects.push_back({i, j, i + 1, j + 1});
            }
        }
        if (unionArea(rects) != 100) fail("grid failed");
    }
    // 一個大矩形 + 許多小矩形
    {
        std::vector<Rect> rects;
        rects.push_back({0, 0, 100, 100});
        for (int i = 0; i < 30; ++i) {
            rects.push_back({i * 3, i * 3, i * 3 + 2, i * 3 + 2});
        }
        if (unionArea(rects) != 10000) fail("big rect failed");
    }
}

// ---------------------------------------------------------------- IO 模式

std::optional<long long> parseInt(const std::string& tok) {
    if (tok.empty()) return std::nullopt;
    std::size_t i = 0;
    if (tok[0] == '+' || tok[0] == '-') i = 1;
    if (i == tok.size()) return std::nullopt;
    for (std::size_t j = i; j < tok.size(); ++j) {
        if (tok[j] < '0' || tok[j] > '9') return std::nullopt;
    }
    // 用 __int128 累加避免溢出；超過 long long 視為非法
    __int128 val = 0;
    bool neg = (tok[0] == '-');
    for (std::size_t j = i; j < tok.size(); ++j) {
        val = val * 10 + (tok[j] - '0');
        if (val > (static_cast<__int128>(1) << 100)) return std::nullopt;
    }
    if (neg) val = -val;
    if (val > static_cast<__int128>(PARSE_LIMIT)) return std::nullopt;
    if (val < -static_cast<__int128>(PARSE_LIMIT)) return std::nullopt;
    return static_cast<long long>(val);
}

void runIo(const std::string& raw) {
    std::vector<std::string> toks;
    {
        std::size_t i = 0;
        while (i < raw.size()) {
            while (i < raw.size() && static_cast<unsigned char>(raw[i]) <= ' ') ++i;
            if (i >= raw.size()) break;
            std::size_t j = i;
            while (j < raw.size() && static_cast<unsigned char>(raw[j]) > ' ') ++j;
            toks.push_back(raw.substr(i, j - i));
            i = j;
        }
    }
    std::size_t pos = 0;
    auto nxt = [&toks, &pos]() -> const std::string* {
        if (pos >= toks.size()) return nullptr;
        return &toks[pos++];
    };
    auto bad = []() { std::cout << "error\n"; };

    const std::string* pt = nxt();
    std::optional<long long> nOpt = std::nullopt;
    if (pt) nOpt = parseInt(*pt);
    if (!nOpt || *nOpt < 1 || *nOpt > NMAX) {
        bad();
        return;
    }
    int n = static_cast<int>(*nOpt);
    std::vector<Rect> rects;
    for (int i = 0; i < n; ++i) {
        long long vals[4];
        for (int j = 0; j < 4; ++j) {
            const std::string* t = nxt();
            std::optional<long long> v = std::nullopt;
            if (t) v = parseInt(*t);
            if (!v || *v < -COORD_MAX || *v > COORD_MAX) {
                bad();
                return;
            }
            vals[j] = *v;
        }
        Rect r{vals[0], vals[1], vals[2], vals[3]};
        if (r.x1 >= r.x2 || r.y1 >= r.y2) {
            bad();
            return;
        }
        rects.push_back(r);
    }
    std::cout << unionArea(rects) << "\n";
}

}  // namespace

int main() {
#ifdef _WIN32
    _setmode(_fileno(stdout), _O_BINARY);
#endif
    std::string raw((std::istreambuf_iterator<char>(std::cin)),
                    std::istreambuf_iterator<char>());
    bool blank = true;
    for (std::size_t i = 0; i < raw.size(); ++i) {
        if (static_cast<unsigned char>(raw[i]) > ' ') {
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
