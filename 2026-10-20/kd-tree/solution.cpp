// K-D Tree：k 維空間的二分搜尋樹，把「最近鄰」從 O(n) 降到期望 O(log n)
//
// 題意：
//     給平面上 n 個點，建一棵 **2-d tree**，支援五種空間查詢，並與 O(n) 暴力掃描
//     逐項對拍，確保每一次剪枝都沒剪掉真正的答案。
//
//     操作碼（每行 op 及其參數）：
//         1  qx qy                 最近鄰：輸出「距離平方 點編號」（平手取編號最小）
//         2  qx qy k               k 近鄰：輸出前 k 個距離平方（升序；不足補 -1）
//         3  x1 y1 x2 y2           矩形範圍計數（邊界含在內）
//         4  qx qy r               半徑計數：距離平方 ≤ r² 的點數
//         5  qx qy x1 y1 x2 y2    帶矩形約束的最近鄰（只考慮落在矩形內的點）
//
// 思路：
//     ### 1. 為什麼二分空間也能「二分搜尋」
//     BST 是對**一維**鍵值二分；k-d tree 把它推廣成輪流對**各個維度**二分。
//     本題在 2D：每次建節點時，看當前這群點在哪個維度上**散得最開**（spread 最大），
//     就沿那個維度排序取中位數切開——這樣切出來的「格子」最接近正方形，
//     剪枝效果最好（固定輪流切 x/y 也可以，但遇到狹長分佈會退化）。
//
//     每個節點額外存**這棵子樹所有點的包圍盒**（bounding box, 4 個數字）。
//     有了包圍盒，就能算出「查詢點到這整棵子樹的**最小可能距離**」，
//     這正是剪枝的依據。
//
//     ### 2. 最近鄰的剪枝（本題核心）
//     遞迴時先走「離查詢點比較近」的那一邊，拿到一個候選答案 best 之後，
//     再決定另一邊要不要走：
//
//         若 minDist(查詢點, 另一邊的包圍盒) > best  →  整棵子樹不可能更近，剪掉
//
//     注意必須是 `>` 才剪（不是 `>=`）：平手時我們要取**編號最小**的點，
//     所以距離相等的子樹仍必須進去看一眼。這是兩語言要一致的第一個細節。
//
//     ### 3. k 近鄰
//     維護一個大小 k 的有序緩衝區，第 k 名就是當前的「門檻」。
//     剪枝條件一樣用 `>`：包圍盒最小距離等於門檻時仍要進去，
//     因為可能存在「距離一樣但編號更小」的點，會把目前的第 k 名擠掉。
//
//     ### 4. 兩語言一致性的三個細節
//     - 切分維度平手時固定選 x（`spread_x >= spread_y` → axis 0）。
//     - 排序鍵必須是 `(該維度座標, 點編號)` 這種**全序**，不能只比座標，
//       否則座標相同的點順序在兩邊可能不同，中位數就換人，整棵樹的形狀跟著不同。
//     - 矩形查詢先把 x1>x2、y1>y2 的輸入正規化（min/max），
//       半徑 r 先夾到非負，避免兩個語言對「退化矩形」有不同解讀。
//
//     ### 5. 複雜度：是「期望」不是「最壞」
//     k-d tree 的查詢在點均勻分佈時是 O(log n)；但遇到點全落在圓周上、
//     或查詢點離所有點都很近的最壞情況，仍會退化成 O(n)。
//     要嚴格的最壞對數複雜度得用 R-tree / BVH / Voronoi，
//     不過實務上 k-d tree 因為常數極小、又好寫，仍是第一選擇。
//
// 應用場景：
//     - 地理資訊：外送 / 打車配對「最近的司機」、附近餐廳、電子圍欄
//     - 機器學習：kNN 分類與迴歸、DBSCAN / 密度聚類的鄰域查詢加速
//     - 電腦圖學：光線追蹤的加速結構（k-d tree 是 BVH 的老前輩）、碰撞檢測
//     - 電腦視覺：SIFT / ORB 特徵點的最近鄰匹配（FLANN 底層就有 k-d tree）
//     - 天文學：星表交叉比對（cross-match），把兩份上百萬筆的星表對起來
//     - 物理模擬：SPH / 分子動力學的鄰居清單（neighbor list）
//
// 複雜度（n 個點、2 維）：
//     建樹                 O(n log² n)（每層排序；用 nth_element 可降到 O(n log n)）
//     最近鄰查詢           平均 O(log n)，最壞 O(n)
//     k 近鄰查詢           平均 O(k log n)
//     矩形 / 半徑計數       平均 O(√n + 輸出數)，最壞 O(n)
//     空間                 O(n)

#include <algorithm>
#include <cerrno>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#ifdef _WIN32
#include <fcntl.h>
#include <io.h>
#endif

namespace {

const int NMAX = 400;               // IO 模式下的點數上限
const int QMAX = 200;               // IO 模式下的查詢數上限
const int KMAX = 8;                 // k 近鄰的 k 上限
const long long CMAX = 1000000;     // 座標絕對值上限
const long long RMAX = 3000000;     // 半徑上限

const long long INF = (1LL << 62);

typedef std::pair<long long, long long> Point;

// ---------------------------------------------------------------- k-d tree

struct KDNode {
    int idx;
    int left;                        // 節點在 pool 裡的索引，-1 表示空
    int right;
    long long minx, maxx, miny, maxy;
};

class KDTree {
public:
    explicit KDTree(const std::vector<Point>& pts) : pts_(pts) {
        if (!pts_.empty()) {
            std::vector<int> ids(pts_.size());
            for (size_t i = 0; i < pts_.size(); ++i) ids[i] = static_cast<int>(i);
            root_ = build(ids);
        } else {
            root_ = -1;
        }
    }

    // -- 最近鄰 --
    std::pair<long long, long long> nearest(long long qx, long long qy) const {
        bd_ = INF;
        bi_ = -1;
        nearestRec(root_, qx, qy);
        if (bi_ < 0) return {-1, -1};
        return {bd_, bi_};
    }

    // -- k 近鄰 --
    std::vector<long long> knn(long long qx, long long qy, int k) const {
        buf_.clear();
        k_ = k;
        knnRec(root_, qx, qy);
        std::sort(buf_.begin(), buf_.end());
        std::vector<long long> out;
        for (const auto& p : buf_) out.push_back(p.first);
        return out;
    }

    // -- 矩形範圍計數 --
    long long rangeCount(long long x1, long long y1, long long x2, long long y2) const {
        if (x1 > x2) std::swap(x1, x2);
        if (y1 > y2) std::swap(y1, y2);
        rx1_ = x1; ry1_ = y1; rx2_ = x2; ry2_ = y2;
        cnt_ = 0;
        rangeRec(root_);
        return cnt_;
    }

    // -- 半徑計數 --
    long long radiusCount(long long qx, long long qy, long long r) const {
        if (r < 0) r = 0;
        qx_ = qx; qy_ = qy; r2_ = r * r;
        cnt_ = 0;
        radiusRec(root_);
        return cnt_;
    }

    // -- 帶矩形約束的最近鄰 --
    std::pair<long long, long long> nearestInRect(long long qx, long long qy,
                                                  long long x1, long long y1,
                                                  long long x2, long long y2) const {
        if (x1 > x2) std::swap(x1, x2);
        if (y1 > y2) std::swap(y1, y2);
        rx1_ = x1; ry1_ = y1; rx2_ = x2; ry2_ = y2;
        bd_ = INF;
        bi_ = -1;
        nearestRectRec(root_, qx, qy);
        if (bi_ < 0) return {-1, -1};
        return {bd_, bi_};
    }

private:
    int build(std::vector<int> ids) {
        if (ids.empty()) return -1;
        long long minx = pts_[ids[0]].first, maxx = minx;
        long long miny = pts_[ids[0]].second, maxy = miny;
        for (int id : ids) {
            minx = std::min(minx, pts_[id].first);
            maxx = std::max(maxx, pts_[id].first);
            miny = std::min(miny, pts_[id].second);
            maxy = std::max(maxy, pts_[id].second);
        }
        // 散得最開的維度優先切；平手時選 x（axis 0）
        int axis = ((maxx - minx) >= (maxy - miny)) ? 0 : 1;
        // 排序鍵必須是全序（座標, 編號）
        std::sort(ids.begin(), ids.end(), [&](int a, int b) {
            long long va = axis == 0 ? pts_[a].first : pts_[a].second;
            long long vb = axis == 0 ? pts_[b].first : pts_[b].second;
            if (va != vb) return va < vb;
            return a < b;
        });
        size_t mid = ids.size() / 2;
        int id = pool_.size();
        pool_.push_back(KDNode());
        std::vector<int> leftIds(ids.begin(), ids.begin() + mid);
        std::vector<int> rightIds(ids.begin() + mid + 1, ids.end());
        int lchild = build(leftIds);
        int rchild = build(rightIds);
        // 注意：build 會讓 pool_ 重新配置，所以回填欄位時要重新取指標
        pool_[id].idx = ids[mid];
        pool_[id].left = lchild;
        pool_[id].right = rchild;
        pool_[id].minx = minx;
        pool_[id].maxx = maxx;
        pool_[id].miny = miny;
        pool_[id].maxy = maxy;
        return id;
    }

    // 包圍盒到點的最小距離平方
    long long boxDist(int t, long long qx, long long qy) const {
        if (t < 0) return INF;
        const KDNode& nd = pool_[t];
        long long dx = 0;
        if (qx < nd.minx) dx = nd.minx - qx;
        else if (qx > nd.maxx) dx = qx - nd.maxx;
        long long dy = 0;
        if (qy < nd.miny) dy = nd.miny - qy;
        else if (qy > nd.maxy) dy = qy - nd.maxy;
        return dx * dx + dy * dy;
    }

    // 包圍盒是否與矩形有交集
    bool boxHit(int t, long long x1, long long y1, long long x2, long long y2) const {
        if (t < 0) return false;
        const KDNode& nd = pool_[t];
        return !(nd.maxx < x1 || nd.minx > x2 || nd.maxy < y1 || nd.miny > y2);
    }

    long long dist2(int i, long long qx, long long qy) const {
        long long dx = pts_[i].first - qx;
        long long dy = pts_[i].second - qy;
        return dx * dx + dy * dy;
    }

    void nearestRec(int t, long long qx, long long qy) const {
        if (t < 0) return;
        long long d = dist2(pool_[t].idx, qx, qy);
        if (d < bd_ || (d == bd_ && pool_[t].idx < bi_)) {
            bd_ = d;
            bi_ = pool_[t].idx;
        }
        long long dl = boxDist(pool_[t].left, qx, qy);
        long long dr = boxDist(pool_[t].right, qx, qy);
        int first = pool_[t].left, second = pool_[t].right;
        long long fd = dl, sd = dr;
        if (dl > dr) {
            std::swap(first, second);
            std::swap(fd, sd);
        }
        // 用 `<=` 才剪（不是 `<`）：平手時還要進去找編號更小的點
        if (fd <= bd_) nearestRec(first, qx, qy);
        if (sd <= bd_) nearestRec(second, qx, qy);
    }

    void knnRec(int t, long long qx, long long qy) const {
        if (t < 0) return;
        long long d = dist2(pool_[t].idx, qx, qy);
        buf_.push_back({d, pool_[t].idx});
        std::sort(buf_.begin(), buf_.end());
        if (static_cast<int>(buf_.size()) > k_) buf_.resize(k_);
        long long thresh = (static_cast<int>(buf_.size()) >= k_) ? buf_.back().first : INF;
        long long dl = boxDist(pool_[t].left, qx, qy);
        long long dr = boxDist(pool_[t].right, qx, qy);
        int first = pool_[t].left, second = pool_[t].right;
        long long fd = dl, sd = dr;
        if (dl > dr) {
            std::swap(first, second);
            std::swap(fd, sd);
        }
        // 等於門檻時仍要進去：可能有「距離一樣但編號更小」的點擠掉第 k 名
        if (fd <= thresh) {
            knnRec(first, qx, qy);
            thresh = (static_cast<int>(buf_.size()) >= k_) ? buf_.back().first : INF;
        }
        if (sd <= thresh) knnRec(second, qx, qy);
    }

    void rangeRec(int t) const {
        if (t < 0) return;
        if (!boxHit(t, rx1_, ry1_, rx2_, ry2_)) return;
        long long px = pts_[pool_[t].idx].first;
        long long py = pts_[pool_[t].idx].second;
        if (rx1_ <= px && px <= rx2_ && ry1_ <= py && py <= ry2_) ++cnt_;
        rangeRec(pool_[t].left);
        rangeRec(pool_[t].right);
    }

    void radiusRec(int t) const {
        if (t < 0) return;
        if (boxDist(t, qx_, qy_) > r2_) return;
        if (dist2(pool_[t].idx, qx_, qy_) <= r2_) ++cnt_;
        radiusRec(pool_[t].left);
        radiusRec(pool_[t].right);
    }

    void nearestRectRec(int t, long long qx, long long qy) const {
        if (t < 0) return;
        if (!boxHit(t, rx1_, ry1_, rx2_, ry2_)) return;
        long long px = pts_[pool_[t].idx].first;
        long long py = pts_[pool_[t].idx].second;
        if (rx1_ <= px && px <= rx2_ && ry1_ <= py && py <= ry2_) {
            long long d = dist2(pool_[t].idx, qx, qy);
            if (d < bd_ || (d == bd_ && pool_[t].idx < bi_)) {
                bd_ = d;
                bi_ = pool_[t].idx;
            }
        }
        long long dl = boxDist(pool_[t].left, qx, qy);
        long long dr = boxDist(pool_[t].right, qx, qy);
        int first = pool_[t].left, second = pool_[t].right;
        long long fd = dl, sd = dr;
        if (dl > dr) {
            std::swap(first, second);
            std::swap(fd, sd);
        }
        if (fd <= bd_) nearestRectRec(first, qx, qy);
        if (sd <= bd_) nearestRectRec(second, qx, qy);
    }

    std::vector<Point> pts_;
    std::vector<KDNode> pool_;
    int root_ = -1;
    // 遞迴過程中的可變狀態（對應 Python 版本的實例欄位）
    mutable long long bd_ = INF;
    mutable long long bi_ = -1;
    mutable std::vector<std::pair<long long, long long>> buf_;
    mutable int k_ = 1;
    mutable long long cnt_ = 0;
    mutable long long rx1_ = 0, ry1_ = 0, rx2_ = 0, ry2_ = 0;
    mutable long long qx_ = 0, qy_ = 0, r2_ = 0;
};

// ---------------------------------------------------------------- 暴力對拍

std::pair<long long, long long> bfNearest(const std::vector<Point>& pts, long long qx, long long qy) {
    std::pair<long long, long long> best = {INF, -1};
    for (size_t i = 0; i < pts.size(); ++i) {
        long long dx = pts[i].first - qx;
        long long dy = pts[i].second - qy;
        std::pair<long long, long long> cand = {dx * dx + dy * dy, static_cast<long long>(i)};
        if (cand < best) best = cand;
    }
    if (best.second < 0) return {-1, -1};
    return best;
}

std::vector<long long> bfKnn(const std::vector<Point>& pts, long long qx, long long qy, int k) {
    std::vector<std::pair<long long, long long>> ds;
    for (size_t i = 0; i < pts.size(); ++i) {
        long long dx = pts[i].first - qx;
        long long dy = pts[i].second - qy;
        ds.push_back({dx * dx + dy * dy, static_cast<long long>(i)});
    }
    std::sort(ds.begin(), ds.end());
    std::vector<long long> out;
    for (int i = 0; i < k && i < static_cast<int>(ds.size()); ++i) out.push_back(ds[i].first);
    return out;
}

long long bfRange(const std::vector<Point>& pts, long long x1, long long y1, long long x2, long long y2) {
    if (x1 > x2) std::swap(x1, x2);
    if (y1 > y2) std::swap(y1, y2);
    long long c = 0;
    for (const Point& p : pts) {
        if (x1 <= p.first && p.first <= x2 && y1 <= p.second && p.second <= y2) ++c;
    }
    return c;
}

long long bfRadius(const std::vector<Point>& pts, long long qx, long long qy, long long r) {
    if (r < 0) r = 0;
    long long r2 = r * r;
    long long c = 0;
    for (const Point& p : pts) {
        long long dx = p.first - qx;
        long long dy = p.second - qy;
        if (dx * dx + dy * dy <= r2) ++c;
    }
    return c;
}

std::pair<long long, long long> bfNearestRect(const std::vector<Point>& pts, long long qx, long long qy,
                                              long long x1, long long y1, long long x2, long long y2) {
    if (x1 > x2) std::swap(x1, x2);
    if (y1 > y2) std::swap(y1, y2);
    std::pair<long long, long long> best = {INF, -1};
    for (size_t i = 0; i < pts.size(); ++i) {
        long long px = pts[i].first, py = pts[i].second;
        if (x1 <= px && px <= x2 && y1 <= py && py <= y2) {
            long long dx = px - qx;
            long long dy = py - qy;
            std::pair<long long, long long> cand = {dx * dx + dy * dy, static_cast<long long>(i)};
            if (cand < best) best = cand;
        }
    }
    if (best.second < 0) return {-1, -1};
    return best;
}

// ---------------------------------------------------------------- 主流程

std::string run(const std::vector<Point>& pts, const std::vector<std::vector<long long>>& qs) {
    KDTree tree(pts);
    std::ostringstream oss;
    int ok = 1;

    for (const auto& q : qs) {
        long long op = q[0];
        if (op == 1) {
            long long qx = q[1], qy = q[2];
            std::pair<long long, long long> got = tree.nearest(qx, qy);
            std::pair<long long, long long> want = bfNearest(pts, qx, qy);
            oss << got.first << " " << got.second << "\n";
            if (got != want) ok = 0;
        } else if (op == 2) {
            long long qx = q[1], qy = q[2];
            int k = static_cast<int>(q[3]);
            std::vector<long long> got = tree.knn(qx, qy, k);
            std::vector<long long> want = bfKnn(pts, qx, qy, k);
            std::vector<long long> g = got, w = want;
            while (static_cast<int>(g.size()) < k) g.push_back(-1);
            while (static_cast<int>(w.size()) < k) w.push_back(-1);
            for (int i = 0; i < k; ++i) {
                if (i) oss << " ";
                oss << g[i];
            }
            oss << "\n";
            if (g != w) ok = 0;
        } else if (op == 3) {
            long long x1 = q[1], y1 = q[2], x2 = q[3], y2 = q[4];
            long long got = tree.rangeCount(x1, y1, x2, y2);
            long long want = bfRange(pts, x1, y1, x2, y2);
            oss << got << "\n";
            if (got != want) ok = 0;
        } else if (op == 4) {
            long long qx = q[1], qy = q[2], r = q[3];
            long long got = tree.radiusCount(qx, qy, r);
            long long want = bfRadius(pts, qx, qy, r);
            oss << got << "\n";
            if (got != want) ok = 0;
        } else if (op == 5) {
            long long qx = q[1], qy = q[2];
            long long x1 = q[3], y1 = q[4], x2 = q[5], y2 = q[6];
            std::pair<long long, long long> got = tree.nearestInRect(qx, qy, x1, y1, x2, y2);
            std::pair<long long, long long> want = bfNearestRect(pts, qx, qy, x1, y1, x2, y2);
            oss << got.first << " " << got.second << "\n";
            if (got != want) ok = 0;
        }
    }
    oss << ok << "\n";
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

std::vector<std::string> splitLines(const std::string& s) {
    std::vector<std::string> out;
    std::istringstream iss(s);
    std::string ln;
    while (std::getline(iss, ln)) {
        if (!ln.empty() && ln.back() == '\r') ln.pop_back();
        out.push_back(ln);
    }
    return out;
}

void runTests() {
    // ---- 固定用例：五個點 ----
    {
        std::vector<Point> pts = {{0, 0}, {3, 4}, {10, 10}, {-5, 2}, {1, 1}};
        std::vector<std::vector<long long>> qs = {
            {1, 0, 0}, {1, 3, 4}, {2, 0, 0, 3}, {3, -10, -10, 2, 2},
            {4, 0, 0, 5}, {5, 0, 0, -10, -10, 2, 2}};
        std::vector<std::string> lines = splitLines(run(pts, qs));
        if (lines.size() != 7) { std::cout << "FATAL f1 size\n"; std::abort(); }
        if (lines[0] != "0 0") { std::cout << "FATAL f1.0 [" << lines[0] << "]\n"; std::abort(); }
        if (lines[1] != "0 1") { std::cout << "FATAL f1.1 [" << lines[1] << "]\n"; std::abort(); }
        if (lines[2] != "0 2 25") { std::cout << "FATAL f1.2 [" << lines[2] << "]\n"; std::abort(); }
        if (lines[3] != "3") { std::cout << "FATAL f1.3 [" << lines[3] << "]\n"; std::abort(); }
        if (lines[4] != "3") { std::cout << "FATAL f1.4 [" << lines[4] << "]\n"; std::abort(); }
        if (lines[5] != "0 0") { std::cout << "FATAL f1.5 [" << lines[5] << "]\n"; std::abort(); }
        if (lines[6] != "1") { std::cout << "FATAL f1.6 [" << lines[6] << "]\n"; std::abort(); }
    }

    // ---- 空點集 ----
    {
        std::vector<Point> pts;
        std::vector<std::vector<long long>> qs = {
            {1, 0, 0}, {2, 0, 0, 2}, {3, 0, 0, 5, 5}, {4, 0, 0, 1}, {5, 0, 0, 0, 0, 9, 9}};
        std::vector<std::string> lines = splitLines(run(pts, qs));
        if (lines.size() != 6) { std::cout << "FATAL f2 size\n"; std::abort(); }
        if (lines[0] != "-1 -1") { std::cout << "FATAL f2.0 [" << lines[0] << "]\n"; std::abort(); }
        if (lines[1] != "-1 -1") { std::cout << "FATAL f2.1 [" << lines[1] << "]\n"; std::abort(); }
        if (lines[2] != "0") { std::cout << "FATAL f2.2 [" << lines[2] << "]\n"; std::abort(); }
        if (lines[3] != "0") { std::cout << "FATAL f2.3 [" << lines[3] << "]\n"; std::abort(); }
        if (lines[4] != "-1 -1") { std::cout << "FATAL f2.4 [" << lines[4] << "]\n"; std::abort(); }
        if (lines[5] != "1") { std::cout << "FATAL f2.5 [" << lines[5] << "]\n"; std::abort(); }
    }

    // ---- 全部點重疊（剪枝的退化情況）----
    {
        std::vector<Point> pts(7, {2, 2});
        std::vector<std::vector<long long>> qs = {{1, 2, 2}, {1, 0, 0}, {2, 2, 2, 4}, {3, 2, 2, 2, 2}};
        std::vector<std::string> lines = splitLines(run(pts, qs));
        if (lines.size() != 5) { std::cout << "FATAL f3 size\n"; std::abort(); }
        if (lines[0] != "0 0") { std::cout << "FATAL f3.0 [" << lines[0] << "]\n"; std::abort(); }
        if (lines[1] != "8 0") { std::cout << "FATAL f3.1 [" << lines[1] << "]\n"; std::abort(); }
        if (lines[2] != "0 0 0 0") { std::cout << "FATAL f3.2 [" << lines[2] << "]\n"; std::abort(); }
        if (lines[3] != "7") { std::cout << "FATAL f3.3 [" << lines[3] << "]\n"; std::abort(); }
        if (lines[4] != "1") { std::cout << "FATAL f3.4 [" << lines[4] << "]\n"; std::abort(); }
    }

    // ---- 隨機對拍（一般分佈）----
    {
        lcgState = 20261020ULL;
        for (int t = 0; t < 300; ++t) {
            int n = static_cast<int>(lcgInt(0, 60));
            std::vector<Point> pts;
            for (int i = 0; i < n; ++i) pts.push_back({lcgInt(-50, 50), lcgInt(-50, 50)});
            std::vector<std::vector<long long>> qs;
            int qn = static_cast<int>(lcgInt(1, 12));
            for (int i = 0; i < qn; ++i) {
                int op = static_cast<int>(lcgInt(1, 5));
                if (op == 1) qs.push_back({1, lcgInt(-60, 60), lcgInt(-60, 60)});
                else if (op == 2) qs.push_back({2, lcgInt(-60, 60), lcgInt(-60, 60), lcgInt(1, 6)});
                else if (op == 3) qs.push_back({3, lcgInt(-60, 60), lcgInt(-60, 60), lcgInt(-60, 60), lcgInt(-60, 60)});
                else if (op == 4) qs.push_back({4, lcgInt(-60, 60), lcgInt(-60, 60), lcgInt(0, 40)});
                else qs.push_back({5, lcgInt(-60, 60), lcgInt(-60, 60), lcgInt(-60, 60), lcgInt(-60, 60), lcgInt(-60, 60), lcgInt(-60, 60)});
            }
            std::string out = run(pts, qs);
            if (out.size() < 2 || out[out.size() - 2] != '1') {
                std::cout << "FATAL fuzz " << t << " [" << out << "]\n";
                std::abort();
            }
        }
    }

    // ---- 隨機對拍（退化分佈：同一條水平線 / 同一條直線 / 圓周上）----
    {
        for (int t = 0; t < 150; ++t) {
            int n = static_cast<int>(lcgInt(1, 50));
            int style = t % 3;
            std::vector<Point> pts;
            if (style == 0) {
                for (int i = 0; i < n; ++i) pts.push_back({lcgInt(-40, 40), 0});
            } else if (style == 1) {
                for (int i = 0; i < n; ++i) {
                    long long x = lcgInt(-40, 40);
                    pts.push_back({x, 3 * x + 1});
                }
            } else {
                for (int i = 0; i < n; ++i) {
                    double a = 2.0 * 3.14159265358979323846 * static_cast<double>(i) / static_cast<double>(n);
                    pts.push_back({static_cast<long long>(1000.0 * std::cos(a)),
                                   static_cast<long long>(1000.0 * std::sin(a))});
                }
            }
            std::vector<std::vector<long long>> qs;
            for (int i = 0; i < 8; ++i) {
                int op = static_cast<int>(lcgInt(1, 5));
                if (op == 1) qs.push_back({1, lcgInt(-1100, 1100), lcgInt(-1100, 1100)});
                else if (op == 2) qs.push_back({2, lcgInt(-1100, 1100), lcgInt(-1100, 1100), lcgInt(1, 5)});
                else if (op == 3) qs.push_back({3, lcgInt(-1100, 1100), lcgInt(-1100, 1100), lcgInt(-1100, 1100), lcgInt(-1100, 1100)});
                else if (op == 4) qs.push_back({4, lcgInt(-1100, 1100), lcgInt(-1100, 1100), lcgInt(0, 600)});
                else qs.push_back({5, lcgInt(-1100, 1100), lcgInt(-1100, 1100), lcgInt(-1100, 1100), lcgInt(-1100, 1100), lcgInt(-1100, 1100), lcgInt(-1100, 1100)});
            }
            std::string out = run(pts, qs);
            if (out.size() < 2 || out[out.size() - 2] != '1') {
                std::cout << "FATAL fuzz2 " << t << " style " << style << " [" << out << "]\n";
                std::abort();
            }
        }
    }

    // ---- 大規模：500 個點、200 次查詢 ----
    {
        lcgState = 777ULL;
        std::vector<Point> pts;
        for (int i = 0; i < 500; ++i) {
            pts.push_back({lcgInt(-CMAX, CMAX), lcgInt(-CMAX, CMAX)});
        }
        std::vector<std::vector<long long>> qs;
        for (int i = 0; i < 100; ++i) qs.push_back({1, lcgInt(-CMAX, CMAX), lcgInt(-CMAX, CMAX)});
        for (int i = 0; i < 100; ++i) qs.push_back({4, lcgInt(-CMAX, CMAX), lcgInt(-CMAX, CMAX), lcgInt(0, 300000)});
        std::string out = run(pts, qs);
        if (out.size() < 2 || out[out.size() - 2] != '1') {
            std::cout << "FATAL big [" << out.substr(out.size() - 40) << "]\n";
            std::abort();
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
    auto cl = [](long long v) -> long long {
        if (v < -CMAX) v = -CMAX;
        if (v > CMAX) v = CMAX;
        return v;
    };

    long long n = nextInt(0);
    if (n < 0) n = 0;
    if (n > NMAX) n = NMAX;
    std::vector<Point> pts;
    for (long long i = 0; i < n; ++i) {
        long long x = cl(nextInt(0));
        long long y = cl(nextInt(0));
        pts.push_back({x, y});
    }
    long long q = nextInt(0);
    if (q < 0) q = 0;
    if (q > QMAX) q = QMAX;
    std::vector<std::vector<long long>> qs;
    for (long long i = 0; i < q; ++i) {
        long long op = nextInt(0);
        if (op == 1) {
            qs.push_back({1, cl(nextInt(0)), cl(nextInt(0))});
        } else if (op == 2) {
            // 順序是 `2 qx qy k`，座標先讀，最後才是 k
            long long qx = cl(nextInt(0));
            long long qy = cl(nextInt(0));
            long long k = nextInt(1);
            if (k < 1) k = 1;
            if (k > KMAX) k = KMAX;
            qs.push_back({2, qx, qy, k});
        } else if (op == 3) {
            qs.push_back({3, cl(nextInt(0)), cl(nextInt(0)), cl(nextInt(0)), cl(nextInt(0))});
        } else if (op == 4) {
            // 順序是 `4 qx qy r`
            long long qx = cl(nextInt(0));
            long long qy = cl(nextInt(0));
            long long r = nextInt(0);
            if (r < 0) r = 0;
            if (r > RMAX) r = RMAX;
            qs.push_back({4, qx, qy, r});
        } else if (op == 5) {
            qs.push_back({5, cl(nextInt(0)), cl(nextInt(0)), cl(nextInt(0)),
                          cl(nextInt(0)), cl(nextInt(0)), cl(nextInt(0))});
        }
    }
    std::cout << run(pts, qs);
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
