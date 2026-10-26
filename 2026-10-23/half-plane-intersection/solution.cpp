// 半平面交（Half-Plane Intersection, HPI）
//
// 題意：
//     給 n 個半平面，第 i 個寫成 a_i·x + b_i·y ≤ c_i。求它們交集的邊界——
//     一個（可能空的、可能退化的）凸區域。為了讓輸出永遠是**有界**多邊形，
//     題目約定把結果裁剪到包圍盒 [-10⁹, 10⁹]² 內輸出。
//
//     IO 格式（有 stdin 時）：
//         第 1 行    n
//         接下來 n 行  a b c          半平面 a·x + b·y ≤ c
//     輸出：
//         第 1 行    k               交集中「頂點」的個數（k = 0 表示空集）
//         接下來 k 行  xn/xd yn/yd   頂點座標，以最簡分數輸出（分母 > 0），
//                                   按逆時針順序、從字典序最小的頂點開始
//
//     例：三個半平面 x ≥ 0、y ≥ 0、x + y ≤ 10 的交集是三角形，
//     輸出 3 與 (0,0)、(10,0)、(0,10)。
//
// 思路：
//     ### 1. 交集長什麼樣
//     半平面交 = 凸多邊形（或空集、或線段 / 點這類退化情形）。
//     有界凸多邊形的每個頂點，都是「兩條約束邊界線的交點」：
//         頂點 p 處至少有兩條約束線同時取等號，且兩條線不平行。
//     反過來也成立：任取兩條不平行的約束線，求交，若交點滿足全部約束，
//     它一定是極值點（頂點）——因為在交點附近，可行域被夾在兩條交叉線
//     形成的楔形裡，交點這個「尖角」抹不平。
//
//     ### 2. 本題的演算法：候選枚舉（O(n³)，精確有理數）
//     把全部約束線（n 個輸入半平面 + 4 條包圍盒線）兩兩求交，
//     用**精確有理數**檢查交點是否滿足所有約束，倖存者就是全部頂點。
//     這個方法直覺、好驗證，且所有座標都是精確分數，兩個語言可以逐位元組一致。
//
//     包圍盒的作用：無界區域（例如只剩一個半平面）被盒子裁成有界多邊形，
//     盒子的四條線也參與「兩兩求交」，於是盒角與「輸入線 ∩ 盒線」都會被枚舉到。
//     輸入範圍保證了：任何有界區域的頂點座標 |value| ≤ 10⁹ < 盒子邊長，
//     盒子不會誤裁有界區域；而無界區域必然被盒子裁到（凸集無界 ⟹ 某方向延伸到無窮）。
//
//     ### 3. 兩個對拍參照（Sutherland–Hodgman 裁剪 + 經典雙端棧）
//     候選枚舉好懂但慢，正式賽事另有兩種主流做法，本題都實作作為內置對拍：
//     - **多邊形裁剪**：以包圍盒為初始凸多邊形，依次被每個半平面裁剪
//       （Sutherland–Hodgman，精確有理數）。直接算出「交集裁到盒子」的結果，
//       天然處理空集 / 點 / 線段 / 多邊形全部退化情形，O(m·k)。
//     - **排序 + 增量雙端棧**：按邊界線極角排序後增量維護，O(m log m)，
//       常數小但平行 / 重複 / 空集等邊界情況極易寫錯，因此加了有效性護欄，
//       只在結果可信時參與比對。
//
//     ### 4. 退化與一致性細節
//     - (0, 0, c)：c ≥ 0 是全平面（跳過）；c < 0 是空集（直接輸出 0 個頂點）。
//     - 同向（法向量同向）的半平面只保留最緊的一個；法向量反向的半平面都要保留。
//     - 頂點排序：取字典序最小的頂點 v0，其餘頂點按極角繞 v0 排序即得逆時針順序
//       （凸多邊形從任一頂點看其餘頂點，張角 < π，純叉積比較就是全序）。
//     - 座標一律以最簡分數輸出，兩個語言都用精確有理數，不經過浮點。
//
//     ### 5. C++ 的精確有理數表示：齊次座標
//     所有頂點都是「兩條整係數約束線的交點」，因此每個頂點唯一地寫成
//         點 = (X/D, Y/D)，  D = det(線₁, 線₂) > 0（符號正規化）
//     X、Y、D 都是 __int128。兩個頂點作差時分子分母有界：
//         |X| ≤ 10¹¹、|D| ≤ 4·10⁴  ⟹  |Δx| 的分子 ≤ 8·10¹⁵
//     繞 v0 的極角比較化成 sign((p−v0)×(q−v0))，其分子
//         |Aₚ·C_q − Cₚ·A_q| ≤ 6.4·10³¹
//     穩穩落在 __int128（≈1.7·10³⁸）內，全程無浮點、無溢位。
//
//     裁剪對拍（Sutherland–Hodgman）另有一個符號化技巧：多邊形的每條邊都
//     貼著某條約束線，於是每個頂點除了座標再記一條「出邊線」編號，
//     裁剪新頂點直接由「舊邊線 ∩ 裁剪線」兩線求交得到，不必走參數插值
//     （那條路的中間值會爆 __int128）。連續重複點去重時保留**後者**，
//     剛好讓「出邊線」標註與實際邊一致；另加共線安全網，一旦標註失效
//     就放棄該次裁剪對拍（回傳 None），絕不輸出錯答案。
//
// 複雜度（m 條約束線，含包圍盒）：
//     候選枚舉        O(m³) 時間、O(m) 空間（m ≤ 64，實測毫秒級）
//     多邊形裁剪      O(m·k) 時間、O(k) 空間（k 為交集頂點數）
//     雙端棧版        O(m log m) 時間
//     輸出規模        O(m)
//
// 應用場景：
//     - 電腦圖學：視錐裁剪（frustum culling）、任意凸多邊形裁剪、CSG 布林運算
//     - 機器人與遊戲：雷達 / 視野範圍（visibility region）、導航可行區域
//     - 機器學習：線性分類器的可行域、支援向量機的間隔區域
//     - 作業研究：線性規劃可行域的外顯表示（半平面交就是 2D 的 H-多面體）
//     - 地理資訊：緩衝區交集、行政區合併的幾何核

#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <deque>
#include <iostream>
#include <iterator>
#include <optional>
#include <random>
#include <string>
#include <vector>

#ifdef _WIN32
#include <fcntl.h>
#include <io.h>
#endif

namespace {

const long long NMAX = 40;             // IO 模式的半平面數上限
const long long AMAX = 100;            // |a|、|b| 上限
const long long CMAX = 10000;          // |c| 上限
const long long CLAMP_ABS = 1000000000000000000LL;  // 通用整數夾取範圍
const long long BOX = 1000000000LL;    // 包圍盒邊長

using i128 = __int128;

// ---------------------------------------------------------------- 測試工具

void fail(const std::string& msg) {
    std::cerr << "FAIL: " << msg << std::endl;
    std::exit(1);
}

// ---------------------------------------------------------------- 基礎幾何

struct Plane {
    long long a, b, c;  // 半平面 a·x + b·y ≤ c
};

// 齊次座標點：(X, Y, D) 表示 (X/D, Y/D)，恆有 D > 0。
struct Pt {
    i128 X, Y, D;
};

i128 gcd128(i128 a, i128 b) {
    if (a < 0) a = -a;
    if (b < 0) b = -b;
    while (b) {
        i128 t = a % b;
        a = b;
        b = t;
    }
    return a;
}

Pt makePt(i128 X, i128 Y, i128 D) {
    if (D < 0) {
        X = -X;
        Y = -Y;
        D = -D;
    }
    return Pt{X, Y, D};
}

// 兩條直線 a·x + b·y = c 的交點；平行（det = 0）回傳 None。
std::optional<Pt> lineIntersect(const Plane& p, const Plane& q) {
    i128 det = (i128)p.a * q.b - (i128)q.a * p.b;
    if (det == 0) return std::nullopt;
    i128 X = (i128)p.c * q.b - (i128)q.c * p.b;
    i128 Y = (i128)p.a * q.c - (i128)q.a * p.c;
    return makePt(X, Y, det);
}

// 點 pt 到半平面 pl 的有號距離（≤ 0 表示在內部，含邊界）。
i128 dist(const Plane& pl, const Pt& pt) {
    return (i128)pl.a * pt.X + (i128)pl.b * pt.Y - (i128)pl.c * pt.D;
}

bool inside(const Plane& pl, const Pt& pt) { return dist(pl, pt) <= 0; }

long long gcd3(long long a, long long b, long long c) {
    return (long long)gcd128(gcd128(a, b), c);
}

bool sameDir(const Plane& p, const Plane& q) {
    // 法向量同向（cross = 0 且 dot > 0）
    return p.a * q.b - q.a * p.b == 0 && p.a * q.a + p.b * q.b > 0;
}

// p·(qn/qd) ≤ r 的精確判斷（qd ≠ 0）
bool ratLE(long long p, long long qn, long long qd, long long r) {
    i128 lhs = (i128)p * qn;
    i128 rhs = (i128)r * qd;
    if (qd > 0) return lhs <= rhs;
    return lhs >= rhs;  // 分母為負，不等號翻轉
}

bool tighter(const Plane& p, const Plane& q) {
    // 同向平行線中，p 的半平面是否比 q 更緊（p 的區域 ⊆ q 的區域）
    if (p.a != 0) return ratLE(p.c, q.a, p.a, q.c);
    return ratLE(p.c, q.b, p.b, q.c);
}

// 點的字典序比較（x 為主，D 恆正）
bool cmpPt(const Pt& p, const Pt& q) {
    i128 lx = p.X * q.D, rx = q.X * p.D;
    if (lx != rx) return lx < rx;
    return p.Y * q.D < q.Y * p.D;
}

bool eqPt(const Pt& p, const Pt& q) {
    return p.X * q.D == q.X * p.D && p.Y * q.D == q.Y * p.D;
}

inline bool operator==(const Pt& p, const Pt& q) { return eqPt(p, q); }

// sign((p − v) × (q − v))：繞 v 的極角比較（凸多邊形上等價於角度大小）。
// 齊次座標下 (p−v) = (Aₚ, Cₚ)/(Dₚ·Dᵥ)，兩個 D 皆正，
// 叉積分子 = Aₚ·C_q − Cₚ·A_q（|·| ≤ 6.4·10³¹，__int128 安全）。
int crossSign(const Pt& p, const Pt& q, const Pt& v) {
    i128 Ap = p.X * v.D - v.X * p.D;
    i128 Cp = p.Y * v.D - v.Y * p.D;
    i128 Aq = q.X * v.D - v.X * q.D;
    i128 Cq = q.Y * v.D - v.Y * q.D;
    i128 t = Ap * Cq - Cp * Aq;
    if (t > 0) return 1;
    if (t < 0) return -1;
    return 0;
}

std::string i128str(i128 v) {
    if (v == 0) return "0";
    bool neg = v < 0;
    unsigned __int128 u =
        neg ? (unsigned __int128)(-v) : (unsigned __int128)v;
    std::string s;
    while (u > 0) {
        s += char('0' + (int)(u % 10));
        u /= 10;
    }
    if (neg) s += '-';
    std::reverse(s.begin(), s.end());
    return s;
}

// 最簡分數輸出：x 與 y 各自約分（與 Python Fraction 一致）
std::string fmtPt(const Pt& p) {
    i128 gx = gcd128(p.X, p.D);
    i128 gy = gcd128(p.Y, p.D);
    return i128str(p.X / gx) + "/" + i128str(p.D / gx) + " " +
           i128str(p.Y / gy) + "/" + i128str(p.D / gy);
}

// ---------------------------------------------------------------- 正規化

struct NormResult {
    std::vector<Plane> planes;
    bool empty;
};

NormResult normalizePlanes(const std::vector<Plane>& planes) {
    bool empty = false;
    std::vector<Plane> cleaned;
    for (const Plane& pl : planes) {
        if (pl.a == 0 && pl.b == 0) {
            if (pl.c < 0) empty = true;
            continue;
        }
        long long g = gcd3(pl.a, pl.b, pl.c);
        if (g > 1) {
            cleaned.push_back(Plane{pl.a / g, pl.b / g, pl.c / g});
        } else {
            cleaned.push_back(pl);
        }
    }
    if (empty) return NormResult{std::vector<Plane>(), true};
    // 包圍盒：把無界區域裁成有界多邊形
    cleaned.push_back(Plane{1, 0, BOX});
    cleaned.push_back(Plane{-1, 0, BOX});
    cleaned.push_back(Plane{0, 1, BOX});
    cleaned.push_back(Plane{0, -1, BOX});
    // 同向去重，保留最緊的
    std::vector<Plane> deduped;
    for (const Plane& pl : cleaned) {
        bool merged = false;
        for (size_t i = 0; i < deduped.size(); ++i) {
            if (sameDir(pl, deduped[i])) {
                if (tighter(pl, deduped[i])) deduped[i] = pl;
                merged = true;
                break;
            }
        }
        if (!merged) deduped.push_back(pl);
    }
    return NormResult{deduped, false};
}

// ---------------------------------------------------------------- 候選枚舉

std::vector<Pt> hpiVertices(const std::vector<Plane>& planes) {
    std::vector<Pt> verts;
    int m = (int)planes.size();
    for (int i = 0; i < m; ++i) {
        for (int j = i + 1; j < m; ++j) {
            std::optional<Pt> pt = lineIntersect(planes[i], planes[j]);
            if (!pt) continue;
            bool ok = true;
            for (int k = 0; k < m; ++k) {
                if (!inside(planes[k], *pt)) {
                    ok = false;
                    break;
                }
            }
            if (ok) verts.push_back(*pt);
        }
    }
    std::sort(verts.begin(), verts.end(), cmpPt);
    verts.erase(std::unique(verts.begin(), verts.end(), eqPt), verts.end());
    return verts;
}

// ---------------------------------------------------------------- 雙端棧版

int angleCmp(const Plane& p, const Plane& q) {
    long long a1 = p.a, b1 = p.b, a2 = q.a, b2 = q.b;
    int h1 = (a1 > 0 || (a1 == 0 && b1 < 0)) ? 0 : 1;
    int h2 = (a2 > 0 || (a2 == 0 && b2 < 0)) ? 0 : 1;
    if (h1 != h2) return h1 < h2 ? -1 : 1;
    long long cr = a1 * b2 - a2 * b1;  // cross(d1, d2) > 0 ⟹ d1 極角更小
    if (cr > 0) return -1;
    if (cr < 0) return 1;
    return 0;
}

// 經典 O(n log n) 半平面交（排序 + 增量雙端棧）；退化時回傳 None。
std::optional<std::vector<Pt>> hpiDeque(const std::vector<Plane>& planes) {
    std::vector<Plane> ps = planes;
    std::stable_sort(ps.begin(), ps.end(),
                     [](const Plane& p, const Plane& q) {
                         return angleCmp(p, q) < 0;
                     });
    std::deque<Plane> dq;
    for (const Plane& q : ps) {
        while (dq.size() >= 2) {
            std::optional<Pt> pt =
                lineIntersect(dq[dq.size() - 2], dq[dq.size() - 1]);
            if (pt && inside(q, *pt)) break;
            dq.pop_back();
        }
        while (dq.size() >= 2) {
            std::optional<Pt> pt = lineIntersect(dq[0], dq[1]);
            if (pt && inside(q, *pt)) break;
            dq.pop_front();
        }
        dq.push_back(q);
    }
    while (dq.size() >= 3) {
        std::optional<Pt> pt =
            lineIntersect(dq[dq.size() - 2], dq[dq.size() - 1]);
        if (pt && inside(dq[0], *pt)) break;
        dq.pop_back();
    }
    while (dq.size() >= 3) {
        std::optional<Pt> pt = lineIntersect(dq[0], dq[1]);
        if (pt && inside(dq[dq.size() - 1], *pt)) break;
        dq.pop_front();
    }
    if (dq.size() < 3) return std::nullopt;
    int k = (int)dq.size();
    std::vector<Pt> vs;
    for (int i = 0; i < k; ++i) {
        std::optional<Pt> pt = lineIntersect(dq[i], dq[(i + 1) % k]);
        if (!pt) return std::nullopt;
        vs.push_back(*pt);
    }
    std::sort(vs.begin(), vs.end(), cmpPt);
    vs.erase(std::unique(vs.begin(), vs.end(), eqPt), vs.end());
    // 有效性護欄：交集為空時經典雙端棧可能吐出幽靈頂點，逐個驗證
    for (const Pt& v : vs) {
        for (const Plane& pl : planes) {
            if (!inside(pl, v)) return std::nullopt;
        }
    }
    return vs;
}

// ---------------------------------------------------------------- 多邊形裁剪

// 裁剪多邊形頂點：座標 + 「出邊」貼著的約束線編號（lines 的下標）。
struct PV {
    Pt p;
    int outLine;
};

bool collinear(const Pt& u, const Pt& v, const Pt& w) {
    // cross((v − u), (w − u)) == 0
    i128 A1 = v.X * u.D - u.X * v.D;
    i128 C1 = v.Y * u.D - u.Y * v.D;
    i128 A2 = w.X * u.D - u.X * w.D;
    i128 C2 = w.Y * u.D - u.Y * w.D;
    return A1 * C2 - C1 * A2 == 0;
}

// 用 lines[k] 裁剪多邊形 poly；不變式被破壞（理論上不發生）時回傳 None。
std::optional<std::vector<PV>> clipOnce(const std::vector<PV>& poly,
                                        const std::vector<Plane>& lines,
                                        int k) {
    std::vector<PV> out;
    int m = (int)poly.size();
    for (int i = 0; i < m; ++i) {
        const PV& U = poly[i];
        const PV& V = poly[(i + 1) % m];
        i128 d1 = dist(lines[k], U.p);
        i128 d2 = dist(lines[k], V.p);
        bool in1 = d1 <= 0;
        bool in2 = d2 <= 0;
        if (in1) out.push_back(U);  // 保留頂點，出邊線不變
        if (in1 != in2) {
            // 邊 (U, V) 貼著 U.outLine；新頂點 = 兩線之交（齊次、有界）
            int edgeLine = U.outLine;
            std::optional<Pt> w = lineIntersect(lines[edgeLine], lines[k]);
            if (!w) return std::nullopt;
            if (!collinear(U.p, V.p, *w)) return std::nullopt;  // 安全網
            // 內→外：下一個輸出必為裁剪點（新邊在 k 上）；
            // 外→內：下一個輸出必為保留點（新邊在舊邊線上）。
            int wOut = in1 ? k : edgeLine;
            out.push_back(PV{*w, wOut});
        }
    }
    // 連續去重，保留**後者**（其出邊線才描述實際的出邊）
    std::vector<PV> ded;
    for (const PV& pv : out) {
        if (!ded.empty() && eqPt(ded.back().p, pv.p)) {
            ded.back() = pv;
        } else {
            ded.push_back(pv);
        }
    }
    return ded;
}

// 以包圍盒為初始多邊形，依次被每個半平面精確裁剪（Sutherland–Hodgman）。
// 回傳排序去重後的頂點列表；不變式被破壞時回傳 None（放棄對拍）。
std::optional<std::vector<Pt>> hpiClip(const std::vector<Plane>& planes) {
    std::vector<Plane> lines;
    for (const Plane& pl : planes) {
        if (pl.a == 0 && pl.b == 0) {
            if (pl.c < 0) return std::vector<Pt>{};  // 空集
            continue;
        }
        lines.push_back(pl);
    }
    int iL = (int)lines.size();
    lines.push_back(Plane{-1, 0, BOX});
    int iR = (int)lines.size();
    lines.push_back(Plane{1, 0, BOX});
    int iD = (int)lines.size();
    lines.push_back(Plane{0, -1, BOX});
    int iU = (int)lines.size();
    lines.push_back(Plane{0, 1, BOX});
    std::vector<PV> poly = {
        PV{makePt(-BOX, -BOX, 1), iD},
        PV{makePt(BOX, -BOX, 1), iR},
        PV{makePt(BOX, BOX, 1), iU},
        PV{makePt(-BOX, BOX, 1), iL},
    };
    int idx = -1;
    for (const Plane& pl : planes) {
        if (pl.a == 0 && pl.b == 0) {
            if (pl.c < 0) return std::vector<Pt>{};
            continue;
        }
        ++idx;
        std::optional<std::vector<PV>> res = clipOnce(poly, lines, idx);
        if (!res) return std::nullopt;
        poly = *res;
        if (poly.empty()) return std::vector<Pt>{};
    }
    std::vector<Pt> vs;
    for (const PV& pv : poly) vs.push_back(pv.p);
    std::sort(vs.begin(), vs.end(), cmpPt);
    vs.erase(std::unique(vs.begin(), vs.end(), eqPt), vs.end());
    return vs;
}

// ---------------------------------------------------------------- 輸出排序

// 凸多邊形頂點：逆時針、從字典序最小點開始。
std::vector<Pt> orderVertices(const std::vector<Pt>& verts) {
    if (verts.size() <= 2) return verts;  // 已排序
    Pt v0 = verts.front();
    std::vector<Pt> rest(verts.begin() + 1, verts.end());
    std::stable_sort(rest.begin(), rest.end(),
                     [&v0](const Pt& p, const Pt& q) {
                         return crossSign(p, q, v0) > 0;  // cross > 0 ⟹ p 在前
                     });
    std::vector<Pt> out;
    out.push_back(v0);
    out.insert(out.end(), rest.begin(), rest.end());
    return out;
}

std::vector<Pt> runVerts(const std::vector<Plane>& planes) {
    NormResult norm = normalizePlanes(planes);
    if (norm.empty) return {};
    std::vector<Pt> verts = hpiVertices(norm.planes);
    return orderVertices(verts);
}

std::string run(const std::vector<Plane>& planes) {
    std::vector<Pt> verts = runVerts(planes);
    std::string out = std::to_string(verts.size()) + "\n";
    for (const Pt& p : verts) out += fmtPt(p) + "\n";
    return out;
}

// ---------------------------------------------------------------- 結構檢查

bool onSegment(const Pt& p, const Pt& q, long long rx, long long ry) {
    // r 是否在線段 pq 上（精確判斷；r 為整數點）
    i128 Aq = q.X * p.D - p.X * q.D;
    i128 Cq = q.Y * p.D - p.Y * q.D;
    i128 Ar = (i128)rx * p.D - p.X;
    i128 Cr = (i128)ry * p.D - p.Y;
    if (Aq * Cr - Cq * Ar != 0) return false;
    // 包圍盒（分母 > 0，可直接交叉相乘）
    bool rxGEmin = !((i128)rx * p.D < p.X) || !((i128)rx * q.D < q.X);
    bool rxLEmax = !((i128)rx * p.D > p.X) || !((i128)rx * q.D > q.X);
    bool ryGEmin = !((i128)ry * p.D < p.Y) || !((i128)ry * q.D < q.Y);
    bool ryLEmax = !((i128)ry * p.D > p.Y) || !((i128)ry * q.D > q.Y);
    return rxGEmin && rxLEmax && ryGEmin && ryLEmax;
}

bool pointInPoly(const std::vector<Pt>& verts, long long rx, long long ry) {
    int k = (int)verts.size();
    if (k == 0) return false;
    if (k == 1) {
        return (i128)rx * verts[0].D == verts[0].X &&
               (i128)ry * verts[0].D == verts[0].Y;
    }
    if (k == 2) return onSegment(verts[0], verts[1], rx, ry);
    for (int i = 0; i < k; ++i) {
        const Pt& A = verts[i];
        const Pt& B = verts[(i + 1) % k];
        i128 A1 = B.X * A.D - A.X * B.D;
        i128 C1 = B.Y * A.D - A.Y * B.D;
        i128 Ar = (i128)rx * A.D - A.X;
        i128 Cr = (i128)ry * A.D - A.Y;
        if (A1 * Cr - C1 * Ar < 0) return false;
    }
    return true;
}

void checkStructure(const std::vector<Plane>& L, const std::vector<Pt>& verts) {
    int k = (int)verts.size();
    if (k < 3) return;
    for (int i = 0; i < k; ++i) {
        const Pt& A = verts[i];
        const Pt& B = verts[(i + 1) % k];
        const Pt& C = verts[(i + 2) % k];
        i128 A1 = B.X * A.D - A.X * B.D;
        i128 C1 = B.Y * A.D - A.Y * B.D;
        i128 A2 = C.X * B.D - B.X * C.D;
        i128 C2 = C.Y * B.D - B.Y * C.D;
        i128 cr = A1 * C2 - C1 * A2;
        if (!(cr > 0)) fail("not strictly convex ccw");
    }
    for (int i = 0; i < k; ++i) {
        const Pt& p = verts[i];
        const Pt& q = verts[(i + 1) % k];
        bool found = false;
        for (const Plane& pl : L) {
            if ((i128)pl.a * p.X + (i128)pl.b * p.Y == (i128)pl.c * p.D &&
                (i128)pl.a * q.X + (i128)pl.b * q.Y == (i128)pl.c * q.D) {
                found = true;
                break;
            }
        }
        if (!found) fail("edge not on a constraint line");
    }
}

i128 ratFloorXY(i128 X, i128 D) {  // floor(X/D)，D > 0
    i128 q = X / D;
    i128 r = X % D;
    if (r != 0 && X < 0) q -= 1;
    return q;
}

i128 ratCeilXY(i128 X, i128 D) {  // ceil(X/D)，D > 0
    i128 q = X / D;
    i128 r = X % D;
    if (r != 0 && X > 0) q += 1;
    return q;
}

void checkMembership(const std::vector<Plane>& orig,
                     const std::vector<Pt>& verts, std::mt19937& rng) {
    if (verts.empty()) return;
    // 頂點的包圍盒（外擴 3）
    Pt loX = verts[0], hiX = verts[0], loY = verts[0], hiY = verts[0];
    for (const Pt& v : verts) {
        if (cmpPt(v, loX)) loX = v;
        if (cmpPt(hiX, v)) hiX = v;
        if (v.Y * loY.D < loY.Y * v.D) loY = v;
        if (v.Y * hiY.D > hiY.Y * v.D) hiY = v;
    }
    long long lx = (long long)ratFloorXY(loX.X, loX.D) - 3;
    long long hx = (long long)ratCeilXY(hiX.X, hiX.D) + 3;
    long long ly = (long long)ratFloorXY(loY.Y, loY.D) - 3;
    long long hy = (long long)ratCeilXY(hiY.Y, hiY.D) + 3;
    std::uniform_int_distribution<long long> dx(lx, hx);
    std::uniform_int_distribution<long long> dy(ly, hy);
    for (int it = 0; it < 120; ++it) {
        long long x = dx(rng);
        long long y = dy(rng);
        bool inpoly = pointInPoly(verts, x, y);
        bool feasible = true;
        for (const Plane& pl : orig) {
            if (pl.a == 0 && pl.b == 0) continue;
            if (!((i128)pl.a * x + (i128)pl.b * y <= pl.c)) {
                feasible = false;
                break;
            }
        }
        if (inpoly != feasible) fail("membership");
    }
}

// ---------------------------------------------------------------- 內置測試

std::vector<Plane> genPlanes(std::mt19937& rng, int style) {
    std::uniform_int_distribution<int> d24(0, 24);
    int n = d24(rng);
    std::vector<Plane> planes;
    if (style == 0) {
        // 保證非空：先挑內點，再環繞它產生半平面
        std::uniform_int_distribution<int> d50(-50, 50);
        std::uniform_int_distribution<int> dSlack(0, 60);
        int x0 = d50(rng), y0 = d50(rng);
        static const int dirs[][2] = {{1, 0},  {0, 1},  {-1, 0}, {0, -1},
                                      {1, 1},  {1, -1}, {-1, 1}, {-1, -1},
                                      {2, 1},  {1, 2}};
        std::uniform_int_distribution<int> dDir(0, 9);
        for (int i = 0; i < n; ++i) {
            int j = dDir(rng);
            int a = dirs[j][0], b = dirs[j][1];
            planes.push_back(Plane{a, b, a * x0 + b * y0 + dSlack(rng)});
        }
    } else if (style == 1) {
        // 純隨機（多數為空集）
        std::uniform_int_distribution<int> d5(-5, 5), d30(-30, 30);
        for (int i = 0; i < n; ++i) {
            int a, b;
            while (true) {
                a = d5(rng);
                b = d5(rng);
                if (a || b) break;
            }
            planes.push_back(Plane{a, b, d30(rng)});
        }
    } else if (style == 2) {
        // 小係數、強制平行 / 同向重複
        static const int dirs[][2] = {{1, 0}, {0, 1},  {-1, 0}, {0, -1},
                                      {1, 1}, {2, 2},  {-1, -1}, {1, 2},
                                      {2, 4}};
        std::uniform_int_distribution<int> dDir(0, 8), d8(-8, 8);
        for (int i = 0; i < n; ++i) {
            int j = dDir(rng);
            planes.push_back(Plane{dirs[j][0], dirs[j][1], d8(rng)});
        }
    } else {
        // 夾帶退化平面 (0, 0, c)
        std::uniform_int_distribution<int> d100(0, 99), d5(-5, 5);
        std::uniform_int_distribution<int> d4(-4, 4), d20(-20, 20);
        for (int i = 0; i < n; ++i) {
            if (d100(rng) < 15) {
                planes.push_back(Plane{0, 0, d5(rng)});
            } else {
                int a = d4(rng), b = d4(rng);
                if (a == 0 && b == 0) b = 1;
                planes.push_back(Plane{a, b, d20(rng)});
            }
        }
    }
    return planes;
}

void checkCase(const std::vector<Plane>& planes, std::mt19937& rng) {
    NormResult norm = normalizePlanes(planes);
    if (norm.empty) {
        if (run(planes) != "0\n") fail("empty case output");
        std::optional<std::vector<Pt>> clip = hpiClip(planes);
        if (!clip || !clip->empty()) fail("empty case clip");
        return;
    }
    std::vector<Pt> ordered = runVerts(planes);
    std::vector<Pt> sortedV = ordered;
    std::sort(sortedV.begin(), sortedV.end(), cmpPt);
    // 參照 1：Sutherland–Hodgman 精確裁剪（空 / 點 / 線段 / 多邊形全覆蓋）
    std::optional<std::vector<Pt>> clip = hpiClip(planes);
    if (!clip) fail("clip invariant broken");
    if (sortedV != *clip) fail("candidate enumeration vs clip");
    checkStructure(norm.planes, ordered);
    checkMembership(planes, ordered, rng);
    // 參照 2：經典雙端棧版（退化 / 空集時回傳 None，直接跳過）
    std::optional<std::vector<Pt>> dq = hpiDeque(norm.planes);
    if (dq && sortedV != *dq) fail("candidate enumeration vs deque");
}

bool samePts(const std::vector<Pt>& got, const std::vector<Pt>& want) {
    if (got.size() != want.size()) return false;
    for (size_t i = 0; i < got.size(); ++i) {
        if (!eqPt(got[i], want[i])) return false;
    }
    return true;
}

Pt ip(long long x, long long y) { return makePt(x, y, 1); }

void runTests() {
    const long long B = BOX;
    std::string boxOut = "4\n-" + std::to_string(B) + "/1 -" + std::to_string(B) +
                         "/1\n" + std::to_string(B) + "/1 -" +
                         std::to_string(B) + "/1\n" + std::to_string(B) +
                         "/1 " + std::to_string(B) + "/1\n-" +
                         std::to_string(B) + "/1 " + std::to_string(B) + "/1\n";

    // ---- 固定用例：三角形 x ≥ 0、y ≥ 0、x + y ≤ 10 ----
    {
        std::vector<Plane> planes = {{-1, 0, 0}, {0, -1, 0}, {1, 1, 10}};
        std::string out = run(planes);
        if (out != "3\n0/1 0/1\n10/1 0/1\n0/1 10/1\n") fail("triangle: " + out);
    }
    // ---- 固定用例：正方形 [0,10]² ----
    {
        std::vector<Plane> planes = {{1, 0, 10}, {-1, 0, 0},
                                     {0, 1, 10}, {0, -1, 0}};
        std::string out = run(planes);
        if (out != "4\n0/1 0/1\n10/1 0/1\n10/1 10/1\n0/1 10/1\n") {
            fail("square: " + out);
        }
    }
    // ---- 固定用例：退化成一個點 ----
    {
        std::vector<Plane> planes = {{1, 0, 0}, {-1, 0, 0},
                                     {0, 1, 0}, {0, -1, 0}};
        if (run(planes) != "1\n0/1 0/1\n") fail("point case");
    }
    // ---- 固定用例：空集（x ≤ 0 且 x ≥ 1）----
    if (run({{1, 0, 0}, {-1, 0, -1}}) != "0\n") fail("empty case");
    // ---- 固定用例：退化平面 ----
    if (run({{0, 0, -1}}) != "0\n") fail("degenerate plane empty");
    if (run({{0, 0, 5}}) != boxOut) fail("degenerate plane full");
    // ---- 固定用例：無界半平面 x ≤ 5 → 被盒子裁成四邊形 ----
    {
        std::vector<Pt> verts = runVerts({{1, 0, 5}});
        std::sort(verts.begin(), verts.end(), cmpPt);
        std::vector<Pt> want = {ip(5, -B), ip(5, B), ip(-B, -B), ip(-B, B)};
        std::sort(want.begin(), want.end(), cmpPt);
        if (!samePts(verts, want)) fail("unbounded x<=5");
    }
    // ---- 固定用例：長條 3 ≤ x ≤ 5 ----
    {
        std::vector<Pt> verts = runVerts({{1, 0, 5}, {-1, 0, -3}});
        std::sort(verts.begin(), verts.end(), cmpPt);
        std::vector<Pt> want = {ip(3, -B), ip(3, B), ip(5, -B), ip(5, B)};
        std::sort(want.begin(), want.end(), cmpPt);
        if (!samePts(verts, want)) fail("strip 3<=x<=5");
    }
    // ---- 空輸入（沒有任何半平面）→ 整個包圍盒 ----
    if (run({}) != boxOut) fail("no planes");

    // ---- 隨機對拍：候選枚舉 vs 多邊形裁剪 + 雙端棧 + 結構 / 成員檢查 ----
    std::mt19937 rng(20261023);
    for (int t = 0; t < 220; ++t) {
        checkCase(genPlanes(rng, t % 4), rng);
    }
    // ---- 隨機對拍（較大規模，40 個平面）----
    for (int t = 0; t < 30; ++t) {
        std::mt19937 rng2(90000 + t);
        std::vector<Plane> planes = genPlanes(rng2, t % 4);
        while ((int)planes.size() < 40) {
            std::vector<Plane> more = genPlanes(rng2, t % 4);
            planes.insert(planes.end(), more.begin(), more.end());
        }
        planes.resize(40);
        checkCase(planes, rng2);
    }
}

// ---------------------------------------------------------------- IO

bool isSpaceByte(unsigned char u) {
    // 與 Python str.split()/strip() 在 latin-1 下的空白集一致
    return u == ' ' || u == '\t' || u == '\n' || u == '\r' || u == '\v' ||
           u == '\f' || (u >= 0x1c && u <= 0x1f) || u == 0x85 || u == 0xa0;
}

bool isIntToken(const std::string& t) {
    if (t.empty()) return false;
    size_t i = 0;
    if (t[0] == '+' || t[0] == '-') i = 1;
    if (i == t.size()) return false;
    for (; i < t.size(); ++i) {
        if (t[i] < '0' || t[i] > '9') return false;
    }
    return true;
}

long long parseClamped(const std::string& t) {
    // t 已通過 isIntToken；溢位一律夾到 ±CLAMP_ABS
    bool neg = false;
    size_t i = 0;
    if (t[0] == '+' || t[0] == '-') {
        neg = (t[0] == '-');
        i = 1;
    }
    unsigned long long v = 0;
    for (; i < t.size(); ++i) {
        unsigned d = (unsigned)(t[i] - '0');
        if (v > (18446744073709551615ULL - d) / 10ULL) return neg ? -CLAMP_ABS : CLAMP_ABS;
        v = v * 10ULL + d;
        if (v > (unsigned long long)CLAMP_ABS) return neg ? -CLAMP_ABS : CLAMP_ABS;
    }
    long long val = (long long)v;
    return neg ? -val : val;
}

void runIo(const std::string& raw) {
    std::vector<std::string> toks;
    std::string cur;
    for (char ch : raw) {
        if (isSpaceByte((unsigned char)ch)) {
            if (!cur.empty()) {
                toks.push_back(cur);
                cur.clear();
            }
        } else {
            cur += ch;
        }
    }
    if (!cur.empty()) toks.push_back(cur);
    size_t idx = 0;
    auto nextInt = [&](long long dflt) -> long long {
        if (idx >= toks.size()) return dflt;
        std::string t = toks[idx++];
        if (!isIntToken(t)) return dflt;
        return parseClamped(t);
    };
    auto cl = [](long long v, long long lo, long long hi) -> long long {
        if (v < lo) return lo;
        if (v > hi) return hi;
        return v;
    };

    long long n = cl(nextInt(0), 0, NMAX);
    std::vector<Plane> planes;
    for (long long i = 0; i < n; ++i) {
        long long a = cl(nextInt(0), -AMAX, AMAX);
        long long b = cl(nextInt(0), -AMAX, AMAX);
        long long c = cl(nextInt(0), -CMAX, CMAX);
        planes.push_back(Plane{a, b, c});
    }
    std::cout << run(planes);
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
        if (!isSpaceByte((unsigned char)raw[i])) {
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
