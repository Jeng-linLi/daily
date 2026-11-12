// SMAWK 演算法：完全單調矩陣（totally monotone matrix）逐列最小值搜尋。
//
// 命令（IO 模式，從 stdin 讀入）：
//   minima n m x0 ... x_{n-1} y0 ... y_{m-1}
//       -> 矩陣 A[i][j] = (x[i] - y[j])^2，x、y 遞增 => 反 Monge => 完全單調
//       -> 每列輸出一行「最小值 argmin」
//
// stdin 無內容時跑內建斷言測試，輸出 all tests passed。
//
// 核心要點：
//   * Reduce：用奇數列維護候選欄 stack；頂欄不優於新欄（>=，平手丟舊欄）
//     時彈棧，stack 滿了丟新欄。被丟欄對所有剩餘列都被某存活欄支配。
//   * 遞迴解奇數列；偶數列 argmin 被相鄰奇數列 argmin 夾住，視窗內
//     以「相鄰欄比較」二分（<= 往左收，> 往右收）。
//   * 複雜度 O(n + m) 次矩陣取值。

#include <cstdio>
#include <cstdlib>
#include <cctype>
#include <algorithm>
#include <string>
#include <vector>
#include <map>
#include <functional>

using namespace std;

typedef long long ll;

static vector<ll> gX, gY;

// 延遲取值：A(i, j) = (x[i] - y[j])^2
static ll matA(int i, int j) {
    ll d = gX[i] - gY[j];
    return d * d;
}

// SMAWK：回傳 rows 中每列的 argmin（限於 cols 中）
static map<int, int> smawk(const vector<int>& rows, const vector<int>& cols) {
    map<int, int> ans;
    if (rows.size() == 1) {
        int r = rows[0];
        int best = cols[0];
        for (int c : cols) {
            if (matA(r, c) < matA(r, best)) best = c;
        }
        ans[r] = best;
        return ans;
    }

    // --- Reduce：奇數列篩選候選欄 ---
    vector<int> odd;
    for (size_t t = 1; t < rows.size(); t += 2) odd.push_back(rows[t]);
    vector<int> stack;
    for (int c : cols) {
        while (!stack.empty() &&
               matA(odd[stack.size() - 1], stack.back()) >= matA(odd[stack.size() - 1], c)) {
            stack.pop_back();
        }
        if ((int)stack.size() < (int)odd.size()) stack.push_back(c);
    }
    map<int, int> oddAns = smawk(odd, stack);

    // --- Interpolate：偶數列在相鄰奇數列答案之間二分 ---
    map<int, int> colPos;
    for (int t = 0; t < (int)cols.size(); ++t) colPos[cols[t]] = t;
    ans = oddAns;
    for (int idx = 0; idx < (int)rows.size(); ++idx) {
        if (idx % 2 == 1) continue;
        int r = rows[idx];
        int lo = (idx > 0) ? colPos[oddAns[rows[idx - 1]]] : 0;
        int hi = (idx + 1 < (int)rows.size()) ? colPos[oddAns[rows[idx + 1]]]
                                              : (int)cols.size() - 1;
        while (lo < hi) {
            int mid = lo + (hi - lo) / 2;
            if (matA(r, cols[mid]) <= matA(r, cols[mid + 1])) hi = mid;
            else lo = mid + 1;
        }
        ans[r] = cols[lo];
    }
    return ans;
}

static vector<pair<ll, int>> nearestSquares(const vector<ll>& x, const vector<ll>& y) {
    vector<pair<ll, int>> out;
    if (x.empty() || y.empty()) return out;
    // 先把 y 去重成嚴格遞增：y 有重複時 (x-y)^2 依索引看在下降段會出現
    // 平台（16,16,16,9,9,4,4,1），相鄰比較二分會停在平台邊緣而漏掉最小值。
    // 去重後每列嚴格單峰（至多一組對稱平手），二分才成立。
    vector<ll> yd;
    map<ll, int> firstIdx;
    for (int j = 0; j < (int)y.size(); ++j) {
        if (!firstIdx.count(y[j])) {
            firstIdx[y[j]] = j;
            yd.push_back(y[j]);
        }
    }
    gX = x;
    gY = yd;
    vector<int> rows(x.size()), cols(yd.size());
    for (int i = 0; i < (int)x.size(); ++i) rows[i] = i;
    for (int j = 0; j < (int)yd.size(); ++j) cols[j] = j;
    map<int, int> argmins = smawk(rows, cols);
    for (int i = 0; i < (int)x.size(); ++i) {
        int j = argmins[i];
        ll d = x[i] - yd[j];
        out.push_back(make_pair(d * d, firstIdx[yd[j]]));
    }
    return out;
}

// O(nm) 暴力（取最左 argmin），僅供內建測試對拍
static vector<pair<ll, int>> brute(const vector<ll>& x, const vector<ll>& y) {
    vector<pair<ll, int>> out;
    for (ll xi : x) {
        int bestJ = 0;
        ll bd = xi - y[0];
        ll bestV = bd * bd;
        for (int j = 1; j < (int)y.size(); ++j) {
            ll d = xi - y[j];
            ll v = d * d;
            if (v < bestV) {
                bestV = v;
                bestJ = j;
            }
        }
        out.push_back(make_pair(bestV, bestJ));
    }
    return out;
}

// 嚴格整數解析：^[+-]?[0-9]+$，非法 token 視為輸入結束
static bool tryParse(const string& s, ll& v) {
    if (s.empty()) return false;
    size_t p = 0;
    if (s[0] == '+' || s[0] == '-') p = 1;
    if (p >= s.size()) return false;
    for (size_t t = p; t < s.size(); ++t) {
        if (s[t] < '0' || s[t] > '9') return false;
    }
    v = strtoll(s.c_str(), NULL, 10);
    return true;
}

#ifdef _WIN32
#include <io.h>
#include <fcntl.h>
#endif

int main() {
#ifdef _WIN32
    _setmode(_fileno(stdout), _O_BINARY);
#endif
    // 讀全部 stdin
    string data, line;
    char buf[65536];
    size_t got;
    while ((got = fread(buf, 1, sizeof(buf), stdin)) > 0) data.append(buf, got);

    bool hasInput = false;
    for (char c : data) {
        if (!isspace((unsigned char)c)) {
            hasInput = true;
            break;
        }
    }

    if (hasInput) {
        vector<string> toks;
        {
            string cur;
            for (char c : data) {
                if (isspace((unsigned char)c)) {
                    if (!cur.empty()) {
                        toks.push_back(cur);
                        cur.clear();
                    }
                } else {
                    cur += c;
                }
            }
            if (!cur.empty()) toks.push_back(cur);
        }
        size_t pos = 0;
        auto nxt = [&](string& outTok) -> bool {
            if (pos >= toks.size()) return false;
            outTok = toks[pos++];
            return true;
        };
        string cmd;
        while (nxt(cmd)) {
            if (cmd != "minima") break;
            string t;
            ll n, m;
            if (!nxt(t) || !tryParse(t, n)) break;
            if (!nxt(t) || !tryParse(t, m)) break;
            vector<ll> x, y;
            bool ok = true;
            for (ll i = 0; i < n && ok; ++i) {
                ll v;
                if (nxt(t) && tryParse(t, v)) x.push_back(v);
                else ok = false;
            }
            for (ll j = 0; j < m && ok; ++j) {
                ll v;
                if (nxt(t) && tryParse(t, v)) y.push_back(v);
                else ok = false;
            }
            if (!ok) break;
            vector<pair<ll, int>> res = nearestSquares(x, y);
            for (size_t i = 0; i < res.size(); ++i) {
                printf("%lld %d\n", res[i].first, res[i].second);
            }
        }
        fflush(stdout);
        return 0;
    }

    // ---- 內建斷言測試 ----
    {
        vector<ll> x = {1, 4, 9}, y = {2, 7};
        vector<pair<ll, int>> res = nearestSquares(x, y);
        // 期望值以實跑為準：x=4 最近 y=2、x=9 最近 y=7
        if (!(res.size() == 3 && res[0].first == 1 && res[0].second == 0 &&
              res[1].first == 4 && res[1].second == 0 &&
              res[2].first == 4 && res[2].second == 1)) return 1;
        if (!nearestSquares({}, {1, 2}).empty()) return 1;
        if (!nearestSquares({1}, {}).empty()) return 1;
    }

    {
        unsigned int st = 20261111u;
        auto rnd = [&](unsigned int mm) -> unsigned int {
            st = st * 1103515245u + 12345u;
            return (st >> 1) % mm;
        };
        for (int tc = 0; tc < 300; ++tc) {
            int n = (int)rnd(30) + 1;
            int m = (int)rnd(30) + 1;
            vector<ll> x, y;
            for (int i = 0; i < n; ++i) x.push_back((ll)rnd(10) + rnd(15));
            for (int j = 0; j < m; ++j) y.push_back((ll)rnd(10) + rnd(15));
            sort(x.begin(), x.end());
            sort(y.begin(), y.end());
            vector<pair<ll, int>> got = nearestSquares(x, y);
            vector<pair<ll, int>> want = brute(x, y);
            if (got.size() != want.size()) return 1;
            for (int i = 0; i < n; ++i) {
                if (got[i].first != want[i].first) return 1;  // 最小值必須一致
                ll d = x[i] - y[got[i].second];
                if (d * d != got[i].first) return 1;          // 位置自洽
            }
        }
        for (int tc = 0; tc < 20; ++tc) {
            int n = (int)rnd(300) + 1;
            int m = (int)rnd(300) + 1;
            vector<ll> x, y;
            for (int i = 0; i < n; ++i) x.push_back((ll)rnd(1000000u));
            for (int j = 0; j < m; ++j) y.push_back((ll)rnd(1000000u));
            sort(x.begin(), x.end());
            sort(y.begin(), y.end());
            vector<pair<ll, int>> got = nearestSquares(x, y);
            vector<pair<ll, int>> want = brute(x, y);
            if (got.size() != want.size()) return 1;
            for (int i = 0; i < n; ++i) {
                if (got[i].first != want[i].first) return 1;  // 最小值必須一致
                ll d = x[i] - y[got[i].second];
                if (d * d != got[i].first) return 1;          // 位置自洽
            }
        }
    }

    printf("all tests passed\n");
    return 0;
}
