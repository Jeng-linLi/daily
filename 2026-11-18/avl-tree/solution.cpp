// AVL 樹（自平衡二元搜尋樹）與順序統計。
//
// 命令（IO 模式，stdin 第一個 token 為操作數 n，之後 n 個操作）：
//   i x  插入 x          -> 輸出插入後元素總數
//   e x  刪除一個 x      -> 刪除後元素總數
//   f x  是否存在 x      -> 1 / 0
//   r x  嚴格小於 x 的個數
//   k x  第 x 小（1-based），越界輸出 -1
//   p x  嚴格小於 x 的最大者，無則 -1
//   s x  嚴格大於 x 的最小者，無則 -1
//   o    中序輸出全部元素（含重複），空則輸出 -
//
// stdin 無內容時跑內建對拍（與排序 vector + lower_bound 比較），輸出 all tests passed。

#include <cassert>
#include <cstdio>
#include <string>
#include <vector>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <cctype>

#ifdef _WIN32
#include <fcntl.h>
#include <io.h>
#endif

using namespace std;

static const long long CAP = 1000000000000000LL;

static bool isInt(const string &s) {
    if (s.empty()) return false;
    size_t i = 0;
    if (s[0] == '+' || s[0] == '-') i = 1;
    if (i == s.size()) return false;
    for (; i < s.size(); i++) if (!isdigit((unsigned char)s[i])) return false;
    return true;
}

static long long parseCap(const string &s) {
    long long v = 0;
    bool neg = false;
    size_t i = 0;
    if (s[0] == '+' || s[0] == '-') { neg = (s[0] == '-'); i = 1; }
    for (; i < s.size(); i++) v = v * 10 + (s[i] - '0');
    if (v > CAP) v = CAP;
    return neg ? -v : v;
}

// ------------------------------- AVL -------------------------------

struct AVL {
    vector<long long> key, cnt, sz;
    vector<int> h, lc, rc;

    AVL() {
        key.push_back(0); cnt.push_back(0); sz.push_back(0);
        h.push_back(0); lc.push_back(0); rc.push_back(0);
    }

    int newNode(long long k) {
        key.push_back(k); cnt.push_back(1); sz.push_back(1);
        h.push_back(1); lc.push_back(0); rc.push_back(0);
        return (int)key.size() - 1;
    }

    void pull(int v) {
        if (v == 0) return;
        h[v] = 1 + max(h[lc[v]], h[rc[v]]);
        sz[v] = cnt[v] + sz[lc[v]] + sz[rc[v]];
    }

    int rotRight(int y) {
        int x = lc[y];
        lc[y] = rc[x];
        rc[x] = y;
        pull(y); pull(x);
        return x;
    }

    int rotLeft(int x) {
        int y = rc[x];
        rc[x] = lc[y];
        lc[y] = x;
        pull(x); pull(y);
        return y;
    }

    int rebalance(int v) {
        if (v == 0) return 0;
        pull(v);
        int bf = h[lc[v]] - h[rc[v]];
        if (bf > 1) {
            if (h[lc[lc[v]]] < h[rc[lc[v]]]) lc[v] = rotLeft(lc[v]);
            return rotRight(v);
        }
        if (bf < -1) {
            if (h[rc[rc[v]]] < h[lc[rc[v]]]) rc[v] = rotRight(rc[v]);
            return rotLeft(v);
        }
        return v;
    }

    int insert(int v, long long k) {
        if (v == 0) return newNode(k);
        if (k == key[v]) cnt[v]++;
        else if (k < key[v]) lc[v] = insert(lc[v], k);
        else rc[v] = insert(rc[v], k);
        return rebalance(v);
    }

    // 移除 v 子樹的最小節點，回傳 (新根, key, cnt)。
    pair<int, pair<long long, long long> > popMin(int v) {
        if (lc[v] == 0) return {rc[v], {key[v], cnt[v]}};
        auto r = popMin(lc[v]);
        lc[v] = r.first;
        return {rebalance(v), r.second};
    }

    int erase(int v, long long k) {
        if (v == 0) return 0;
        if (k < key[v]) lc[v] = erase(lc[v], k);
        else if (k > key[v]) rc[v] = erase(rc[v], k);
        else {
            if (cnt[v] > 1) cnt[v]--;
            else if (lc[v] == 0 || rc[v] == 0) v = lc[v] ? lc[v] : rc[v];
            else {
                auto r = popMin(rc[v]);
                rc[v] = r.first;
                key[v] = r.second.first;
                cnt[v] = r.second.second;
            }
        }
        return rebalance(v);
    }

    bool find(int v, long long k) const {
        while (v) {
            if (k == key[v]) return true;
            v = (k < key[v]) ? lc[v] : rc[v];
        }
        return false;
    }

    long long rankLess(int v, long long k) const {
        long long ans = 0;
        while (v) {
            if (k <= key[v]) v = lc[v];
            else { ans += sz[lc[v]] + cnt[v]; v = rc[v]; }
        }
        return ans;
    }

    long long kth(int v, long long k) const {
        if (k < 1 || k > sz[v]) return -1;
        while (v) {
            long long ls = sz[lc[v]];
            if (k <= ls) v = lc[v];
            else if (k > ls + cnt[v]) { k -= ls + cnt[v]; v = rc[v]; }
            else return key[v];
        }
        return -1;
    }

    long long pred(int v, long long k) const {
        long long ans = -1;
        bool has = false;
        while (v) {
            if (key[v] < k) { ans = key[v]; has = true; v = rc[v]; }
            else v = lc[v];
        }
        (void)has;
        return ans;
    }

    long long succ(int v, long long k) const {
        long long ans = -1;
        bool has = false;
        while (v) {
            if (key[v] > k) { ans = key[v]; has = true; v = lc[v]; }
            else v = rc[v];
        }
        (void)has;
        return ans;
    }

    void inorder(int v, vector<long long> &out) const {
        if (v == 0) return;
        inorder(lc[v], out);
        for (long long i = 0; i < cnt[v]; i++) out.push_back(key[v]);
        inorder(rc[v], out);
    }

    int check(int v, long long lo, long long hi, bool hasLo, bool hasHi) {
        if (v == 0) return 0;
        if (hasLo) assert(key[v] > lo);
        if (hasHi) assert(key[v] < hi);
        int lh = check(lc[v], lo, key[v], hasLo, true);
        int rh = check(rc[v], key[v], hi, true, hasHi);
        assert(abs(lh - rh) <= 1);
        assert(h[v] == 1 + max(lh, rh));
        assert(sz[v] == cnt[v] + sz[lc[v]] + sz[rc[v]]);
        return h[v];
    }
};

// ------------------------- 輸出：一邊從樹算，一邊從排序陣列算 -------------------------

// 只從 AVL 樹計算輸出（IO 模式與對拍都用這個）。
static string treeOut(AVL &t, int root, const string &op, long long arg) {
    if (op == "i") return to_string(t.sz[root]);
    if (op == "e") return to_string(t.sz[root]);
    if (op == "f") return t.find(root, arg) ? "1" : "0";
    if (op == "r") return to_string(t.rankLess(root, arg));
    if (op == "k") return to_string(t.kth(root, arg));
    if (op == "p") return to_string(t.pred(root, arg));
    if (op == "s") return to_string(t.succ(root, arg));
    vector<long long> seq;
    t.inorder(root, seq);
    ostringstream oss;
    for (size_t i = 0; i < seq.size(); i++) { if (i) oss << ' '; oss << seq[i]; }
    string s = oss.str();
    return s.empty() ? "-" : s;
}

// 獨立參考：從排序 vector 計算輸出，同時更新它。
static string refOut(vector<long long> &ref, const string &op, long long arg) {
    vector<long long>::iterator it;
    if (op == "i") {
        it = lower_bound(ref.begin(), ref.end(), arg);
        ref.insert(it, arg);
        return to_string(ref.size());
    }
    if (op == "e") {
        it = lower_bound(ref.begin(), ref.end(), arg);
        if (it != ref.end() && *it == arg) ref.erase(it);
        return to_string(ref.size());
    }
    if (op == "f") {
        it = lower_bound(ref.begin(), ref.end(), arg);
        return (it != ref.end() && *it == arg) ? "1" : "0";
    }
    if (op == "r") return to_string(lower_bound(ref.begin(), ref.end(), arg) - ref.begin());
    if (op == "k") {
        if (arg < 1 || arg > (long long)ref.size()) return "-1";
        return to_string(ref[(size_t)(arg - 1)]);
    }
    if (op == "p") {
        it = lower_bound(ref.begin(), ref.end(), arg);
        return (it == ref.begin()) ? "-1" : to_string(*(it - 1));
    }
    if (op == "s") {
        it = upper_bound(ref.begin(), ref.end(), arg);
        return (it == ref.end()) ? "-1" : to_string(*it);
    }
    ostringstream oss;
    for (size_t i = 0; i < ref.size(); i++) { if (i) oss << ' '; oss << ref[i]; }
    string s = oss.str();
    return s.empty() ? "-" : s;
}

// 對拍用：同時套用 op 到樹與排序陣列，回傳 (樹輸出, 參考輸出)。
static pair<string, string> applyOp(AVL &t, int &root, vector<long long> &ref,
                                    const string &op, long long arg) {
    if (op == "i") root = t.insert(root, arg);
    else if (op == "e") root = t.erase(root, arg);
    return {treeOut(t, root, op, arg), refOut(ref, op, arg)};
}

// ------------------------------- IO 模式 -------------------------------

static void run_io(const vector<string> &toks) {
    if (toks.empty() || !isInt(toks[0])) return;
    long long n = parseCap(toks[0]);
    if (n < 0) n = 0;
    if (n > 200000) n = 200000;
    AVL t;
    int root = 0;
    string out;
    size_t pos = 1;
    for (long long q = 0; q < n; q++) {
        if (pos >= toks.size()) break;
        string op = toks[pos++];
        if (op != "i" && op != "e" && op != "f" && op != "r" &&
            op != "k" && op != "p" && op != "s" && op != "o") break;
        long long arg = 0;
        if (op != "o") {
            if (pos >= toks.size() || !isInt(toks[pos])) break;
            arg = parseCap(toks[pos++]);
        }
        string line;
        if (op == "i") root = t.insert(root, arg);
        else if (op == "e") root = t.erase(root, arg);
        line = treeOut(t, root, op, arg);
        out += line + "\n";
    }
    fputs(out.c_str(), stdout);
}

// ------------------------------ 內建測試 ------------------------------

static unsigned long long g_state = 20261118ULL;
static int rnd(int m) {
    g_state = (g_state * 1103515245ULL + 12345ULL) & 0x7FFFFFFFULL;
    return (int)(g_state % (unsigned)m);
}

static void run_tests() {
    // 固定情境（輸出已實跑確認）
    {
        AVL t;
        int root = 0;
        vector<long long> ref;
        vector<pair<string, long long> > ops = {
            {"i",5},{"i",3},{"i",7},{"i",3},{"o",0},
            {"f",3},{"f",4},{"r",6},{"k",2},{"k",5},
            {"p",6},{"s",6},{"e",3},{"o",0},{"e",3},{"o",0}
        };
        vector<string> want = {
            "1","2","3","4","3 3 5 7","1","0","3","3","-1",
            "5","7","3","3 5 7","2","5 7"
        };
        for (size_t i = 0; i < ops.size(); i++) {
            pair<string, string> r = applyOp(t, root, ref, ops[i].first, ops[i].second);
            assert(r.first == want[i]);        // 樹輸出 == 固定期望
            assert(r.second == want[i]);       // 參考輸出 == 固定期望
        }
    }

    // 隨機對拍：AVL vs 排序 vector，並檢查不變量
    for (int tc = 0; tc < 160; tc++) {
        AVL t;
        int root = 0;
        vector<long long> ref;
        for (int step = 0; step < 120; step++) {
            static const char ops[] = "iefrkpso";
            string op(1, ops[rnd(8)]);
            long long arg = rnd(41) - 20;
            pair<string, string> r = applyOp(t, root, ref, op, arg);
            assert(r.first == r.second);       // 樹 vs 排序陣列逐操作對拍
            t.check(root, 0, 0, false, false);
        }
    }

    // 單調插入 / 刪除
    {
        AVL t;
        int root = 0;
        vector<long long> ref;
        for (int i = 0; i < 400; i++) {
            pair<string, string> r = applyOp(t, root, ref, "i", i);
            assert(r.first == to_string(i + 1));
            assert(r.first == r.second);
        }
        for (int i = 0; i < 400; i++) {
            pair<string, string> r = applyOp(t, root, ref, "e", i);
            assert(r.first == to_string(399 - i));
            assert(r.first == r.second);
        }
    }
}

int main() {
#ifdef _WIN32
    _setmode(_fileno(stdout), _O_BINARY);
#endif
    ostringstream buf;
    buf << cin.rdbuf();
    string input = buf.str();
    if (!input.empty() && input.find_first_not_of(" \t\r\n") != string::npos) {
        istringstream iss(input);
        vector<string> toks;
        string x;
        while (iss >> x) toks.push_back(x);
        run_io(toks);
    } else {
        run_tests();
        fputs("all tests passed\n", stdout);
    }
    return 0;
}
