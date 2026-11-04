// 可持久化 01-Trie（Persistent Binary Trie）：區間最大異或
//
// 題意：
//     給定長度 n 的數列 a[1..n]（0 <= a_i < 2^20）與 q 筆查詢：
//         l r x —— 求 max{ a_i XOR x : l <= i <= r }
//     輸入：第一行 n q，第二行 n 個整數，接下來 q 行每行 l r x（1 起算）。
//     輸出：每筆查詢一行答案。
//
// 思路：
//     版本 i 的 Trie = 依次插入 a_1..a_i 的 01-Trie（版本 0 為空樹）。
//     插入時只複製被修改路徑上的節點，其餘孩子共享舊版本指針，
//     每個版本只多 O(B) 個節點（B = 20）。節點維護 cnt = 子樹內插入過的數的個數。
//     查詢 [l, r] 時在新版本 root[r] 與舊版本 root[l-1] 上同時向下走：
//     想走的方向 d（與 x 該位相反）若 cntR[d] - cntOld[d] > 0，代表區間內
//     存在帶該前綴的數 → 走 d 並把答案該位置 1；否則走 x 該位的方向。
//     每筆查詢 O(B)。
//
// 複雜度：
//     建樹 O(nB) 節點 / 時間；查詢每筆 O(B)；空間 O(nB)。
//
// 應用場景：
//     區間異或查詢（遊戲數值對比 / 訊號指紋比對）、離線日誌時間視窗分析、
//     資料庫位圖索引的範圍 + 異或組合查詢。

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

#ifdef _WIN32
#include <io.h>
#include <fcntl.h>
#endif

typedef long long ll;

const int BITS = 20;
const int VAL_MAX = 1 << BITS;   // 1048576
const int NMAX = 50000;
const int QMAX = 50000;

struct PersistentTrie {
    // 節點 0 是空節點：ch = {0,0}、cnt = 0
    std::vector<int> ch0, ch1, cnt;
    std::vector<int> roots;

    PersistentTrie() {
        ch0.push_back(0); ch1.push_back(0); cnt.push_back(0);
        roots.push_back(0);
    }
    int new_node() {
        ch0.push_back(0); ch1.push_back(0); cnt.push_back(0);
        return (int)ch0.size() - 1;
    }
    // 在版本 prev 基礎上插入 val，回傳新版本根
    int insert(int prev, int val) {
        int cur = new_node();
        int root = cur;
        int p = prev;
        cnt[cur] = cnt[p] + 1;
        for (int b = BITS - 1; b >= 0; --b) {
            int d = (val >> b) & 1;
            // 另一側孩子直接共享舊版本
            if (d == 0) ch1[cur] = ch1[p];
            else        ch0[cur] = ch0[p];
            int q = d == 0 ? ch0[p] : ch1[p];
            int nxt = new_node();
            if (d == 0) ch0[cur] = nxt; else ch1[cur] = nxt;
            cur = nxt;
            p = q;
            cnt[cur] = cnt[p] + 1;
        }
        return root;
    }
    void build(const std::vector<int> &vals) {
        for (size_t i = 0; i < vals.size(); ++i)
            roots.push_back(insert(roots.back(), vals[i]));
    }
    // max{ a_i XOR x : l <= i <= r }，l / r 為 1 起算
    ll query_max_xor(int l, int r, int x) const {
        int u = roots[r];
        int v = roots[l - 1];
        ll ans = 0;
        for (int b = BITS - 1; b >= 0; --b) {
            int d = (x >> b) & 1;
            int want = d ^ 1;
            int cu = d == 0 ? ch0[u] : ch1[u];
            int cv = d == 0 ? ch0[v] : ch1[v];
            int cr = want == 0 ? ch0[u] : ch1[u];
            int co = want == 0 ? ch0[v] : ch1[v];
            if (cnt[cr] - cnt[co] > 0) {
                ans |= 1LL << b;
                u = cr; v = co;
            } else {
                u = cu; v = cv;
            }
        }
        return ans;
    }
};

// 參考實作：O(r-l+1) 枚舉
static ll brute_force(const std::vector<int> &a, int l, int r, int x) {
    ll best = -1;
    for (int i = l - 1; i < r; ++i) {
        ll cur = (ll)(a[i] ^ x);
        if (cur > best) best = cur;
    }
    return best;
}

// 自製 LCG（與 Python 版同一序列）
static unsigned int lcg_state = 20261103u;
static unsigned int lcg_next() {
    lcg_state = (unsigned int)((ll)lcg_state * 1103515245LL + 12345LL) & 0x7FFFFFFFu;
    return lcg_state;
}

static void emit(const std::string &s) { fwrite(s.data(), 1, s.size(), stdout); }

static bool run_one_test(const std::vector<int> &vals, int n,
                         const std::vector<int> &qs /* 3 個一組 */) {
    PersistentTrie trie;
    trie.build(vals);
    for (size_t i = 0; i + 3 <= qs.size(); i += 3) {
        int l = qs[i], r = qs[i + 1], x = qs[i + 2];
        ll got = trie.query_max_xor(l, r, x);
        ll want = brute_force(vals, l, r, x);
        if (got != want) return false;
    }
    return true;
}

static void run_tests() {
    // ---- 固定用例 ----
    {
        PersistentTrie t; std::vector<int> a; a.push_back(1); a.push_back(2); a.push_back(3);
        t.build(a);
        if (t.query_max_xor(1, 3, 0) != 3) { emit("FAIL 1\n"); return; }
        if (t.query_max_xor(1, 3, 1) != 3) { emit("FAIL 2\n"); return; }
        if (t.query_max_xor(2, 2, 5) != (2 ^ 5)) { emit("FAIL 3\n"); return; }
        if (t.query_max_xor(3, 3, 7) != (3 ^ 7)) { emit("FAIL 4\n"); return; }
    }
    {
        PersistentTrie t; std::vector<int> a(3, 0);
        t.build(a);
        if (t.query_max_xor(1, 3, 0) != 0) { emit("FAIL 5\n"); return; }
        if (t.query_max_xor(1, 3, 7) != 7) { emit("FAIL 6\n"); return; }
    }
    {
        PersistentTrie t; std::vector<int> a; a.push_back(10); a.push_back(10);
        a.push_back(10); a.push_back(3);
        t.build(a);
        if (t.query_max_xor(1, 3, 0) != 10) { emit("FAIL 7\n"); return; }
        if (t.query_max_xor(1, 4, 9) != 10) { emit("FAIL 8\n"); return; }
    }

    // ---- 與暴力對拍 ----
    for (int it = 0; it < 200; ++it) {
        int n = 1 + (int)(lcg_next() % 12u);
        std::vector<int> vals;
        for (int i = 0; i < n; ++i) vals.push_back((int)(lcg_next() % 64u));
        std::vector<int> qs;
        for (int j = 0; j < 10; ++j) {
            int l = 1 + (int)(lcg_next() % (unsigned)n);
            int r = l + (int)(lcg_next() % (unsigned)(n - l + 1));
            int x = (int)(lcg_next() % 64u);
            qs.push_back(l); qs.push_back(r); qs.push_back(x);
        }
        if (!run_one_test(vals, n, qs)) { emit("FAIL random\n"); return; }
    }
    // 大值域（< 2^20）
    for (int it = 0; it < 60; ++it) {
        int n = 1 + (int)(lcg_next() % 60u);
        std::vector<int> vals;
        for (int i = 0; i < n; ++i) vals.push_back((int)(lcg_next() % (unsigned)VAL_MAX));
        std::vector<int> qs;
        for (int j = 0; j < 8; ++j) {
            int l = 1 + (int)(lcg_next() % (unsigned)n);
            int r = l + (int)(lcg_next() % (unsigned)(n - l + 1));
            int x = (int)(lcg_next() % (unsigned)VAL_MAX);
            qs.push_back(l); qs.push_back(r); qs.push_back(x);
        }
        if (!run_one_test(vals, n, qs)) { emit("FAIL random big\n"); return; }
    }
}

// ---- IO 模式 ----
// 嚴格整數解析：^[+-]?[0-9]+$，超出 ±10^18 回 false
static bool parse_ll(const char *s, ll &out) {
    if (!s || !*s) return false;
    int i = 0;
    if (s[0] == '+' || s[0] == '-') i = 1;
    if (!s[i]) return false;
    for (int j = i; s[j]; ++j)
        if (s[j] < '0' || s[j] > '9') return false;
    if (strlen(s) - i > 19) return false;
    long long v = 0;
    for (int j = i; s[j]; ++j) v = v * 10 + (s[j] - '0');
    if (v > 1000000000000000000LL) return false;
    out = (s[0] == '-') ? -v : v;
    return true;
}

int main() {
#ifdef _WIN32
    _setmode(_fileno(stdout), _O_BINARY);
#endif
    std::string raw;
    {
        char buf[1 << 16];
        size_t n;
        while ((n = fread(buf, 1, sizeof(buf), stdin)) > 0) raw.append(buf, n);
    }
    std::vector<std::string> toks;
    {
        std::string cur;
        for (size_t i = 0; i <= raw.size(); ++i) {
            char ch = (i < raw.size()) ? raw[i] : ' ';
            if (ch == ' ' || ch == '\n' || ch == '\r' || ch == '\t') {
                if (!cur.empty()) { toks.push_back(cur); cur.clear(); }
            } else {
                cur += ch;
            }
        }
    }
    if (toks.empty()) {
        run_tests();
        emit("all tests passed\n");
        return 0;
    }

    size_t pos = 0;
    ll n_ll, q_ll;
    if (pos >= toks.size() || !parse_ll(toks[pos++].c_str(), n_ll) ||
        n_ll < 1 || n_ll > NMAX) { emit("error\n"); return 0; }
    if (pos >= toks.size() || !parse_ll(toks[pos++].c_str(), q_ll) ||
        q_ll < 1 || q_ll > QMAX) { emit("error\n"); return 0; }
    int n = (int)n_ll, q = (int)q_ll;

    std::vector<int> a;
    for (int i = 0; i < n; ++i) {
        ll v;
        if (pos >= toks.size() || !parse_ll(toks[pos++].c_str(), v) ||
            v < 0 || v >= VAL_MAX) { emit("error\n"); return 0; }
        a.push_back((int)v);
    }

    PersistentTrie trie;
    trie.build(a);

    std::string out;
    char buf[32];
    for (int i = 0; i < q; ++i) {
        ll l, r, x;
        if (pos >= toks.size() || !parse_ll(toks[pos++].c_str(), l) ||
            pos >= toks.size() || !parse_ll(toks[pos++].c_str(), r) ||
            pos >= toks.size() || !parse_ll(toks[pos++].c_str(), x)) {
            emit("error\n"); return 0;
        }
        if (!(1 <= l && l <= r && r <= n) || !(0 <= x && x < VAL_MAX)) {
            emit("error\n"); return 0;
        }
        sprintf(buf, "%lld\n", trie.query_max_xor((int)l, (int)r, (int)x));
        out += buf;
    }
    emit(out);
    return 0;
}
