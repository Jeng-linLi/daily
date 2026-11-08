// 博弈論與 Sprague-Grundy 函數（mex DP / 組合遊戲和）
//
// 題意：多筆遊戲，每筆有 k「堆」，每堆屬於兩種類型之一：
//   type 0 取石遊戲：一次可從該堆取 1..3 顆（取完即消失），g(x) = x mod 4；
//   type 1 拆分遊戲：一次可把大小 x 的堆拆成 i 與 x-i（1 <= i < x）兩堆，
//                    g(x) = mex{ g(i) xor g(x-i) }（組合遊戲和定理）。
// 組合遊戲和的整體 SG = 各堆 SG 的 XOR；XOR != 0 ⟺ 先手必勝。
//
// IO 模式（stdin 有內容時）：
//   t
//   每筆：k，接著 k 行「type x」（type ∈ {0,1}，0 <= x <= 3000）
//   每筆輸出兩行：第一行 k 個 SG 值；第二行 first 或 second。
// 無 stdin 時跑內建斷言（與極小極大記憶化搜索對拍）並輸出 all tests passed。
#include <cstdio>
#include <cctype>
#include <cstring>
#include <sstream>
#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <set>
#include <algorithm>
#include <cassert>

using namespace std;
typedef long long ll;

static const int LIMIT = 3000;
static int gt[LIMIT + 1], gs[LIMIT + 1];

static int mex(const set<int>& s) {
    int m = 0;
    while (s.count(m)) m++;
    return m;
}

static void buildTables() {
    gt[0] = gs[0] = 0;
    for (int x = 1; x <= LIMIT; x++) {
        set<int> s1;
        for (int r = 1; r <= 3; r++)
            if (x - r >= 0) s1.insert(gt[x - r]);
        gt[x] = mex(s1);
        set<int> s2;
        for (int i = 1; i < x; i++) s2.insert(gs[i] ^ gs[x - i]);
        gs[x] = mex(s2);
    }
}

static inline int sg(int t, int x) { return t == 0 ? gt[x] : gs[x]; }

// ---- 極小極大暴力（局面 = 排序後的 (type,x) 序列，x=0 的堆丟棄）----
typedef vector<pair<int, int> > State;

static map<State, char> winMemo;

static State normalize(const State& st) {
    State r;
    for (size_t i = 0; i < st.size(); i++)
        if (st[i].second > 0) r.push_back(st[i]);
    sort(r.begin(), r.end());
    return r;
}

static bool win(const State& raw) {
    State key = normalize(raw);
    if (key.empty()) return false;  // 無子可走：當前玩家輸
    map<State, char>::iterator it = winMemo.find(key);
    if (it != winMemo.end()) return it->second != 0;
    bool res = false;
    for (size_t idx = 0; idx < key.size() && !res; idx++) {
        int t = key[idx].first, x = key[idx].second;
        State base;
        for (size_t j = 0; j < key.size(); j++)
            if (j != idx) base.push_back(key[j]);
        if (t == 0) {
            for (int r = 1; r <= 3 && !res; r++) {
                if (x - r < 0) continue;
                State nxt = base;
                if (x - r > 0) nxt.push_back(make_pair(0, x - r));
                if (!win(nxt)) res = true;
            }
        } else if (x >= 2) {
            for (int i = 1; i < x && !res; i++) {
                State nxt = base;
                nxt.push_back(make_pair(1, i));
                nxt.push_back(make_pair(1, x - i));
                if (!win(nxt)) res = true;
            }
        }
    }
    winMemo[key] = res ? 1 : 0;
    return res;
}

static void runAllTests() {
    // 1) 取石閉式解：g(x) = x mod 4
    for (int x = 0; x <= LIMIT; x++) assert(gt[x] == x % 4);
    // 2) split 的 SG 交替性（實跑規律）：偶數 1、奇數 0
    for (int x = 1; x <= LIMIT; x++) assert(gs[x] == (x % 2 == 0 ? 1 : 0));
    // 3) 混合局面：XOR != 0 ⟺ 先手必勝（與極小極大暴力對拍）
    unsigned long long seed = 20261108ULL;
    auto rnd = [&](int n) -> int {
        seed = seed * 1103515245ULL + 12345ULL;
        seed &= 0x7fffffffULL;
        return (int)(seed % (unsigned)n);
    };
    for (int tc = 0; tc < 400; tc++) {
        int k = 1 + rnd(4);
        State st;
        int xorAll = 0;
        for (int i = 0; i < k; i++) {
            int t = rnd(2), x = rnd(11);
            st.push_back(make_pair(t, x));
            xorAll ^= sg(t, x);
        }
        assert((xorAll != 0) == win(st));
    }
}

int main() {
    buildTables();
    stringstream buf;
    buf << cin.rdbuf();
    string input = buf.str();
    bool anyTok = false;
    for (size_t i = 0; i < input.size(); i++)
        if (!isspace((unsigned char)input[i])) { anyTok = true; break; }
    if (!anyTok) { runAllTests(); printf("all tests passed\n"); return 0; }
    istringstream it(input);
    long long t;
    if (!(it >> t)) { return 0; }
    ostringstream oss;
    for (int tc = 0; tc < t; tc++) {
        long long k;
        it >> k;
        vector<int> vals;
        int xorAll = 0;
        for (long long i = 0; i < k; i++) {
            long long tp, x;
            it >> tp >> x;
            if ((tp != 0 && tp != 1) || x < 0 || x > LIMIT) { tp = 0; x = 0; }
            int g = sg((int)tp, (int)x);
            vals.push_back(g);
            xorAll ^= g;
        }
        for (size_t i = 0; i < vals.size(); i++) {
            if (i) oss << ' ';
            oss << vals[i];
        }
        oss << '\n' << (xorAll != 0 ? "first" : "second") << '\n';
    }
    fputs(oss.str().c_str(), stdout);
    return 0;
}
