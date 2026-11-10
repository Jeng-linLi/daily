// Stern-Brocot 樹與 Farey 數列（C++17）
//
// 三種命令（IO 模式，從 stdin 讀入）：
//   farey n        -> 輸出 F_n 的全部既約分數（升冪，空格分隔）
//   path p q       -> 輸出 Stern-Brocot 樹中從 1/1 走到 p/q 的壓縮路徑（如 L3 R5）
//   cmp a b c d    -> 精確比較 a/b 與 c/d，輸出 -1 / 0 / 1
//
// stdin 無內容時跑內建斷言測試，輸出 all tests passed。
//
// 核心要點：
//   * Farey 相鄰項遞推：當前項 a/b、前一項 c/d，則後繼 = (k*a - c)/(k*b - d)，
//     其中 k = floor((n + d) / b)。O(|F_n|) 而非 O(n^2 log n) 排序。
//   * Farey 鄰域性質：相鄰兩項 a/b < c/d 滿足 b*c - a*d = 1。
//   * |F_n| = 1 + sum_{m=1..n} phi(m)。
//   * Stern-Brocot walk：區間 [a/b, c/d] 從哨兵 0/1 與 1/0 出發，
//     mediant = (a+c)/(b+d)；p/q 偏小走 L、偏大走 R。每步 mediant
//     分母嚴格遞增，步數上界為 p+q（約簡後）。全程整數交叉相乘，精確無誤差。

#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>
#include <utility>
#include <numeric>
#include <algorithm>
#include <iostream>

using namespace std;
typedef long long ll;
typedef pair<ll, ll> pll;

vector<pll> farey(ll n) {
    vector<pll> res;
    if (n <= 0) return res;
    // 初始兩項 0/1, 1/n；c/d 為前一項
    res.push_back({0, 1});
    res.push_back({1, n});
    ll c = 0, d = 1, a = 1, b = n;
    while (!(a == 1 && b == 1)) {
        ll k = (n + d) / b;
        ll na = k * a - c, nb = k * b - d;
        c = a; d = b;
        a = na; b = nb;
        res.push_back({a, b});
    }
    return res;
}

vector<char> sternBrocotPath(ll p, ll q) {
    // p=0 或 q=0（哨兵）以及 p/q = 1/1（根本身）回傳空路徑
    vector<char> steps;
    if (p <= 0 || q <= 0) return steps;
    ll g = gcd(p, q);
    p /= g; q /= g;
    if (p == q) return steps;
    ll a = 0, b = 1, c = 1, d = 0;  // 區間 [a/b, c/d]，c/d = 1/0 = 無窮大
    while (true) {
        ll mn = a + c, md = b + d;  // mediant
        ll lhs = p * md, rhs = q * mn;
        if (lhs == rhs) break;
        if (lhs < rhs) {  // p/q < mediant -> 左子樹
            steps.push_back('L');
            c = mn; d = md;
        } else {          // p/q > mediant -> 右子樹
            steps.push_back('R');
            a = mn; b = md;
        }
    }
    return steps;
}

string rle(const vector<char>& steps) {
    string out;
    size_t i = 0;
    while (i < steps.size()) {
        size_t j = i;
        while (j < steps.size() && steps[j] == steps[i]) ++j;
        out += steps[i];
        out += to_string(j - i);
        i = j;
        if (i < steps.size()) out += ' ';
    }
    return out;
}

int cmpFrac(ll a, ll b, ll c, ll d) {
    // 精確比較 a/b 與 c/d：交叉相乘 a*d vs b*c
    ll lhs = a * d, rhs = b * c;
    if (lhs < rhs) return -1;
    if (lhs > rhs) return 1;
    return 0;
}

ll phiBrute(ll m) {
    ll cnt = 0;
    for (ll x = 1; x <= m; ++x)
        if (gcd(x, m) == 1) ++cnt;
    return cnt;
}

// 自製 LCG（與倉庫其他題目一致的線性同餘產生器）
static ll g_state = 20261109;
ll rnd(ll m) {
    g_state = (g_state * 1103515245 + 12345) & 0x7FFFFFFF;
    return g_state % m;
}

void runTests() {
    // 1) Farey 與暴力枚舉排序一致
    for (ll n = 1; n <= 40; ++n) {
        vector<pll> brute;
        for (ll b = 1; b <= n; ++b)
            for (ll a = 0; a <= b; ++a)
                if (gcd(a, b) == 1) brute.push_back({a, b});
        sort(brute.begin(), brute.end(),
             [](const pll& x, const pll& y) {
                 return x.first * y.second < y.first * x.second;
             });
        vector<pll> f = farey(n);
        assert(f == brute);
    }

    // 2) 鄰域性質：相鄰兩項 b*c - a*d = 1
    for (ll n = 1; n <= 40; ++n) {
        vector<pll> f = farey(n);
        for (size_t i = 0; i + 1 < f.size(); ++i)
            assert(f[i].second * f[i + 1].first - f[i].first * f[i + 1].second == 1);
    }

    // 3) |F_n| = 1 + sum(phi(m))
    for (ll n = 1; n <= 60; ++n) {
        ll total = 1;
        for (ll m = 1; m <= n; ++m) total += phiBrute(m);
        assert((ll)farey(n).size() == total);
    }

    // 4) path 重放：終點 mediant 恰為 p/q；步數 <= p+q；rle 可無損還原
    for (ll p = 1; p <= 60; ++p) {
        for (ll q = 1; q <= 60; ++q) {
            ll g = gcd(p, q);
            ll pp = p / g, qq = q / g;
            vector<char> steps = sternBrocotPath(p, q);
            ll a = 0, b = 1, c = 1, d = 0;
            for (char s : steps) {
                ll mn = a + c, md = b + d;
                if (s == 'L') { c = mn; d = md; }
                else { a = mn; b = md; }
            }
            assert((a + c) * qq == (b + d) * pp);
            assert((ll)steps.size() <= pp + qq);
            // rle 解碼還原
            vector<char> dec;
            string enc = rle(steps);
            size_t i = 0;
            while (i < enc.size()) {
                char ch = enc[i++];
                ll cnt = 0;
                while (i < enc.size() && enc[i] >= '0' && enc[i] <= '9')
                    cnt = cnt * 10 + (enc[i++] - '0');
                dec.insert(dec.end(), (size_t)cnt, ch);
            }
            assert(dec == steps);
        }
    }

    // 5) 對稱性：p/q 與 q/p 的路徑 L/R 互換
    for (ll p = 1; p <= 40; ++p) {
        for (ll q = 1; q <= 40; ++q) {
            vector<char> sp = sternBrocotPath(p, q);
            vector<char> sq = sternBrocotPath(q, p);
            assert(sp.size() == sq.size());
            for (size_t i = 0; i < sp.size(); ++i)
                assert(sq[i] == (sp[i] == 'L' ? 'R' : 'L'));
        }
    }

    // 6) cmp：整數交叉相乘為真值，浮點僅確認方向（差至少 1/2500）
    for (int t = 0; t < 400; ++t) {
        ll b = rnd(51), d = rnd(51);
        ll a = (b == 0) ? rnd(50) + 1 : rnd(51);
        ll c = (d == 0) ? rnd(50) + 1 : rnd(51);
        int got = cmpFrac(a, b, c, d);
        int exp;
        if (a * d == b * c) exp = 0;
        else {
            double f1 = (b == 0) ? 1e300 : (double)a / b;
            double f2 = (d == 0) ? 1e300 : (double)c / d;
            exp = (f1 < f2) ? -1 : 1;
        }
        assert(got == exp);
    }

    // 7) 固定輸出用例（實跑結果作斷言，防回歸）
    {
        vector<pll> f5 = farey(5);
        string s;
        for (size_t i = 0; i < f5.size(); ++i) {
            if (i) s += ' ';
            s += to_string(f5[i].first) + "/" + to_string(f5[i].second);
        }
        assert(s == "0/1 1/5 1/4 1/3 2/5 1/2 3/5 2/3 3/4 4/5 1/1");
        assert(rle(sternBrocotPath(5, 3)) == "R1 L1 R1");
        assert(rle(sternBrocotPath(3, 5)) == "L1 R1 L1");
        assert(rle(sternBrocotPath(1, 1)) == "");
        assert(rle(sternBrocotPath(0, 7)) == "");
        assert(cmpFrac(2, 4, 1, 2) == 0);
        assert(cmpFrac(22, 7, 3, 1) == 1);
        assert(cmpFrac(1, 0, 999999, 1) == 1);
    }
}

int main() {
    // 判斷 stdin 是否有內容：有 -> IO 模式；無 -> 內建斷言測試
    string cmd;
    if (cin >> cmd) {
        // 已讀到第一個命令，直接進入處理循環
        bool any = true;
        while (any) {
            if (cmd == "farey") {
                ll n;
                if (!(cin >> n)) break;
                vector<pll> f = farey(n);
                for (size_t i = 0; i < f.size(); ++i) {
                    if (i) cout << ' ';
                    cout << f[i].first << '/' << f[i].second;
                }
                cout << '\n';
            } else if (cmd == "path") {
                ll p, q;
                if (!(cin >> p >> q)) break;
                cout << rle(sternBrocotPath(p, q)) << '\n';
            } else if (cmd == "cmp") {
                ll a, b, c, d;
                if (!(cin >> a >> b >> c >> d)) break;
                cout << cmpFrac(a, b, c, d) << '\n';
            } else {
                break;  // 未知命令
            }
            if (!(cin >> cmd)) break;
        }
    } else {
        runTests();
        printf("all tests passed\n");
    }
    return 0;
}
