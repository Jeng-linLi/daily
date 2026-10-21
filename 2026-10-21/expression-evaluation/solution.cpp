// 表達式求值：調度場算法（Shunting-yard）、三種表示法互轉、精確有理數運算
//
// 題意：
//     給定若干變數的值，再給 q 個詢問。每個詢問是一種「表達式表示法」的運算：
//
//         mode 1  <中綴>      求值（走調度場算法）
//         mode 2  <中綴>      轉成後綴（逆波蘭式 RPN）
//         mode 3  <中綴>      轉成前綴（波蘭式 PN）
//         mode 4  <後綴>      後綴求值
//         mode 5  <前綴>      前綴求值
//         mode 6  <中綴>      輸出「完全括號化」的中綴（表達式樹的形狀）
//
//     運算子：`+` `-` `*` `/` `^`（右結合）以及括號；另支援**一元負號**，
//     在後綴／前綴裡以 `~` 表示。運算元是**非負整數**或**單一小寫字母變數**。
//
//     所有運算都以**精確有理數**進行（不做浮點），輸出最簡分數 `p/q`（q = 1 時只印 p）。
//     任何錯誤（語法錯、除以零、未定義變數、指數非負整數以外、溢位）一律輸出 `error`。
//
// 思路：
//     ### 1. 調度場算法（Shunting-yard）
//     一邊掃描中綴 token、一邊用一個**運算子棧**決定輸出順序：
//
//         運算元  → 直接輸出
//         `(`     → 推入棧
//         `)`     → 一直彈出到遇到 `(`（沒遇到就是括號不匹配）
//         運算子 o1 → 先把棧頂「優先序更高、或同階但左結合」的 o2 彈出，再推入 o1
//
//     掃完後把棧彈空；若中途看到 `(` 就是括號沒閉合。
//
//     ### 2. 一元負號的優先序（最容易寫錯的地方）
//     數學與 Python 的慣例是 `-a^b == -(a^b)`，也就是 `^` 綁得比一元負號**緊**。
//     所以優先序定成：
//
//         `^` = 4（右結合）  >  `~` = 3（一元）  >  `*` `/` = 2  >  `+` `-` = 1
//
//     而且**一元（前綴）運算子推入時一律不彈棧** —— 它只作用於緊接在後面的運算元，
//     彈出別人反而會把 `a ^ ~ b` 算成 `(a^b) 取負`。
//
//     ### 3. 為什麼要有第二套實作
//     本題除了調度場，還另外寫了一套**遞歸下降 + 優先序爬升**（precedence climbing）
//     的剖析器。兩套演算法思路完全不同，卻必須得到同一個值 ——
//     這是抓「結合性寫反」「一元優先序錯」這類 bug 最有效的手段。
//     再加上「中綴求值 == 後綴求值 == 前綴求值 == 括號化後再求值」四路交叉驗證。
//
//     ### 4. 精確有理數，不用浮點
//     分數以 `(num, den)` 表示，強制 `den > 0`、每次運算後用 gcd 約分。
//     這樣 Python 與 C++ 的輸出才可能**逐位元組一致**（浮點做不到）。
//     為了避免 C++ 的 int64 溢位，統一訂一個上限 `LIMIT = 10^9`：
//     約分後的分子或分母一旦超過就視為 `error`（兩個語言同一條規則）。
//     指數運算限制指數必須是 `[0, 60]` 內的整數，`0^0` 定義為 1。
//
// 應用場景：
//     - 試算表 / 計算機：Excel 公式、SQL 計算欄位、金融報價的運算式引擎
//     - 編譯器前端：中綴 → AST → 三地址碼，本題的表達式樹就是 AST 的雛形
//     - 設定檔與規則引擎：促銷折扣公式、風控規則（`amount > 1000 and ...`）
//     - 符號計算：CAS（如 SymPy）底層同樣走精確有理數，避免 0.1+0.2 != 0.3
//     - 早期 Fortran / BASIC 編譯器、HP 計算機的 RPN 輸入法
//
// 複雜度（表達式長度 m，變數數固定）：
//     剖析 / 求值 / 三種表示法互轉        O(m)
//     空間                              O(m)

#include <algorithm>
#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <map>
#include <numeric>
#include <random>
#include <sstream>
#include <string>
#include <vector>

#ifdef _WIN32
#include <fcntl.h>
#include <io.h>
#endif

using namespace std;
using ll = long long;

static const ll LIMIT = 1000000000LL;
static const ll VMAX = 1000000000LL;
static const int QMAX = 200;
static const int TMAX = 60;
static const int EMAX = 60;

struct EvalError {};

static int precOf(const string& t) {
    if (t == "^") return 4;
    if (t == "~") return 3;
    if (t == "*" || t == "/") return 2;
    if (t == "+" || t == "-") return 1;
    return 0;
}

static bool isRightAssoc(const string& t) { return t == "^" || t == "~"; }
static bool isBinary(const string& t) {
    return t == "^" || t == "*" || t == "/" || t == "+" || t == "-";
}

static bool isNumberTok(const string& s) {
    if (s.empty()) return false;
    for (char c : s) if (c < '0' || c > '9') return false;
    return true;
}

static bool isVarTok(const string& s) {
    return s.size() == 1 && s[0] >= 'a' && s[0] <= 'z';
}

// ------------------------------------------------------------- 有理數運算

struct Rat {
    ll n, d;
    Rat(ll nn = 0, ll dd = 1) : n(nn), d(dd) {}
};

static Rat rMake(ll n, ll d) {
    if (d == 0) throw EvalError();
    if (d < 0) { n = -n; d = -d; }
    if (n == 0) return Rat(0, 1);
    ll g = std::gcd(n > 0 ? n : -n, d);
    n /= g;
    d /= g;
    if (n > LIMIT || n < -LIMIT || d > LIMIT) throw EvalError();
    return Rat(n, d);
}

static Rat rAdd(const Rat& a, const Rat& b) { return rMake(a.n * b.d + b.n * a.d, a.d * b.d); }
static Rat rSub(const Rat& a, const Rat& b) { return rMake(a.n * b.d - b.n * a.d, a.d * b.d); }
static Rat rMul(const Rat& a, const Rat& b) { return rMake(a.n * b.n, a.d * b.d); }
static Rat rNeg(const Rat& a) { return Rat(-a.n, a.d); }

static Rat rDiv(const Rat& a, const Rat& b) {
    if (b.n == 0) throw EvalError();
    return rMake(a.n * b.d, a.d * b.n);
}

static Rat rPow(const Rat& a, const Rat& b) {
    if (b.d != 1) throw EvalError();
    ll e = b.n;
    if (e < 0 || e > EMAX) throw EvalError();
    Rat res(1, 1), base = a;
    while (e > 0) {
        if (e & 1) res = rMul(res, base);
        e >>= 1;
        if (e) base = rMul(base, base);
    }
    return res;
}

static string rStr(const Rat& a) {
    if (a.d == 1) return to_string(a.n);
    return to_string(a.n) + "/" + to_string(a.d);
}

// ------------------------------------------------------------------ 表達式樹

struct Node {
    string tok;
    int left = -1, right = -1;
};

static vector<Node> pool;

static int newNode(const string& tok, int l = -1, int r = -1) {
    Node nd;
    nd.tok = tok;
    nd.left = l;
    nd.right = r;
    pool.push_back(nd);
    return (int)pool.size() - 1;
}

static void goPost(int x, vector<string>& out) {
    if (pool[x].left >= 0) goPost(pool[x].left, out);
    if (pool[x].right >= 0) goPost(pool[x].right, out);
    out.push_back(pool[x].tok);
}

static void goPre(int x, vector<string>& out) {
    out.push_back(pool[x].tok);
    if (pool[x].left >= 0) goPre(pool[x].left, out);
    if (pool[x].right >= 0) goPre(pool[x].right, out);
}

static void goParen(int x, vector<string>& out) {
    if (pool[x].left < 0 && pool[x].right < 0) {
        out.push_back(pool[x].tok);
    } else if (pool[x].tok == "~") {
        out.push_back("(");
        out.push_back("~");
        goParen(pool[x].left, out);
        out.push_back(")");
    } else {
        out.push_back("(");
        goParen(pool[x].left, out);
        out.push_back(pool[x].tok);
        goParen(pool[x].right, out);
        out.push_back(")");
    }
}

// ------------------------------------------------------------------ 剖析

static vector<string> toUnary(const vector<string>& toks) {
    vector<string> out;
    bool hasPrev = false;
    string prev;
    for (const string& t : toks) {
        bool unaryPos = !hasPrev || prev == "(" || precOf(prev) > 0;
        if (t == "-" && unaryPos) {
            out.push_back("~");
        } else if (t == "+" && unaryPos) {
            // 一元正號是恆等運算，直接省略
        } else {
            out.push_back(t);
        }
        prev = t;
        hasPrev = true;
    }
    return out;
}

static vector<string> shuntingYard(const vector<string>& toks) {
    vector<string> out, st;
    for (const string& t : toks) {
        if (isNumberTok(t) || isVarTok(t)) {
            out.push_back(t);
        } else if (t == "(") {
            st.push_back(t);
        } else if (t == ")") {
            bool found = false;
            while (!st.empty()) {
                string top = st.back();
                st.pop_back();
                if (top == "(") { found = true; break; }
                out.push_back(top);
            }
            if (!found) throw EvalError();
        } else if (isBinary(t)) {
            while (!st.empty() && st.back() != "(") {
                string o2 = st.back();
                bool pop;
                if (isRightAssoc(t)) pop = precOf(o2) > precOf(t);
                else pop = precOf(o2) >= precOf(t);
                if (!pop) break;
                out.push_back(st.back());
                st.pop_back();
            }
            st.push_back(t);
        } else if (t == "~") {
            // 前綴運算子：只作用於緊接在後的運算元，推入時不彈棧
            st.push_back(t);
        } else {
            throw EvalError();
        }
    }
    while (!st.empty()) {
        string top = st.back();
        st.pop_back();
        if (top == "(") throw EvalError();
        out.push_back(top);
    }
    if (out.empty()) throw EvalError();
    return out;
}

static int buildTree(const vector<string>& post) {
    vector<int> st;
    for (const string& t : post) {
        if (t == "~") {
            if (st.empty()) throw EvalError();
            int c = st.back();
            st.pop_back();
            st.push_back(newNode("~", c, -1));
        } else if (isBinary(t)) {
            if (st.size() < 2) throw EvalError();
            int r = st.back(); st.pop_back();
            int l = st.back(); st.pop_back();
            st.push_back(newNode(t, l, r));
        } else if (isNumberTok(t) || isVarTok(t)) {
            st.push_back(newNode(t));
        } else {
            throw EvalError();
        }
    }
    if (st.size() != 1) throw EvalError();
    return st[0];
}

// ------------------------------------------------------------------ 求值

using Env = map<string, Rat>;

static Rat leafValue(const string& t, const Env& env) {
    if (isNumberTok(t)) {
        if (t.size() > 18) throw EvalError();
        ll v = strtoll(t.c_str(), nullptr, 10);
        if (v > VMAX) throw EvalError();
        return Rat(v, 1);
    }
    if (isVarTok(t)) {
        auto it = env.find(t);
        if (it == env.end()) throw EvalError();
        return it->second;
    }
    throw EvalError();
}

static Rat evalNode(int x, const Env& env) {
    const Node& nd = pool[x];
    if (nd.left < 0 && nd.right < 0) return leafValue(nd.tok, env);
    if (nd.tok == "~") return rNeg(evalNode(nd.left, env));
    Rat a = evalNode(nd.left, env);
    Rat b = evalNode(nd.right, env);
    if (nd.tok == "+") return rAdd(a, b);
    if (nd.tok == "-") return rSub(a, b);
    if (nd.tok == "*") return rMul(a, b);
    if (nd.tok == "/") return rDiv(a, b);
    if (nd.tok == "^") return rPow(a, b);
    throw EvalError();
}

static Rat evalPostfix(const vector<string>& post, const Env& env) {
    vector<Rat> st;
    for (const string& t : post) {
        if (t == "~") {
            if (st.empty()) throw EvalError();
            Rat v = st.back(); st.pop_back();
            st.push_back(rNeg(v));
        } else if (isBinary(t)) {
            if (st.size() < 2) throw EvalError();
            Rat b = st.back(); st.pop_back();
            Rat a = st.back(); st.pop_back();
            if (t == "+") st.push_back(rAdd(a, b));
            else if (t == "-") st.push_back(rSub(a, b));
            else if (t == "*") st.push_back(rMul(a, b));
            else if (t == "/") st.push_back(rDiv(a, b));
            else st.push_back(rPow(a, b));
        } else {
            st.push_back(leafValue(t, env));
        }
    }
    if (st.size() != 1) throw EvalError();
    return st[0];
}

static Rat evalPrefix(const vector<string>& pre, const Env& env) {
    // 從右往左掃：遇到運算元推棧，遇到運算子取棧頂
    vector<Rat> st;
    for (int i = (int)pre.size() - 1; i >= 0; --i) {
        const string& t = pre[i];
        if (t == "~") {
            if (st.empty()) throw EvalError();
            Rat v = st.back(); st.pop_back();
            st.push_back(rNeg(v));
        } else if (isBinary(t)) {
            if (st.size() < 2) throw EvalError();
            Rat a = st.back(); st.pop_back();   // 左運算元離運算子較近，先被掃到
            Rat b = st.back(); st.pop_back();
            if (t == "+") st.push_back(rAdd(a, b));
            else if (t == "-") st.push_back(rSub(a, b));
            else if (t == "*") st.push_back(rMul(a, b));
            else if (t == "/") st.push_back(rDiv(a, b));
            else st.push_back(rPow(a, b));
        } else {
            st.push_back(leafValue(t, env));
        }
    }
    if (st.size() != 1) throw EvalError();
    return st[0];
}

// ------------------------------------------- 第二套實作：遞歸下降 + 優先序爬升

class Parser {
public:
    const vector<string>& toks;
    const Env& env;
    int i;

    Parser(const vector<string>& t, const Env& e) : toks(t), env(e), i(0) {}

    bool hasTok() const { return i < (int)toks.size(); }
    string peek() const { return hasTok() ? toks[i] : string(); }

    Rat primary() {
        if (!hasTok()) throw EvalError();
        string t = toks[i++];
        if (t == "~") return rNeg(expr(3));
        if (t == "(") {
            Rat v = expr(1);
            if (!hasTok() || toks[i] != ")") throw EvalError();
            ++i;
            return v;
        }
        return leafValue(t, env);
    }

    Rat expr(int minPrec) {
        Rat lhs = primary();
        while (true) {
            if (!hasTok()) break;
            string t = toks[i];
            if (!isBinary(t)) break;
            int p = precOf(t);
            if (p < minPrec) break;
            ++i;
            int nxt = isRightAssoc(t) ? p : p + 1;
            Rat rhs = expr(nxt);
            if (t == "+") lhs = rAdd(lhs, rhs);
            else if (t == "-") lhs = rSub(lhs, rhs);
            else if (t == "*") lhs = rMul(lhs, rhs);
            else if (t == "/") lhs = rDiv(lhs, rhs);
            else lhs = rPow(lhs, rhs);
        }
        return lhs;
    }

    Rat parseAll() {
        Rat v = expr(1);
        if (i != (int)toks.size()) throw EvalError();
        return v;
    }
};

// ------------------------------------------------------------------ 驅動

static Rat infixValue(const vector<string>& toks, const Env& env, int* rootOut) {
    vector<string> post = shuntingYard(toUnary(toks));
    int root = buildTree(post);
    if (rootOut) *rootOut = root;
    return evalNode(root, env);
}

static string run(const Env& env, const vector<vector<string>>& queries) {
    ostringstream oss;
    for (const auto& q : queries) {
        if (q.empty()) { oss << "error\n"; continue; }
        const string& modeTok = q[0];
        bool modeOk = false;
        int m = 0;
        if (!modeTok.empty()) {
            bool allDigit = true;
            size_t st = (modeTok[0] == '+' || modeTok[0] == '-') ? 1 : 0;
            if (st >= modeTok.size()) allDigit = false;
            for (size_t k = st; k < modeTok.size(); ++k)
                if (modeTok[k] < '0' || modeTok[k] > '9') allDigit = false;
            if (allDigit) { m = (int)strtol(modeTok.c_str(), nullptr, 10); modeOk = true; }
        }
        if (!modeOk) { oss << "error\n"; continue; }
        vector<string> toks(q.begin() + 1, q.end());
        if ((int)toks.size() > TMAX) toks.resize(TMAX);

        try {
            if (m == 1) {
                int root = -1;
                Rat v1 = infixValue(toks, env, &root);
                vector<string> unary1 = toUnary(toks);
                Parser p(unary1, env);
                Rat v2 = p.parseAll();
                if (v1.n != v2.n || v1.d != v2.d) {
                    fprintf(stderr, "infix mismatch\n");
                    abort();
                }
                oss << rStr(v1) << "\n";
            } else if (m == 2) {
                int root = -1;
                Rat v = infixValue(toks, env, &root);
                vector<string> post;
                goPost(root, post);
                Rat w = evalPostfix(post, env);
                if (w.n != v.n || w.d != v.d) { fprintf(stderr, "postfix mismatch\n"); abort(); }
                for (size_t k = 0; k < post.size(); ++k) oss << (k ? " " : "") << post[k];
                oss << "\n";
            } else if (m == 3) {
                int root = -1;
                Rat v = infixValue(toks, env, &root);
                vector<string> pre;
                goPre(root, pre);
                Rat w = evalPrefix(pre, env);
                if (w.n != v.n || w.d != v.d) { fprintf(stderr, "prefix mismatch\n"); abort(); }
                for (size_t k = 0; k < pre.size(); ++k) oss << (k ? " " : "") << pre[k];
                oss << "\n";
            } else if (m == 4) {
                oss << rStr(evalPostfix(toks, env)) << "\n";
            } else if (m == 5) {
                oss << rStr(evalPrefix(toks, env)) << "\n";
            } else if (m == 6) {
                int root = -1;
                Rat v = infixValue(toks, env, &root);
                vector<string> par;
                goParen(root, par);
                int r2 = -1;
                Rat w = infixValue(par, env, &r2);
                if (w.n != v.n || w.d != v.d) { fprintf(stderr, "paren mismatch\n"); abort(); }
                for (size_t k = 0; k < par.size(); ++k) oss << (k ? " " : "") << par[k];
                oss << "\n";
            } else {
                oss << "error\n";
            }
        } catch (EvalError&) {
            oss << "error\n";
        }
    }
    return oss.str();
}

static void check(const Env& env, const vector<vector<string>>& queries) {
    run(env, queries);
}

static Env makeEnv(const vector<pair<string, ll>>& ps) {
    Env e;
    for (auto& p : ps) e[p.first] = Rat(p.second, 1);
    return e;
}

static void runTests() {
    Env e = makeEnv({{"a", 6}, {"b", 3}, {"c", 0}});

    check(e, {{"1", "a", "+", "b", "*", "2"}});
    check(e, {{"2", "a", "+", "b", "*", "2"}, {"3", "a", "+", "b", "*", "2"}});
    check(e, {{"6", "(", "a", "+", "b", ")", "*", "2"}});

    check(e, {{"1", "2", "^", "3", "^", "2"}, {"2", "2", "^", "3", "^", "2"}});
    check(e, {{"1", "~", "b", "^", "2"}, {"2", "-", "b", "^", "2"}});
    check(e, {{"1", "a", "^", "~", "b"}, {"3", "a", "^", "-", "b"}});

    check(e, {{"4", "a", "b", "2", "*", "+"}, {"5", "+", "a", "*", "b", "2"}});

    check(e, {{"1", "b", "/", "a"}, {"1", "1", "/", "3", "+", "1", "/", "6"}});

    check(e, {{"1", "a", "/", "c"}, {"1", "z"}, {"1", "(", "a"}, {"1", "a", ")"},
               {"1", "a", "^", "~", "1"}, {"1", "a", "^", "1", "/", "2"},
               {"1"}, {"1", "a", "+"}});

    check(e, {{"1", "1000000000", "*", "1000000000"}});

    // 隨機對拍：隨機 token 序列，六種模式全問一遍（含大量非法輸入）
    mt19937 rng(20261021);
    auto ri = [&](int a, int b) -> int { return (int)(rng() % (unsigned)(b - a + 1)) + a; };
    vector<Env> envs;
    envs.push_back(makeEnv({{"a", 1}, {"b", 2}, {"c", 3}}));
    envs.push_back(makeEnv({{"a", -7}, {"b", 5}, {"c", 11}}));
    envs.push_back(makeEnv({{"a", 0}, {"b", 1}, {"c", 2}}));

    static const char* OPS[] = {"+", "-", "*", "/", "^"};

    for (int t = 0; t < 300; ++t) {
        const Env& env = envs[ri(0, (int)envs.size() - 1)];
        int len = ri(1, 12);
        vector<string> toks;
        for (int i = 0; i < len; ++i) {
            int r = ri(0, 9);
            if (r < 3) toks.push_back(to_string(ri(0, 9)));
            else if (r < 6) { static const char* VS[] = {"a", "b", "c"}; toks.push_back(VS[ri(0, 2)]); }
            else if (r < 8) toks.push_back(OPS[ri(0, 4)]);
            else if (r == 8) toks.push_back("(");
            else toks.push_back(")");
        }
        for (int m = 1; m <= 6; ++m) {
            vector<string> one;
            one.push_back(to_string(m));
            for (auto& s : toks) one.push_back(s);
            check(env, {one});
        }
    }

    // 隨機對拍：合法後綴 / 前綴序列
    for (int t = 0; t < 200; ++t) {
        const Env& env = envs[ri(0, (int)envs.size() - 1)];
        vector<string> post;
        int nleaf = 0;
        int k = ri(1, 5);
        for (int i = 0; i < k; ++i) {
            if (ri(0, 1) == 0) post.push_back(to_string(ri(0, 9)));
            else { static const char* VS[] = {"a", "b", "c"}; post.push_back(VS[ri(0, 2)]); }
            ++nleaf;
            if (ri(0, 9) < 6 && nleaf >= 2) { post.push_back(OPS[ri(0, 4)]); --nleaf; }
            if (ri(0, 9) < 2) post.push_back("~");
        }
        vector<string> q4;
        q4.push_back("4");
        for (auto& s : post) q4.push_back(s);
        check(env, {q4});
        vector<string> q5;
        q5.push_back("5");
        for (int i = (int)post.size() - 1; i >= 0; --i) q5.push_back(post[i]);
        check(env, {q5});
    }
}

// -- IO 解析 --

static bool isIntTokStrict(const string& s) {
    if (s.empty()) return false;
    size_t i = 0;
    if (s[0] == '+' || s[0] == '-') i = 1;
    if (i >= s.size()) return false;
    for (size_t k = i; k < s.size(); ++k)
        if (s[k] < '0' || s[k] > '9') return false;
    return true;
}

int main() {
#ifdef _WIN32
    _setmode(_fileno(stdout), _O_BINARY);
#endif
    string raw, line;
    ostringstream buf;
    while (getline(cin, line)) buf << line << "\n";
    raw = buf.str();

    bool hasInput = false;
    for (char c : raw) if (!isspace((unsigned char)c)) { hasInput = true; break; }
    if (!hasInput) {
        runTests();
        fputs("all tests passed\n", stdout);
        return 0;
    }

    vector<string> lines;
    {
        string cur;
        for (char c : raw) {
            if (c == '\n') { lines.push_back(cur); cur.clear(); }
            else cur.push_back(c);
        }
        if (!cur.empty()) lines.push_back(cur);
    }
    size_t pos = 0;
    auto nextLine = [&]() -> string {
        if (pos >= lines.size()) return string();
        return lines[pos++];
    };
    auto trim = [](const string& s) -> string {
        size_t a = 0, b = s.size();
        while (a < b && isspace((unsigned char)s[a])) ++a;
        while (b > a && isspace((unsigned char)s[b - 1])) --b;
        return s.substr(a, b - a);
    };
    auto nextInt = [&](ll defVal) -> ll {
        string t = trim(nextLine());
        if (!isIntTokStrict(t)) return defVal;
        ll v = strtoll(t.c_str(), nullptr, 10);
        if (v > VMAX) v = VMAX;
        if (v < -VMAX) v = -VMAX;
        return v;
    };

    ll ml = nextInt(0);
    if (ml < 0) ml = 0;
    if (ml > 26) ml = 26;
    Env env;
    for (ll i = 0; i < ml; ++i) {
        string s = nextLine();
        vector<string> parts;
        {
            string cur;
            for (char c : s) {
                if (isspace((unsigned char)c)) { if (!cur.empty()) { parts.push_back(cur); cur.clear(); } }
                else cur.push_back(c);
            }
            if (!cur.empty()) parts.push_back(cur);
        }
        if (parts.size() < 2) continue;
        if (!isVarTok(parts[0])) continue;
        if (!isIntTokStrict(parts[1])) continue;
        ll v = strtoll(parts[1].c_str(), nullptr, 10);
        if (v > VMAX) v = VMAX;
        if (v < -VMAX) v = -VMAX;
        env[parts[0]] = Rat(v, 1);
    }

    ll ql = nextInt(0);
    if (ql < 0) ql = 0;
    if (ql > QMAX) ql = QMAX;
    vector<vector<string>> queries;
    for (ll i = 0; i < ql; ++i) {
        string s = nextLine();
        vector<string> parts;
        string cur;
        for (char c : s) {
            if (isspace((unsigned char)c)) { if (!cur.empty()) { parts.push_back(cur); cur.clear(); } }
            else cur.push_back(c);
        }
        if (!cur.empty()) parts.push_back(cur);
        queries.push_back(parts);
    }

    string out = run(env, queries);
    fwrite(out.data(), 1, out.size(), stdout);
    return 0;
}
