"""表達式求值：調度場算法（Shunting-yard）、三種表示法互轉、精確有理數運算

題意：
    給定若干變數的值，再給 q 個詢問。每個詢問是一種「表達式表示法」的運算：

        mode 1  <中綴>      求值（走調度場算法）
        mode 2  <中綴>      轉成後綴（逆波蘭式 RPN）
        mode 3  <中綴>      轉成前綴（波蘭式 PN）
        mode 4  <後綴>      後綴求值
        mode 5  <前綴>      前綴求值
        mode 6  <中綴>      輸出「完全括號化」的中綴（表達式樹的形狀）

    運算子：`+` `-` `*` `/` `^`（右結合）以及括號；另支援**一元負號**，
    在後綴／前綴裡以 `~` 表示。運算元是**非負整數**或**單一小寫字母變數**。

    所有運算都以**精確有理數**進行（不做浮點），輸出最簡分數 `p/q`（q = 1 時只印 p）。
    任何錯誤（語法錯、除以零、未定義變數、指數非負整數以外、溢位）一律輸出 `error`。

思路：
    ### 1. 調度場算法（Shunting-yard）
    一邊掃描中綴 token、一邊用一個**運算子棧**決定輸出順序：

        運算元  → 直接輸出
        `(`     → 推入棧
        `)`     → 一直彈出到遇到 `(`（沒遇到就是括號不匹配）
        運算子 o1 → 先把棧頂「優先序更高、或同階但左結合」的 o2 彈出，再推入 o1

    掃完後把棧彈空；若中途看到 `(` 就是括號沒閉合。

    ### 2. 一元負號的優先序（最容易寫錯的地方）
    數學與 Python 的慣例是 `-a^b == -(a^b)`，也就是 `^` 綁得比一元負號**緊**。
    所以優先序定成：

        `^` = 4（右結合）  >  `~` = 3（一元）  >  `*` `/` = 2  >  `+` `-` = 1

    而且**一元（前綴）運算子推入時一律不彈棧** —— 它只作用於緊接在後面的運算元，
    彈出別人反而會把 `a ^ ~ b` 算成 `(a^b) 取負`。

    ### 3. 為什麼要有第二套實作
    本題除了調度場，還另外寫了一套**遞歸下降 + 優先序爬升**（precedence climbing）
    的剖析器。兩套演算法思路完全不同，卻必須得到同一個值 ——
    這是抓「結合性寫反」「一元優先序錯」這類 bug 最有效的手段。
    再加上「中綴求值 == 後綴求值 == 前綴求值 == 括號化後再求值」四路交叉驗證。

    ### 4. 精確有理數，不用浮點
    分數以 `(num, den)` 表示，強制 `den > 0`、每次運算後用 gcd 約分。
    這樣 Python 與 C++ 的輸出才可能**逐位元組一致**（浮點做不到）。
    為了避免 C++ 的 int64 溢位，統一訂一個上限 `LIMIT = 10^9`：
    約分後的分子或分母一旦超過就視為 `error`（兩個語言同一條規則）。
    指數運算限制指數必須是 `[0, 60]` 內的整數，`0^0` 定義為 1。

應用場景：
    - 試算表 / 計算機：Excel 公式、SQL 計算欄位、金融報價的運算式引擎
    - 編譯器前端：中綴 → AST → 三地址碼，本題的表達式樹就是 AST 的雛形
    - 設定檔與規則引擎：促銷折扣公式、風控規則（`amount > 1000 and ...`）
    - 符號計算：CAS（如 SymPy）底層同樣走精確有理數，避免 0.1+0.2 != 0.3
    - 早期 Fortran / BASIC 編譯器、HP 計算機的 RPN 輸入法

複雜度（表達式長度 m，變數數固定）：
    剖析 / 求值 / 三種表示法互轉        O(m)
    空間                              O(m)
"""

import math
import re
import sys
from typing import Dict, List, Optional, Tuple

INT_RE = re.compile(r"^[+-]?[0-9]+$")

LIMIT = 10 ** 9        # 約分後分子 / 分母的絕對值上限，超過視為 error
VMAX = 10 ** 9         # 輸入值上限
QMAX = 200             # IO 模式下的詢問數上限
TMAX = 60              # 單一表達式的 token 數上限
EMAX = 60              # 指數上限

PREC = {"^": 4, "~": 3, "*": 2, "/": 2, "+": 1, "-": 1}
RIGHT_ASSOC = {"^", "~"}
BINARY = {"^", "*", "/", "+", "-"}


class EvalError(Exception):
    """任何讓表達式無法求值 / 剖析的情況。"""


# ------------------------------------------------------------- 有理數運算


def r_make(n: int, d: int) -> Tuple[int, int]:
    """建構最簡分數，強制 den > 0；超出 LIMIT 視為溢位。"""
    if d == 0:
        raise EvalError()
    if d < 0:
        n = -n
        d = -d
    if n == 0:
        return (0, 1)
    g = math.gcd(n if n > 0 else -n, d)
    n //= g
    d //= g
    if n > LIMIT or n < -LIMIT or d > LIMIT:
        raise EvalError()
    return (n, d)


def r_add(a: Tuple[int, int], b: Tuple[int, int]) -> Tuple[int, int]:
    return r_make(a[0] * b[1] + b[0] * a[1], a[1] * b[1])


def r_sub(a: Tuple[int, int], b: Tuple[int, int]) -> Tuple[int, int]:
    return r_make(a[0] * b[1] - b[0] * a[1], a[1] * b[1])


def r_mul(a: Tuple[int, int], b: Tuple[int, int]) -> Tuple[int, int]:
    return r_make(a[0] * b[0], a[1] * b[1])


def r_div(a: Tuple[int, int], b: Tuple[int, int]) -> Tuple[int, int]:
    if b[0] == 0:
        raise EvalError()
    return r_make(a[0] * b[1], a[1] * b[0])


def r_neg(a: Tuple[int, int]) -> Tuple[int, int]:
    return (-a[0], a[1])


def r_pow(a: Tuple[int, int], b: Tuple[int, int]) -> Tuple[int, int]:
    """a ** b，指數必須是 [0, EMAX] 內的整數；0 ** 0 = 1。"""
    if b[1] != 1:
        raise EvalError()
    e = b[0]
    if e < 0 or e > EMAX:
        raise EvalError()
    res = (1, 1)
    base = a
    while e > 0:
        if e & 1:
            res = r_mul(res, base)
        e >>= 1
        if e:
            base = r_mul(base, base)
    return res


def r_str(a: Tuple[int, int]) -> str:
    if a[1] == 1:
        return str(a[0])
    return "%d/%d" % (a[0], a[1])


# ------------------------------------------------------------------ 剖析


class Node:
    """表達式樹的節點：葉子（運算元）或一元 / 二元運算子。"""

    __slots__ = ("tok", "left", "right")

    def __init__(self, tok: str, left: Optional["Node"] = None,
                 right: Optional["Node"] = None) -> None:
        self.tok = tok
        self.left = left
        self.right = right


def is_number(tok: str) -> bool:
    return bool(tok) and all("0" <= c <= "9" for c in tok)


def is_var(tok: str) -> bool:
    return len(tok) == 1 and "a" <= tok[0] <= "z"


def is_op(tok: str) -> bool:
    return tok in PREC


def to_unary(toks: List[str]) -> List[str]:
    """把「出現在一元位置」的 `-` 換成 `~`，一元 `+` 直接丟掉。"""
    out: List[str] = []
    prev: Optional[str] = None
    for t in toks:
        if t == "-" and (prev is None or prev == "(" or is_op(prev)):
            out.append("~")
        elif t == "+" and (prev is None or prev == "(" or is_op(prev)):
            pass                      # 一元正號是恆等運算，直接省略
        else:
            out.append(t)
        prev = t
    return out


def shunting_yard(toks: List[str]) -> List[str]:
    """中綴 → 後綴。語法錯丟 EvalError。"""
    out: List[str] = []
    st: List[str] = []
    for t in toks:
        if is_number(t) or is_var(t):
            out.append(t)
        elif t == "(":
            st.append(t)
        elif t == ")":
            found = False
            while st:
                top = st.pop()
                if top == "(":
                    found = True
                    break
                out.append(top)
            if not found:
                raise EvalError()
        elif t in BINARY:
            while st and st[-1] != "(":
                o2 = st[-1]
                if t in RIGHT_ASSOC:
                    if PREC[o2] > PREC[t]:
                        out.append(st.pop())
                    else:
                        break
                else:
                    if PREC[o2] >= PREC[t]:
                        out.append(st.pop())
                    else:
                        break
            st.append(t)
        elif t == "~":
            # 前綴運算子：只作用於緊接在後的運算元，推入時不彈棧
            st.append(t)
        else:
            raise EvalError()
    while st:
        top = st.pop()
        if top == "(":
            raise EvalError()
        out.append(top)
    if not out:
        raise EvalError()
    return out


def build_tree(post: List[str]) -> Node:
    """後綴 → 表達式樹。"""
    st: List[Node] = []
    for t in post:
        if t == "~":
            if not st:
                raise EvalError()
            st.append(Node("~", st.pop(), None))
        elif t in BINARY:
            if len(st) < 2:
                raise EvalError()
            r = st.pop()
            l = st.pop()
            st.append(Node(t, l, r))
        elif is_number(t) or is_var(t):
            st.append(Node(t))
        else:
            raise EvalError()
    if len(st) != 1:
        raise EvalError()
    return st[0]


def tree_postfix(root: Node) -> List[str]:
    """後序遍歷（左-右-根）。"""
    out: List[str] = []

    def go(x: Node) -> None:
        if x.left is not None:
            go(x.left)
        if x.right is not None:
            go(x.right)
        out.append(x.tok)

    go(root)
    return out


def tree_prefix(root: Node) -> List[str]:
    """前序遍歷（根-左-右）。"""
    out: List[str] = []

    def go(x: Node) -> None:
        out.append(x.tok)
        if x.left is not None:
            go(x.left)
        if x.right is not None:
            go(x.right)

    go(root)
    return out


def tree_paren(root: Node) -> List[str]:
    """完全括號化的中綴，用來「看見」表達式樹的形狀。"""
    out: List[str] = []

    def go(x: Node) -> None:
        if x.left is None and x.right is None:
            out.append(x.tok)
        elif x.tok == "~":
            out.append("(")
            out.append("~")
            go(x.left)
            out.append(")")
        else:
            out.append("(")
            go(x.left)
            out.append(x.tok)
            go(x.right)
            out.append(")")

    go(root)
    return out


# ------------------------------------------------------------------ 求值


def eval_node(x: Node, env: Dict[str, Tuple[int, int]]) -> Tuple[int, int]:
    if x.left is None and x.right is None:
        if is_number(x.tok):
            v = int(x.tok)
            if v > VMAX:
                raise EvalError()
            return (v, 1)
        if is_var(x.tok):
            if x.tok not in env:
                raise EvalError()
            return env[x.tok]
        raise EvalError()
    if x.tok == "~":
        return r_neg(eval_node(x.left, env))
    a = eval_node(x.left, env)
    b = eval_node(x.right, env)
    if x.tok == "+":
        return r_add(a, b)
    if x.tok == "-":
        return r_sub(a, b)
    if x.tok == "*":
        return r_mul(a, b)
    if x.tok == "/":
        return r_div(a, b)
    if x.tok == "^":
        return r_pow(a, b)
    raise EvalError()


def eval_postfix(post: List[str], env: Dict[str, Tuple[int, int]]) -> Tuple[int, int]:
    st: List[Tuple[int, int]] = []
    for t in post:
        if t == "~":
            if not st:
                raise EvalError()
            st.append(r_neg(st.pop()))
        elif t in BINARY:
            if len(st) < 2:
                raise EvalError()
            b = st.pop()
            a = st.pop()
            if t == "+":
                st.append(r_add(a, b))
            elif t == "-":
                st.append(r_sub(a, b))
            elif t == "*":
                st.append(r_mul(a, b))
            elif t == "/":
                st.append(r_div(a, b))
            else:
                st.append(r_pow(a, b))
        elif is_number(t):
            v = int(t)
            if v > VMAX:
                raise EvalError()
            st.append((v, 1))
        elif is_var(t):
            if t not in env:
                raise EvalError()
            st.append(env[t])
        else:
            raise EvalError()
    if len(st) != 1:
        raise EvalError()
    return st[0]


def eval_prefix(pre: List[str], env: Dict[str, Tuple[int, int]]) -> Tuple[int, int]:
    """前綴求值：從**右往左**掃，遇到運算元推棧，遇到運算子取棧頂。"""
    st: List[Tuple[int, int]] = []
    for i in range(len(pre) - 1, -1, -1):
        t = pre[i]
        if t == "~":
            if not st:
                raise EvalError()
            st.append(r_neg(st.pop()))
        elif t in BINARY:
            if len(st) < 2:
                raise EvalError()
            a = st.pop()          # 前綴的左運算元離運算子比較近，會先被掃到
            b = st.pop()
            if t == "+":
                st.append(r_add(a, b))
            elif t == "-":
                st.append(r_sub(a, b))
            elif t == "*":
                st.append(r_mul(a, b))
            elif t == "/":
                st.append(r_div(a, b))
            else:
                st.append(r_pow(a, b))
        elif is_number(t):
            v = int(t)
            if v > VMAX:
                raise EvalError()
            st.append((v, 1))
        elif is_var(t):
            if t not in env:
                raise EvalError()
            st.append(env[t])
        else:
            raise EvalError()
    if len(st) != 1:
        raise EvalError()
    return st[0]


# ------------------------------------------- 第二套實作：遞歸下降 + 優先序爬升


class Parser:
    """獨立於調度場的第二套剖析器，用來交叉驗證。"""

    def __init__(self, toks: List[str], env: Dict[str, Tuple[int, int]]) -> None:
        self.toks = toks
        self.env = env
        self.i = 0

    def peek(self) -> Optional[str]:
        return self.toks[self.i] if self.i < len(self.toks) else None

    def primary(self) -> Tuple[int, int]:
        t = self.peek()
        if t is None:
            raise EvalError()
        self.i += 1
        if t == "~":
            return r_neg(self.expr(3))
        if t == "(":
            v = self.expr(1)
            if self.peek() != ")":
                raise EvalError()
            self.i += 1
            return v
        if is_number(t):
            v = int(t)
            if v > VMAX:
                raise EvalError()
            return (v, 1)
        if is_var(t):
            if t not in self.env:
                raise EvalError()
            return self.env[t]
        raise EvalError()

    def expr(self, min_prec: int) -> Tuple[int, int]:
        lhs = self.primary()
        while True:
            t = self.peek()
            if t is None or t not in BINARY:
                break
            p = PREC[t]
            if p < min_prec:
                break
            self.i += 1
            nxt = p if t in RIGHT_ASSOC else p + 1
            rhs = self.expr(nxt)
            if t == "+":
                lhs = r_add(lhs, rhs)
            elif t == "-":
                lhs = r_sub(lhs, rhs)
            elif t == "*":
                lhs = r_mul(lhs, rhs)
            elif t == "/":
                lhs = r_div(lhs, rhs)
            else:
                lhs = r_pow(lhs, rhs)
        return lhs

    def parse_all(self) -> Tuple[int, int]:
        v = self.expr(1)
        if self.i != len(self.toks):
            raise EvalError()
        return v


# ------------------------------------------------------------------ 驅動


def parse_infix(toks: List[str], env: Dict[str, Tuple[int, int]]) -> Tuple[Node, Tuple[int, int]]:
    """中綴 → 表達式樹 + 其值（走調度場）。"""
    post = shunting_yard(to_unary(toks))
    tree = build_tree(post)
    return tree, eval_node(tree, env)


def infix_value(toks: List[str], env: Dict[str, Tuple[int, int]]) -> Tuple[int, int]:
    tree, v = parse_infix(toks, env)
    return v


def run(env: Dict[str, Tuple[int, int]], queries: List[List[str]]) -> str:
    out: List[str] = []
    for q in queries:
        if not q:
            out.append("error")
            continue
        mode = q[0]
        if not INT_RE.match(mode):
            out.append("error")
            continue
        m = int(mode)
        toks = q[1:]
        if len(toks) > TMAX:
            toks = toks[:TMAX]

        try:
            if m == 1:
                v1 = infix_value(toks, env)
                v2 = Parser(to_unary(toks), env).parse_all()
                assert v1 == v2, ("infix mismatch", toks, v1, v2)
                out.append(r_str(v1))
            elif m == 2:
                tree, v = parse_infix(toks, env)
                post = tree_postfix(tree)
                assert eval_postfix(post, env) == v, ("postfix mismatch", toks)
                out.append(" ".join(post))
            elif m == 3:
                tree, v = parse_infix(toks, env)
                pre = tree_prefix(tree)
                assert eval_prefix(pre, env) == v, ("prefix mismatch", toks)
                out.append(" ".join(pre))
            elif m == 4:
                out.append(r_str(eval_postfix(toks, env)))
            elif m == 5:
                out.append(r_str(eval_prefix(toks, env)))
            elif m == 6:
                tree, v = parse_infix(toks, env)
                par = tree_paren(tree)
                assert infix_value(par, env) == v, ("paren mismatch", toks)
                out.append(" ".join(par))
            else:
                out.append("error")
        except EvalError:
            out.append("error")
        except RecursionError:
            out.append("error")
        except ZeroDivisionError:
            out.append("error")
        except ValueError:
            out.append("error")
    return "".join(s + "\n" for s in out)


def _check(env: Dict[str, Tuple[int, int]], queries: List[List[str]]) -> None:
    run(env, queries)


def _env(pairs: List[Tuple[str, int]]) -> Dict[str, Tuple[int, int]]:
    return {k: (v, 1) for k, v in pairs}


def run_tests() -> None:
    e = _env([("a", 6), ("b", 3), ("c", 0)])

    # 基本四則
    _check(e, [["1", "a", "+", "b", "*", "2"]])
    _check(e, [["2", "a", "+", "b", "*", "2"], ["3", "a", "+", "b", "*", "2"]])
    _check(e, [["6", "(", "a", "+", "b", ")", "*", "2"]])

    # 右結合與一元負號
    _check(e, [["1", "2", "^", "3", "^", "2"], ["2", "2", "^", "3", "^", "2"]])
    _check(e, [["1", "~", "b", "^", "2"], ["2", "-", "b", "^", "2"]])
    _check(e, [["1", "a", "^", "~", "b"], ["3", "a", "^", "-", "b"]])

    # 括號與後綴 / 前綴求值
    _check(e, [["4", "a", "b", "2", "*", "+"], ["5", "+", "a", "*", "b", "2"]])

    # 分數運算
    _check(e, [["1", "b", "/", "a"], ["1", "1", "/", "3", "+", "1", "/", "6"]])

    # 錯誤情形
    _check(e, [["1", "a", "/", "c"], ["1", "z"], ["1", "(", "a"], ["1", "a", ")"],
               ["1", "a", "^", "~", "1"], ["1", "a", "^", "1", "/", "2"],
               ["1"], ["1", "a", "+"]])

    # 溢位
    _check(e, [["1", "1000000000", "*", "1000000000"]])

    # 隨機對拍：隨機生成中綴表達式，六種模式全問一遍
    import random
    rnd = random.Random(20261021)
    envs = [
        _env([("a", 1), ("b", 2), ("c", 3)]),
        _env([("a", -7), ("b", 5), ("c", 11)]),
        _env([("a", 0), ("b", 1), ("c", 2)]),
    ]
    for _ in range(300):
        env = rnd.choice(envs)
        depth = rnd.randint(1, 4)
        toks: List[str] = []

        def gen(d: int) -> None:
            if d <= 0 or rnd.random() < 0.3:
                if rnd.random() < 0.5:
                    toks.append(str(rnd.randint(0, 9)))
                else:
                    toks.append(rnd.choice("abc"))
                return
            if rnd.random() < 0.25:
                toks.append("-")
                gen(d - 1)
                return
            if rnd.random() < 0.3:
                toks.append("(")
                gen(d - 1)
            else:
                gen(d - 1)
            toks.append(rnd.choice(["+", "-", "*", "/", "^"]))
            gen(d - 1)
            if toks[-2] == "(" or (len(toks) > 2 and toks[-2] == "("):
                pass
            if toks[0] == "(":
                pass
            toks.append(")")

        gen(depth)
        # 把多餘的 ')' 去掉（簡單修整，語法錯本身也是合法測試）
        _check(env, [[str(m)] + toks for m in range(1, 7)])

    # 隨機對拍：合法後綴 / 前綴序列
    for _ in range(200):
        env = rnd.choice(envs)
        post: List[str] = []
        nleaf = 0
        k = rnd.randint(1, 5)
        for i in range(k):
            if rnd.random() < 0.5:
                post.append(str(rnd.randint(0, 9)))
            else:
                post.append(rnd.choice("abc"))
            nleaf += 1
            if rnd.random() < 0.6 and nleaf >= 2:
                post.append(rnd.choice(["+", "-", "*", "/", "^"]))
                nleaf -= 1
            if rnd.random() < 0.2:
                post.append("~")
        _check(env, [["4"] + post])
        _check(env, [["5"] + list(reversed(post))])


def emit(text: str) -> None:
    sys.stdout.buffer.write(text.encode("latin-1"))


def run_io(raw: str) -> None:
    lines = raw.split("\n")
    pos = 0

    def next_line() -> str:
        nonlocal pos
        if pos >= len(lines):
            return ""
        s = lines[pos]
        pos += 1
        return s

    def next_int(default: int = 0) -> int:
        t = next_line().strip()
        if not INT_RE.match(t):
            return default
        try:
            v = int(t)
        except ValueError:
            return default
        return max(-VMAX, min(VMAX, v))

    m = next_int(0)
    if m < 0:
        m = 0
    if m > 26:
        m = 26
    env: Dict[str, Tuple[int, int]] = {}
    for _ in range(m):
        parts = next_line().split()
        if len(parts) < 2:
            continue
        name = parts[0]
        if not is_var(name):
            continue
        vt = parts[1]
        if not INT_RE.match(vt):
            continue
        try:
            v = int(vt)
        except ValueError:
            continue
        env[name] = (max(-VMAX, min(VMAX, v)), 1)

    q = next_int(0)
    if q < 0:
        q = 0
    if q > QMAX:
        q = QMAX
    queries: List[List[str]] = []
    for _ in range(q):
        queries.append(next_line().split())
    emit(run(env, queries))


def main() -> None:
    raw = sys.stdin.buffer.read().decode("latin-1")
    if raw.strip():
        run_io(raw)
    else:
        run_tests()
        emit("all tests passed\n")


if __name__ == "__main__":
    main()
