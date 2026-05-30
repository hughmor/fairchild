//! Behavioural-source expression grammar.
//!
//! Used by SPICE B-elements (`Bname n+ n- V=<expr>` or `I=<expr>`) and — in a
//! future commit — `.measure` and `.param` evaluation.  The expression is
//! parsed once into an `Expr` AST that can be evaluated repeatedly against a
//! `Context` providing node voltages, branch currents, and simulation time.
//!
//! Grammar (Pratt-style precedence climbing):
//!   expr   := ternary
//!   ternary := orExpr ['?' expr ':' expr]
//!   orExpr  := andExpr {'||' andExpr}
//!   andExpr := cmpExpr {'&&' cmpExpr}
//!   cmpExpr := sum {('<'|'>'|'<='|'>='|'=='|'!=') sum}
//!   sum    := product {('+'|'-') product}
//!   product := power {('*'|'/') power}
//!   power  := unary {'^' unary}     (right associative; handled iteratively)
//!   unary  := ['+'|'-'|'!'] atom
//!   atom   := number | name | V(name [, name]) | I(name) | TIME |
//!             name '(' [args] ')' | '(' expr ')'

use std::fmt;

#[derive(Debug, Clone, PartialEq)]
pub enum BinOp {
    Add,
    Sub,
    Mul,
    Div,
    Pow,
    Lt,
    Le,
    Gt,
    Ge,
    Eq,
    Ne,
    And,
    Or,
}

#[derive(Debug, Clone)]
pub enum Expr {
    /// Numeric literal.
    Num(f64),
    /// Single-node voltage reference, `V(node)`.
    NodeV(String),
    /// Differential voltage, `V(node1, node2)` = V(node1) − V(node2).
    NodeDiffV(String, String),
    /// Branch current through a voltage source, `I(Vname)`.
    BranchI(String),
    /// Simulation time `TIME`.
    Time,
    /// Unary `-`.
    Neg(Box<Expr>),
    /// Logical not `!`.
    Not(Box<Expr>),
    /// Binary operator.
    Bin(BinOp, Box<Expr>, Box<Expr>),
    /// `cond ? a : b`.
    If(Box<Expr>, Box<Expr>, Box<Expr>),
    /// Named function call: `sin(x)`, `pow(a, b)`, `if(c, a, b)`, etc.
    Call(String, Vec<Expr>),
    /// Bare scalar variable, resolved via [`EvalContext::variable`]. Used by
    /// device-internal constitutive maps (e.g. a `.model` expression
    /// `dneff_dV = "-3.1e-5*V"` where `V` is the device's bias) — not by
    /// B-source expressions, which reference nodes via `V(...)`.
    Var(String),
}

/// Things the evaluator needs from the surrounding circuit.
pub trait EvalContext {
    fn node_voltage(&self, node: &str) -> f64;
    fn branch_current(&self, vsrc: &str) -> f64;
    fn time(&self) -> f64;
    /// Resolve a bare scalar variable ([`Expr::Var`]). Default 0.0 — only
    /// contexts that evaluate device constitutive maps (over `V`, `T`, `lambda`,
    /// …) need to override it; B-source contexts never produce `Var`.
    fn variable(&self, _name: &str) -> f64 {
        0.0
    }
}

#[derive(Debug)]
pub enum ExprError {
    UnexpectedToken(String, usize),
    UnexpectedEof,
    BadNumber(String),
    UnknownFunction(String),
}

impl fmt::Display for ExprError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            ExprError::UnexpectedToken(t, i) => write!(f, "unexpected token '{t}' at offset {i}"),
            ExprError::UnexpectedEof => write!(f, "unexpected end of expression"),
            ExprError::BadNumber(s) => write!(f, "invalid number '{s}'"),
            ExprError::UnknownFunction(s) => write!(f, "unknown function '{s}'"),
        }
    }
}
impl std::error::Error for ExprError {}

impl Expr {
    /// Parse an expression string into an AST.
    pub fn parse(input: &str) -> Result<Expr, ExprError> {
        let tokens = tokenize(input)?;
        let mut p = Parser { tokens, pos: 0 };
        let e = p.parse_expr()?;
        if p.pos != p.tokens.len() {
            return Err(ExprError::UnexpectedToken(
                format!("{:?}", &p.tokens[p.pos]),
                p.pos,
            ));
        }
        Ok(e)
    }

    /// Evaluate this expression against `ctx`.
    pub fn eval<C: EvalContext>(&self, ctx: &C) -> f64 {
        match self {
            Expr::Num(v) => *v,
            Expr::NodeV(n) => ctx.node_voltage(n),
            Expr::NodeDiffV(a, b) => ctx.node_voltage(a) - ctx.node_voltage(b),
            Expr::BranchI(n) => ctx.branch_current(n),
            Expr::Var(n) => ctx.variable(n),
            Expr::Time => ctx.time(),
            Expr::Neg(e) => -e.eval(ctx),
            Expr::Not(e) => {
                if e.eval(ctx) != 0.0 {
                    0.0
                } else {
                    1.0
                }
            }
            Expr::Bin(op, a, b) => {
                let av = a.eval(ctx);
                let bv = b.eval(ctx);
                match op {
                    BinOp::Add => av + bv,
                    BinOp::Sub => av - bv,
                    BinOp::Mul => av * bv,
                    BinOp::Div => av / bv,
                    BinOp::Pow => av.powf(bv),
                    BinOp::Lt => (av < bv) as i64 as f64,
                    BinOp::Le => (av <= bv) as i64 as f64,
                    BinOp::Gt => (av > bv) as i64 as f64,
                    BinOp::Ge => (av >= bv) as i64 as f64,
                    BinOp::Eq => (av == bv) as i64 as f64,
                    BinOp::Ne => (av != bv) as i64 as f64,
                    BinOp::And => {
                        if av != 0.0 && bv != 0.0 {
                            1.0
                        } else {
                            0.0
                        }
                    }
                    BinOp::Or => {
                        if av != 0.0 || bv != 0.0 {
                            1.0
                        } else {
                            0.0
                        }
                    }
                }
            }
            Expr::If(c, a, b) => {
                if c.eval(ctx) != 0.0 {
                    a.eval(ctx)
                } else {
                    b.eval(ctx)
                }
            }
            Expr::Call(name, args) => {
                let vs: Vec<f64> = args.iter().map(|a| a.eval(ctx)).collect();
                eval_fn(name, &vs).unwrap_or(0.0)
            }
        }
    }

    /// Walk the AST collecting every `V(node)` / `V(n1,n2)` / `I(vsrc)`
    /// reference (in any nesting depth) into the provided vectors.
    pub fn collect_refs(&self, v_nodes: &mut Vec<String>, i_srcs: &mut Vec<String>) {
        match self {
            Expr::NodeV(n) => v_nodes.push(n.clone()),
            Expr::NodeDiffV(a, b) => {
                v_nodes.push(a.clone());
                v_nodes.push(b.clone());
            }
            Expr::BranchI(s) => i_srcs.push(s.clone()),
            Expr::Num(_) | Expr::Time | Expr::Var(_) => {}
            Expr::Neg(e) | Expr::Not(e) => e.collect_refs(v_nodes, i_srcs),
            Expr::Bin(_, a, b) => {
                a.collect_refs(v_nodes, i_srcs);
                b.collect_refs(v_nodes, i_srcs);
            }
            Expr::If(c, a, b) => {
                c.collect_refs(v_nodes, i_srcs);
                a.collect_refs(v_nodes, i_srcs);
                b.collect_refs(v_nodes, i_srcs);
            }
            Expr::Call(_, args) => {
                for a in args {
                    a.collect_refs(v_nodes, i_srcs);
                }
            }
        }
    }
}

fn eval_fn(name: &str, args: &[f64]) -> Option<f64> {
    let one = |f: fn(f64) -> f64| args.first().copied().map(f);
    match name.to_lowercase().as_str() {
        "sin" => one(f64::sin),
        "cos" => one(f64::cos),
        "tan" => one(f64::tan),
        "asin" => one(f64::asin),
        "acos" => one(f64::acos),
        "atan" => one(f64::atan),
        "sinh" => one(f64::sinh),
        "cosh" => one(f64::cosh),
        "tanh" => one(f64::tanh),
        "exp" => one(f64::exp),
        "log" | "ln" => one(f64::ln),
        "log10" => one(f64::log10),
        "sqrt" => one(f64::sqrt),
        "abs" => one(f64::abs),
        "sgn" => one(|x| {
            if x > 0.0 {
                1.0
            } else if x < 0.0 {
                -1.0
            } else {
                0.0
            }
        }),
        "ceil" => one(f64::ceil),
        "floor" => one(f64::floor),
        "min" => args.first().and_then(|a| args.get(1).map(|b| a.min(*b))),
        "max" => args.first().and_then(|a| args.get(1).map(|b| a.max(*b))),
        "pow" => args.first().and_then(|a| args.get(1).map(|b| a.powf(*b))),
        "atan2" => args.first().and_then(|a| args.get(1).map(|b| a.atan2(*b))),
        "if" => {
            // if(cond, then, else)
            if args.len() >= 3 {
                Some(if args[0] != 0.0 { args[1] } else { args[2] })
            } else {
                None
            }
        }
        _ => None,
    }
}

// ─── tokenizer ───────────────────────────────────────────────────────────────

#[derive(Debug, Clone, PartialEq)]
enum Tok {
    Num(f64),
    Ident(String),
    LParen,
    RParen,
    Comma,
    Plus,
    Minus,
    Star,
    Slash,
    Caret,
    Question,
    Colon,
    Lt,
    Le,
    Gt,
    Ge,
    Eq,
    Ne,
    And,
    Or,
    Bang,
}

fn tokenize(input: &str) -> Result<Vec<Tok>, ExprError> {
    let bytes = input.as_bytes();
    let mut i = 0usize;
    let mut out = Vec::new();
    while i < bytes.len() {
        let c = bytes[i] as char;
        if c.is_whitespace() {
            i += 1;
            continue;
        }
        // Two-char operators first.
        if i + 1 < bytes.len() {
            let s = &input[i..i + 2];
            match s {
                "<=" => {
                    out.push(Tok::Le);
                    i += 2;
                    continue;
                }
                ">=" => {
                    out.push(Tok::Ge);
                    i += 2;
                    continue;
                }
                "==" => {
                    out.push(Tok::Eq);
                    i += 2;
                    continue;
                }
                "!=" => {
                    out.push(Tok::Ne);
                    i += 2;
                    continue;
                }
                "&&" => {
                    out.push(Tok::And);
                    i += 2;
                    continue;
                }
                "||" => {
                    out.push(Tok::Or);
                    i += 2;
                    continue;
                }
                _ => {}
            }
        }
        match c {
            '(' => {
                out.push(Tok::LParen);
                i += 1;
            }
            ')' => {
                out.push(Tok::RParen);
                i += 1;
            }
            ',' => {
                out.push(Tok::Comma);
                i += 1;
            }
            '+' => {
                out.push(Tok::Plus);
                i += 1;
            }
            '-' => {
                out.push(Tok::Minus);
                i += 1;
            }
            '*' => {
                out.push(Tok::Star);
                i += 1;
            }
            '/' => {
                out.push(Tok::Slash);
                i += 1;
            }
            '^' => {
                out.push(Tok::Caret);
                i += 1;
            }
            '?' => {
                out.push(Tok::Question);
                i += 1;
            }
            ':' => {
                out.push(Tok::Colon);
                i += 1;
            }
            '<' => {
                out.push(Tok::Lt);
                i += 1;
            }
            '>' => {
                out.push(Tok::Gt);
                i += 1;
            }
            '=' => {
                out.push(Tok::Eq);
                i += 1;
            } // tolerate single `=` as ==
            '!' => {
                out.push(Tok::Bang);
                i += 1;
            }
            c if c.is_ascii_digit() || c == '.' => {
                let start = i;
                while i < bytes.len() {
                    let b = bytes[i] as char;
                    if b.is_ascii_digit() || b == '.' {
                        i += 1;
                    } else if b == 'e' || b == 'E' {
                        i += 1;
                        if i < bytes.len() && (bytes[i] as char == '+' || bytes[i] as char == '-') {
                            i += 1;
                        }
                    } else {
                        break;
                    }
                }
                // SPICE-suffix support (k, meg, m, u, n, p, f, g, t):
                let num_end = i;
                let suffix_start = i;
                while i < bytes.len() && (bytes[i] as char).is_ascii_alphabetic() {
                    i += 1;
                }
                let num_str = &input[start..num_end];
                let suffix = &input[suffix_start..i].to_lowercase();
                let mult = match suffix.as_str() {
                    "meg" => 1e6,
                    "k" => 1e3,
                    "m" => 1e-3,
                    "u" => 1e-6,
                    "n" => 1e-9,
                    "p" => 1e-12,
                    "f" => 1e-15,
                    "g" => 1e9,
                    "t" => 1e12,
                    "" => 1.0,
                    _ => {
                        // The suffix isn't recognised; back up so it lexes as
                        // an identifier instead.
                        i = num_end;
                        1.0
                    }
                };
                let n: f64 = num_str
                    .parse()
                    .map_err(|_| ExprError::BadNumber(num_str.into()))?;
                out.push(Tok::Num(n * mult));
            }
            c if c.is_ascii_alphabetic() || c == '_' => {
                let start = i;
                while i < bytes.len() {
                    let b = bytes[i] as char;
                    if b.is_ascii_alphanumeric() || b == '_' {
                        i += 1;
                    } else {
                        break;
                    }
                }
                out.push(Tok::Ident(input[start..i].to_string()));
            }
            other => return Err(ExprError::UnexpectedToken(other.to_string(), i)),
        }
    }
    Ok(out)
}

// ─── Pratt parser ────────────────────────────────────────────────────────────

struct Parser {
    tokens: Vec<Tok>,
    pos: usize,
}

impl Parser {
    fn peek(&self) -> Option<&Tok> {
        self.tokens.get(self.pos)
    }
    fn bump(&mut self) -> Option<Tok> {
        let t = self.tokens.get(self.pos).cloned();
        if t.is_some() {
            self.pos += 1;
        }
        t
    }
    fn expect(&mut self, want: &Tok) -> Result<(), ExprError> {
        match self.bump() {
            Some(t) if &t == want => Ok(()),
            Some(t) => Err(ExprError::UnexpectedToken(format!("{t:?}"), self.pos)),
            None => Err(ExprError::UnexpectedEof),
        }
    }

    fn parse_expr(&mut self) -> Result<Expr, ExprError> {
        let cond = self.parse_or()?;
        if matches!(self.peek(), Some(Tok::Question)) {
            self.bump();
            let a = self.parse_expr()?;
            self.expect(&Tok::Colon)?;
            let b = self.parse_expr()?;
            Ok(Expr::If(Box::new(cond), Box::new(a), Box::new(b)))
        } else {
            Ok(cond)
        }
    }

    fn parse_or(&mut self) -> Result<Expr, ExprError> {
        let mut e = self.parse_and()?;
        while matches!(self.peek(), Some(Tok::Or)) {
            self.bump();
            let r = self.parse_and()?;
            e = Expr::Bin(BinOp::Or, Box::new(e), Box::new(r));
        }
        Ok(e)
    }

    fn parse_and(&mut self) -> Result<Expr, ExprError> {
        let mut e = self.parse_cmp()?;
        while matches!(self.peek(), Some(Tok::And)) {
            self.bump();
            let r = self.parse_cmp()?;
            e = Expr::Bin(BinOp::And, Box::new(e), Box::new(r));
        }
        Ok(e)
    }

    fn parse_cmp(&mut self) -> Result<Expr, ExprError> {
        let mut e = self.parse_sum()?;
        while let Some(op) = match self.peek() {
            Some(Tok::Lt) => Some(BinOp::Lt),
            Some(Tok::Le) => Some(BinOp::Le),
            Some(Tok::Gt) => Some(BinOp::Gt),
            Some(Tok::Ge) => Some(BinOp::Ge),
            Some(Tok::Eq) => Some(BinOp::Eq),
            Some(Tok::Ne) => Some(BinOp::Ne),
            _ => None,
        } {
            self.bump();
            let r = self.parse_sum()?;
            e = Expr::Bin(op, Box::new(e), Box::new(r));
        }
        Ok(e)
    }

    fn parse_sum(&mut self) -> Result<Expr, ExprError> {
        let mut e = self.parse_product()?;
        loop {
            let op = match self.peek() {
                Some(Tok::Plus) => BinOp::Add,
                Some(Tok::Minus) => BinOp::Sub,
                _ => break,
            };
            self.bump();
            let r = self.parse_product()?;
            e = Expr::Bin(op, Box::new(e), Box::new(r));
        }
        Ok(e)
    }

    fn parse_product(&mut self) -> Result<Expr, ExprError> {
        let mut e = self.parse_power()?;
        loop {
            let op = match self.peek() {
                Some(Tok::Star) => BinOp::Mul,
                Some(Tok::Slash) => BinOp::Div,
                _ => break,
            };
            self.bump();
            let r = self.parse_power()?;
            e = Expr::Bin(op, Box::new(e), Box::new(r));
        }
        Ok(e)
    }

    fn parse_power(&mut self) -> Result<Expr, ExprError> {
        let base = self.parse_unary()?;
        if matches!(self.peek(), Some(Tok::Caret)) {
            self.bump();
            // right-associative
            let exp = self.parse_power()?;
            Ok(Expr::Bin(BinOp::Pow, Box::new(base), Box::new(exp)))
        } else {
            Ok(base)
        }
    }

    fn parse_unary(&mut self) -> Result<Expr, ExprError> {
        match self.peek() {
            Some(Tok::Plus) => {
                self.bump();
                self.parse_unary()
            }
            Some(Tok::Minus) => {
                self.bump();
                let e = self.parse_unary()?;
                Ok(Expr::Neg(Box::new(e)))
            }
            Some(Tok::Bang) => {
                self.bump();
                let e = self.parse_unary()?;
                Ok(Expr::Not(Box::new(e)))
            }
            _ => self.parse_atom(),
        }
    }

    fn parse_atom(&mut self) -> Result<Expr, ExprError> {
        match self.bump() {
            Some(Tok::Num(n)) => Ok(Expr::Num(n)),
            Some(Tok::LParen) => {
                let e = self.parse_expr()?;
                self.expect(&Tok::RParen)?;
                Ok(e)
            }
            Some(Tok::Ident(name)) => {
                let name_lc = name.to_lowercase();
                if matches!(self.peek(), Some(Tok::LParen)) {
                    self.bump();
                    // Special: V(...) and I(...)
                    if name_lc == "v" {
                        // V(node) or V(n1, n2)
                        let n1 = self.parse_ident_or_str()?;
                        if matches!(self.peek(), Some(Tok::Comma)) {
                            self.bump();
                            let n2 = self.parse_ident_or_str()?;
                            self.expect(&Tok::RParen)?;
                            return Ok(Expr::NodeDiffV(n1, n2));
                        }
                        self.expect(&Tok::RParen)?;
                        return Ok(Expr::NodeV(n1));
                    }
                    if name_lc == "i" {
                        let n = self.parse_ident_or_str()?;
                        self.expect(&Tok::RParen)?;
                        return Ok(Expr::BranchI(n));
                    }
                    // Generic function call.
                    let mut args = Vec::new();
                    if !matches!(self.peek(), Some(Tok::RParen)) {
                        loop {
                            args.push(self.parse_expr()?);
                            if matches!(self.peek(), Some(Tok::Comma)) {
                                self.bump();
                            } else {
                                break;
                            }
                        }
                    }
                    self.expect(&Tok::RParen)?;
                    Ok(Expr::Call(name_lc, args))
                } else if name_lc == "time" {
                    Ok(Expr::Time)
                } else {
                    // Bare identifier → a scalar variable resolved by the
                    // EvalContext (e.g. `V`, `T`, `lambda` in a device
                    // constitutive map). B-source contexts leave
                    // `EvalContext::variable` at its 0.0 default, so a stray
                    // bare ident there reads as zero rather than failing — node
                    // references must still use the explicit `V(...)` form.
                    Ok(Expr::Var(name_lc))
                }
            }
            Some(t) => Err(ExprError::UnexpectedToken(format!("{t:?}"), self.pos)),
            None => Err(ExprError::UnexpectedEof),
        }
    }

    fn parse_ident_or_str(&mut self) -> Result<String, ExprError> {
        match self.bump() {
            Some(Tok::Ident(n)) => Ok(n),
            Some(t) => Err(ExprError::UnexpectedToken(format!("{t:?}"), self.pos)),
            None => Err(ExprError::UnexpectedEof),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    struct StubCtx {
        v: std::collections::HashMap<&'static str, f64>,
        i: std::collections::HashMap<&'static str, f64>,
        t: f64,
    }
    impl EvalContext for StubCtx {
        fn node_voltage(&self, node: &str) -> f64 {
            *self.v.get(node).unwrap_or(&0.0)
        }
        fn branch_current(&self, vsrc: &str) -> f64 {
            *self.i.get(vsrc).unwrap_or(&0.0)
        }
        fn time(&self) -> f64 {
            self.t
        }
    }

    fn ctx() -> StubCtx {
        let mut v = std::collections::HashMap::new();
        v.insert("n1", 2.0);
        v.insert("n2", 0.5);
        let mut i = std::collections::HashMap::new();
        i.insert("v1", 1e-3);
        StubCtx { v, i, t: 1e-6 }
    }

    #[test]
    fn arithmetic_and_precedence() {
        let e = Expr::parse("2 + 3 * 4").unwrap();
        assert_eq!(e.eval(&ctx()), 14.0);
        let e = Expr::parse("(2 + 3) * 4").unwrap();
        assert_eq!(e.eval(&ctx()), 20.0);
        let e = Expr::parse("2 ^ 3 ^ 2").unwrap(); // right-assoc
        assert_eq!(e.eval(&ctx()), 512.0);
    }

    /// Bare scalar variables resolve via `EvalContext::variable` — the path a
    /// device constitutive map (`dneff_dV = "-3.1e-5*V - 1.2e-5*V*V"`) uses.
    #[test]
    fn bare_variables_resolve_via_context() {
        struct VarCtx {
            vars: std::collections::HashMap<&'static str, f64>,
        }
        impl EvalContext for VarCtx {
            fn node_voltage(&self, _: &str) -> f64 {
                0.0
            }
            fn branch_current(&self, _: &str) -> f64 {
                0.0
            }
            fn time(&self) -> f64 {
                0.0
            }
            fn variable(&self, name: &str) -> f64 {
                *self.vars.get(name).unwrap_or(&0.0)
            }
        }
        let mut vars = std::collections::HashMap::new();
        vars.insert("v", 2.0);
        let ctx = VarCtx { vars };
        // SPICE suffix (0.1m) + bare variable V, both in one constitutive map.
        let e = Expr::parse("-3.1e-5*V - 1.2e-5*V*V").unwrap();
        let expected = -3.1e-5 * 2.0 - 1.2e-5 * 2.0 * 2.0;
        assert!((e.eval(&ctx) - expected).abs() < 1e-18, "got {}", e.eval(&ctx));
        // Unknown variable → 0 via the default-backed map.
        let e2 = Expr::parse("foo + 1").unwrap();
        assert_eq!(e2.eval(&ctx), 1.0);
    }

    #[test]
    fn node_refs() {
        let e = Expr::parse("V(n1) * 2 + V(n2)").unwrap();
        assert_eq!(e.eval(&ctx()), 2.0 * 2.0 + 0.5);
        let e = Expr::parse("V(n1, n2)").unwrap();
        assert_eq!(e.eval(&ctx()), 1.5);
    }

    #[test]
    fn branch_current_and_time() {
        let e = Expr::parse("I(v1) * 1k + TIME").unwrap();
        let val = e.eval(&ctx());
        assert!((val - (1e-3 * 1e3 + 1e-6)).abs() < 1e-12);
    }

    #[test]
    fn functions() {
        let e = Expr::parse("sqrt(V(n1)*V(n1)+V(n2)*V(n2))").unwrap();
        let expected = (4.0 + 0.25_f64).sqrt();
        assert!((e.eval(&ctx()) - expected).abs() < 1e-12);
    }

    #[test]
    fn comparison_and_ternary() {
        let e = Expr::parse("V(n1) > 1 ? 5 : 10").unwrap();
        assert_eq!(e.eval(&ctx()), 5.0);
        let e = Expr::parse("V(n2) > 1 ? 5 : 10").unwrap();
        assert_eq!(e.eval(&ctx()), 10.0);
    }

    #[test]
    fn if_function() {
        let e = Expr::parse("if(V(n1) > 1, 5, 10)").unwrap();
        assert_eq!(e.eval(&ctx()), 5.0);
    }

    #[test]
    fn collect_refs_works() {
        let e = Expr::parse("V(out) + I(vsrc1) * V(in, gnd)").unwrap();
        let mut vs = Vec::new();
        let mut is = Vec::new();
        e.collect_refs(&mut vs, &mut is);
        assert_eq!(
            vs,
            vec!["out".to_string(), "in".to_string(), "gnd".to_string()]
        );
        assert_eq!(is, vec!["vsrc1".to_string()]);
    }
}
