//! Where one `key=value` ends and the next begins.
//!
//! Both dialects have their own tokeniser, because their bracket and comment
//! rules differ. Neither may have its own idea of what *opens an assignment* —
//! that is one concept, and two copies of it is how a value gets cut in half.
//!
//! # Why a token boundary is not enough
//!
//! A tokeniser that keeps bracketed runs whole still splits a value at a space
//! that sits at bracket depth zero, and a real expression has plenty:
//!
//! ```text
//! rwire=(extr==1) ? (1e-4) : ((5.3585/(w*nf) - 0.000194)*(nf==1))
//! rsw1 = r1*(ev==0)+ 0*(ev==1)+ 0*(ev==2)
//! ```
//!
//! The parens balance after `(extr==1)`, so the space before `?` is depth zero and
//! the value ends there. Everything after it becomes a stray token, and a
//! `.subckt` header reads a stray token as a port — which is how a three-port
//! foundry cell came to declare eleven (#105).
//!
//! Chasing the operators is whack-a-mole: `?`, `:`, a trailing `+`, a bare `-`.
//! The rule that holds is the other way round — **a value runs until something
//! that can only be the next assignment**, and nothing else ends it.

/// Does this token begin a new `key=value`?
///
/// True when the token carries an `=` at bracket depth zero whose left side is a
/// plausible parameter name. Depth matters: `(nf==1)` is a comparison *inside* a
/// value and must not look like the start of the next one, and it is the single
/// most common way a real card would break this rule.
///
/// A leading `=` is not an opener either — `name =value` puts the `=` at the head
/// of its token, and the name is the token before it, so the caller glues rather
/// than starting something new.
pub(crate) fn opens_assignment(tok: &str) -> bool {
    let mut depth = 0i32;
    let mut span: Option<char> = None;
    for (i, c) in tok.char_indices() {
        match span {
            Some(open) => {
                let closes = match open {
                    '{' => c == '}',
                    q => c == q,
                };
                if closes {
                    span = None;
                }
            }
            None => match c {
                '{' | '\'' | '"' => span = Some(c),
                '(' | '[' => depth += 1,
                ')' | ']' => depth -= 1,
                '=' if depth == 0 => {
                    // `a==b` is a comparison, not an assignment, and `<=`/`>=`/`!=`
                    // are operators. Only a bare `=` after a name opens one.
                    let before = &tok[..i];
                    let after = tok[i + 1..].chars().next();
                    if after == Some('=') {
                        return false;
                    }
                    if before.ends_with(['=', '<', '>', '!']) {
                        return false;
                    }
                    return is_name(before);
                }
                _ => {}
            },
        }
    }
    false
}

/// A plausible parameter name: non-empty, and made of the characters a card uses.
///
/// Deliberately permissive about what a name may contain and strict about what it
/// may *start* with, because the job here is to tell a name from an expression,
/// not to validate an identifier — that happens later, with a better error.
fn is_name(s: &str) -> bool {
    let s = s.trim();
    let mut chars = s.chars();
    match chars.next() {
        Some(c) if c.is_ascii_alphabetic() || c == '_' => {}
        _ => return false,
    }
    chars.all(|c| c.is_ascii_alphanumeric() || c == '_' || c == '.' || c == '$')
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn an_assignment_opens_and_a_comparison_does_not() {
        for yes in [
            "w=10u",
            "rwire=(a)",
            "a.b=1",
            "x_1=2",
            "p=",
            "scale={a + b}",
        ] {
            assert!(opens_assignment(yes), "{yes} opens an assignment");
        }
        for no in [
            "(nf==1)", // a comparison inside brackets
            "nf==1",   // a comparison at depth zero
            "a>=1",    // an operator
            "a<=1",
            "a!=1",
            "?", // a ternary's arms
            ":",
            "(1e-4)",
            "0*(ev==1)+", // a continued sum
            "=4",         // `name =value`: the name is the token before
            "1e-4=2",     // not a name on the left
            "",
        ] {
            assert!(!opens_assignment(no), "{no:?} does not open an assignment");
        }
    }

    /// The case that motivated this: a ternary whose condition closes its bracket.
    #[test]
    fn a_ternarys_arms_do_not_open_assignments() {
        let toks = [
            "rwire=(extr==1)",
            "?",
            "(1e-4)",
            ":",
            "((5.3/(w*nf))*(nf==1))",
        ];
        let opens: Vec<bool> = toks.iter().map(|t| opens_assignment(t)).collect();
        assert_eq!(
            opens,
            vec![true, false, false, false, false],
            "only the first token opens an assignment: {toks:?}"
        );
    }
}
