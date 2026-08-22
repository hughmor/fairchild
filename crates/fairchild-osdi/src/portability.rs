//! Verilog-A constructs the compiler cannot be trusted with on this platform.
//!
//! # The one entry today
//!
//! On macOS/aarch64, **any** formatted-output or severity task in a model takes
//! the process down during `setup_instance`. OpenVAF declares `snprintf` as a
//! variadic with three fixed parameters and then lowers the call with all of
//! them on the stack. AArch64 Apple's variadic convention passes the *fixed*
//! arguments in x0–x2 and only the tail on the stack, so `snprintf` reads
//! whatever is in x0 — the OSDI handle — as its destination buffer:
//!
//! ```text
//!   7ac: stp  xzr, xzr, [sp]     ; dest, size  → stack
//!   7b4: add  x25, x25, #0xb60   ; format string
//!   7b8: stp  x25, x25, [sp, #0x10]
//!   7bc: bl   snprintf           ; x0/x1/x2 never loaded
//! ```
//!
//! On x86-64 SysV the two conventions coincide for these arguments, which is why
//! it never shows up upstream. Confirmed by probe on openvaf 23.5.0: `$strobe`,
//! `$display`, `$write`, `$monitor`, `$debug`, `$info`, `$warning`, `$error` and
//! `$fatal` each segfault, and the same model with the call removed runs and
//! gives the right answer. It is an upstream codegen bug, not ours (see #42);
//! the entry below comes out when upstream fixes it, and until then this is the
//! difference between "BSIM4 crashes" and "BSIM4 runs".
//!
//! # Why remove rather than refuse
//!
//! Refusing to compile would be the usual answer in this codebase, and here it
//! is the wrong one: the alternative is not a wrong number, it is a `SIGSEGV`
//! and no simulation at all. A model whose diagnostics are removed still solves
//! the circuit it describes — the constitutive equations are untouched.
//!
//! What *is* lost is the model talking to us, and that includes `$error` and
//! `$fatal`: a foundry card validating its own parameters can no longer say a
//! value is out of range. So every removal is reported, per task per file,
//! saying exactly that. Losing a diagnostic loudly is the trade; losing one
//! quietly would not be.
//!
//! Note the model's *other* channel is unaffected: `OsdiInitInfo` errors are
//! read back and surfaced (`device.rs`), so a range violation reported through
//! the OSDI interface rather than through `$error` still reaches the user.
//!
//! # Adding an entry
//!
//! One [`Rule`] per (platform, construct) pair. Both the `os`/`arch` match and
//! the stripping are plain functions over explicit strings, so a rule for a
//! platform you are not on is still testable on the platform you are on — which
//! is the only way the macOS rule gets covered by CI at all.

use fairchild_core::warn_user;

/// Skip every transformation: compile the source exactly as written.
///
/// For checking whether an upstream fix has landed, and for the case where
/// somebody would rather have the crash than the silence.
pub const KEEP_ENV: &str = "FAIRCHILD_VA_KEEP_UNSUPPORTED";

/// One platform, and the system tasks its toolchain cannot compile correctly.
pub struct Rule {
    /// `std::env::consts::OS` value this applies to.
    pub os: &'static str,
    /// `std::env::consts::ARCH` value, or `None` for every architecture.
    pub arch: Option<&'static str>,
    /// System tasks to remove, written as they appear in source.
    pub tasks: &'static [&'static str],
    /// What the user loses, in a sentence that finishes "…, so <this>".
    pub cost: &'static str,
}

/// Every rule. One entry, and it should stay that way.
pub const RULES: &[Rule] = &[Rule {
    os: "macos",
    arch: Some("aarch64"),
    // Every formatted-output and severity task, because they share the one
    // broken lowering — probed individually rather than assumed.
    tasks: &[
        "$strobe",
        "$display",
        "$write",
        "$monitor",
        "$debug",
        "$info",
        "$warning",
        "$error",
        "$fatal",
        "$fstrobe",
        "$fdisplay",
        "$fwrite",
    ],
    cost: "any condition the model would have reported — including a parameter \
           it considers out of range — is now invisible",
}];

/// The rules that apply to `(os, arch)`.
pub fn rules_for<'a>(os: &'a str, arch: &'a str) -> impl Iterator<Item = &'static Rule> + 'a {
    RULES
        .iter()
        .filter(move |r| r.os == os && r.arch.is_none_or(|a| a == arch))
}

/// What [`sanitize`] took out: a task, and how many calls to it went.
#[derive(Debug, PartialEq, Eq)]
pub struct Removed {
    pub task: &'static str,
    pub calls: usize,
    pub cost: &'static str,
}

/// Whether the caller should transform at all.
///
/// Read here rather than inside [`sanitize`] so the transformation itself stays a
/// pure function of its arguments: an env-var read inside it made every other
/// test in this module depend on the one that set the variable, and they run in
/// parallel — a flaky test about a guard against silent wrongness would be a
/// poor joke.
pub fn keep_as_written() -> bool {
    std::env::var_os(KEEP_ENV).is_some()
}

/// Remove every call this platform's toolchain cannot compile.
///
/// Returns the text to compile and what was taken out of it; an empty second
/// element means the source is untouched and should be compiled as written.
/// Pure: the caller decides whether to apply it (see [`keep_as_written`]).
pub fn sanitize(text: &str, os: &str, arch: &str) -> (String, Vec<Removed>) {
    let mut out = text.to_string();
    let mut removed = Vec::new();
    for rule in rules_for(os, arch) {
        for task in rule.tasks {
            let (next, calls) = strip_calls(&out, task);
            if calls > 0 {
                removed.push(Removed {
                    task,
                    calls,
                    cost: rule.cost,
                });
            }
            out = next;
        }
    }
    (out, removed)
}

/// Report a sanitised source on stderr, once per task.
pub fn report(source: &std::path::Path, removed: &[Removed]) {
    for r in removed {
        warn_user!(
            "removed {} call(s) to {} from '{}' before compiling: the Verilog-A \
             compiler miscompiles it on {}-{} and the model would crash the run \
             instead of printing. The circuit is unchanged, but {}. Set {}=1 to \
             compile the source as written.",
            r.calls,
            r.task,
            source.display(),
            std::env::consts::OS,
            std::env::consts::ARCH,
            r.cost,
            KEEP_ENV
        );
    }
}

/// Delete every `task(...)` statement from `text`, replacing each with the null
/// statement `;`.
///
/// `;` rather than nothing, because a call is a *statement* and can be the whole
/// body of a branch: deleting the text of `if (x) $strobe("a"); else y;` outright
/// leaves `if (x) else y;`, which does not parse. A null statement is legal
/// Verilog-A exactly where the call was (checked against the compiler, not
/// assumed).
///
/// String literals, line comments and block comments are skipped, so a `$display`
/// inside a comment or a format string is left alone.
fn strip_calls(text: &str, task: &str) -> (String, usize) {
    let bytes = text.as_bytes();
    let mut out = String::with_capacity(text.len());
    let mut i = 0usize;
    let mut calls = 0usize;
    while i < bytes.len() {
        // Comments and strings are copied through verbatim; a `$display` inside
        // one is text, not a call.
        if let Some(end) = skip_span(bytes, i) {
            out.push_str(&text[i..end]);
            i = end;
            continue;
        }
        if bytes[i] == b'$' && text[i..].starts_with(task) {
            // A task name must not be a prefix of a longer identifier:
            // `$displayfoo` is not `$display`.
            let after = i + task.len();
            let boundary = bytes
                .get(after)
                .is_none_or(|c| !(c.is_ascii_alphanumeric() || *c == b'_' || *c == b'$'));
            if boundary {
                if let Some(open) = next_paren(bytes, after) {
                    if let Some(close) = matching_paren(bytes, open) {
                        // Swallow a trailing `;` if the statement has one, so a
                        // single `;` replaces the whole thing rather than two.
                        let mut j = close + 1;
                        while j < bytes.len() && (bytes[j] as char).is_whitespace() {
                            j += 1;
                        }
                        if bytes.get(j) == Some(&b';') {
                            j += 1;
                        }
                        out.push(';');
                        i = j;
                        calls += 1;
                        continue;
                    }
                }
            }
        }
        // Not a call: copy one whole char (not one byte — the source may hold
        // UTF-8 in a comment or a string).
        let ch = text[i..].chars().next().expect("in bounds");
        out.push(ch);
        i += ch.len_utf8();
    }
    (out, calls)
}

/// If `i` opens a comment or a string literal, the index just past its end.
fn skip_span(b: &[u8], i: usize) -> Option<usize> {
    match (b[i], b.get(i + 1)) {
        (b'/', Some(b'/')) => Some(
            b[i..]
                .iter()
                .position(|c| *c == b'\n')
                .map_or(b.len(), |n| i + n),
        ),
        (b'/', Some(b'*')) => {
            let mut j = i + 2;
            while j + 1 < b.len() && !(b[j] == b'*' && b[j + 1] == b'/') {
                j += 1;
            }
            Some((j + 2).min(b.len()))
        }
        (b'"', _) => {
            let mut j = i + 1;
            while j < b.len() {
                match b[j] {
                    b'\\' => j += 2,
                    b'"' => return Some(j + 1),
                    _ => j += 1,
                }
            }
            Some(b.len())
        }
        _ => None,
    }
}

/// The `(` opening a call, if only whitespace separates it from `from`.
fn next_paren(b: &[u8], from: usize) -> Option<usize> {
    let mut j = from;
    while j < b.len() && (b[j] as char).is_whitespace() {
        j += 1;
    }
    (b.get(j) == Some(&b'(')).then_some(j)
}

/// The `)` matching the `(` at `open`, counting nesting and ignoring parens
/// inside strings and comments.
fn matching_paren(b: &[u8], open: usize) -> Option<usize> {
    let mut depth = 0usize;
    let mut i = open;
    while i < b.len() {
        if let Some(end) = skip_span(b, i) {
            i = end;
            continue;
        }
        match b[i] {
            b'(' => depth += 1,
            b')' => {
                depth -= 1;
                if depth == 0 {
                    return Some(i);
                }
            }
            _ => {}
        }
        i += 1;
    }
    None
}

#[cfg(test)]
mod tests {
    use super::*;

    const MAC: (&str, &str) = ("macos", "aarch64");

    /// The platform match is by string, so the macOS rule is testable from
    /// Linux — otherwise CI would never run any of this.
    #[test]
    fn a_rule_applies_only_to_its_platform() {
        let src = "analog begin $strobe(\"hi\"); end\n";
        let (out, removed) = sanitize(src, MAC.0, MAC.1);
        assert_eq!(removed.len(), 1);
        assert_eq!(out, "analog begin ; end\n");

        for (os, arch) in [
            ("linux", "x86_64"),
            ("macos", "x86_64"),
            ("windows", "x86_64"),
        ] {
            let (out, removed) = sanitize(src, os, arch);
            assert!(removed.is_empty(), "{os}-{arch} should be untouched");
            assert_eq!(out, src);
        }
    }

    /// A call is a statement, and can be a whole branch body — so what replaces
    /// it has to be a statement too.
    #[test]
    fn a_stripped_call_leaves_a_null_statement() {
        let (out, _) = sanitize(
            "if (x > 0) $error(\"bad %g\", x); else y = 1;\n",
            MAC.0,
            MAC.1,
        );
        assert_eq!(out, "if (x > 0) ; else y = 1;\n");
    }

    /// Nested parens, commas, and a format string containing both.
    #[test]
    fn the_whole_argument_list_goes_however_nested() {
        let (out, removed) = sanitize(
            "$display(\"a) ; b(\", foo(bar(1), 2), \"c\");\nI(a) <+ V(a);\n",
            MAC.0,
            MAC.1,
        );
        assert_eq!(out, ";\nI(a) <+ V(a);\n");
        assert_eq!(removed[0].calls, 1);
    }

    /// Text that only looks like a call: inside a comment, inside a string, or
    /// part of a longer name.
    #[test]
    fn only_real_calls_are_removed() {
        let src = "// $strobe(\"in a comment\");\n\
                   /* $display(\"in a block\"); */\n\
                   x = \"$write(1)\";\n\
                   $displayfoo(1);\n\
                   $monitorbar = 2;\n";
        let (out, removed) = sanitize(src, MAC.0, MAC.1);
        assert!(removed.is_empty(), "removed {removed:?} from {src}");
        assert_eq!(out, src);
    }

    /// Several calls to several tasks, counted per task, because the warning is
    /// per task and a count is what tells the user how much they lost.
    #[test]
    fn every_call_is_counted_per_task() {
        let src = "$strobe(\"a\"); $strobe(\"b\");\n$warning(\"c\");\n";
        let (out, removed) = sanitize(src, MAC.0, MAC.1);
        assert_eq!(out, "; ;\n;\n");
        let strobe = removed.iter().find(|r| r.task == "$strobe").unwrap();
        assert_eq!(strobe.calls, 2);
        assert_eq!(
            removed.iter().find(|r| r.task == "$warning").unwrap().calls,
            1
        );
    }

    /// The escape hatch, for checking whether upstream has fixed it. Tests the
    /// switch alone — `sanitize` does not read it, so nothing else here can be
    /// disturbed by this test setting a process-wide variable.
    #[test]
    fn the_keep_env_is_read() {
        assert!(!keep_as_written(), "{KEEP_ENV} was already set");
        std::env::set_var(KEEP_ENV, "1");
        let seen = keep_as_written();
        std::env::remove_var(KEEP_ENV);
        assert!(seen);
    }
}
