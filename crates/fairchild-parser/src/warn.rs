//! One switch for every user-facing warning the library emits.
//!
//! `--quiet` promised to "suppress all warning messages" and suppressed only the
//! CLI's own: the twenty-odd warnings raised inside the parser and the solver
//! printed regardless. Threading a flag through every function that can warn —
//! `collect_defs`, `SimOptions::from_netlist`, the device registry, the sanity
//! check, the Newton loop — would touch signatures that have nothing else to do
//! with output, so the switch is process-wide instead, the way a log level is.
//!
//! Frontends set it once at startup; nothing in the library reads it except
//! [`warn_user!`]. Warnings are cosmetic, which is what makes one global
//! acceptable here: a Python program running two `Circuit`s on two threads
//! shares the setting, and the cost of that is a missing line of stderr.

use std::sync::atomic::{AtomicBool, Ordering};

static QUIET: AtomicBool = AtomicBool::new(false);

/// Silence every [`warn_user!`] in the library. Errors are unaffected — a
/// warning says the run continued, and an error is still returned to the caller.
pub fn set_quiet(quiet: bool) {
    QUIET.store(quiet, Ordering::Relaxed);
}

/// Whether warnings are currently suppressed.
pub fn quiet() -> bool {
    QUIET.load(Ordering::Relaxed)
}

/// Print one warning to stderr unless [`set_quiet`] turned them off.
///
/// Takes the message *without* the `warning: ` prefix — the macro writes it, so
/// every warning in the tree reads the same way and there is one place to change
/// if that ever needs a stream other than stderr.
#[macro_export]
macro_rules! warn_user {
    ($($arg:tt)*) => {
        if !$crate::warn::quiet() {
            eprintln!("warning: {}", format_args!($($arg)*));
        }
    };
}
