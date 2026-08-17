//! The `fairchild` binary. Everything is in the library so that the Python
//! wheel can expose the same command without a second copy of it.

fn main() {
    std::process::exit(fairchild_cli::run(std::env::args_os()));
}
