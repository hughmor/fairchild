fn main() {
    // On macOS, Python extension modules must allow undefined symbols at link
    // time (they are resolved when the .so is imported by the Python runtime).
    if std::env::var("CARGO_CFG_TARGET_OS").as_deref() == Ok("macos") {
        println!("cargo:rustc-link-arg=-undefined");
        println!("cargo:rustc-link-arg=dynamic_lookup");
    }
}
