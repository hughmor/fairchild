use std::path::{Path, PathBuf};

use thiserror::Error;

#[derive(Debug, Error)]
pub enum OsdiError {
    /// The library loaded, but speaks an interface this runtime cannot read.
    ///
    /// Worth naming the compiler: when fairchild drove the compile itself, the
    /// user's next question is which binary produced this, and `openvaf` 23.x
    /// emits an older OSDI than `openvaf-r` does.
    #[error(
        "OSDI version {major}.{minor} not supported (need 0.4+) — this library was built by a \
         compiler emitting an older interface; OpenVAF-Reloaded (`openvaf-r`) emits 0.4"
    )]
    Version { major: u32, minor: u32 },

    #[error("OSDI_DESCRIPTOR_SIZE={got} is smaller than OsdiDescriptor ({expected} bytes)")]
    DescriptorSizeMismatch { expected: usize, got: usize },

    #[error("dlopen failed: {0}")]
    DlOpen(String),

    #[error("symbol '{symbol}' not found: {detail}")]
    Symbol {
        symbol: &'static str,
        detail: String,
    },

    /// Any of the above, with the file it happened to. Added by
    /// [`OsdiError::with_context`] at the load site, where the path is known.
    #[error("cannot load OSDI library '{}': {detail}", path.display())]
    Load { path: PathBuf, detail: String },

    // ── driving the Verilog-A compile ─────────────────────────────────────
    /// No compiler. Names every binary looked for and both ways to point at
    /// one — never a fallback: a deck missing a device is a wrong circuit, not
    /// a degraded one.
    #[error(
        "no Verilog-A compiler found (tried: {}). Install OpenVAF-Reloaded and put it on PATH, \
         or name it with --openvaf <path> or FAIRCHILD_OPENVAF=<path>. To use only pre-compiled \
         artefacts instead, pass --no-va-compile and name them with .osdi",
        tried.join(", ")
    )]
    CompilerNotFound { tried: Vec<String> },

    #[error("Verilog-A compiler '{}' could not be run: {detail}", path.display())]
    CompilerFailed { path: PathBuf, detail: String },

    #[error("compiling '{}' failed:\n{stderr}", path.display())]
    CompileFailed { path: PathBuf, stderr: String },

    #[error("Verilog-A source '{}' does not exist", path.display())]
    VaSourceMissing { path: PathBuf },

    #[error("--va-include '{}' is not a directory", path.display())]
    IncludeDirMissing { path: PathBuf },

    /// `--no-va-compile` and a `.va` source. Refused whatever the cache holds:
    /// a flag whose effect depends on an invisible directory is not the
    /// reproducible offline route it claims to be.
    #[error(
        "'{}' is Verilog-A source and --no-va-compile is set. Compile it yourself and name the \
         result with .osdi, or drop --no-va-compile",
        path.display()
    )]
    CompileDisabled { path: PathBuf },

    #[error("Verilog-A cache '{}': {detail}", path.display())]
    CacheDir { path: PathBuf, detail: String },
}

impl OsdiError {
    /// Attach the file this failure was about.
    ///
    /// A bare "dlopen failed: …" is useless in a deck naming six libraries.
    pub fn with_context(self, path: &Path) -> Self {
        match self {
            // Already carries a path, or is about a different file entirely.
            Self::Load { .. }
            | Self::CompilerNotFound { .. }
            | Self::CompilerFailed { .. }
            | Self::CompileFailed { .. }
            | Self::VaSourceMissing { .. }
            | Self::IncludeDirMissing { .. }
            | Self::CompileDisabled { .. }
            | Self::CacheDir { .. } => self,
            other => Self::Load {
                path: path.to_path_buf(),
                detail: other.to_string(),
            },
        }
    }
}
