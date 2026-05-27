use thiserror::Error;

#[derive(Debug, Error)]
pub enum OsdiError {
    #[error("OSDI version {major}.{minor} not supported (need 0.4+)")]
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
}
