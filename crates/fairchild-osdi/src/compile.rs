//! Driving the Verilog-A compile ourselves: `.va` source in, loaded device out.
//!
//! A user hands us a deck plus Verilog-A sources. Everything between that and a
//! stampable device — invoking OpenVAF, choosing an output path, deciding
//! whether an existing artefact is still good — lives here, so the CLI and the
//! Python binding cannot disagree about it.
//!
//! # Why a subprocess and not a linked-in compiler
//!
//! OpenVAF-Reloaded is GPL-3.0 and this tree is Apache-2.0, so linking it would
//! force the licence. That is the reason it *cannot* be linked; it is not the
//! reason it *should not* be. Upstream publishes no library crate — only the
//! driver binary — so linking means depending on the internals of a compiler
//! workspace that never offered an API, plus LLVM 18+ as a hard build
//! dependency for everyone who builds fairchild, wheels included. And the
//! `.osdi` file is the interop artefact: ngspice and VACASK read the same one,
//! so producing it keeps a user's compiled PDK portable. An in-process JIT
//! would make it fairchild-only, and upstream has no JIT anyway.
//!
//! # Staleness
//!
//! A stale `.osdi` is a silently wrong device, which is the worst outcome this
//! codebase admits. So the cache key is not the source file's mtime, and not a
//! hash of the source alone: it is a hash of the compiler's own
//! `--print-expansion` output — the fully preprocessed source, every
//! `` `include `` resolved and every macro expanded, as the compiler sees it —
//! together with the compiler's version. Nothing this crate can misread about
//! Verilog-A include semantics enters the key, because this crate does not read
//! Verilog-A. The whole include closure is covered by construction.
//!
//! The cost is one preprocessor run per source per process, even on a cache
//! hit. That buys "cannot be stale", and it is the cheap half of a compile.

use std::collections::{BTreeMap, BTreeSet};
use std::hash::{Hash, Hasher};
use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::Arc;

use fairchild_core::device_registry::DeviceRegistry;

use crate::error::OsdiError;
use crate::loader::OsdiLibrary;

/// Binaries tried, in order, when no compiler is named explicitly.
///
/// `openvaf-r` first: the reloaded driver is the one that emits OSDI 0.4, which
/// is the only version this runtime loads.
const CANDIDATES: [&str; 2] = ["openvaf-r", "openvaf"];

/// How to turn `.va` into `.osdi`, as the frontend was told to.
#[derive(Debug, Clone, Default)]
pub struct VaOptions {
    /// `--openvaf <path>`: use exactly this binary, no PATH search.
    pub compiler: Option<PathBuf>,
    /// `--va-include <dir>`, in order: passed through as `-I`, searched after
    /// the source's own directory. A PDK depends on the order.
    pub include_dirs: Vec<PathBuf>,
    /// Where artefacts land. `None` picks the per-user cache.
    pub cache_dir: Option<PathBuf>,
    /// `--no-va-compile`: refuse a `.va` source rather than compile it.
    pub no_compile: bool,
    /// `--emit-generated <dir>`: where to write the expanded per-N source of a
    /// bundle-dialect model. `None` keeps it beside the artefact cache.
    pub generated_dir: Option<PathBuf>,
}

impl VaOptions {
    /// Options from the environment, for a frontend with no flags of its own.
    ///
    /// `FAIRCHILD_OPENVAF` names the compiler, `FAIRCHILD_VA_CACHE` the cache
    /// directory. Both exist so the Python binding — which takes no command
    /// line — can still be pointed at a toolchain.
    pub fn from_env() -> Self {
        Self {
            compiler: std::env::var_os("FAIRCHILD_OPENVAF").map(PathBuf::from),
            include_dirs: Vec::new(),
            cache_dir: std::env::var_os("FAIRCHILD_VA_CACHE").map(PathBuf::from),
            generated_dir: std::env::var_os("FAIRCHILD_VA_GENERATED").map(PathBuf::from),
            no_compile: false,
        }
    }

    /// Merge environment defaults under explicit flags: a flag wins, and the
    /// environment fills what the flag left empty.
    pub fn or_env(mut self) -> Self {
        let env = Self::from_env();
        self.compiler = self.compiler.or(env.compiler);
        self.cache_dir = self.cache_dir.or(env.cache_dir);
        self
    }

    /// Where generated per-N sources are written. Beside the artefact cache by
    /// default; `--emit-generated` points it somewhere the author can read.
    pub fn generated_dir(&self) -> PathBuf {
        self.generated_dir
            .clone()
            .unwrap_or_else(|| self.cache_dir().join("generated"))
    }

    fn cache_dir(&self) -> PathBuf {
        if let Some(dir) = &self.cache_dir {
            return dir.clone();
        }
        // No `dirs` crate for one path. The temp-dir fallback is last so a
        // container with no HOME still gets a cache rather than an error.
        if let Some(xdg) = std::env::var_os("XDG_CACHE_HOME") {
            return PathBuf::from(xdg).join("fairchild").join("va");
        }
        if let Some(home) = std::env::var_os("HOME") {
            return PathBuf::from(home)
                .join(".cache")
                .join("fairchild")
                .join("va");
        }
        std::env::temp_dir().join("fairchild-va")
    }
}

/// A located compiler, and the version string that goes into every cache key.
#[derive(Debug, Clone)]
pub struct VaCompiler {
    pub path: PathBuf,
    pub version: String,
}

impl VaCompiler {
    /// Find the compiler, or say precisely what was looked for and how to fix it.
    ///
    /// Never falls back to "skip the model": a deck missing a device is a wrong
    /// circuit, not a degraded one.
    pub fn find(opts: &VaOptions) -> Result<Self, OsdiError> {
        if let Some(path) = &opts.compiler {
            // Named explicitly, so a miss is that name's fault and nothing
            // else's — and *why* it missed is the whole message. A binary that
            // is there but will not run (not executable, not a binary, being
            // rewritten) reported "not found", which names the wrong problem
            // and sends the reader looking for a file that is already there.
            return match version_of(path) {
                Ok(version) => Ok(Self {
                    path: path.clone(),
                    version,
                }),
                Err(VersionError::NotFound) => Err(OsdiError::CompilerNotFound {
                    tried: vec![path.display().to_string()],
                }),
                Err(VersionError::WouldNotRun(detail)) => Err(OsdiError::CompilerFailed {
                    path: path.clone(),
                    detail,
                }),
            };
        }
        for name in CANDIDATES {
            let path = PathBuf::from(name);
            match version_of(&path) {
                Ok(version) => return Ok(Self { path, version }),
                // Searching a list of names: one that is not installed is the
                // normal case and the search goes on. One that *is* installed
                // and will not run is not — stopping there beats reporting that
                // none of four names exists when the first one does.
                Err(VersionError::NotFound) => continue,
                Err(VersionError::WouldNotRun(detail)) => {
                    return Err(OsdiError::CompilerFailed { path, detail })
                }
            }
        }
        Err(OsdiError::CompilerNotFound {
            tried: CANDIDATES.iter().map(|s| s.to_string()).collect(),
        })
    }

    fn command(&self, opts: &VaOptions, src: &Path) -> Command {
        let mut cmd = Command::new(&self.path);
        for dir in &opts.include_dirs {
            cmd.arg("-I").arg(dir);
        }
        // The source's own directory, last, so an explicit `-I` wins. OpenVAF
        // resolves a quoted `include` against the source anyway; naming it
        // keeps an angle-bracketed one working too.
        if let Some(parent) = src.parent().filter(|p| !p.as_os_str().is_empty()) {
            cmd.arg("-I").arg(parent);
        }
        cmd.arg(src);
        cmd
    }
}

/// Why asking a candidate for its version did not produce one.
enum VersionError {
    /// Nothing of that name to run — the ordinary miss while searching.
    NotFound,
    /// It is there and it did not run, or ran and failed. The detail is the
    /// operating system's own words, which are the only ones that name the
    /// actual obstacle: a directory, a missing execute bit, a script with no
    /// interpreter, a file another process is still writing.
    WouldNotRun(String),
}

/// Ask the compiler its version.
fn version_of(path: &Path) -> Result<String, VersionError> {
    let out = match Command::new(path).arg("--version").output() {
        Ok(out) => out,
        Err(e) if e.kind() == std::io::ErrorKind::NotFound => return Err(VersionError::NotFound),
        Err(e) => return Err(VersionError::WouldNotRun(e.to_string())),
    };
    if !out.status.success() {
        let stderr = String::from_utf8_lossy(&out.stderr).trim().to_string();
        let detail = if stderr.is_empty() {
            format!("`--version` exited with {}", out.status)
        } else {
            format!("`--version` exited with {}: {stderr}", out.status)
        };
        return Err(VersionError::WouldNotRun(detail));
    }
    Ok(String::from_utf8_lossy(&out.stdout).trim().to_string())
}

/// Compile `src` to an `.osdi` in the cache, or return the cached artefact.
///
/// There is deliberately no in-process memo on top of the on-disk cache. One
/// was tried, keyed by source path so that a corner sweep would not re-run the
/// preprocessor once per corner — and it served a stale artefact the moment a
/// long-lived process (a notebook holding the Python binding) edited a `.va`
/// between runs. The saving was one preprocessor run; the cost was the exact
/// failure this whole design exists to prevent. If per-corner expansion ever
/// shows up in a profile, the fix is to compile once before the corner loop,
/// not to key a cache on something weaker than content.
pub fn compile(compiler: &VaCompiler, src: &Path, opts: &VaOptions) -> Result<PathBuf, OsdiError> {
    // Before anything: `--no-va-compile` refuses a `.va` source outright, cache
    // or no cache. Serving one from the cache would make the flag's behaviour
    // depend on a directory the user cannot see, which is the opposite of the
    // reproducible offline route it exists to provide.
    if opts.no_compile {
        return Err(OsdiError::CompileDisabled {
            path: src.to_path_buf(),
        });
    }
    if !src.is_file() {
        return Err(OsdiError::VaSourceMissing {
            path: src.to_path_buf(),
        });
    }
    // OpenVAF rejects a `-I` that is not a directory before it reads a line of
    // source, and its message names only the path. Ours names the flag too.
    for dir in &opts.include_dirs {
        if !dir.is_dir() {
            return Err(OsdiError::IncludeDirMissing {
                path: dir.to_path_buf(),
            });
        }
    }

    let key = cache_key(compiler, src, opts)?;
    let stem = src
        .file_stem()
        .map(|s| s.to_string_lossy().into_owned())
        .unwrap_or_else(|| "model".into());
    let dir = opts.cache_dir();
    let out = dir.join(format!("{stem}-{key:016x}.osdi"));
    if out.is_file() {
        return Ok(out);
    }
    std::fs::create_dir_all(&dir).map_err(|e| OsdiError::CacheDir {
        path: dir.clone(),
        detail: e.to_string(),
    })?;

    // Compile to a private name and rename into place. Two fairchild runs
    // sharing a cache would otherwise be able to hand each other a half-written
    // library, which dlopen reports as a corrupt file at best.
    // Unique per *call*, not just per process. Two threads compiling the same
    // source — two devices from one deck, two tests in one binary — otherwise
    // pick the same temp name, both run the compiler into it, and the second
    // rename fails with ENOENT because the first already moved the file away.
    // The final path is shared on purpose: same key, same bytes, last writer
    // wins.
    static SEQ: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(0);
    let seq = SEQ.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
    let tmp = dir.join(format!(
        ".{stem}-{key:016x}.{}.{seq}.tmp",
        std::process::id()
    ));
    // OpenVAF's `--output` parser requires the parent to exist but not the
    // file; it also refuses to overwrite anything that is not a plain file.
    let status = compiler
        .command(opts, src)
        .arg("-o")
        .arg(&tmp)
        .output()
        .map_err(|e| OsdiError::CompilerFailed {
            path: compiler.path.clone(),
            detail: e.to_string(),
        })?;
    if !status.status.success() {
        let _ = std::fs::remove_file(&tmp);
        // The compiler's own diagnostics, verbatim: they carry file, line and
        // column into the Verilog-A source, which nothing here can improve on.
        return Err(OsdiError::CompileFailed {
            path: src.to_path_buf(),
            stderr: format!(
                "{}{}",
                String::from_utf8_lossy(&status.stderr),
                String::from_utf8_lossy(&status.stdout)
            )
            .trim()
            .to_string(),
        });
    }
    std::fs::rename(&tmp, &out).map_err(|e| OsdiError::CacheDir {
        path: out.clone(),
        detail: format!("cannot install compiled model: {e}"),
    })?;
    Ok(out)
}

/// Hash of everything that can change what the compile produces.
fn cache_key(compiler: &VaCompiler, src: &Path, opts: &VaOptions) -> Result<u64, OsdiError> {
    let expanded = compiler
        .command(opts, src)
        .arg("--print-expansion")
        .output()
        .map_err(|e| OsdiError::CompilerFailed {
            path: compiler.path.clone(),
            detail: e.to_string(),
        })?;
    if !expanded.status.success() {
        return Err(OsdiError::CompileFailed {
            path: src.to_path_buf(),
            stderr: String::from_utf8_lossy(&expanded.stderr).trim().to_string(),
        });
    }

    // SipHash, via the std default. Not stable across Rust releases — which
    // costs a recompile after a toolchain bump and can never serve a stale
    // artefact, the only direction that matters here.
    let mut h = std::collections::hash_map::DefaultHasher::new();
    compiler.version.hash(&mut h);
    // `-I` order is part of the answer: two directories holding the same
    // filename resolve differently, and a PDK does that on purpose.
    for dir in &opts.include_dirs {
        dir.hash(&mut h);
    }
    expanded.stdout.hash(&mut h);
    Ok(h.finish())
}

/// Everything a deck asked to be loaded: `.va` sources compiled first, then
/// `.osdi` artefacts, each registered into `registry`.
///
/// Returns the artefact paths in load order, which is what a frontend echoes
/// back to the user. Relative paths resolve against `base_dir` — the deck's own
/// directory — except where the parser already made them absolute while
/// splicing an include.
///
/// The caller still owns `.model` cards: they alias descriptors these libraries
/// define, so `DeviceRegistry::register_loaded_model_cards` has to run after
/// this, not inside it.
pub fn load_libraries(
    osdi_paths: &[String],
    va_sources: &[String],
    base_dir: Option<&Path>,
    opts: &VaOptions,
    registry: &mut DeviceRegistry,
) -> Result<Vec<PathBuf>, OsdiError> {
    load_libraries_with_widths(
        osdi_paths,
        va_sources,
        base_dir,
        opts,
        registry,
        &BTreeMap::new(),
        3,
    )
}

/// [`load_libraries`], told which terminal counts the deck instantiates each
/// model at.
///
/// A source written against the bundle-port dialect has no fixed port count — it
/// is generated for the width the deck declared — so loading one needs to know
/// the widths in advance. `wanted` maps a model name to the flattened terminal
/// counts seen on its `X` lines, which the first of the two load passes has
/// already worked out. Ordinary Verilog-A ignores all of this.
pub fn load_libraries_with_widths(
    osdi_paths: &[String],
    va_sources: &[String],
    base_dir: Option<&Path>,
    opts: &VaOptions,
    registry: &mut DeviceRegistry,
    wanted: &BTreeMap<String, BTreeSet<usize>>,
    wpc: usize,
) -> Result<Vec<PathBuf>, OsdiError> {
    let mut loaded = Vec::with_capacity(osdi_paths.len() + va_sources.len());

    // Told not to compile, and asked to: say so about the deck's own first
    // source. Better than "no compiler found" for someone who wanted neither.
    if opts.no_compile {
        if let Some(first) = va_sources.first() {
            return Err(OsdiError::CompileDisabled {
                path: resolve(first, base_dir),
            });
        }
    }
    // One lookup for the whole deck, and only when something needs it: a deck
    // of plain `.osdi` lines must not require a compiler on PATH.
    let compiler = if va_sources.is_empty() {
        None
    } else {
        Some(VaCompiler::find(opts)?)
    };

    for src in va_sources {
        let src = resolve(src, base_dir);
        let compiler = compiler
            .as_ref()
            .expect("sources present, compiling allowed");

        // Bundle-port dialect? Then this source is a template, and one artefact
        // per requested channel count comes out of it.
        let text = std::fs::read_to_string(&src).map_err(|e| OsdiError::CompileFailed {
            path: src.clone(),
            stderr: format!("cannot read: {e}"),
        })?;
        if let Some(m) = crate::dialect::scan(&text)? {
            let widths: Vec<usize> = wanted
                .get(&m.name)
                .into_iter()
                .flatten()
                .filter_map(|&flat| m.channels_for(flat, wpc))
                .collect();
            if widths.is_empty() {
                // Nothing in the deck instantiates it. Generating a default N
                // would be a guess, and a guessed width is a wrong device, so
                // register the shape and let the arity oracle refuse anything
                // that does not fit.
                registry.declare_arity(
                    m.name.clone(),
                    fairchild_core::ArityDecl::Bundle {
                        scalars: m.scalars.len(),
                        per_channel: m.bundles.len() * wpc,
                    },
                );
                continue;
            }
            // The generated source lives in the cache, not beside the author's
            // file, so its relative `include`s would no longer resolve — the
            // compiler searches the source's own directory, and that is now the
            // cache. Put the original's directory on the include path so the
            // author's `include "optical.vams"` means what they wrote.
            let mut gen_opts = opts.clone();
            if let Some(dir) = src.parent() {
                gen_opts.include_dirs.push(dir.to_path_buf());
            }
            for n in widths {
                let expanded = crate::dialect::expand(&text, &m, n, wpc)?;
                let gen_path = generated_path(&src, n, opts)?;
                std::fs::write(&gen_path, &expanded).map_err(|e| OsdiError::CompileFailed {
                    path: gen_path.clone(),
                    stderr: format!("cannot write generated source: {e}"),
                })?;
                let compiled = compile(compiler, &gen_path, &gen_opts)?;
                load_one(&compiled, registry)?;
                loaded.push(compiled);
            }
            registry.declare_arity(
                m.name.clone(),
                fairchild_core::ArityDecl::Bundle {
                    scalars: m.scalars.len(),
                    per_channel: m.bundles.len() * wpc,
                },
            );
            continue;
        }

        let compiled = compile(compiler, &src, opts)?;
        load_one(&compiled, registry)?;
        loaded.push(compiled);
    }
    for path in osdi_paths {
        let path = resolve(path, base_dir);
        load_one(&path, registry)?;
        loaded.push(path);
    }
    Ok(loaded)
}

/// Where a generated per-N source lands: beside the cache, named for its width,
/// and kept rather than written to a temporary. An author debugging a model that
/// works at N=1 and misbehaves at N=8 needs to be able to read what was actually
/// compiled.
fn generated_path(src: &Path, n: usize, opts: &VaOptions) -> Result<PathBuf, OsdiError> {
    let dir = opts.generated_dir();
    std::fs::create_dir_all(&dir).map_err(|e| OsdiError::CacheDir {
        path: dir.clone(),
        detail: e.to_string(),
    })?;
    let stem = src.file_stem().map(|s| s.to_string_lossy().to_string());
    Ok(dir.join(format!(
        "{}.n{n}.va",
        stem.unwrap_or_else(|| "model".into())
    )))
}

fn resolve(path: &str, base_dir: Option<&Path>) -> PathBuf {
    let p = Path::new(path);
    match base_dir {
        Some(dir) if p.is_relative() => dir.join(p),
        _ => p.to_path_buf(),
    }
}

fn load_one(path: &Path, registry: &mut DeviceRegistry) -> Result<(), OsdiError> {
    // SAFETY: the path came from the deck; `OsdiLibrary::open` validates the
    // OSDI version and descriptor layout before anything is called through.
    let lib = unsafe { OsdiLibrary::open(path) }.map_err(|e| e.with_context(path))?;
    Arc::new(lib).register_into(registry);
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn write(path: &Path, text: &str) {
        std::fs::create_dir_all(path.parent().unwrap()).unwrap();
        std::fs::write(path, text).unwrap();
    }

    /// Write an executable stub, and wait until the OS will actually run it.
    ///
    /// Writing an executable and exec'ing it from a multithreaded process is a
    /// race on Linux: exec returns `ETXTBSY` while any process holds a write
    /// descriptor to that inode, and a sibling test that forks while this thread
    /// is still writing hands its child exactly such a descriptor — briefly,
    /// even though Rust opens files close-on-exec. It cost a CI job: the same
    /// commit passed ubuntu in one run and failed it in the other, inside
    /// `VaCompiler::find`, on a stub that had just been chmod'ed. The old error
    /// path then reported it as "compiler not found", which is why that is fixed
    /// too.
    ///
    /// One successful exec closes the question — afterwards no descriptor to
    /// this file is open anywhere, so every later run of it is safe. So the wait
    /// lives here, once, rather than as a retry inside `version_of`: production
    /// code that silently re-runs a compiler is a worse thing to own than a test
    /// that waits for its own fixture.
    fn write_stub(path: &Path, script: &str) {
        write(path, script);
        set_exec(path);
        let deadline = std::time::Instant::now() + std::time::Duration::from_secs(10);
        loop {
            match Command::new(path).arg("--version").output() {
                // Only that it *ran* matters here, not what it said.
                Ok(_) => return,
                Err(e)
                    if e.kind() == std::io::ErrorKind::ExecutableFileBusy
                        && std::time::Instant::now() < deadline =>
                {
                    std::thread::sleep(std::time::Duration::from_millis(5));
                }
                Err(e) => panic!("stub '{}' will not run: {e}", path.display()),
            }
        }
    }

    fn scratch(tag: &str) -> PathBuf {
        // Per-process, so parallel `cargo test` runs cannot delete each other's.
        let dir = std::env::temp_dir().join(format!("fc_va_{tag}_{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        dir
    }

    /// The whole point of the design: an include the top file does not itself
    /// contain still changes the key, because the key is the compiler's
    /// expansion. Uses a stub compiler so it runs without OpenVAF.
    #[test]
    fn key_follows_the_include_closure() {
        let dir = scratch("closure");
        let top = dir.join("top.va");
        let inc = dir.join("inc.va");
        write(&top, "`include \"inc.va\"\nmodule m(a); end\n");
        write(&inc, "// version one\n");

        // A stub that does what `--print-expansion` does: concatenate the
        // closure. Real OpenVAF is not needed to prove the key tracks it.
        let stub = dir.join("stub-openvaf");
        write_stub(
            &stub,
            "#!/bin/sh\n\
             for a in \"$@\"; do case \"$a\" in --version) echo 'stub 1'; exit 0;; esac; done\n\
             for a in \"$@\"; do case \"$a\" in *.va) cat \"$a\";; esac; done\n\
             for a in \"$@\"; do case \"$a\" in *.va) d=$(dirname \"$a\");; esac; done\n\
             cat \"$d\"/inc.va\n",
        );

        let opts = VaOptions {
            compiler: Some(stub.clone()),
            ..Default::default()
        };
        let compiler = VaCompiler::find(&opts).expect("stub answers --version");

        let before = cache_key(&compiler, &top, &opts).unwrap();
        write(&inc, "// version two, materially different\n");
        let after = cache_key(&compiler, &top, &opts).unwrap();

        assert_ne!(
            before, after,
            "an edited include must invalidate the cache: a stale .osdi is a \
             silently wrong device"
        );
    }

    /// A compiler upgrade invalidates too — the OSDI it emits is an ABI.
    #[test]
    fn key_follows_the_compiler_version() {
        let dir = scratch("version");
        let top = dir.join("top.va");
        write(&top, "module m(a); end\n");
        let stub = dir.join("stub-openvaf");
        write_stub(&stub, "#!/bin/sh\ncat \"$2\" 2>/dev/null; exit 0\n");

        let opts = VaOptions {
            compiler: Some(stub),
            ..Default::default()
        };
        let mut compiler = VaCompiler::find(&opts).unwrap();
        let before = cache_key(&compiler, &top, &opts).unwrap();
        compiler.version = "some later build".into();
        let after = cache_key(&compiler, &top, &opts).unwrap();
        assert_ne!(before, after);
    }

    /// A missing compiler names what was looked for and both ways to point at
    /// one. The alternative — skipping the model — is a wrong circuit.
    #[test]
    fn missing_compiler_says_how_to_fix_it() {
        let opts = VaOptions {
            compiler: Some(PathBuf::from("/nonexistent/openvaf-r")),
            ..Default::default()
        };
        let msg = VaCompiler::find(&opts).unwrap_err().to_string();
        assert!(msg.contains("/nonexistent/openvaf-r"), "{msg}");
        assert!(msg.contains("--openvaf"), "{msg}");
        assert!(msg.contains("FAIRCHILD_OPENVAF"), "{msg}");
    }

    /// A compiler that is *there* and will not run must not be reported as
    /// missing. "Not found" sends the reader looking for a file that is already
    /// on disk; what they need is the operating system's reason — no execute
    /// bit, not a binary, a directory, or another process still writing it.
    #[test]
    fn a_compiler_that_will_not_run_is_not_reported_as_missing() {
        let dir = scratch("wontrun");

        // There, but not executable.
        let unrunnable = dir.join("openvaf-r");
        write(&unrunnable, "not a program\n");
        let msg = VaCompiler::find(&VaOptions {
            compiler: Some(unrunnable.clone()),
            ..Default::default()
        })
        .unwrap_err()
        .to_string();
        assert!(msg.contains("could not be run"), "{msg}");
        assert!(msg.contains("openvaf-r"), "{msg}");
        assert!(
            !msg.contains("not found") && !msg.contains("No Verilog-A compiler"),
            "a file that exists must not be reported as missing: {msg}"
        );

        // There, runs, and fails — the version is part of every cache key, so a
        // compiler that will not say what it is cannot be used silently either.
        let angry = dir.join("angry-openvaf");
        write_stub(&angry, "#!/bin/sh\necho 'boom' >&2\nexit 3\n");
        let msg = VaCompiler::find(&VaOptions {
            compiler: Some(angry),
            ..Default::default()
        })
        .unwrap_err()
        .to_string();
        assert!(msg.contains("could not be run"), "{msg}");
        assert!(msg.contains("boom"), "the compiler's own words: {msg}");
    }

    /// `--no-va-compile` refuses; it does not quietly load nothing.
    #[test]
    fn no_compile_refuses_rather_than_skips() {
        let dir = scratch("nocompile");
        let top = dir.join("top.va");
        write(&top, "module m(a); end\n");
        let stub = dir.join("stub-openvaf");
        write_stub(&stub, "#!/bin/sh\ncat \"$2\" 2>/dev/null; exit 0\n");

        let opts = VaOptions {
            compiler: Some(stub),
            cache_dir: Some(dir.join("cache")),
            no_compile: true,
            ..Default::default()
        };
        let compiler = VaCompiler::find(&opts).unwrap();
        let err = compile(&compiler, &top, &opts).unwrap_err().to_string();
        assert!(err.contains("no-va-compile"), "{err}");
    }

    /// A `.va` path that is not there fails by name, before any compiler runs.
    #[test]
    fn missing_source_is_named() {
        let opts = VaOptions::default();
        let compiler = VaCompiler {
            path: PathBuf::from("/nonexistent/openvaf-r"),
            version: "stub".into(),
        };
        let err = compile(&compiler, Path::new("/nope/absent.va"), &opts)
            .unwrap_err()
            .to_string();
        assert!(err.contains("absent.va"), "{err}");
    }

    fn set_exec(path: &Path) {
        #[cfg(unix)]
        {
            use std::os::unix::fs::PermissionsExt;
            std::fs::set_permissions(path, std::fs::Permissions::from_mode(0o755)).unwrap();
        }
    }
}
