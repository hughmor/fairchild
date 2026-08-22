"""Refuse to import a compiled extension that is older than its sources.

`cargo build` does not rebuild `fairchild*.so` — only `maturin develop` does — so
in a source checkout the extension can be months older than the Rust it claims to
be. It still imports, and every layer above it assumes it is current. That has
cost real time in both directions:

* the bindings and the CLI once disagreed by a factor of 116 on a noise figure,
  because the extension predated the fix and was reporting confident,
  self-consistent, wrong numbers;
* a Verilog-A runtime test failed on clean master for a whole session against a
  stale mock library, and passed the moment the workspace was rebuilt.

Nothing connected either artefact's freshness to the source it was built from, so
it loaded successfully and the layer above believed it. A hard error naming the
rebuild is strictly better than a wrong number.

Run this file directly for its self-check:

    python3 python/fairchild/_freshness.py
"""

from __future__ import annotations

import datetime
import os
import pathlib

ALLOW_STALE = "FAIRCHILD_ALLOW_STALE"

# Crates that are *not* linked into the extension. Editing the CLI or the C ABI
# cannot change what `import fairchild` does, and flagging it would train people
# to set ALLOW_STALE permanently — which is the failure mode this guards against,
# one level up.
NOT_IN_EXTENSION = ("fairchild-cli", "fairchild-c")

_SOURCE_GLOBS = (
    "crates/**/*.rs",
    "crates/**/Cargo.toml",
    "Cargo.toml",
    "Cargo.lock",
)


def repo_root(package_dir: pathlib.Path) -> pathlib.Path | None:
    """The checkout `package_dir` lives in, or `None` for an installed wheel.

    An installed wheel has no Rust beside it and nothing to be stale against, so
    the check stays out of its way entirely: a false alarm in a released package
    would be worse than the bug this exists to catch.
    """
    root = package_dir.resolve().parent.parent
    if (root / "crates").is_dir() and (root / "Cargo.toml").is_file():
        return root
    return None


def newest_source(root: pathlib.Path) -> tuple[float, pathlib.Path] | None:
    """The most recently modified file the extension is built from.

    Deliberately mtimes rather than a content hash: the question is "was this
    built after what it was built from", and mtime answers it in a few
    milliseconds over ~150 files. A hash answers a question nobody asked — did
    the *content* change — at a cost paid on every import.

    Note a `git checkout` rewrites files and so moves their mtimes. That is not a
    false positive: after switching branches the extension really does not match
    the checkout any more.
    """
    newest: tuple[float, pathlib.Path] | None = None
    for pattern in _SOURCE_GLOBS:
        for path in root.glob(pattern):
            posix = path.as_posix()
            if "/target/" in posix or any(f"/{c}/" in posix for c in NOT_IN_EXTENSION):
                continue
            try:
                mtime = path.stat().st_mtime
            except OSError:
                continue
            if newest is None or mtime > newest[0]:
                newest = (mtime, path)
    return newest


def staleness_error(package_dir: pathlib.Path, ext_suffix: str) -> str | None:
    """The message to refuse with, or `None` if the extension may be imported.

    Returns rather than raises so the self-check below can assert on both
    outcomes without catching exceptions — and so the caller decides whether a
    stale artefact is fatal.
    """
    root = repo_root(package_dir)
    if root is None or os.environ.get(ALLOW_STALE):
        return None

    # The extension *this interpreter* would import, named for its own ABI. That
    # naming is what makes the stale case so easy to hit: two interpreters in one
    # checkout have two `.so` files, and rebuilding for one leaves the other
    # exactly as old as it was.
    ext = package_dir / f"fairchild{ext_suffix}"
    if not ext.is_file():
        return None  # nothing built here; the import itself will say so

    newest = newest_source(root)
    if newest is None:
        return None
    source_mtime, culprit = newest
    built = ext.stat().st_mtime
    if built >= source_mtime:
        return None

    def when(t: float) -> str:
        return datetime.datetime.fromtimestamp(t).strftime("%Y-%m-%d %H:%M")

    return (
        f"the compiled fairchild extension is older than the sources it was built "
        f"from, so importing it would run code that is not in this checkout:\n"
        f"  {ext.name}: {when(built)}\n"
        f"  {culprit.relative_to(root)}: {when(source_mtime)}\n"
        f"Rebuild it:\n"
        f"  maturin develop --release\n"
        f"`cargo build` does not do this — it builds the crate as a plain library "
        f"and leaves the extension module alone. Set {ALLOW_STALE}=1 to import it "
        f"anyway, but a stale extension has already produced confidently wrong "
        f"numbers here once (a noise figure off by 116x)."
    )


def _self_check() -> None:
    """Both directions, on a directory laid out like a checkout.

    Synthetic rather than driving a real build: the thing under test is the
    comparison, and a real `maturin develop` would take minutes to tell us less.
    The negative direction matters as much as the positive one — a guard that
    fires on a fresh extension gets switched off, and then it is worse than none.
    """
    import shutil
    import tempfile
    import time

    suffix = ".cpython-test.so"
    tmp = pathlib.Path(tempfile.mkdtemp(prefix="fairchild_freshness_"))
    try:
        pkg = tmp / "python" / "fairchild"
        crate = tmp / "crates" / "fairchild-core" / "src"
        crate.mkdir(parents=True)
        pkg.mkdir(parents=True)
        (tmp / "Cargo.toml").write_text("[workspace]\n")
        source = crate / "lib.rs"
        source.write_text("// pretend\n")
        ext = pkg / f"fairchild{suffix}"
        ext.write_bytes(b"\x00")

        # Fresh: built after the source.
        os.utime(source, (time.time() - 60, time.time() - 60))
        os.utime(ext, (time.time(), time.time()))
        assert staleness_error(pkg, suffix) is None, "flagged a fresh extension"

        # Stale: the source moved and nothing rebuilt.
        os.utime(source, (time.time() + 60, time.time() + 60))
        msg = staleness_error(pkg, suffix)
        assert msg is not None, "a stale extension was accepted"
        assert "maturin develop" in msg, msg
        assert "lib.rs" in msg, "the message must name what moved: " + msg

        # The escape hatch works, and only when asked for.
        os.environ[ALLOW_STALE] = "1"
        assert staleness_error(pkg, suffix) is None, "ALLOW_STALE was ignored"
        del os.environ[ALLOW_STALE]

        # A crate the extension does not link cannot make it stale.
        cli = tmp / "crates" / "fairchild-cli" / "src"
        cli.mkdir(parents=True)
        cli_src = cli / "main.rs"
        cli_src.write_text("// pretend\n")
        os.utime(source, (time.time() - 60, time.time() - 60))
        os.utime(cli_src, (time.time() + 60, time.time() + 60))
        assert staleness_error(pkg, suffix) is None, "the CLI is not in the extension"

        # An installed wheel has no checkout around it and must be left alone.
        lonely = tmp / "site-packages" / "fairchild"
        lonely.mkdir(parents=True)
        (lonely / f"fairchild{suffix}").write_bytes(b"\x00")
        assert staleness_error(lonely, suffix) is None, "flagged an installed wheel"
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    print("ok: the freshness guard fires on a stale extension and only on that")


if __name__ == "__main__":
    _self_check()
