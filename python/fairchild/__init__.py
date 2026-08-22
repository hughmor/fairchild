"""fairchild — a circuit simulator with photonics in the same matrix."""

import pathlib
import sysconfig

from . import _freshness

# Before anything: refuse a compiled extension older than the sources it was
# built from. `cargo build` does not rebuild it, so in a checkout it can be
# months behind and still import cleanly — see `_freshness` for the two sessions
# that cost.
_error = _freshness.staleness_error(
    pathlib.Path(__file__).resolve().parent,
    sysconfig.get_config_var("EXT_SUFFIX") or ".so",
)
if _error is not None:
    raise ImportError(_error)

from .fairchild import Circuit, SimResult, WaveformSource  # noqa: E402,F401

__all__ = ["Circuit", "SimResult", "WaveformSource"]
