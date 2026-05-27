"""
demo_circuit.py — Basic fairchild Python API demo using the ring resonator.

NOTE: The ring resonator uses OSDI-compiled photonic models (cw_laser,
directional_coupler, waveguide, photodetector).  Those models must be compiled
(cd legacy/va-models && bash build.sh) and the .osdi paths in the netlist must be valid.
Without compiled models the simulation will raise a RuntimeError for
'unknown model'; the API calls themselves are still illustrative.

Build & install the Python extension first:
    pip install maturin
    maturin develop          # installs fairchild into the current venv

Then run:
    python examples/photonic/demo_circuit.py
"""

import pathlib
import sys

try:
    import fairchild as fc
except ImportError:
    sys.exit(
        "fairchild not installed — run 'maturin develop' from the repo root first."
    )

# Path to the ring resonator netlist (relative to repo root).
NETLIST = pathlib.Path(__file__).parent / "ring_resonator_sweep.sp"

# ---------------------------------------------------------------------------
# 1. Load the netlist
# ---------------------------------------------------------------------------
ckt = fc.Circuit()
ckt.load(str(NETLIST))
print(f"Loaded: {NETLIST.name}")

# ---------------------------------------------------------------------------
# 2. Override the laser wavelength and power programmatically
# ---------------------------------------------------------------------------
ckt.set_param("Xlaser", "wavelength_nm", 1550.0)
ckt.set_param("Xlaser", "power_mW", 1.0)

# ---------------------------------------------------------------------------
# 3. Run a DC operating point
# ---------------------------------------------------------------------------
print("\nRunning DC operating point at 1550.0 nm ...")
try:
    result = ckt.run("op")
    print("Available signals:", result.signals())
    ph_a = result["V(ph_a)"]
    print(f"  V(ph_a) = {ph_a[0]:.6e} V  (photodetector output)")
except RuntimeError as exc:
    print(f"  Simulation error: {exc}")
    print("  (Compile OSDI models with 'cd legacy/va-models && bash build.sh' to run this demo.)")

# ---------------------------------------------------------------------------
# 4. Parametric wavelength sweep
# ---------------------------------------------------------------------------
wavelengths = [1544.0, 1546.0, 1548.0, 1550.0, 1552.0, 1554.0]
print(f"\nSweeping wavelength over {wavelengths} nm ...")
try:
    results = ckt.sweep("Xlaser.wavelength_nm", wavelengths, "op")
    for wl, res in zip(wavelengths, results):
        v = res["V(ph_a)"][0]
        print(f"  wavelength_nm={wl:.1f}  V(ph_a)={v:.6e} V")
except RuntimeError as exc:
    print(f"  Sweep error: {exc}")
    print("  (Compile OSDI models with 'cd legacy/va-models && bash build.sh' to run this demo.)")
