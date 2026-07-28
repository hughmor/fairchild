#!/usr/bin/env python3
"""Run every example netlist and assert the physics.

    ./build.sh && ./check.py

Nothing here is a substitute for the crate test suite — it is the smallest
thing that fails if the Verilog-A path stops working end to end.
"""
import csv
import io
import math
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]


def run(netlist, probe, extra=()):
    exe = next(
        (p for p in (ROOT / "target/release/fairchild", ROOT / "target/debug/fairchild")
         if p.exists()),
        None,
    )
    if exe is None:
        sys.exit("build fairchild first:  cargo build --release -p fairchild-cli")
    out = subprocess.run(
        [str(exe), "-f", str(HERE / netlist), "--probe", probe, "--format", "csv", *extra],
        cwd=HERE, capture_output=True, text=True,
    )
    if out.returncode != 0:
        sys.exit(f"{netlist} failed:\n{out.stderr}")
    rows = list(csv.DictReader(io.StringIO(out.stdout)))
    # `.op` output leads with a non-numeric "analysis" column; drop it.
    cols = {}
    for k in rows[0]:
        try:
            cols[k] = [float(r[k]) for r in rows]
        except ValueError:
            pass
    return cols


if not sorted((HERE / "build").glob("*.osdi")):
    sys.exit(
        "no compiled models in build/ — run ./build.sh first.\n"
        "It needs OpenVAF-Reloaded:  OPENVAF=/path/to/openvaf-r ./build.sh"
    )


def close(a, b, tol, what):
    assert abs(a - b) <= tol, f"{what}: {a:.6g} != {b:.6g} (tol {tol:g})"
    print(f"  ok  {what}: {a:.6g}")


print("rectifier.sp — Verilog-A diode + native R, C, V")
r = run("rectifier.sp", "v(in),v(out)")
vout, vin = r["V(out)"], r["V(in)"]
peak = max(vout)
# Rectified output tracks the source peak less one diode drop, and never
# exceeds it.  Vf lands near 0.9 V at ~0.4 mA through Is=1e-14, N=1.
close(peak, max(vin) - 0.9, 0.15, "peak V(out) = Vpk - Vf")
assert peak <= max(vin), "output exceeded the source — diode is conducting backwards"
# Ripple: the cap discharges through Rload between conduction windows, and the
# floor must stay above a bare-RC decay over one full period (1 ms / 10 ms).
settled = vout[len(vout) // 2:]
floor = min(settled)
assert floor > peak * math.exp(-1e-3 / (1e-6 * 10e3)), "droop faster than Rload*Cload"
assert floor < peak, "no ripple at all — is Cload connected?"
print(f"  ok  ripple {peak - floor:.4g} V on a {peak:.4g} V rail")

print("cmos_inverter.sp — Verilog-A transistors via .model cards")
c = run("cmos_inverter.sp", "v(in),v(out)")
vin_c, vout_c = c["V(in)"], c["V(out)"]
# Edges overshoot slightly — that is real Miller feedthrough through CGDO, so
# bound it rather than demanding the rails exactly.
assert -0.3 < min(vout_c) and max(vout_c) < 3.6, \
    f"output left the rails by more than Miller overshoot: {min(vout_c):.3g}..{max(vout_c):.3g}"
# Inverting, sampled just before each edge so the output has settled.
# PULSE(0 3.3 1n 200p 200p 4n 10n): high 1.2-5.2 ns, low again from 5.4 ns.
def at(t):
    i = min(range(len(c["time"])), key=lambda k: abs(c["time"][k] - t))
    return vin_c[i], vout_c[i]

for t, want_in, want_out in ((5.0e-9, 3.3, 0.0), (9.9e-9, 0.0, 3.3)):
    got_in, got_out = at(t)
    close(got_in, want_in, 1e-9, f"V(in) at {t * 1e9:.1f} ns")
    close(got_out, want_out, 1e-3, f"V(out) at {t * 1e9:.1f} ns (inverted)")
# That the card and instance params actually reach the device is pinned
# properly (and without OpenVAF) by crates/fairchild-osdi/tests/osdi_model_card.rs.
# `--param` cannot help here: it only reaches X/R/C/L elements, not M/Q/D.

print("eam_link.sp — Verilog-A modulator inside a native photonic link")
e = run("eam_link.sp", "v(drv),v(eam_a),v(pd_out)")
va, pd = e["V(eam_a)"], e["V(pd_out)"]
# Unbiased: the modulator only inserts il_dB, and generates no photocurrent,
# so the drive node sits exactly at the source.
close(va[0], 0.0, 1e-9, "V(eam_a) at 0 V drive (no photocurrent unbiased)")
hi, lo = max(pd), min(pd)
er = 10 * math.log10(hi / lo)
assert er > 5.0, f"only {er:.3g} dB of extinction — is the modulator being driven?"
print(f"  ok  extinction ratio {er:.4g} dB")
# The detected level tracks the bias exactly: everything but the electro-
# absorption term is common to both levels, so the ratio is 10^(-ea_dB/10)
# with ea_dB = er_dB*(vr/v_full)^2.  er_dB=10, v_full=2 from the netlist.
i = va.index(min(va))
ea_dB = 10.0 * (min(va) / -2.0) ** 2
close(pd[i] / hi, 10 ** (-ea_dB / 10), 2e-3,
      "P(detected) ratio at peak bias = 10^(-er_dB*(vr/v_full)^2/10)")
# Recovery: back to the unbiased level once the drive returns to 0.
close(pd[-1], pd[0], 1e-3, "V(pd_out) recovers after the pulse")

print("va_link.sp — a Mach-Zehnder built entirely in Verilog-A")
# Sweep just over one FSR. FSR = lambda^2/(n_g*dL) = 1550^2/(4.2*20e3) ~ 28.6 nm.
sweep = []
for wl in range(1540, 1576, 5):
    p = [f"--param=X{e}.wavelength_nm={wl}" for e in ("laser", "wgtop", "wgbot")]
    r = run("va_link.sp", "v(bar_i),v(cross_i)", p)
    sweep.append((wl, r["V(bar_i)"][0], r["V(cross_i)"][0]))

# Unitary coupler + lossy arms: the two ports are complementary about the mean
# arm transmission, at every wavelength. This is the whole-chain power budget.
loss = [10 ** (-2.0 * (l * 1e-4) / 10) for l in (400.0, 420.0)]
total = sum(loss) / 2
for wl, bar, cross in sweep:
    close(bar + cross, total, 5e-4, f"bar+cross at {wl} nm")
# …and it has to actually interfere, or "complementary" is trivially satisfied.
bars = [b for _, b, _ in sweep]
assert max(bars) - min(bars) > 0.9 * total, "no MZI fringe — arms are not interfering"
print(f"  ok  fringe {min(bars):.4g} .. {max(bars):.4g} V over 1540-1575 nm")

print("va_waveguide vs native fc_waveguide")
# The promoted model fixes a factor-2 loss bug the legacy tree still has; this
# is the assertion that keeps it fixed. 1 mm at 3 dB/cm must pass 10^-0.3.
wg = run("wg_compare.sp", "v(va_i),v(native_i)")
va, native = wg["V(va_i)"][0], wg["V(native_i)"][0]
close(va, native, 1e-9, "detected power, Verilog-A vs native")
# V = responsivity * P * R_load, so P/P_in = V / (0.8 A/W * 1k) with P_in = 1 mW.
close(va / (0.8 * 1e3) / 1e-3, 10 ** (-0.3 / 10), 2e-4,
      "1 mm at 3 dB/cm passes 10^(-0.3/10)")

print("\nall checks passed")
