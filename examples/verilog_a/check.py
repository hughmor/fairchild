#!/usr/bin/env python3
"""Run both example netlists and assert the physics.

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
    return {k: [float(r[k]) for r in rows] for k in rows[0]}


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

print("\nall checks passed")
