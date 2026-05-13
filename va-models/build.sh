#!/bin/bash
# Compile all Verilog-A models to OSDI shared libraries.
# Requires OpenVAF-Reloaded (openvaf-r) on PATH or OPENVAF env var.
#
# Usage:
#   cd va-models
#   DYLD_LIBRARY_PATH=/opt/homebrew/opt/llvm@18/lib ./build.sh
#
# To compile a subset, set MODEL_FILTER (space-separated partial names):
#   MODEL_FILTER="mrr_modulator pn_phase" ./build.sh

set -e

OPENVAF="${OPENVAF:-openvaf-r}"
OUTDIR="build"
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

mkdir -p "$OUTDIR"

compile() {
    local src="$1"
    local name
    name="$(basename "$src" .va)"

    # Apply optional name filter
    if [ -n "${MODEL_FILTER:-}" ]; then
        local match=0
        for pat in $MODEL_FILTER; do
            if echo "$name" | grep -q "$pat"; then match=1; break; fi
        done
        [ $match -eq 1 ] || return 0
    fi

    echo "  Compiling $src → $OUTDIR/$name.osdi"
    # -I "$SCRIPT_DIR" ensures `include "disciplines/optical.vams" resolves from va-models/
    "$OPENVAF" -I "$SCRIPT_DIR" "$src" -o "$OUTDIR/$name.osdi"
}

echo "=== Fairchild Verilog-A model build ==="
echo "Compiler: $(${OPENVAF} --version 2>&1 | head -1 || echo '(unknown)')"
echo ""

echo "-- Electronic models --"
compile electronic/diode_shockley.va
compile electronic/nmos_l1.va
compile electronic/pmos_l1.va

echo ""
echo "-- Photonic passive / sources --"
compile photonic/cw_laser.va
compile photonic/waveguide.va
compile photonic/directional_coupler.va

echo ""
echo "-- Photodetector models --"
compile photonic/photodetector.va
compile photonic/photodetector_l2.va

echo ""
echo "-- PN junction phase shifters --"
compile photonic/pn_phase_shifter_l1.va
compile photonic/pn_phase_shifter_l2.va

echo ""
echo "-- Thermo-optic phase shifters --"
compile photonic/thermo_phase_shifter_l1.va
compile photonic/thermo_phase_shifter_l2.va

echo ""
echo "-- MRR modulators (PN junction) --"
compile photonic/mrr_modulator_l1.va
compile photonic/mrr_modulator_l2.va
compile photonic/mrr_modulator_l3.va

echo ""
echo "-- MRR heater-tuned (N-doped) --"
compile photonic/mrr_heater_l1.va
compile photonic/mrr_heater_l2.va

echo ""
echo "Done. Outputs in $SCRIPT_DIR/$OUTDIR/"
