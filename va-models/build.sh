#!/bin/bash
# Compile all Verilog-A models to OSDI shared libraries.
# Requires OpenVAF-Reloaded (openvaf-r) on PATH or OPENVAF env var.
#
# Usage:
#   cd va-models
#   DYLD_LIBRARY_PATH=/opt/homebrew/opt/llvm@18/lib ./build.sh

set -e

OPENVAF="${OPENVAF:-openvaf-r}"
OUTDIR="build"
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

mkdir -p "$OUTDIR"

compile() {
    local src="$1"
    local name
    name="$(basename "$src" .va)"
    echo "Compiling $src..."
    # -I "$SCRIPT_DIR" ensures `include "disciplines/optical.vams" resolves from va-models/
    "$OPENVAF" -I "$SCRIPT_DIR" "$src" -o "$OUTDIR/$name.osdi"
}

# Electronic models
compile electronic/diode_shockley.va
compile electronic/nmos_l1.va
compile electronic/pmos_l1.va

# Photonic models (Phase 2+)
compile photonic/cw_laser.va
compile photonic/waveguide.va
compile photonic/directional_coupler.va
compile photonic/photodetector.va

echo "Done. Outputs in $OUTDIR/"
