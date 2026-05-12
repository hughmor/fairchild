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
mkdir -p "$OUTDIR"

compile() {
    local src="$1"
    local name="$(basename "$src" .va)"
    echo "Compiling $src..."
    "$OPENVAF" "$src" -o "$OUTDIR/$name.osdi"
}

# Electrical models
compile diode_shockley.va
compile nmos_l1.va
compile pmos_l1.va

# Phase 2 photonic models
compile cw_laser.va
compile waveguide.va
compile directional_coupler.va
compile photodetector.va

echo "Done. Outputs in $OUTDIR/"
