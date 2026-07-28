#!/bin/bash
# Compile the example Verilog-A models to OSDI shared libraries.
#
# Needs OpenVAF-Reloaded.  Point OPENVAF at the binary if it is not on PATH:
#   OPENVAF=~/src/OpenVAF-Reloaded/target/release/openvaf-r ./build.sh
#
# On macOS openvaf-r also needs LLVM 18 on the dynamic loader path; the export
# below is a no-op elsewhere.
set -euo pipefail

cd "$(dirname "$0")"

OPENVAF="${OPENVAF:-openvaf-r}"
export DYLD_LIBRARY_PATH="${DYLD_LIBRARY_PATH:-/opt/homebrew/opt/llvm@18/lib}"

mkdir -p build

for src in models/va_diode.va models/va_eam.va; do
    name="$(basename "$src" .va)"
    echo "  $src -> build/$name.osdi"
    # -I models/ so `include "optical.vams"` resolves.
    "$OPENVAF" -I models "$src" -o "build/$name.osdi"
done

echo "Done. Run ./check.py to verify."
