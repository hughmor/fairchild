#!/usr/bin/env bash
# Regenerate the KiCad IPC protobuf bindings from a KiCad source checkout.
#
# The bindings MUST come from the same commit as the KiCad you run, or the
# schematic messages drift silently. The PyPI `kicad-python` package is not a
# substitute: as of 0.7.1 its schematic wrappers import names its own generated
# protos don't define, so `import kipy.schematic` raises ImportError.
#
#   ./regen_protos.sh [/path/to/kicad/src]     # default: ~/Local/src/kicad
set -euo pipefail

SRC="${1:-$HOME/Local/src/kicad}/api/proto"
OUT="$(cd "$(dirname "$0")" && pwd)/_proto"
PY="$(cd "$(dirname "$0")/../../.." && pwd)/.venv/bin/python"

[ -d "$SRC" ] || { echo "no proto dir at $SRC" >&2; exit 1; }

rm -rf "$OUT" && mkdir -p "$OUT"
(cd "$SRC" && find . -name '*.proto' | sed 's|^\./||') | \
    xargs "$PY" -m grpc_tools.protoc -I "$SRC" --python_out="$OUT" --pyi_out="$OUT"

# protoc emits root-relative imports ("from common.types import base_types_pb2"),
# which would need _proto/ on sys.path and would squat the top-level names
# `common`, `board`, `schematic`. Rewrite them to this package instead.
PKG=fairchild.kicad._proto
find "$OUT" -name '*_pb2.py' -o -name '*_pb2.pyi' | while read -r f; do
    /usr/bin/sed -i '' -E \
        -e "s/^from (common|board|schematic)([._a-zA-Z0-9]*) import /from $PKG.\1\2 import /" \
        -e "s/^import (common|board|schematic)([._a-zA-Z0-9]*) as /import $PKG.\1\2 as /" \
        "$f"
done
find "$OUT" -type d -exec touch {}/__init__.py \;

echo "generated $(find "$OUT" -name '*_pb2.py' | wc -l | tr -d ' ') modules from $SRC"
git -C "$(dirname "$SRC")/.." describe --tags 2>/dev/null | sed 's/^/kicad source: /'
