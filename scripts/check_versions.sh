#!/usr/bin/env bash
#
# Assert that every version string in the workspace agrees, and — when given an
# argument — that a release tag agrees with them too.
#
#   scripts/check_versions.sh          # internal consistency
#   scripts/check_versions.sh v0.3.1   # ...and that the tag matches
#
# Why this exists: publishing to crates.io requires a `version` on every
# internal dependency (it replaces `path` in the published manifest), and Cargo
# offers no `version.workspace = true` inside [workspace.dependencies]. So the
# version is written once in [workspace.package] and copied once per internal
# dependency, and nothing but this script notices when a copy is missed.
#
# The failure it prevents is not a build error. A stale copy publishes crate X
# depending on an *older* X-dep that is still on crates.io, so it resolves, it
# compiles, and the user gets a mismatched pair. v0.3.0 was already tagged once
# against a tree that said 0.2.0; that is the same class of mistake.

set -euo pipefail

root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
manifest="$root/Cargo.toml"

want="$(sed -n '/^\[workspace\.package\]/,/^\[[^w]/ s/^version *= *"\([^"]*\)".*/\1/p' \
  "$manifest" | head -1)"
if [ -z "$want" ]; then
  echo "check_versions: no version found in [workspace.package] of $manifest" >&2
  exit 1
fi

fail=0

# Every `fairchild-*` line in [workspace.dependencies] must carry exactly $want.
deps="$(sed -n '/^\[workspace\.dependencies\]/,/^\[/p' "$manifest" | grep '^fairchild-' || true)"
if [ -z "$deps" ]; then
  echo "check_versions: no fairchild-* entries in [workspace.dependencies] —" \
    "either the table moved or the internal deps stopped being shared" >&2
  exit 1
fi

while IFS= read -r line; do
  name="$(printf '%s' "$line" | sed 's/ *=.*//')"
  got="$(printf '%s' "$line" | sed -n 's/.*version *= *"\([^"]*\)".*/\1/p')"
  if [ -z "$got" ]; then
    echo "FAIL $name has no version — cargo publish will refuse it" >&2
    fail=1
  elif [ "$got" != "$want" ]; then
    echo "FAIL $name is \"$got\", [workspace.package] is \"$want\"" >&2
    fail=1
  fi
done <<<"$deps"

# No crate may pin its own version: they all inherit, and one that stopped
# inheriting would drift silently past the check above.
while IFS= read -r f; do
  if grep -qE '^version *= *"' "$f"; then
    echo "FAIL ${f#"$root/"} sets its own version; use version.workspace = true" >&2
    fail=1
  fi
done < <(find "$root/crates" -mindepth 2 -maxdepth 2 -name Cargo.toml)

# Optional: a release tag, with or without the leading `v`.
if [ $# -ge 1 ]; then
  tag="${1#v}"
  if [ "$tag" != "$want" ]; then
    echo "FAIL tag \"$1\" does not match workspace version \"$want\"" >&2
    fail=1
  fi
fi

if [ "$fail" -ne 0 ]; then
  echo "check_versions: FAILED (workspace version is $want)" >&2
  exit 1
fi

echo "check_versions: ok, everything says $want${1:+ (tag $1)}"
