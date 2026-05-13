---
name: feedback-examples
description: How to handle example scripts — plot output directory and conventions
metadata:
  type: feedback
---

All example scripts that generate plots must save them into `docs/plots/`.

**Why:** User preference stated 2026-05-12 — keeps generated assets in one findable place.

**How to apply:** Any Python (or other) script created as an example or utility that calls matplotlib savefig (or equivalent) should use `docs/plots/<name>.png` as the output path, not `examples/` or `/tmp/`.
