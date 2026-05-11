---
name: Testing philosophy - ngspice as ground truth
description: Every new simulator capability must ship with an ngspice golden comparison test
type: feedback
---

Every new simulation capability must have a corresponding ngspice comparison test in the same commit. Do not consider a feature complete without it.

**Why:** User explicitly stressed this as a core workflow requirement. We also caught a real bug (current source sign convention) immediately via the golden tests on Day 1.

**How to apply:**
- Golden netlists go in `tests/golden/*.sp`
- Comparison tests go in `crates/fairchild-core/tests/ngspice_golden.rs`
- Use the `golden_test!` macro defined in that file
- Tolerance: `max(abs_floor, REL_TOL * |expected|)` where REL_TOL=1e-5 (10 ppm). Never use a pure absolute tolerance on quantities that may be large.
- ngspice is found via `find_ngspice()` which checks PATH and `/opt/homebrew/bin/ngspice`
