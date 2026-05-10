---
name: fairchild project context
description: Context for the electronic-photonic co-simulation suite project (fairchild)
type: project
---

Open-source EO co-simulation suite targeting SiPh. Rust core, PyO3 Python bindings, CLI.

**Why:** No open-source time-domain EO co-simulator exists. Cadence Spectre Photonics is the only working implementation; it's proprietary and extremely expensive.

**Key architectural decisions:**
- Verilog-A compiler: OpenVAF-Reloaded (arpadbuermen fork on Codeberg) → OSDI v0.4 shared libs
- Solver: Port VACASK (C++, Codeberg: arpadbuermen/VACASK, FOSDEM 2025) architecture to Rust rather than deriving from scratch
- Simulation mode: time-domain DAE first-class; S-matrix augmentation for passive components without time-domain models
- Differentiable: design adjoint sensitivity from the start (not bolted on)
- Photonic signal: complex envelope (SVEA), two nodes per optical port (amplitude + phase), following MIT Sorace-Agaskar 2015 convention
- Netlist formats: SPICE, Verilog-AMS, YAML/TOML
- PDK targets: SiEPIC_EBeam_PDK (primary open), GF45SPCLO (stretch), custom .va models

**Application domain:** Silicon photonics (SiPh) ICs.

**Timeline:** 6-12 months to working prototype.

**Phase 0 falsifiability gate (Weeks 1-6):** CMOS inverter transient via BSIM4 OSDI model.

**Key open questions to validate immediately:**
- Does OpenVAF-Reloaded need modification to handle a custom optical discipline? (Week 1-2 task)
- Can OSDI v0.4 carry optical port semantics without compiler changes?

**Plan document:** `/Users/hugh/Local/src/veriloga/PLAN.md`

**How to apply:** Use this context when implementing any part of the simulator. The plan doc is the authoritative reference.
