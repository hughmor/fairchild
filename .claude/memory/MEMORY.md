# Memory index

- [project context](project_fairchild.md) — phase status table, crate structure, key architectural decisions, build commands
- [Phase 2: Photonic Discipline](phase_2.md) — optical discipline plan, SVEA signal model, keikawa model porting, carrier freq design
- [Phase 3: Python Bindings](phase_3.md) — PyO3 API design, numpy waveform bridge, maturin packaging
- [Phase 4: Differentiable Simulation](phase_4.md) — adjoint method, OSDI parameter Jacobian, gradient API design
- [Update after commits](feedback_update_after_commit.md) — update README, user guide, phase files, and status table after every feature commit
- [Testing philosophy](feedback_testing.md) — every new capability needs an ngspice golden test in the same commit; golden_test! macro, 10 ppm tolerance
- [Commit workflow](feedback_workflow.md) — bite-sized commits per logical unit; commit cadence and ordering preferences
