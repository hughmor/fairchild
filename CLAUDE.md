# fairchild — Claude Instructions

## Memory

At the start of each session, read `.claude/memory/project_fairchild.md` for project context, phase status, key architectural decisions, and build commands.

When implementing a feature for a specific phase, also read the relevant phase file:
- Phase 2 (photonic discipline): `.claude/memory/phase_2.md`
- Phase 3 (Python bindings): `.claude/memory/phase_3.md`
- Phase 4 (differentiable sim): `.claude/memory/phase_4.md`
- Phases 5+: `PLAN.md` Parts 3 and 5–7

**After completing and committing any feature**: update the status table in `project_fairchild.md` and mark completed steps in the relevant phase file. This is not optional — it keeps context accurate for future sessions.

## Architecture reference

`PLAN.md` contains the research landscape (Part 1), full architecture diagrams (Part 2), technical risks (Part 4), and differentiator analysis (Part 5). Read it when making architectural decisions or when the user asks about design rationale. Don't read it for routine feature work.
