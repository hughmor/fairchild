---
name: Update docs and phase files after each feature commit
description: After committing any feature, update README, user guide, project_fairchild.md, and the relevant phase file
type: feedback
---
After completing and committing any feature for fairchild:

1. **`README.md`** — update the Features table if new elements/analyses/capabilities were added; keep the status line and roadmap current.
2. **`docs/user-guide.md`** — add or update the relevant section (waveform syntax, analysis directives, solver theory, CLI flags, etc.) so the guide stays in sync with the implementation.
3. **`.claude/memory/project_fairchild.md`** — update the Status table (mark items ✅ done, update "next" pointer).
4. **The relevant phase file** (`phase_2.md`, `phase_3.md`, `phase_4.md`) — mark completed steps, note what remains, update any implementation details discovered during the work.

**Why:** Without this, docs drift stale and memory drifts stale. Both cause wasted time — users hit undocumented behaviour, and future sessions re-derive context that's already known.

**How to apply:** Do items 1–4 in the same session as the feature, before ending. Commit docs separately from the feature code if convenient, but don't skip them.
