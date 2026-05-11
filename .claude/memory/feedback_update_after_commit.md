---
name: Update phase files after each feature commit
description: After committing any feature, update project_fairchild.md status table and the relevant phase file
type: feedback
originSessionId: 8fabe026-50fc-44f1-81d6-cbca5aa71816
---
After completing and committing any feature for fairchild, always update memory:

1. **`project_fairchild.md`** — update the Status table (mark items ✅ done, update "next" pointer)
2. **The relevant phase file** (`phase_2.md`, `phase_3.md`, `phase_4.md`) — mark completed steps, note what remains, update any implementation details discovered during the work

**Why:** Without this, the memory drifts stale and future sessions have to re-read PLAN.md (693 lines) or the full codebase to reconstruct what's been done. The phase files are the single source of truth for roadmap state; PLAN.md is the research/architecture reference.

**How to apply:** After the git commit succeeds, write the memory updates before ending the session. This is not optional — it's what makes the context system useful across sessions.
