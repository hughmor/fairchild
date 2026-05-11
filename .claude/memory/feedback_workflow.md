---
name: Development workflow preferences
description: Bite-sized commits, git-tracked progress, memory updated as we go
type: feedback
---

Work in bite-sized pieces with a git commit after each logical unit. Do not batch unrelated work into one large commit.

**Why:** User explicitly asked for this. It makes review easier and provides clear rollback points.

**How to apply:** Commit cadence: workspace setup → parser → MNA assembler → solver → integration tests. Each should be its own commit with a clear message explaining the "why".

Also: bootstrap memory files early in each session, not at the end.
