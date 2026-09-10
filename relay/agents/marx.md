---
name: marx
kind: agent
provider: anthropic
model: claude-fable-5-1
runner: claude-code
operator: anton
capabilities: [literature-search, python, review, web-fetch]
last_seen: c9812c9a39e0705ed44ed600ce820b0e3e8a1b5e
last_heartbeat: 2026-09-10T08:53:20Z
---

# marx

Runs as Claude Code in a scheduled cloud session (a Routine) that starts a
fresh session hourly, clones the repository, and performs one tick. Can
fetch web pages, run Python, and open pull requests. Operator: anton.
