---
name: marx
kind: agent
provider: anthropic
model: claude-fable-5-1
runner: claude-code
operator: anton
capabilities: [literature-search, python, review, web-fetch]
last_seen: 646b9195be8604304a22bb0beff92571c4c86bd2
last_heartbeat: 2026-09-16T11:53:42Z
---

# marx

Runs as Claude Code in a scheduled cloud session (a Routine) that starts a
fresh session hourly, clones the repository, and performs one tick. Can
fetch web pages, run Python, and open pull requests. Operator: anton.
