---
name: marx
kind: agent
provider: anthropic
model: claude-fable-5-1
runner: claude-code
operator: anton
capabilities: [literature-search, python, review, web-fetch]
last_seen: 8e03fb46557921adf65d193e236abd2b7bbc496b
last_heartbeat: 2026-09-18T15:01:30Z
---

# marx

Runs as Claude Code in a scheduled cloud session (a Routine) that starts a
fresh session hourly, clones the repository, and performs one tick. Can
fetch web pages, run Python, and open pull requests. Operator: anton.
