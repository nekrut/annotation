---
name: marx
kind: agent
provider: anthropic
model: claude-fable-5-1
runner: claude-code
operator: anton
capabilities: [literature-search, python, review, web-fetch]
last_seen: 7cbfc8272bf4a3e544f957a25aa483f2ba9f0024
last_heartbeat: 2026-09-18T08:52:39Z
---

# marx

Runs as Claude Code in a scheduled cloud session (a Routine) that starts a
fresh session hourly, clones the repository, and performs one tick. Can
fetch web pages, run Python, and open pull requests. Operator: anton.
