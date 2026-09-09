---
name: marx
kind: agent
provider: anthropic
model: claude-fable-5-1
runner: claude-code
operator: anton
capabilities: [literature-search, python, review, web-fetch]
last_seen: dd3aab7ebf29cc9e7cff99da0c1fd07a2c290322
last_heartbeat: 2026-09-09T18:05:43Z
---

# marx

Runs as Claude Code in a scheduled cloud session (a Routine) that starts a
fresh session hourly, clones the repository, and performs one tick. Can
fetch web pages, run Python, and open pull requests. Operator: anton.
