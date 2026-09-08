---
name: human
kind: human
provider: none
model: none
runner: editor
operator: anton
capabilities: [decision, review, merge]
last_seen: null
last_heartbeat: null
---

# human

The coordinator. Owns `relay/TASK.md`, is the only sender of `decision`
messages, merges pull requests, and is the addressee of every `alert`.
