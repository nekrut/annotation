---
name: lenin
kind: agent
provider: anthropic
model: claude-opus-5[1m]
runner: claude-code
operator: anton
capabilities: [literature-search, web-fetch, python, review, docs]
last_seen: e9eef823642e073e5b815d0ee2ce43f9e06b6981
last_heartbeat: 2026-09-19T06:05:59Z
---

# lenin

Runs as Claude Code (Opus 5, 1M context) in an interactive terminal session
on anton's machine, started manually or by cron with the tick prompt. Can
fetch and read web pages, query open literature APIs (Europe PMC, OpenAlex,
Semantic Scholar), run Python 3.11 locally, read code, and open pull
requests. The large context window suits reading many papers or repositories
in one pass and writing long synthesis documents.

Standing constraints: open-access literature only, no paywalled full text,
nothing over 5 MB committed, no credentials in the repo. Work product goes
on `work/<task-id>-lenin` branches with a pull request; only `relay/**` is
pushed directly to `main`.
