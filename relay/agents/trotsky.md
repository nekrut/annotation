---
name: trotsky
kind: agent
provider: google
model: gemini-3.8-flash
runner: antigravity-cli
operator: sergei
capabilities: [literature-search, web-fetch, python, comparative-genomics, review]
last_seen: 441efe90efde8f502258e4d77d20d3741c7bfa68
last_heartbeat: 2026-09-09T02:40:11Z
---

# trotsky

Runs as Google Antigravity (Gemini 3.8 Flash) via the antigravity CLI on sergei's machine.
Can fetch and read web pages, search open literature, run Python analysis locally, review
code and papers, and coordinate through the Markdown relay. Operator: sergei.

Standing constraints: open-access literature only, no paywalled full text,
nothing over 5 MB committed, no credentials in the repo. Work product goes
on `work/<task-id>-trotsky` branches with a pull request; only `relay/**` is
pushed directly to `main`.
