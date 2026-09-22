---
id: T-human-001
title: Smoke test the relay: every agent registers and says hello
status: done
owner: null
created_by: human
created: 2026-09-08T17:30:00Z
lease_until: null
depends_on: []
touches: [relay/agents]
pr: null
---

## Goal

Prove that every agent can complete one full tick: pull, validate, read the
inbox, write a message, update its roster file, push.

## Definition of done

Each agent named in the roster has a `relay/agents/<name>.md` with
`last_heartbeat` set, and has sent one `note` to `all` with the title
`hello from <name>` describing its runner and capabilities. The coordinator
then moves this task to `done`.

This task is intentionally shared: it is the one task that is not claimed,
because every agent has to do its own part. Do not claim it.

## Log

- 2026-09-08 human: created.
- 2026-09-22 human: -> done.
