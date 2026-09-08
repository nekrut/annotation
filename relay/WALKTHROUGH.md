# Walkthrough: one working day on the relay

A concrete trace of the protocol with the four agents `marx`, `engels`,
`lenin`, and `stalin`, plus the coordinator `human`. Nothing here is special
to these names; swap in any roster.

Each agent runs on its own schedule, on its own machine, in its own vendor's
tool. They never talk directly. Everything below is a file landing on `main`.

## 09:00 Kickoff

`human` edits `relay/TASK.md` with the real goal, creates three tasks, and
posts a decision:

```
relay/tasks/T-human-002.md   "Ingest raw input files"          touches: [src/ingest.py]
relay/tasks/T-human-003.md   "Define the output schema"         touches: [docs/schema.md]
relay/tasks/T-human-004.md   "Write the validator"  depends_on: [T-human-003]  touches: [src/validate.py]
relay/messages/20260909T090000Z-human-0002.md   decision: "Charter v1 is final; claim away"
```

One commit, `relay: charter v1 and first tasks`, pushed to `main`.

## 09:15 marx ticks

marx's runner (a Routine, cron job, or GitHub Action) starts it with the
tick prompt. marx pulls, validates, and reads its inbox:

```
$ python3 relay/bin/relay.py inbox --for marx
20260908T173000Z-human-0001  [decision]  Relay protocol v0.1 is live
20260909T090000Z-human-0002  [decision]  Charter v1 is final; claim away
```

It owns nothing, so it picks the highest open task with no unmet
dependencies and no `touches` overlap: `T-human-002`.

```
$ python3 relay/bin/relay.py claim T-human-002 --as marx --push
claimed T-human-002 as marx until 2026-09-09T11:15:00Z
```

marx creates branch `work/T-human-002-marx`, writes the first version of
`src/ingest.py`, pushes the branch, and ends its tick with one note and its
cursor update:

```
$ python3 relay/bin/relay.py new --from marx --to all --type note \
    --title "T-human-002: ingest skeleton on work/T-human-002-marx" --task T-human-002
$ python3 relay/bin/relay.py heartbeat --as marx
$ python3 relay/bin/relay.py inbox --for marx --mark
$ git add relay && git commit -m "relay: marx tick" && git push origin main
```

## 09:20 engels and lenin tick at almost the same time

Both see `T-human-003` open. Both run `claim --push`.

- engels's push lands first.
- lenin's push is rejected. Its `pull --rebase` conflicts on the `status:`
  and `owner:` lines of `T-human-003.md`. The helper aborts the rebase,
  resets to `origin/main`, and prints `lost the race for T-human-003`.
- lenin re-reads the board. `T-human-004` depends on `T-human-003`, which is
  not done, so nothing is claimable. lenin posts a proposal instead:

```
type: proposal   to: [human]
title: "Proposed task: test fixtures for ingest"
```

and ends its tick. No human intervened; git decided.

## 10:05 stalin has a question for marx

stalin is reviewing the charter and needs to know marx's file layout before
it can start on anything. It posts:

```
type: question   to: [marx]   task: T-human-002
title: "Where do parsed records go: one file per input, or one table?"
```

marx is not running right now. The question waits on `main` as a file.

## 10:15 marx ticks again

marx's inbox now contains stalin's question. Answering it is the first thing
marx does, before touching its own task:

```
type: answer   to: [stalin]   reply_to: 20260909T100500Z-stalin-0001
title: "One table, Parquet, one row per record"
```

Then marx renews its lease (`release T-human-002 --as marx --status
in_progress` bumps `lease_until`), does another bounded chunk of work on its
branch, and ends the tick.

## 11:40 engels finishes and asks for review

engels opens a pull request from `work/T-human-003-engels` and moves the
task:

```
$ python3 relay/bin/relay.py release T-human-003 --as engels --status review \
    --pr https://github.com/nekrut/annotation/pull/7
```

`T-human-004` now has its dependency in `review`, not `done`, so it is still
not claimable. That is deliberate: the schema might change in review.

## 12:10 lenin reviews

lenin's tick sees the task in `review` and a PR it did not write. It reads
the diff and posts:

```
type: review   to: [engels, human]   task: T-human-003
title: "Schema: two required fields missing, otherwise good"
```

Reviews are advisory. engels fixes the fields on its next tick and pushes to
the same branch. The task stays in `review`.

## 14:00 human merges

human reads the review thread, merges PR #7, and moves the task:

```
$ python3 relay/bin/relay.py release T-human-003 --as human --status done
```

Only the coordinator (or whoever the charter delegates) moves tasks to
`done`. On its next tick, stalin sees `T-human-004` unblocked and claims it.

## 16:30 something goes wrong

marx's runner crashes mid-tick at 16:30. Its lease on `T-human-002` runs
out at 17:15. At 17:20 engels ticks, sees the expired lease in `validate`
output, and posts a note before taking over:

```
type: note   to: [all]   task: T-human-002
title: "marx's lease on T-human-002 expired; claiming"
```

```
$ python3 relay/bin/relay.py claim T-human-002 --as engels --push
claimed T-human-002 as engels (previous lease by marx expired)
```

engels checks out `work/T-human-002-marx`, reads marx's log entries in the
task file and its notes in the relay, and continues from there. When marx's
runner comes back, its next tick shows it owns nothing, and it moves on.

## What the human sees

At any time, in any clone:

```
$ python3 relay/bin/relay.py status
TASKS
  T-human-002      in_progress  engels     Ingest raw input files
  T-human-004      claimed      stalin     Write the validator
  T-human-003      done         -          Define the output schema
AGENTS
  engels   agent  codex        heartbeat=2026-09-09T17:20:00Z  unread=0
  lenin    agent  gemini-cli   heartbeat=2026-09-09T12:10:00Z  unread=2
  marx     agent  claude-code  heartbeat=2026-09-09T16:30:00Z  unread=4
  stalin   agent  script       heartbeat=2026-09-09T14:05:00Z  unread=1
```

And `git log --oneline -- relay/` is the complete, ordered, attributable
history of who said and did what.

## What each operator has to set up, once

| agent   | runner (fill in `relay/agents/<name>.md`) | wakeup |
|---------|--------------------------------------------|--------|
| marx    | e.g. Claude Code                           | scheduled Routine on `main` |
| engels  | e.g. Codex CLI                             | cron on a machine with the CLI |
| lenin   | e.g. Gemini CLI                            | cron, or GitHub Action on push to `relay/**` |
| stalin  | e.g. a script calling any model API        | cron |

Every runner gets the same prompt: "You are `<name>`. Read
`relay/PROTOCOL.md` and `relay/TASK.md`, then perform one tick as described
in section 6." Every runner needs a full clone and push rights to `main`
(for `relay/**`) and to `work/*` branches.
