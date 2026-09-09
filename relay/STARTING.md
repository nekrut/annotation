# Starting the team

A checklist for the coordinator. Everything the agents need is already in
the repository; this is what you do.

## Once

1. Merge this branch into `main`, or make it `main`. The agents' prompts
   assume `main`.
2. Give each runner push access to the repository. Every agent pushes
   relay commits to `main` and work branches under `work/*`.
3. Decide which vendor runs which name. It does not matter for the
   protocol. Open each CLI in a fresh clone of the repository.
4. Optional but recommended: branch protection on `main` requiring the
   `relay-validate` check, and pull request review for paths outside `relay/`.

## Per agent, first run

Open the agent's CLI (Claude Code, Codex, Gemini CLI, or a script around
any model API) inside its clone and paste the **First run** block from
`relay/prompts/<name>.md`. The agent registers itself, says hello, claims a
review slot, works for a bounded time, and pushes.

Start the four agents a few minutes apart so their first claims do not all
collide. Collisions are handled correctly, but staggering wastes fewer runs.

## Per agent, every later run

Paste the **Every later run** block from the same file, or automate it:

| runner        | automation                                                       |
|---------------|------------------------------------------------------------------|
| Claude Code   | a scheduled Routine with that block as its prompt                |
| Codex, Gemini CLI, others | cron on a machine with the CLI, running it non-interactively with that block |
| a script      | cron: pull, build the prompt from `relay.py inbox`, call the model, apply edits, push |

Every 30 to 60 minutes per agent is a good starting cadence for review
work. The lease length is two hours, so an agent that ticks hourly never
loses its task.

## What you do while they run

- `python3 relay/bin/relay.py status` in any up-to-date clone shows the
  board and who is alive.
- `python3 relay/bin/relay.py inbox --for human` shows questions and alerts
  addressed to you. Answer with `relay.py new --from human --type answer`
  or `--type decision`.
- When a review task reaches `review`, read the artifact and run
  `python3 relay/bin/relay.py release T-human-00N --as human --status done`.
- Merge pull requests for Phase 2 and 3 tasks after at least one agent
  review, then mark them `done` the same way.
- Edit `relay/TASK.md` when the plan changes, and announce it with a
  `decision`.

## Running unattended

`relay/bin/tick.sh <name>` runs one tick for a local agent without a human:
it reads the agent's runner from its roster file, launches the command you
configured for that runner with the recurring prompt, refuses to start if
the previous tick is still running, kills a tick that exceeds 50 minutes,
and logs to `~/relay-logs/<name>-<timestamp>.log`.

Each CLI has a flag for running without interactive approval. Look it up in
that CLI's `--help` and set the full command once per runner in the crontab
environment. The variable name is `RELAY_CMD_` plus the runner with hyphens
replaced by underscores:

```
RELAY_CMD_claude_code='claude -p ...'
RELAY_CMD_codex='codex exec ...'
RELAY_CMD_gemini_cli='gemini -p ...'

5  * * * * /path/to/lenin/relay/bin/tick.sh   lenin
20 * * * * /path/to/engels/relay/bin/tick.sh  engels
35 * * * * /path/to/stalin/relay/bin/tick.sh  stalin
50 * * * * /path/to/trotsky/relay/bin/tick.sh trotsky
```

The minutes are staggered so local agents never push at once, and marx's
cloud Routine runs at seven past. Before relying on cron, run
`relay/bin/tick.sh lenin` once by hand with the variables set and read the
log; if the CLI stops waiting for approval, the command is wrong.

Do not run cron ticks and hand-pasted ticks in the same checkout at the same
time. The script's lock only knows about its own runs. Let any interactive
session finish before the first cron minute for that agent.

Agents often leave uncommitted files in their own `relay/artifacts/<task>/`
between ticks. That is fine: the script pulls with `--autostash`, and the
agent's next tick commits them.

The machine must stay awake. If it sleeps, ticks stop, leases expire after
two hours, and work pauses until the next tick; nothing is lost, but
nothing happens either.

Health check: `python3 relay/bin/relay.py status` shows every agent's last
heartbeat. An agent whose heartbeat is more than two hours old has stopped
ticking; look in `~/relay-logs/` for its last log.
