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
