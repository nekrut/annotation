# annotation

An experiment in multi-agent research: five autonomous AI agents from
different vendors, coordinated by one human, working on a single scientific
problem and communicating only through Markdown files committed to this
repository.

## The problem

Eukaryotic gene prediction is still slow, expensive, and clade-specific.
NCBI's EGAPx needs 32 CPUs and 256 GB of RAM and excludes fungi, protists,
and nematodes; on Galaxy, 1,409 EGAPx runs consumed about 416,000 CPU-hours
with a 45% failure rate. Yet a single comparative signal, the KA/KS ratio of
human-mouse alignment windows, already separates coding from non-coding
sequence with a few percent error (Nekrutenko, Makova, Li 2002).

The goal is a small, fast, geometrically principled model, in the spirit of
HyphAeon (about 2M parameters, alignment plus tree as coordinates), that
predicts gene structure accurately in any eukaryotic genome without
per-clade retraining. Full charter: [`relay/TASK.md`](relay/TASK.md).

## The team

| name      | runner                          | role                          |
|-----------|---------------------------------|-------------------------------|
| `human`   | Anton, with an editor           | coordinator; owns the charter, issues decisions, merges |
| `marx`    | Claude Code, scheduled cloud session (hourly Routine) | agent |
| `lenin`   | Claude Code, local CLI          | agent                         |
| `engels`  | OpenAI Codex, local CLI         | agent                         |
| `stalin`  | OpenAI Codex, local CLI         | agent                         |
| `trotsky` | Google Antigravity, local CLI   | agent                         |

Names are stable identities; the model and runner behind each name are
recorded in `relay/agents/<name>.md` and can change.

## How the collaboration works

Git is the message bus. There is no server and no shared runtime. Each agent
wakes up on its own schedule, pulls `main`, reads its inbox, does a bounded
amount of work, writes messages and task updates as new Markdown files, and
pushes. Concurrency is safe because every file has exactly one writer, and
a race for a task is settled by whose push git accepts first.

- [`relay/PROTOCOL.md`](relay/PROTOCOL.md): the contract every agent follows.
- [`relay/TASK.md`](relay/TASK.md): the charter, owned by the coordinator.
- [`relay/tasks/`](relay/tasks/): one file per unit of work, with owner and lease.
- [`relay/messages/`](relay/messages/): the conversation, one file per message.
- [`relay/agents/`](relay/agents/): the roster; each agent's own cursor and heartbeat.
- [`relay/artifacts/`](relay/artifacts/): review documents and other large outputs.
- [`relay/prompts/`](relay/prompts/): the instructions each agent was started with.
- [`relay/WALKTHROUGH.md`](relay/WALKTHROUGH.md): a worked example of one day.
- [`relay/STARTING.md`](relay/STARTING.md): the coordinator's checklist.
- [`relay/bin/relay.py`](relay/bin/relay.py): helper for validating, reading
  inboxes, posting messages, and claiming tasks. Standard-library Python only.

`AGENTS.md` at the root is the entry point every agent CLI reads first;
`CLAUDE.md`, `GEMINI.md`, and `.github/copilot-instructions.md` point to it.

## Plan

1. **Review** (now). Each agent independently reviews the gene prediction
   literature and the software that is still alive, blind to the others.
2. **Synthesis, benchmark, data, costs.** One merged review; a cross-clade
   benchmark with held-out species; an inventory of alignments and tracks
   (UCSC, Ensembl, RefSeq, Zoonomia); measured costs of existing tools; the
   KA/KS baseline as the floor to beat.
3. **Design.** Two or three candidate geometric models and a recommendation,
   ending in a coordinator decision.
4. **Prototype.** Defined by a charter revision after the decision.

## Following along

The commit history of `relay/` is the complete, attributable record of who
said and did what. From any clone:

```
git pull --rebase origin main
python3 relay/bin/relay.py status           # board, owners, heartbeats
python3 relay/bin/relay.py inbox --for human  # questions and alerts for the coordinator
```

Every push runs `relay-validate` (GitHub Actions), which rejects malformed
messages and tasks in impossible states.
