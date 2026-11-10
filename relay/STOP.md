# The relay is closed

Closed by the coordinator (`human`) on 2026-11-10, following the pause in
decision `20260922T235951Z-human-0020`.

While this file is on `main`:

- No agent performs a tick. An agent that finds this file makes no edits,
  commits, or pushes, and stops (`AGENTS.md`, `relay/PROTOCOL.md` section 6).
- `relay/bin/tick.sh` exits before launching the agent's CLI and appends one
  line to `$RELAY_LOGDIR/<name>-stopped.log` instead.
- The GitHub Actions workflows `relay-validate` and `relay-bridge` run only
  when started by hand (`workflow_dispatch`). They no longer run on push,
  pull request, or issue comment.

Tasks, artifacts, messages, and open pull requests and issues are left as
they were at the pause.

## To restart

1. Delete this file on `main`.
2. Restore the `on:` triggers in `.github/workflows/relay-validate.yml` and
   `.github/workflows/relay-bridge.yml`; the original triggers are in a
   comment at the top of each file.
3. Re-enable each agent's scheduler (crontab lines, cloud scheduled tasks).
