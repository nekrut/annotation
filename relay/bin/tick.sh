#!/usr/bin/env bash
# Run one relay tick for a local agent, unattended. Usage: tick.sh NAME
#
# Reads the runner from relay/agents/NAME.md and launches the command you
# configured for that runner, non-interactively, with the recurring prompt.
# Uses a lock so ticks never overlap, kills a tick that runs too long, and
# logs each run to $RELAY_LOGDIR (default ~/relay-logs).
#
# You must set one environment variable per runner, naming the CLI and the
# flags that let it run without permission prompts. Look the flag up in the
# CLI's own --help; do not guess. Examples of the shape (not the flags):
#   RELAY_CMD_claude_code='claude -p <unattended flag>'
#   RELAY_CMD_codex='codex exec <unattended flag>'
#   RELAY_CMD_gemini_cli='gemini -p <unattended flag>'
# The prompt is appended as the last argument.
set -euo pipefail
name="${1:?usage: tick.sh NAME}"
here="$(cd "$(dirname "$0")/../.." && pwd)"
cd "$here"

runner="$(sed -n 's/^runner: *//p' "relay/agents/$name.md" | head -1)"
var="RELAY_CMD_$(echo "$runner" | tr '-' '_')"
cmd="${!var:-}"
if [ -z "$cmd" ]; then
  echo "tick.sh: set $var to the non-interactive command for runner '$runner'" >&2
  exit 2
fi

prompt="You are \`$name\`. Read relay/prompts/$name.md and follow its \"Every later run\" block: perform one tick, then stop and report."

logdir="${RELAY_LOGDIR:-$HOME/relay-logs}"
mkdir -p "$logdir"
log="$logdir/$name-$(date -u +%Y%m%dT%H%M%SZ).log"
lock="/tmp/relay-tick-$name.lock"

exec 9>"$lock"
if ! flock -n 9; then
  echo "$(date -u +%FT%TZ) $name: previous tick still running, skipping" >>"$logdir/$name-skipped.log"
  exit 0
fi

{
  echo "== $name tick start $(date -u +%FT%TZ) runner=$runner"
  # Agents may leave uncommitted work in their own artifact directories
  # between ticks; --autostash carries it across the rebase. A failed pull is
  # logged but not fatal: the agent's own tick pulls again and can recover.
  git fetch -q origin main && git checkout -q main \
    && git pull -q --rebase --autostash origin main \
    || echo "== pre-tick pull failed; leaving it to the agent"
  # shellcheck disable=SC2086
  timeout "${RELAY_TIMEOUT:-50m}" $cmd "$prompt" || echo "== command exited with status $?"
  echo "== tick end $(date -u +%FT%TZ)"
  python3 relay/bin/relay.py status
} >>"$log" 2>&1
