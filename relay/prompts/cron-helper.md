# Instructions for a local helper session that sets up cron

Paste everything below the line into a Claude Code session opened in any of
the agent checkouts on the machine that runs the local agents. This session
is a helper for the human, not a relay agent: it must not create or edit
anything under `relay/`, must not claim tasks, and must not push.

---

You are helping the coordinator run four local relay agents unattended via
cron. Read `relay/STARTING.md` (section "Running unattended") and
`relay/bin/tick.sh` first. Do not modify anything under `relay/` and do not
push to the repository. Then do the following, confirming each step's result
before moving on.

1. **Find the checkouts.** Locate the four agent checkouts on this machine
   (directories named `lenin`, `engels`, `stalin`, `trotsky`, each a clone of
   nekrut/annotation). Confirm each is on `main`, has a clean tree, and can
   `git pull --rebase origin main` and `git push --dry-run origin main`
   without prompting for credentials. If a push would prompt, stop and tell
   the user how to fix credentials for that checkout before continuing.

2. **Identify each agent's CLI.** Read `runner:` from each checkout's
   `relay/agents/<name>.md`. For each runner, find the executable on PATH and
   read its `--help` to find (a) the flag for non-interactive or "print"
   mode that takes a prompt as an argument and exits when done, and (b) the
   flag that lets it run without stopping for approval. Quote the exact
   help text you based each choice on. Do not guess flags.

3. **Assemble the commands.** For each runner, write the full command that
   `relay/bin/tick.sh` expects in `RELAY_CMD_<runner with hyphens as
   underscores>`, of the form `<cli> <non-interactive flag> <no-approval
   flag>`; the script appends the prompt as the last argument. Show all of
   them to the user and get an explicit yes before using them, since the
   no-approval flags let the agents run commands unattended.

4. **Test one tick by hand.** With the variables exported, run
   `relay/bin/tick.sh lenin` in the lenin checkout and wait for it. Read the
   log in `~/relay-logs/`. It should show the CLI running, a commit and push
   under `relay/`, and `relay.py status` at the end. If the CLI stopped at a
   prompt, the flag is wrong: go back to step 2. Repeat for each agent once.

5. **Install the crontab.** Add these lines to the user's crontab, with the
   real absolute paths and the variables from step 3, keeping the staggered
   minutes so agents never push at once (marx runs in the cloud at :07):

   ```
   RELAY_CMD_...=...           # one line per runner
   RELAY_LOGDIR=$HOME/relay-logs
   PATH=<the PATH that contains every CLI, python3, git, flock, timeout>
   5  * * * * /abs/path/lenin/relay/bin/tick.sh   lenin
   20 * * * * /abs/path/engels/relay/bin/tick.sh  engels
   35 * * * * /abs/path/stalin/relay/bin/tick.sh  stalin
   50 * * * * /abs/path/trotsky/relay/bin/tick.sh trotsky
   ```

   Show the final crontab with `crontab -l`. On macOS, if `flock` or
   `timeout` is missing, install `coreutils` and `flock` via Homebrew, or
   adjust the script's calls to `gtimeout`, and tell the user what you did.

6. **Keep the machine awake and cron alive.** On macOS, explain that the
   machine must not sleep (System Settings, or a `caffeinate -i` process)
   and that cron needs Full Disk Access for the terminal or cron binary if
   the checkouts are under a protected folder. On Linux, confirm the cron
   service is running. Any environment the CLIs need (API keys, config
   files in the home directory) must be reachable from cron's minimal
   environment; check by running one tick from `env -i` with only the
   crontab variables set.

7. **Verify after the first scheduled hour.** Tell the user to run
   `python3 relay/bin/relay.py status` after the next full hour and confirm
   every local agent's heartbeat has advanced, and to look in
   `~/relay-logs/` for any `<name>-skipped.log` (overlap) or a log ending
   without a push.

Report at the end: the crontab installed, the exact commands per runner,
which ticks you tested by hand and their outcomes, and anything you could
not verify.
