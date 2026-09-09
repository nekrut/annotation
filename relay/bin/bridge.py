#!/usr/bin/env python3
"""GitHub Issues bridge for the coordinator.

Run by .github/workflows/relay-bridge.yml. Standard library only.

  bridge.py sync-inbox       on push to main: open one issue per message that
                             needs the coordinator (question/alert/proposal
                             addressed to human or all) and per task in
                             `review`, if no issue exists yet.
  bridge.py handle-comment   on issue_comment: turn a collaborator's comment on
                             a bridge issue into a relay message or a task
                             release, commit it to main as `human`, and reply.

Markers in issue bodies tie issues to relay objects:
  <!-- relay-msg: <message id> -->     <!-- relay-task: <task id> -->
"""
import datetime as dt
import json
import os
import re
import subprocess
import sys
import urllib.error
import urllib.request
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import relay  # noqa: E402

REPO = os.environ.get("GITHUB_REPOSITORY", "")
TOKEN = os.environ.get("GITHUB_TOKEN", "")
API = "https://api.github.com"
LABEL = "relay-inbox"
COORDINATOR = "human"
NEEDS_HUMAN = {"question", "alert", "proposal"}
MARK_MSG = re.compile(r"<!-- relay-msg: (\S+) -->")
MARK_TASK = re.compile(r"<!-- relay-task: (\S+) -->")


# ---------------------------------------------------------------- github api

def gh(method, path, body=None):
    req = urllib.request.Request(API + path, method=method)
    req.add_header("Authorization", "Bearer " + TOKEN)
    req.add_header("Accept", "application/vnd.github+json")
    data = None
    if body is not None:
        data = json.dumps(body).encode()
        req.add_header("Content-Type", "application/json")
    try:
        with urllib.request.urlopen(req, data) as r:
            return json.loads(r.read() or b"null")
    except urllib.error.HTTPError as e:
        detail = e.read().decode(errors="replace")
        raise SystemExit("GitHub API %s %s -> %s %s" % (method, path, e.code, detail[:300]))


def ensure_label():
    try:
        gh("POST", "/repos/%s/labels" % REPO, {"name": LABEL, "color": "0e8a16",
                                                "description": "relay: needs the coordinator"})
    except SystemExit as e:
        if "422" not in str(e):
            raise


def open_bridge_issues():
    """Return {marker: issue} for all open issues carrying the label."""
    out = {}
    page = 1
    while True:
        items = gh("GET", "/repos/%s/issues?state=all&labels=%s&per_page=100&page=%d" % (REPO, LABEL, page))
        if not items:
            return out
        for it in items:
            if "pull_request" in it:
                continue
            body = it.get("body") or ""
            for rx in (MARK_MSG, MARK_TASK):
                m = rx.search(body)
                if m:
                    out[m.group(1)] = it
        page += 1


# ---------------------------------------------------------------- sync-inbox

def sync_inbox():
    ensure_label()
    existing = open_bridge_issues()
    messages = relay.load_dir(relay.MESSAGES)
    tasks = relay.load_dir(relay.TASKS)
    created = 0

    for stem, (fm, body, _) in messages.items():
        if fm.get("from") == COORDINATOR or fm.get("type") not in NEEDS_HUMAN:
            continue
        to = fm.get("to") or []
        if COORDINATOR not in to and "all" not in to:
            continue
        if stem in existing:
            continue
        kind = fm["type"]
        title = "[relay] %s from %s: %s" % (kind, fm["from"], fm["title"])
        reply_kind = "answer" if kind == "question" else "decision"
        text = (
            "<!-- relay-msg: %s -->\n"
            "**%s** from `%s` to `%s`%s\n\n"
            "Message: `relay/messages/%s.md`\n\n---\n\n%s\n\n---\n"
            "**To reply:** comment on this issue. Your comment is committed to the relay as a "
            "`%s` from `human` addressed to `%s`, and this issue is closed. "
            "Only collaborators' comments are used."
        ) % (stem, kind, fm["from"], ", ".join(to),
             (" about `%s`" % fm["task"]) if fm.get("task") else "",
             stem, body.strip(), reply_kind, fm["from"])
        gh("POST", "/repos/%s/issues" % REPO, {"title": title[:250], "body": text, "labels": [LABEL]})
        created += 1

    for tid, (fm, body, _) in tasks.items():
        if fm.get("status") != "review" or tid in existing:
            continue
        title = "[relay] review requested: %s (%s)" % (tid, fm["title"])
        arts = relay.ROOT / "artifacts" / tid
        art_list = ""
        if arts.is_dir():
            art_list = "\n".join("- `relay/artifacts/%s/%s`" % (tid, p.name) for p in sorted(arts.iterdir()))
        text = (
            "<!-- relay-task: %s -->\n"
            "Task `%s` was moved to `review` by `%s`.\n\n"
            "Task file: `relay/tasks/%s.md`%s\n\nArtifacts:\n%s\n\n---\n"
            "**To accept:** comment `done`. The task is marked done and this issue closes.\n"
            "**To send it back:** comment anything else. Your comment is posted as a "
            "`review` message to `%s` and the task returns to `in_progress`."
        ) % (tid, tid, fm.get("owner"), tid,
             ("\nPull request: %s" % fm["pr"]) if fm.get("pr") else "",
             art_list or "- (none under relay/artifacts)", fm.get("owner"))
        gh("POST", "/repos/%s/issues" % REPO, {"title": title[:250], "body": text, "labels": [LABEL]})
        created += 1

    print("sync-inbox: %d issue(s) created" % created)


# ---------------------------------------------------------------- comments

def write_message(to, mtype, title, body, reply_to=None, task=None):
    t = relay.now()
    mid = "%s-%s-%04d" % (relay.stamp(t), COORDINATOR, relay._next_seq(COORDINATOR))
    fm = ["---", "id: " + mid, "from: " + COORDINATOR, "to: [" + to + "]", "type: " + mtype,
          "title: " + title.replace("\n", " ")[:200], "task: " + (task or "null"), "thread: null",
          "reply_to: " + (reply_to or "null"), "created: " + relay.iso(t), "---", ""]
    path = relay.MESSAGES / (mid + ".md")
    path.write_text("\n".join(fm) + body.strip() + "\n", encoding="utf-8")
    return mid


def git(*cmd, check=True):
    r = subprocess.run(["git", *cmd], cwd=relay.ROOT.parent, text=True, capture_output=True)
    if check and r.returncode != 0:
        raise SystemExit("git %s failed: %s" % (" ".join(cmd), r.stderr.strip()))
    return r


def commit_and_push(message):
    git("config", "user.name", "human (via GitHub)")
    git("config", "user.email", "human@relay.local")
    git("add", "relay")
    git("commit", "-m", message)
    for attempt in range(5):
        if git("push", "origin", "HEAD:main", check=False).returncode == 0:
            return
        git("pull", "--rebase", "origin", "main")
    raise SystemExit("could not push after 5 attempts")


def handle_comment():
    event = json.loads(Path(os.environ["GITHUB_EVENT_PATH"]).read_text())
    issue = event["issue"]
    comment = event["comment"]
    if "pull_request" in issue:
        return
    if LABEL not in {l["name"] for l in issue.get("labels", [])}:
        return
    if comment["user"]["type"] == "Bot":
        return
    if comment.get("author_association") not in ("OWNER", "MEMBER", "COLLABORATOR"):
        print("ignoring comment from non-collaborator")
        return
    text = (comment.get("body") or "").strip()
    if not text:
        return
    number = issue["number"]
    ibody = issue.get("body") or ""
    if relay.validate(None) != 0:
        raise SystemExit("relay does not validate; refusing to write")

    m = MARK_MSG.search(ibody)
    if m:
        mid = m.group(1)
        messages = relay.load_dir(relay.MESSAGES)
        if mid not in messages:
            raise SystemExit("message %s not found" % mid)
        fm, _, _ = messages[mid]
        mtype = "answer" if fm["type"] == "question" else "decision"
        new_id = write_message(fm["from"], mtype, "Re: " + fm["title"], text,
                               reply_to=mid, task=fm.get("task"))
        if relay.validate(None) != 0:
            raise SystemExit("new message failed validation")
        commit_and_push("relay: human replies to %s via GitHub issue #%d" % (fm["from"], number))
        gh("POST", "/repos/%s/issues/%d/comments" % (REPO, number),
           {"body": "Posted to the relay as `%s` (`relay/messages/%s.md`). `%s` sees it on its next tick." % (mtype, new_id, fm["from"])})
        gh("PATCH", "/repos/%s/issues/%d" % (REPO, number), {"state": "closed"})
        return

    m = MARK_TASK.search(ibody)
    if m:
        tid = m.group(1)
        tasks = relay.load_dir(relay.TASKS)
        if tid not in tasks:
            raise SystemExit("task %s not found" % tid)
        fm, _, path = tasks[tid]
        owner = fm.get("owner")
        if text.lower().split()[0].rstrip(".!") == "done":
            t = path.read_text(encoding="utf-8")
            t = relay.set_field(t, "status", "done")
            t = relay.set_field(t, "owner", None)
            t = relay.set_field(t, "lease_until", None)
            path.write_text(relay._log_line(t, "human: -> done (accepted via GitHub issue #%d)." % number), encoding="utf-8")
            if relay.validate(None) != 0:
                raise SystemExit("task failed validation")
            commit_and_push("relay: human accepts %s via GitHub issue #%d" % (tid, number))
            gh("POST", "/repos/%s/issues/%d/comments" % (REPO, number), {"body": "`%s` marked done." % tid})
            gh("PATCH", "/repos/%s/issues/%d" % (REPO, number), {"state": "closed"})
        else:
            new_id = write_message(owner or "all", "review", "Review of %s" % tid, text, task=tid)
            t = path.read_text(encoding="utf-8")
            if fm.get("status") == "review":
                t = relay.set_field(t, "status", "in_progress")
                t = relay.set_field(t, "lease_until", relay.iso(relay.now() + dt.timedelta(hours=2)))
            path.write_text(relay._log_line(t, "human: review posted (%s); back to in_progress." % new_id), encoding="utf-8")
            if relay.validate(None) != 0:
                raise SystemExit("review failed validation")
            commit_and_push("relay: human reviews %s via GitHub issue #%d" % (tid, number))
            gh("POST", "/repos/%s/issues/%d/comments" % (REPO, number),
               {"body": "Posted as a `review` message (`relay/messages/%s.md`); `%s` is back in `in_progress`. This issue stays open; comment `done` when the next revision is acceptable." % (new_id, owner)})
        return
    print("issue #%d carries no relay marker; nothing to do" % number)


if __name__ == "__main__":
    cmd = sys.argv[1] if len(sys.argv) > 1 else ""
    if not REPO or not TOKEN:
        raise SystemExit("GITHUB_REPOSITORY and GITHUB_TOKEN are required")
    {"sync-inbox": sync_inbox, "handle-comment": handle_comment}.get(cmd, lambda: sys.exit(__doc__))()
