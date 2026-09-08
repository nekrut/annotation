#!/usr/bin/env python3
"""Helper for the Markdown relay protocol (relay/PROTOCOL.md).

Standard library only. Python 3.8+.

    relay.py validate
    relay.py status
    relay.py inbox --for NAME [--mark]
    relay.py new --from NAME --to a,b --type TYPE --title "..." [--task T] [--thread ID] [--reply-to ID] [--body-file F]
    relay.py task new --by NAME --title "..." [--touches p1,p2] [--depends-on T1,T2]
    relay.py claim TASK --as NAME [--hours 2] [--push]
    relay.py release TASK --as NAME --status STATUS [--pr URL]
    relay.py heartbeat --as NAME
"""
import argparse
import datetime as dt
import os
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
AGENTS = ROOT / "agents"
TASKS = ROOT / "tasks"
MESSAGES = ROOT / "messages"
TEMPLATES = ROOT / "templates"

NAME_RE = re.compile(r"^[a-z][a-z0-9-]{0,23}$")
MSG_FILE_RE = re.compile(r"^(\d{8}T\d{6}Z)-([a-z][a-z0-9-]{0,23})-([a-z0-9]{4,})$")
TASK_ID_RE = re.compile(r"^T-([a-z][a-z0-9-]{0,23})-(\d{3,})$")
STAMP_RE = re.compile(r"^\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}Z$")

MSG_TYPES = {"note", "question", "answer", "proposal", "review", "handoff", "decision", "alert"}
TASK_STATUSES = {"open", "claimed", "in_progress", "review", "done", "blocked", "dropped"}
OWNED_STATUSES = {"claimed", "in_progress", "review"}
MSG_REQUIRED = ["id", "from", "to", "type", "title", "created"]
TASK_REQUIRED = ["id", "title", "status", "owner", "created_by", "created", "lease_until", "depends_on", "touches", "pr"]
AGENT_REQUIRED = ["name", "kind", "provider", "model", "runner", "operator", "capabilities", "last_seen", "last_heartbeat"]


# ---------------------------------------------------------------- front matter

def _scalar(v):
    v = v.strip()
    if v in ("", "null", "~"):
        return None
    if len(v) >= 2 and v[0] == v[-1] and v[0] in "\"'":
        return v[1:-1]
    return v


def parse_front_matter(text):
    """Return (dict, body). Supports `key: value`, `key: [a, b]`, and block lists."""
    lines = text.splitlines()
    if not lines or lines[0].strip() != "---":
        raise ValueError("missing front matter (file must start with ---)")
    fm, i, key = {}, 1, None
    while i < len(lines):
        line = lines[i]
        if line.strip() == "---":
            return fm, "\n".join(lines[i + 1:])
        if line.startswith("  - ") and key is not None and isinstance(fm.get(key), list):
            fm[key].append(_scalar(line[4:]))
        elif ":" in line and not line.startswith(" "):
            key, _, raw = line.partition(":")
            key, raw = key.strip(), raw.strip()
            if raw.startswith("[") and raw.endswith("]"):
                inner = raw[1:-1].strip()
                fm[key] = [_scalar(x) for x in inner.split(",")] if inner else []
            elif raw == "":
                fm[key] = []  # block list follows (or empty)
            else:
                fm[key] = _scalar(raw)
        else:
            raise ValueError("cannot parse front matter line %d: %r" % (i + 1, line))
        i += 1
    raise ValueError("unterminated front matter")


def set_field(text, key, value):
    """Rewrite one scalar `key:` line inside the front matter block."""
    rendered = "null" if value is None else str(value)
    head, sep, rest = text.partition("\n---")  # first line is '---'
    if not sep:
        raise ValueError("no front matter")
    pattern = re.compile(r"^%s:.*$" % re.escape(key), re.M)
    if not pattern.search(head):
        raise ValueError("field %r not present" % key)
    return pattern.sub("%s: %s" % (key, rendered), head, count=1) + sep + rest


def now():
    return dt.datetime.now(dt.timezone.utc).replace(microsecond=0)


def iso(t):
    return t.strftime("%Y-%m-%dT%H:%M:%SZ")


def stamp(t):
    return t.strftime("%Y%m%dT%H%M%SZ")


def parse_iso(s):
    return dt.datetime.strptime(s, "%Y-%m-%dT%H:%M:%SZ").replace(tzinfo=dt.timezone.utc)


def load_dir(directory, suffix=".md"):
    out = {}
    for p in sorted(directory.glob("*" + suffix)):
        if p.name.startswith("_"):
            continue
        fm, body = parse_front_matter(p.read_text(encoding="utf-8"))
        out[p.stem] = (fm, body, p)
    return out


# ---------------------------------------------------------------- validate

def validate(_args):
    errors, warnings = [], []
    err = lambda p, m: errors.append("%s: %s" % (p.relative_to(ROOT.parent), m))
    warn = lambda p, m: warnings.append("%s: %s" % (p.relative_to(ROOT.parent), m))

    def load(directory):
        out = {}
        for p in sorted(directory.glob("*.md")):
            try:
                out[p.stem] = parse_front_matter(p.read_text(encoding="utf-8")) + (p,)
            except ValueError as e:
                err(p, str(e))
        return out

    agents, tasks, messages = load(AGENTS), load(TASKS), load(MESSAGES)
    names = set()

    for stem, (fm, _, p) in agents.items():
        for k in AGENT_REQUIRED:
            if k not in fm:
                err(p, "missing field %r" % k)
        if not NAME_RE.match(stem):
            err(p, "bad agent name %r" % stem)
        if fm.get("name") != stem:
            err(p, "name %r does not match file name" % fm.get("name"))
        if fm.get("kind") not in ("agent", "human"):
            err(p, "kind must be agent or human")
        seen = fm.get("last_seen")
        if seen and not re.match(r"^[0-9a-f]{7,40}$", seen):
            err(p, "last_seen must be a git commit SHA (set it with `relay.py inbox --for %s --mark`)" % stem)
        elif seen and _git("cat-file", "-e", seen + "^{commit}", check=False).returncode != 0:
            warn(p, "last_seen %s is not in this clone's history (shallow clone?)" % seen)
        names.add(stem)
    if "human" not in names and not any(fm.get("kind") == "human" for fm, _, _ in agents.values()):
        warn(AGENTS, "no coordinator (kind: human) in roster")

    t = now()
    for stem, (fm, body, p) in tasks.items():
        for k in TASK_REQUIRED:
            if k not in fm:
                err(p, "missing field %r" % k)
        m = TASK_ID_RE.match(stem)
        if not m:
            err(p, "bad task id %r (expected T-<name>-<nnn>)" % stem)
        elif m.group(1) not in names:
            err(p, "creator %r in id is not in the roster" % m.group(1))
        if fm.get("id") != stem:
            err(p, "id %r does not match file name" % fm.get("id"))
        if fm.get("status") not in TASK_STATUSES:
            err(p, "bad status %r" % fm.get("status"))
        if fm.get("created_by") not in names:
            err(p, "created_by %r not in roster" % fm.get("created_by"))
        owner = fm.get("owner")
        if fm.get("status") in OWNED_STATUSES:
            if owner not in names:
                err(p, "status %s requires an owner from the roster" % fm["status"])
            lease = fm.get("lease_until")
            if not lease:
                err(p, "status %s requires lease_until" % fm["status"])
            elif not STAMP_RE.match(lease):
                err(p, "lease_until must look like 2026-09-08T18:00:00Z")
            elif parse_iso(lease) < t and fm["status"] != "review":
                warn(p, "lease expired at %s; task may be treated as open" % lease)
        elif owner is not None:
            err(p, "status %s must have owner: null" % fm.get("status"))
        for dep in fm.get("depends_on") or []:
            if dep not in tasks:
                err(p, "depends_on %r does not exist" % dep)
        for section in ("## Goal", "## Definition of done", "## Log"):
            if section not in body:
                warn(p, "missing section %r" % section)

    t = now()
    for stem, (fm, body, p) in messages.items():
        for k in MSG_REQUIRED:
            if k not in fm:
                err(p, "missing field %r" % k)
        m = MSG_FILE_RE.match(stem)
        if not m:
            err(p, "bad file name (expected <YYYYMMDDTHHMMSSZ>-<from>-<seq>.md)")
            continue
        if fm.get("id") != stem:
            err(p, "id %r does not match file name" % fm.get("id"))
        if fm.get("from") != m.group(2):
            err(p, "from %r does not match file name" % fm.get("from"))
        if fm.get("from") not in names:
            err(p, "sender %r not in roster" % fm.get("from"))
        created = fm.get("created")
        if not created or not STAMP_RE.match(created):
            err(p, "created must look like 2026-09-08T18:00:00Z")
        elif stamp(parse_iso(created)) != m.group(1):
            err(p, "created %s does not match file name timestamp" % created)
        elif parse_iso(created) > t + dt.timedelta(minutes=10):
            warn(p, "created %s is in the future; check the sender's clock" % created)
        to = fm.get("to")
        if not isinstance(to, list) or not to:
            err(p, "to must be a non-empty list")
        else:
            for name in to:
                if name != "all" and name not in names:
                    err(p, "recipient %r not in roster" % name)
        mtype = fm.get("type")
        if mtype not in MSG_TYPES:
            err(p, "bad type %r" % mtype)
        if mtype == "decision" and agents.get(fm.get("from"), ({},))[0].get("kind") != "human":
            err(p, "decision messages may only come from the coordinator (kind: human)")
        if mtype == "answer" and not fm.get("reply_to"):
            err(p, "answer requires reply_to")
        if mtype == "review" and not fm.get("task"):
            err(p, "review requires task")
        for ref in ("thread", "reply_to"):
            if fm.get(ref) and fm[ref] not in messages:
                err(p, "%s %r does not exist" % (ref, fm[ref]))
        if fm.get("task") and fm["task"] not in tasks:
            err(p, "task %r does not exist" % fm["task"])
        if len(body.split()) > 2000:
            warn(p, "body is over 2000 words; move content to artifacts/")

    for w in warnings:
        print("warning:", w)
    for e in errors:
        print("error:", e)
    print("%d agents, %d tasks, %d messages; %d errors, %d warnings"
          % (len(agents), len(tasks), len(messages), len(errors), len(warnings)))
    return 1 if errors else 0


# ---------------------------------------------------------------- status

def status(_args):
    agents, tasks, messages = load_dir(AGENTS), load_dir(TASKS), load_dir(MESSAGES)
    t = now()
    print("TASKS")
    order = ["in_progress", "claimed", "review", "open", "blocked", "done", "dropped"]
    for stem, (fm, _, _) in sorted(tasks.items(), key=lambda kv: (order.index(kv[1][0]["status"]) if kv[1][0]["status"] in order else 99, kv[0])):
        lease = fm.get("lease_until")
        flag = " (LEASE EXPIRED)" if lease and fm["status"] in OWNED_STATUSES and parse_iso(lease) < t else ""
        print("  %-16s %-12s %-10s %s%s" % (stem, fm["status"], fm.get("owner") or "-", fm["title"], flag))
    print("AGENTS")
    for stem, (fm, _, _) in agents.items():
        unread = len(_inbox_for(stem, fm, messages))
        print("  %-12s %-6s %-12s heartbeat=%-20s unread=%d" % (stem, fm.get("kind"), fm.get("runner"), fm.get("last_heartbeat") or "-", unread))
    print("MESSAGES: %d total, newest %s" % (len(messages), max(messages) if messages else "-"))
    return 0


# ---------------------------------------------------------------- inbox

def _arrived_since(sha):
    """Stems of message files added to history after commit `sha`, or None if git can't tell."""
    r = _git("log", "--format=", "--name-only", "--diff-filter=A", "%s..HEAD" % sha, "--", "relay/messages", check=False)
    if r.returncode != 0:
        return None
    return {Path(line).stem for line in r.stdout.splitlines() if line.strip()}


def _inbox_for(name, agent_fm, messages):
    last = agent_fm.get("last_seen")
    since = _arrived_since(last) if last else None
    if last and since is None:
        print("warning: last_seen %s is not in this clone's history; treating everything as unread" % last, file=sys.stderr)
    out = []
    for stem, (fm, body, p) in messages.items():
        if fm.get("from") == name:
            continue
        if since is not None and stem not in since:
            continue
        to = fm.get("to") or []
        if "all" in to or name in to:
            out.append((stem, fm, body, p))
    return out


def inbox(args):
    name = args.for_
    agent_path = AGENTS / (name + ".md")
    if not agent_path.exists():
        sys.exit("no roster file for %r" % name)
    agent_fm, _ = parse_front_matter(agent_path.read_text(encoding="utf-8"))
    messages = load_dir(MESSAGES)
    items = _inbox_for(name, agent_fm, messages)
    if not items:
        print("inbox for %s: empty" % name)
    for stem, fm, body, p in items:
        print("=" * 72)
        print("%s  [%s]  from %s  to %s" % (stem, fm["type"], fm["from"], ",".join(fm["to"])))
        print("title: %s" % fm["title"])
        for k in ("task", "thread", "reply_to"):
            if fm.get(k):
                print("%s: %s" % (k, fm[k]))
        print()
        print(body.strip())
    if args.mark:
        head = _git("rev-parse", "HEAD").stdout.strip()
        agent_path.write_text(set_field(agent_path.read_text(encoding="utf-8"), "last_seen", head), encoding="utf-8")
        print("=" * 72)
        print("last_seen -> %s" % head)
    return 0


# ---------------------------------------------------------------- new message

def _next_seq(name):
    n = 0
    for p in MESSAGES.glob("*-%s-*.md" % name):
        m = MSG_FILE_RE.match(p.stem)
        if m and m.group(2) == name and m.group(3).isdigit():
            n = max(n, int(m.group(3)))
    return n + 1


def new_message(args):
    name = args.from_
    if not (AGENTS / (name + ".md")).exists():
        sys.exit("no roster file for %r" % name)
    if args.type not in MSG_TYPES:
        sys.exit("type must be one of %s" % ", ".join(sorted(MSG_TYPES)))
    t = now()
    mid = "%s-%s-%04d" % (stamp(t), name, _next_seq(name))
    body = Path(args.body_file).read_text(encoding="utf-8") if args.body_file else (sys.stdin.read() if not sys.stdin.isatty() else "")
    to = [x.strip() for x in args.to.split(",") if x.strip()]
    fm = [
        "---",
        "id: " + mid,
        "from: " + name,
        "to: [" + ", ".join(to) + "]",
        "type: " + args.type,
        "title: " + args.title,
        "task: " + (args.task or "null"),
        "thread: " + (args.thread or "null"),
        "reply_to: " + (args.reply_to or "null"),
        "created: " + iso(t),
        "---",
        "",
    ]
    path = MESSAGES / (mid + ".md")
    path.write_text("\n".join(fm) + body.rstrip() + "\n", encoding="utf-8")
    print(path.relative_to(ROOT.parent))
    return 0


# ---------------------------------------------------------------- tasks

def task_new(args):
    name = args.by
    if not (AGENTS / (name + ".md")).exists():
        sys.exit("no roster file for %r" % name)
    n = 0
    for p in TASKS.glob("T-%s-*.md" % name):
        m = TASK_ID_RE.match(p.stem)
        if m and m.group(1) == name:
            n = max(n, int(m.group(2)))
    tid = "T-%s-%03d" % (name, n + 1)
    text = (TEMPLATES / "task.md").read_text(encoding="utf-8")
    text = text.replace("T-NAME-000", tid).replace("TITLE", args.title).replace("NAME", name)
    text = set_field(text, "created", iso(now()))
    if args.touches:
        text = set_field(text, "touches", "[" + ", ".join(x.strip() for x in args.touches.split(",")) + "]")
    if args.depends_on:
        text = set_field(text, "depends_on", "[" + ", ".join(x.strip() for x in args.depends_on.split(",")) + "]")
    path = TASKS / (tid + ".md")
    path.write_text(text, encoding="utf-8")
    print(path.relative_to(ROOT.parent))
    return 0


def _git(*cmd, check=True):
    r = subprocess.run(["git", *cmd], cwd=ROOT.parent, text=True, capture_output=True)
    if check and r.returncode != 0:
        sys.exit("git %s failed:\n%s" % (" ".join(cmd), r.stderr.strip()))
    return r


def _log_line(text, line):
    return text.rstrip("\n") + "\n- %s %s\n" % (now().strftime("%Y-%m-%d"), line)


def claim(args):
    name, tid = args.as_, args.task
    path = TASKS / (tid + ".md")
    if not path.exists():
        sys.exit("no such task %r" % tid)
    if not (AGENTS / (name + ".md")).exists():
        sys.exit("no roster file for %r" % name)

    def attempt():
        text = path.read_text(encoding="utf-8")
        fm, _ = parse_front_matter(text)
        t = now()
        expired = fm["status"] in OWNED_STATUSES and fm.get("lease_until") and parse_iso(fm["lease_until"]) < t
        if fm["status"] != "open" and not expired:
            return None, "task %s is %s (owner %s); not claimable" % (tid, fm["status"], fm.get("owner"))
        text = set_field(text, "status", "claimed")
        text = set_field(text, "owner", name)
        text = set_field(text, "lease_until", iso(t + dt.timedelta(hours=args.hours)))
        note = "%s: claimed%s." % (name, " (previous lease by %s expired)" % fm.get("owner") if expired else "")
        path.write_text(_log_line(text, note), encoding="utf-8")
        return True, "claimed %s as %s until %s" % (tid, name, iso(t + dt.timedelta(hours=args.hours)))

    if not args.push:
        ok, msg = attempt()
        print(msg)
        return 0 if ok else 1

    if _git("status", "--porcelain", "--untracked-files=no").stdout.strip():
        sys.exit("claim --push needs a clean working tree (commit or stash first)")
    _git("fetch", "origin", "main")
    if _git("rev-list", "origin/main..HEAD").stdout.strip():
        sys.exit("claim --push needs no unpushed commits (push them first)")
    for round_ in range(4):
        _git("pull", "--rebase", "origin", "main")
        ok, msg = attempt()
        if not ok:
            print(msg)
            return 1
        _git("add", str(path))
        _git("commit", "-m", "relay: %s claims %s" % (name, tid))
        pushed = _git("push", "origin", "main", check=False)
        if pushed.returncode == 0:
            print(msg)
            return 0
        print("push rejected; rebasing (round %d)" % (round_ + 1))
        rebased = _git("pull", "--rebase", "origin", "main", check=False)
        if rebased.returncode != 0:
            _git("rebase", "--abort", check=False)
            _git("reset", "--hard", "origin/main")
            print("lost the race for %s; another agent claimed it" % tid)
            return 1
        # rebase applied cleanly: drop our commit and re-check the file's current state
        _git("reset", "--hard", "origin/main")
    print("gave up after 4 rounds")
    return 1


def release(args):
    name, tid = args.as_, args.task
    path = TASKS / (tid + ".md")
    if not path.exists():
        sys.exit("no such task %r" % tid)
    text = path.read_text(encoding="utf-8")
    fm, _ = parse_front_matter(text)
    if fm.get("owner") != name and name != "human":
        sys.exit("%s does not own %s (owner is %s)" % (name, tid, fm.get("owner")))
    if args.status not in TASK_STATUSES:
        sys.exit("status must be one of %s" % ", ".join(sorted(TASK_STATUSES)))
    text = set_field(text, "status", args.status)
    if args.status in OWNED_STATUSES:
        text = set_field(text, "lease_until", iso(now() + dt.timedelta(hours=args.hours)))
    else:
        text = set_field(text, "owner", None)
        text = set_field(text, "lease_until", None)
    if args.pr:
        text = set_field(text, "pr", args.pr)
    path.write_text(_log_line(text, "%s: -> %s%s." % (name, args.status, " (%s)" % args.pr if args.pr else "")), encoding="utf-8")
    print("%s -> %s" % (tid, args.status))
    return 0


def heartbeat(args):
    path = AGENTS / (args.as_ + ".md")
    if not path.exists():
        sys.exit("no roster file for %r" % args.as_)
    path.write_text(set_field(path.read_text(encoding="utf-8"), "last_heartbeat", iso(now())), encoding="utf-8")
    print("heartbeat %s" % args.as_)
    return 0


# ---------------------------------------------------------------- cli

def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("validate").set_defaults(fn=validate)
    sub.add_parser("status").set_defaults(fn=status)

    p = sub.add_parser("inbox")
    p.add_argument("--for", dest="for_", required=True)
    p.add_argument("--mark", action="store_true", help="advance last_seen to the newest listed message")
    p.set_defaults(fn=inbox)

    p = sub.add_parser("new", help="create a message (body from --body-file or stdin)")
    p.add_argument("--from", dest="from_", required=True)
    p.add_argument("--to", required=True, help="comma-separated names, or all")
    p.add_argument("--type", required=True)
    p.add_argument("--title", required=True)
    p.add_argument("--task")
    p.add_argument("--thread")
    p.add_argument("--reply-to", dest="reply_to")
    p.add_argument("--body-file", dest="body_file")
    p.set_defaults(fn=new_message)

    p = sub.add_parser("task")
    tsub = p.add_subparsers(dest="tcmd", required=True)
    q = tsub.add_parser("new")
    q.add_argument("--by", required=True)
    q.add_argument("--title", required=True)
    q.add_argument("--touches")
    q.add_argument("--depends-on", dest="depends_on")
    q.set_defaults(fn=task_new)

    p = sub.add_parser("claim")
    p.add_argument("task")
    p.add_argument("--as", dest="as_", required=True)
    p.add_argument("--hours", type=float, default=2)
    p.add_argument("--push", action="store_true", help="commit and push, retrying on rejection")
    p.set_defaults(fn=claim)

    p = sub.add_parser("release")
    p.add_argument("task")
    p.add_argument("--as", dest="as_", required=True)
    p.add_argument("--status", required=True)
    p.add_argument("--pr")
    p.add_argument("--hours", type=float, default=2)
    p.set_defaults(fn=release)

    p = sub.add_parser("heartbeat")
    p.add_argument("--as", dest="as_", required=True)
    p.set_defaults(fn=heartbeat)

    args = ap.parse_args(argv)
    return args.fn(args)


if __name__ == "__main__":
    sys.exit(main())
