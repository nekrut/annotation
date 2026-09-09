#!/usr/bin/env python3
"""Refresh public GitHub metadata; standard library, compatible with Python 3.11.

Does not install or execute reviewed software. GitHub's commits endpoint counts
commits reachable from the default branch, including merges, by committer date.
"""
import csv
import argparse
import datetime as dt
import hashlib
import json
import pathlib
import re
import time
import urllib.error
import urllib.parse
import urllib.request

OUT = pathlib.Path(__file__).resolve().parent
REPOS = [
    "Gaius-Augustus/Tiberius", "ncbi/egapx", "usadellab/Helixer",
    "Gaius-Augustus/Augustus", "Gaius-Augustus/BRAKER",
    "gatech-genemark/GeneMark-ETP", "KorfLab/SNAP", "Yandell-Lab/maker",
    "Gaius-Augustus/clamsa", "mlin/PhyloCSF", "nekrut/axomeme",
    "nekrut/scalingPaper", "hexagonbio/geneML",
]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", action="append", help="Repository to inspect; repeatable")
    parser.add_argument("--append", action="store_true", help="Retain previous repositories and request provenance")
    args = parser.parse_args()
    now = dt.datetime.now(dt.timezone.utc).replace(microsecond=0)
    since = now.replace(year=now.year - 1)
    iso = lambda x: x.isoformat().replace("+00:00", "Z")
    audit = {"accessed_utc": iso(now), "since": iso(since), "until": iso(now),
             "commit_definition": "Default-branch reachable commits including merges; GitHub since/until filter.",
             "requests": [], "repositories": []}
    previous = None
    if args.append and (OUT / "repo-snapshot.json").exists():
        previous = json.loads((OUT / "repo-snapshot.json").read_text())
        audit["requests"] = previous["requests"]
        audit["initial_accessed_utc"] = previous.get("initial_accessed_utc", previous["accessed_utc"])

    def get(path):
        url = "https://api.github.com/" + path
        req = urllib.request.Request(url, headers={"User-Agent": "engels-relay-literature-review",
                                                  "Accept": "application/vnd.github+json"})
        with urllib.request.urlopen(req, timeout=30) as resp:
            raw = resp.read()
            audit["requests"].append({"url": url, "sha256": hashlib.sha256(raw).hexdigest(),
                                      "accessed_utc": iso(dt.datetime.now(dt.timezone.utc)),
                                      "etag": resp.headers.get("ETag"),
                                      "rate_limit_remaining": resp.headers.get("X-RateLimit-Remaining")})
            data, links = json.loads(raw), resp.headers.get("Link", "")
        time.sleep(0.3)
        return data, links

    rows = list(previous["repositories"]) if previous else []
    for repo in (args.repo or REPOS):
        row = {"url": "https://github.com/" + repo, "readme_install_fresh_env": "not attempted",
               "notes": "Metadata only; installation and inference untested.", "snapshot_utc": iso(now)}
        try:
            meta, _ = get("repos/" + repo)
            row.update(url=meta["html_url"], default_branch=meta["default_branch"],
                       stars=meta["stargazers_count"], language=meta.get("language") or "unknown",
                       licence=(meta.get("license") or {}).get("spdx_id", "unknown"),
                       archived=meta["archived"])
            # GitHub open_issues_count includes PRs. Subtract all open PRs.
            prs, page = [], 1
            while True:
                part, links = get(f"repos/{repo}/pulls?state=open&per_page=100&page={page}")
                prs.extend(part)
                if 'rel="next"' not in links:
                    break
                page += 1
            row["open_issues"] = meta["open_issues_count"] - len(prs)
            params = urllib.parse.urlencode({"since": iso(since), "until": iso(now), "per_page": 1})
            commits, links = get(f"repos/{repo}/commits?{params}")
            last = re.search(r'<([^>]+)>; rel="last"', links)
            count = int(urllib.parse.parse_qs(urllib.parse.urlparse(last.group(1)).query)["page"][0]) if last else len(commits)
            row["commits_last_12_months"] = count
            if not commits:
                commits, _ = get(f"repos/{repo}/commits?per_page=1")
            row["last_commit_date"] = commits[0]["commit"]["committer"]["date"]
            row["default_branch_sha"] = commits[0]["sha"]
            row["count_since_utc"], row["count_until_utc"] = iso(since), iso(now)
        except (urllib.error.HTTPError, urllib.error.URLError, ValueError, KeyError) as exc:
            row["notes"] = "Metadata incomplete: " + str(exc)
        rows = [old for old in rows if old["url"].lower() != row["url"].lower()]
        rows.append(row)
        audit["repositories"] = rows
        (OUT / "repo-snapshot.json").write_text(json.dumps(audit, indent=2) + "\n")
        print(json.dumps(row), flush=True)
    check = OUT / "tiberius-launcher-check.json"
    if check.exists():
        result = json.loads(check.read_text())
        for row in rows:
            if row["url"] == result["repository"] and row.get("default_branch_sha") == result["revision"]:
                row["readme_install_fresh_env"] = result["result"]
                row["notes"] = "Fresh isolated venv; launcher only. See tiberius-launcher-check.json."
    annotations = OUT / "repo-audit.json"
    if annotations.exists():
        for annotation in json.loads(annotations.read_text()):
            for row in rows:
                if row["url"] == annotation["url"] and row.get("default_branch_sha") == annotation["revision"]:
                    row.setdefault("github_detected_licence", row.get("licence", "unknown"))
                    row["licence"] = annotation["licence"]
                    row["notes"] = row["notes"].split(" License audit:")[0] + " License audit: " + annotation["source"]
    audit["repositories"] = rows
    (OUT / "repo-snapshot.json").write_text(json.dumps(audit, indent=2) + "\n")
    fields = ["url", "last_commit_date", "commits_last_12_months", "open_issues", "stars", "language",
              "licence", "readme_install_fresh_env", "notes", "snapshot_utc", "default_branch",
              "default_branch_sha", "count_since_utc", "count_until_utc", "archived", "github_detected_licence"]
    with (OUT / "repos.tsv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
