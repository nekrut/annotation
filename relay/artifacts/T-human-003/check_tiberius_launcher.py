#!/usr/bin/env python3
"""Small isolated launcher check; does not install TensorFlow or predict genes."""
import datetime as dt
import json
import os
import pathlib
import platform
import subprocess
import sys
import tempfile
import time

OUT = pathlib.Path(__file__).resolve().parent
SHA = "e73844bfc7170665bc6632f1806a05dfe59fae8a"


def main():
    report = {"started_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
              "repository": "https://github.com/Gaius-Augustus/Tiberius",
              "revision": SHA, "python": sys.version, "platform": platform.platform(),
              "scope": "Fresh venv without system site packages: pip install . and --list_cfg only; no model weights, TensorFlow, containers, or inference.",
              "commands": []}
    with tempfile.TemporaryDirectory(prefix=".tiberius-check-", dir=OUT) as name:
        base = pathlib.Path(name)
        env = dict(os.environ)
        env.update(TMPDIR=name, PYTHONNOUSERSITE="1", PIP_NO_CACHE_DIR="1")
        env.pop("PYTHONPATH", None)
        env.pop("PYTHONHOME", None)

        def run(cmd, cwd=base, timeout=180):
            start = time.monotonic()
            r = subprocess.run(cmd, cwd=cwd, env=env, capture_output=True, text=True, timeout=timeout)
            report["commands"].append({"argv": [str(x).replace(name, "<TEMP>") for x in cmd],
                                       "cwd": str(cwd).replace(name, "<TEMP>"), "exit_code": r.returncode,
                                       "elapsed_seconds": round(time.monotonic() - start, 3),
                                       "stdout": r.stdout.replace(name, "<TEMP>"),
                                       "stderr": r.stderr.replace(name, "<TEMP>")})
            print(cmd[0], cmd[1], "exit", r.returncode, flush=True)
            if r.returncode:
                raise RuntimeError("Command failed; see recorded output")
            return r.stdout

        try:
            repo = base / "Tiberius"
            run(["git", "clone", "--depth", "1", "--no-tags", report["repository"], str(repo)])
            actual = run(["git", "rev-parse", "HEAD"], cwd=repo).strip()
            if actual != SHA:
                raise RuntimeError("Repository head changed; inspect and pin before rerunning")
            venv = base / "venv"
            run([sys.executable, "-m", "venv", str(venv)])
            py = str(venv / "bin" / "python")
            run([py, "-m", "pip", "--isolated", "--disable-pip-version-check", "install",
                 "--no-cache-dir", "--index-url", "https://pypi.org/simple", "."], cwd=repo)
            run([py, "tiberius.py", "--list_cfg"], cwd=repo, timeout=30)
            run([py, "-m", "pip", "--isolated", "--disable-pip-version-check", "freeze"])
            report["result"] = "PASS: launcher install and model-configuration listing; full inference untested"
        except (RuntimeError, subprocess.TimeoutExpired) as exc:
            report["result"] = "INCOMPLETE: " + str(exc)
    report["finished_utc"] = dt.datetime.now(dt.timezone.utc).isoformat()
    (OUT / "tiberius-launcher-check.json").write_text(json.dumps(report, indent=2) + "\n")
    print(report["result"], flush=True)


if __name__ == "__main__":
    main()
