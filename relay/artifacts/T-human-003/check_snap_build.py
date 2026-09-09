#!/usr/bin/env python3
"""Check pinned SNAP README commands in a fresh temporary source/build directory.

Uses the host's existing GNU build tools, not a fresh OS/container. Installs no
system packages. Records command results and hashes, then removes the checkout.
Python 3.11 standard library only; no genome-scale accuracy benchmark is run.
"""
import datetime as dt
import hashlib
import json
from pathlib import Path
import platform
import subprocess
import tempfile
import time

OUT = Path(__file__).resolve().parent
REPOSITORY = "https://github.com/KorfLab/SNAP"
REVISION = "4ad1e957cd8e68b63857cc1cb3380d39a7b518b1"


def main():
    report = {
        "repository": REPOSITORY, "revision": REVISION,
        "checked_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "environment": platform.platform(),
        "isolation": "Fresh temporary source/build directory; host compiler and libc; no system installation",
        "readme": f"{REPOSITORY}/blob/{REVISION}/README.md",
        "commands": [], "input_hashes": {},
    }
    with tempfile.TemporaryDirectory(prefix="snap-check-", dir=OUT) as scratch:
        root = Path(scratch)
        env = {"PATH": "/usr/bin:/bin", "LC_ALL": "C"}

        def run(argv, cwd=root, timeout=120):
            started = time.monotonic()
            try:
                p = subprocess.run(argv, cwd=cwd, env=env, capture_output=True,
                                   timeout=timeout)
                code, stdout, stderr = p.returncode, p.stdout, p.stderr
            except subprocess.TimeoutExpired as exc:
                code, stdout, stderr = "timeout", exc.stdout or b"", exc.stderr or b""
            entry = {
                "argv": argv, "exit_code": code,
                "wall_seconds": round(time.monotonic() - started, 4),
                "stdout_bytes": len(stdout),
                "stdout_sha256": hashlib.sha256(stdout).hexdigest(),
                "stdout_excerpt": stdout.decode(errors="replace")[:2000].replace(scratch, "<temporary>"),
                "stderr_excerpt": stderr.decode(errors="replace")[-4000:].replace(scratch, "<temporary>"),
            }
            report["commands"].append(entry)
            return code, stdout

        for argv in [["gcc", "--version"], ["make", "--version"], ["ldd", "--version"]]:
            run(argv)
        code, _ = run(["git", "clone", "--quiet", "--no-checkout", "--filter=blob:none", REPOSITORY + ".git", "src"])
        src = root / "src"
        if code == 0:
            code, _ = run(["git", "checkout", "--quiet", REVISION], cwd=src)
        if code == 0:
            code, sha = run(["git", "rev-parse", "HEAD"], cwd=src)
            assert sha.decode().strip() == REVISION
            for name in ["README.md", "Makefile", "DNA/thale.dna.gz", "DNA/worm.dna.gz"]:
                path = src / name
                if path.exists():
                    report["input_hashes"][name] = hashlib.sha256(path.read_bytes()).hexdigest()
            code, _ = run(["make"], cwd=src)
            report["build_exit_code"] = code
            if code == 0:
                outcomes = []
                for model, dna in [("thale", "thale"), ("worm", "worm")]:
                    code, output = run(["./snap", f"HMM/{model}", f"DNA/{dna}.dna.gz"], cwd=src, timeout=60)
                    outcomes.append(code == 0 and bool(output.strip()))
                report["readme_examples_passed"] = all(outcomes)
        report["result"] = (
            "PASS: fresh source build and both README example predictions; genome accuracy untested"
            if report.get("readme_examples_passed") else
            "FAIL: see command results; no source or command repair applied"
        )
    report["temporary_checkout_removed"] = True
    (OUT / "snap-build-check.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({k: report[k] for k in ["result", "checked_utc", "temporary_checkout_removed"]}))


if __name__ == "__main__":
    main()
