Body of an `alert` to `human` asking for cluster time. Copy the headings;
fill in every one. `gagarin` runs exactly what is written here and nothing
more, so a request that leaves a field blank is sent back with a `question`.

## Task

T-xxx-nnn (the task whose deliverable needs this run).

## What to run

- Command or job script: (a script committed on your work branch, with the
  branch name and commit, or a one-line command).
- Container or environment: (Singularity image URL and digest, conda env
  file on the branch, or "none").
- Inputs: one line per file: URL or path, size, sha256.
- Outputs to return: paths, expected sizes. Only files under 5 MB come back
  through the relay; name a scoring or summary command for anything larger.
- Scoring or summary command: (what produces the numbers you need).

## Resources

- CPUs:
- Memory (GB):
- GPUs (0, 1, or 2; 24 GB each):
- Wall clock (HH:MM):
- Scratch storage (GB):
- Estimated total: CPU-hours, GPU-hours.

## Why the laptop budget is not enough

One or two sentences.

## Where the results go

The file in your deliverable that the numbers feed (for example
`docs/cost-baseline/measured.tsv`).
