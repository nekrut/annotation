# Shared task charter (EXAMPLE)

This is a filled-in example of `relay/TASK.md`. The project here is invented
to show the level of detail that works. Replace every part with your own.

## Goal

Produce a consistent structural annotation for the 12 bacterial genome
assemblies in `data/assemblies/`: gene models in GFF3, one file per genome,
plus a single summary table comparing gene counts, coding density, and
rRNA/tRNA counts across all 12. The result should be reproducible from a
single command and pass the checks in `src/validate.py`. Finished means: all
12 GFF3 files and the summary table exist in `results/`, the validator
passes, and the README describes how to rerun everything.

## Scope

- In scope: running and comparing annotation tools, writing the pipeline
  glue, writing the validator, writing the summary and the README.
- Out of scope: functional annotation (no GO terms, no pathway assignment),
  improving the assemblies themselves, anything requiring a GPU.

## Deliverables

1. `results/gff/<assembly>.gff3` for each of the 12 assemblies, GFF3 1.26,
   sorted, with `gene`, `mRNA`/`CDS`, `tRNA`, and `rRNA` features.
2. `results/summary.tsv` with one row per assembly and the columns listed in
   `docs/schema.md`.
3. `src/` containing the pipeline, runnable as `make all` from a clean clone.
4. `src/validate.py` that exits non-zero if any deliverable is missing or
   malformed. This is what "done" means for the whole project.
5. `README.md` explaining how to rerun and how long it takes.

## Constraints

- Python 3.11, standard library plus what is listed in `requirements.txt`.
  Adding a dependency needs a `proposal` and a `decision`.
- Annotation tools must be available via conda; record exact versions in
  `environment.yml`.
- Input data under `data/` is read-only. Never modify or delete it.
- Runs must complete on a laptop with 16 GB of RAM. If a step needs more,
  post an `alert`.
- No network calls at run time except downloading tool databases in a
  clearly separated `make databases` step.

## Priorities

1. A working end-to-end pipeline on one assembly, even if crude.
2. The validator, so everyone agrees on what correct output looks like.
3. All 12 assemblies.
4. The summary table.
5. README and cleanup.

When choosing between open tasks, prefer the one that unblocks the most
other tasks, then the one highest on this list.

## Review and acceptance

A task moves to `done` when its pull request is merged by the coordinator,
after at least one `review` message from an agent other than the owner.
The coordinator may waive the review for documentation-only tasks.

## How to split the work

Anyone may create tasks. A good task is finishable in one to three ticks,
names the paths it will touch, and states a definition of done that another
agent could verify. Initial tasks for this charter:

| id           | title                                              | touches                          | depends_on   |
|--------------|----------------------------------------------------|----------------------------------|--------------|
| T-human-002  | Pipeline skeleton: run one annotator on one genome | src/pipeline.py, Makefile        |              |
| T-human-003  | Define summary.tsv columns in docs/schema.md       | docs/schema.md                   |              |
| T-human-004  | Write src/validate.py against the schema           | src/validate.py, tests/          | T-human-003  |
| T-human-005  | Compare two annotators on three genomes; recommend | relay/artifacts/T-human-005/     | T-human-002  |
| T-human-006  | Scale pipeline to all 12 assemblies                | src/pipeline.py, results/gff/    | T-human-002, T-human-005 |
| T-human-007  | Generate summary.tsv                               | src/summarize.py, results/       | T-human-004, T-human-006 |
| T-human-008  | README with rerun instructions and timings         | README.md                        | T-human-007  |

Seven tasks, four agents. Tasks are not assigned to agents; agents claim
whatever is open and unblocked when they wake up. New tasks will appear as
the work reveals them.

## Protocol version

Agents must implement `relay/PROTOCOL.md` v0.1.
