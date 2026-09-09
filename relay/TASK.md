# Shared task charter: fast, accurate, species-independent eukaryotic gene prediction

Owned by the human coordinator (`human`). Agents read this on every tick;
only the coordinator edits it. Changes are announced with a `decision`
message. Version: charter v1, 2026-09-09.

## Goal

Build a small, fast, geometrically principled model that predicts protein
coding gene structures (exons, introns, splice sites, start and stop codons,
and eventually UTRs) in any eukaryotic genome, with accuracy at or above the
best current tools and at a small fraction of their compute, and that
generalizes across clades with very different intron length distributions,
genome sizes, and GC content without per-clade retraining.

The project has two halves. First, a rigorous evaluation of what exists:
methods, data, and their real costs. Second, a design and prototype of the
new model, informed by that evaluation. Finished for this charter means the
Phase 1 to Phase 3 deliverables below are merged and the coordinator has
issued a `decision` selecting a model design to prototype. Prototyping is
Phase 4 and will get its own charter revision.

## Background the agents should start from

- Gene prediction is not intrinsically hard. A single comparative signal,
  the KA/KS ratio of human-mouse alignment windows, separates coding from
  non-coding with about 9.5% false negatives and 2 to 3% false positives
  (Nekrutenko, Makova, Li 2002, PMC155263). Simple signals with the right
  geometry go a long way.
- Current production tools are slow, brittle, and clade-limited. NCBI's
  EGAPx chains miniprot/ProSplign, STAR, minimap2, and Gnomon; it needs
  32 CPUs and 256 GB RAM, and declares fungi, protists, and nematodes out
  of scope. On Galaxy, 1,409 EGAPx jobs over 290 genomes consumed about
  416,000 CPU-hours with a 45% failure rate and 24% CPU efficiency
  (nekrut/scalingPaper). Tiberius (CNN + LSTM + differentiable HMM) is
  accurate ab initio but ships separate clade models and wants an 8 GB+ GPU.
- HyphAeon (nekrut/axomeme) shows the design pattern we want to transfer:
  about 2M parameters, a 2D axial transformer over a codon alignment
  (site axis) and taxa (phylogeny axis), with the tree injected as a metric
  through 4D MDS embeddings and Tree-RoPE, so the network never has to learn
  the phylogeny. It amortizes hours of maximum-likelihood optimization into
  one forward pass. The working hypothesis is that gene structure has the
  same property: once alignment, tree, and codon geometry are given as
  coordinates, the remaining decision surface is low-dimensional.
- Alignments of non-coding sequence matter as much as coding. Introns,
  splice sites, and UTRs have their own conservation signatures. UCSC
  provides whole-genome multiple alignments (multiz and Cactus based),
  conservation tracks (phyloP, phastCons), RNA-seq, regulation, and
  variation tracks for many assemblies (genome.ucsc.edu,
  hgdownload.soe.ucsc.edu). Ensembl, NCBI RefSeq, and Zoonomia are further
  sources of alignments and reference annotations.

## Scope

- In scope: literature and software review of gene prediction (ab initio,
  comparative, evidence-based, deep learning); inventory of usable data;
  a benchmark with held-out species across clades; cost baselines of
  existing tools; simple comparative baselines; design of the new model;
  a decision on what to prototype.
- Out of scope for this charter: functional annotation, non-coding RNA
  genes, organelle genomes, genome assembly, building a production
  pipeline, and training a large model. A prototype is Phase 4, after a
  `decision`.

## Deliverables

Phase 1, review (every agent participates; see tasks T-human-002 to 005):

1. `relay/artifacts/T-human-00N/review.md` per agent: an independent
   review of publications (full text where openly available) and of
   active GitHub repositories on gene prediction, in the structured format
   the task file specifies, plus a one-page opinion on what to build.
2. `docs/review/` (T-human-006): one synthesized review, a curated
   bibliography in `docs/refs/refs.bib`, a method comparison table, and an
   explicit list of where the four reviews disagreed and why.

Phase 2, ground truth and baselines:

3. `docs/benchmark.md` and `benchmark/` (T-human-007): a species panel
   spanning clades and intron length regimes, the reference annotations
   used, metrics, and leakage rules for held-out species.
4. `docs/data-sources.md` and `scripts/data/` (T-human-008): what
   alignments, conservation, expression, and annotation data exist, how to
   fetch a locus window with its tree, sizes, and licences.
5. `docs/cost-baseline.md` (T-human-009): measured or documented runtime,
   memory, and hardware for existing tools, and a target budget for ours.
6. `baselines/kaks/` (T-human-010): the KA/KS style comparative classifier
   run on the benchmark, as the floor every later model must beat.

Phase 3, design:

7. `docs/design/proposal.md` (T-human-011): two or three candidate model
   designs in the HyphAeon spirit, with inputs, geometry, parameter budget,
   training data plan, expected failure modes, and a recommendation. Ends
   with a `proposal` message to `human` requesting a `decision`.

## Constraints

- Literature access: use open sources only (PubMed Central OA, Europe PMC,
  bioRxiv, arXiv, publisher OA pages, Semantic Scholar or OpenAlex APIs).
  Do not scrape paywalled content. Commit DOIs, metadata, abstracts, and
  your notes; commit full text only when its licence allows.
- Data: never commit files over 5 MB. Commit fetch scripts and manifests
  with checksums. Respect UCSC and NCBI usage policies and rate limits.
- Code: Python 3.11, dependencies declared in `pyproject.toml`; PyTorch for
  models. Adding a dependency outside the standard scientific stack needs
  a `proposal`.
- Compute: everything in Phases 1 to 3 must run on a laptop or one
  consumer GPU. If a step needs more, post an `alert` instead of running it.
- Claims: every quantitative statement in a review or proposal carries a
  citation or a link to the artifact that produced it.
- Independence: the four Phase 1 reviews are done blind. Do not read another
  agent's review artifact until your own task is in `review`.

## Priorities

1. Phase 1 reviews. Nothing else is claimable until each agent has its
   review task in `review` or `done`.
2. Synthesis (T-human-006).
3. Data inventory and benchmark design (T-human-008, T-human-007), which
   can proceed in parallel.
4. Cost baseline and KA/KS baseline (T-human-009, T-human-010).
5. Design proposal (T-human-011).

Among open tasks, prefer the one that unblocks the most others.

## Review and acceptance

A task moves to `done` when the coordinator merges its pull request after
at least one `review` message from an agent other than the owner. Review
artifacts under `relay/artifacts/` do not need a pull request; the owner
moves the task to `review` and the coordinator marks it `done`.

## Protocol version

Agents must implement `relay/PROTOCOL.md` v0.1.
