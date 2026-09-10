# Where the five reviews disagreed

Task `T-human-006`. Author: `lenin`, 2026-09-10.

The charter asks for "every point where the reviews reached different
conclusions, what each said, and either a resolution with evidence or an open
question for the coordinator". This is that list.

Two kinds of disagreement are separated deliberately. **Factual conflicts**
(§1–§4) are cases where the reviews cannot all be right and the matter is
decidable; each carries a resolution and the evidence for it. **Design
disagreements** (§5) are cases where the reviews read the same evidence and
drew different conclusions; those are not resolved here, because resolving them
is `T-human-011`'s job, and stating them sharply is the most useful thing this
document can do.

---

## 1. Citation conflicts

**What happened.** Four reviews recorded a DOI for the same tool and one
recorded a different one. This is decidable, so it was decided: every DOI in
the merged bibliography was resolved against Crossref on 2026-09-10 by
[`scripts/review/verify_dois.py`](../../scripts/review/verify_dois.py), and the
full result is in [`doi-verification.tsv`](doi-verification.tsv).

Of 126 merged entries: 110 resolve to a work whose registered title matches;
4 have no DOI by design (repository documentation, cited by URL); 2 are arXiv
DOIs, which Crossref legitimately does not serve because arXiv registers with
DataCite; 1 is a journal supplement Crossref stores without a title. **Nine
resolve to an unrelated paper or do not resolve at all. All nine are from
`trotsky`'s bibliography** (`T-human-012`), and in every case the other reviews
that carried the same work agree on a different DOI that does resolve
correctly.

| Tool | DOI the other reviews recorded | Resolves to | `trotsky`'s DOI | Resolves to |
|---|---|---|---|---|
| KA/KS ratio test | `10.1101/gr.200901` | the KA/KS paper ✔ | `10.1101/gr.190902` | **nothing** |
| GeneMark-ES 2005 | `10.1093/nar/gki937` | GeneMark-ES ✔ | `10.1093/nar/gki987` | "Analysis of repetitive element DNA methylation by MethyLight" |
| N-SCAN 2006 | `10.1089/cmb.2006.13.379` | N-SCAN ✔ | `10.1101/gr.5389206` | **nothing** |
| CONTRAST 2007 | `10.1186/gb-2007-8-12-r269` | CONTRAST ✔ | `10.1101/gr.6734807` | **nothing** |
| GeneMark-EP+ 2020 | `10.1093/nargab/lqaa026` | GeneMark-EP+ ✔ | `10.1186/s12859-020-03703-2` | "SQMtools: automated processing and visual analysis of 'omics data" |
| Helixer (journal) | `10.1038/s41592-025-02939-1` | Helixer ✔ | `10.1093/nar/gkad1003` | "VEuPathDB: the eukaryotic pathogen, vector and host bioinformatics resource centre" |
| Tiberius 2024 | `10.1093/bioinformatics/btae685` | Tiberius ✔ | `10.1101/2024.07.19.604245` | "Diversity and functional specialization of oyster immune cells" |
| BRAKER3 2024 | `10.1101/gr.278090.123` | BRAKER3 ✔ | `10.1093/bioinformatics/btae010` | "FAVA: high-quality functional association networks…" |
| GeneMark-ETP 2024 | `10.1101/gr.278373.123` | GeneMark-ETP ✔ | `10.1093/gigascience/giae016` | "Streamlining remote nanopore data access with slow5curl" |

Author attributions in the same entries are also wrong in ways consistent with
the DOIs: the Tiberius entry lists "Grewe, Felix and Gabriel, Lars" (the paper's
first author is Gabriel), the Helixer 2021 entry lists "Stiehler, Felix and
Gabriel, L", and the Helixer 2023 entry lists "Holst, Felix and Bushmanova, E".

**Resolution.** The merged bibliography keeps all nine entries so the record is
complete, but each now carries a `verified = {DO NOT CITE: …}` field naming what
its DOI actually resolves to, so nothing downstream can cite them by accident.
Every claim in [`README.md`](README.md) that these entries would have supported
is instead sourced to the DOI the other reviews recorded.

**A related list that is not this one.**
[`scripts/review/check_conflicts.py`](../../scripts/review/check_conflicts.py)
groups the five source bibliographies by title and reports every group carrying
more than one DOI. Four of those groups are conflicts of the kind tabulated
above; four are a preprint DOI and its journal DOI for the same work (Helixer,
ANNEVO, Nucleotide Transformer, SegmentNT), where **both identifiers are
valid**. The script now labels the two categories separately, so a later reader
of its output does not attribute a legitimate preprint/journal pair to whichever
review recorded the preprint.

**Not a resolution of the underlying question.** This says the identifiers are
wrong. It does not say the review's *reading* of those papers is wrong — that
would need the numbers checked one by one against the correct sources, which
this tick did not do. See §3.

---

## 2. Repository identity conflicts

Five reviews inventoried overlapping repository sets under five different
schemas; the merged table is [`repos.tsv`](repos.tsv) and it flags 23
repositories where the recorded values differ. The disputed rows were
re-queried from the GitHub API on 2026-09-10
([`verify_repos.py`](../../scripts/review/verify_repos.py),
[`repo-verification.tsv`](repo-verification.tsv)).

**2.1 Three different repositories were called "Helixer".** `lenin`, `engels`
and `stalin` used `usadellab/Helixer`; `marx` used `weberlab-hhu/Helixer`;
`trotsky` used `gglyptodon/helixer`.

*Resolved.* `weberlab-hhu/Helixer` **redirects** to `usadellab/Helixer`, so
those two are the same repository under an old and a new name — 305 stars,
GPL-3.0, last pushed 2026-07-22. `gglyptodon/Helixer` is a **different
repository**: a personal fork with 0 stars, last pushed 2024-03-24.

This changes a conclusion. `trotsky`'s review states that Helixer's "commits
have stalled in the last 12 months", recording 0 commits. That is a true
statement about the fork it queried and a false one about the project: the
canonical repository has 37–49 commits in the trailing 12 months, on which the
other four reviews agree.

**2.2 Two code links do not point at the software they name.**
`washut/contrast`, given as CONTRAST's repository, returns **404**.
`NBISweden/GAAS`, given as MAKER2's repository, exists but is an unrelated
annotation toolkit (Perl, GPL-3.0, 221 stars); MAKER is at `Yandell-Lab/maker`,
which three reviews recorded.

**2.3 `BrentLab/Twinscan` exists.** `engels` recorded TWINSCAN's distribution
as "historical, pending"; `stalin` linked the repository. It is real: C, MIT,
2 stars, last pushed 2023-06-23. Resolved in favour of `stalin`.

**2.4 Licences: the GitHub detector disagrees with the repositories.** For
AUGUSTUS, BRAKER, GALBA, SNAP, egapx, GeneMark-ETP, MAKER, ANNEVO and SpliceAI
the five reviews recorded up to five different values each. The pattern is
consistent: reviews that queried the GitHub API recorded `NOASSERTION`, `None`
or `unknown`; reviews that opened the LICENSE file or README recorded the
actual terms. `engels` states the general rule explicitly — GitHub's licence
detector is not authoritative.

*Resolved in favour of the reviews that read the licence documents*, with two
consequences worth carrying forward: SNAP's current LICENSE says **MIT** while
its 2004 article described GPL, so the release and the current source must be
distinguished (`trotsky` records GPL-2.0, which matches neither); and
`ncbi/egapx` is a **public-domain NCBI notice plus bundled third-party
licences**, which is materially different from the `NOASSERTION` three reviews
recorded from the API.

**2.5 Commit counts differ by method, not by drift.** All five snapshots were
taken on 2026-09-09. `marx` counted with `git clone --shallow-since`; the others
used GitHub REST `since` filters including merge commits. Differences of a few
commits (Helixer 37 vs 49, Augustus 2 vs 5, funannotate 279 vs 480) are
method, and the merged table keeps both values rather than picking one. One
gap is larger than method explains: **Tiberius, 174 commits in four reviews
versus 100 in `trotsky`'s**. Unresolved; low stakes.

---

## 3. Numeric conflicts in the publications tables

These are cases where two reviews report a different number for the same
quantity from the same paper. Where four reviews cite a PMC identifier or a
full-text read and one does not, the resolution below follows the sourced
value; that is a statement about provenance, not a re-reading of the papers.

| Quantity | Sourced value | Conflicting value | Status |
|---|---|---|---|
| N-SCAN, human, gene Sn/precision | 35.6 / 25.1, from [gross2007contrast] Table 1, cited by `engels`, `marx`, `stalin` with PMC2246271 | `trotsky`: gene Sn 0.68, Sp 0.65; exon 0.86/0.84 | Follow the sourced value. The gap is large enough (35.6 vs 68) to change how N-SCAN reads against CONTRAST, which is the crux of §5.3 |
| CONTRAST, human, gene Sn/precision | 58.6 / 35.5, exon 92.8 / 72.5, from Table 1, cited by three reviews | `trotsky`: gene 0.59/0.61, exon 0.89/0.87 | Sensitivity agrees; precision does not. `stalin` adds the important caveat that the paper's "specificity" column is precision, not a true-negative rate, and that incomplete CCDS makes it an underestimate |
| AUGUSTUS 2003, human | exact gene Sn/prec 46%/45%, exon 80%/80% (`engels`, from the author manuscript, Tables 1 and 3) | `trotsky`: exon Sn 0.86, Sp 0.69 | Follow the sourced value |
| Tiberius, mammals | exon/gene F1 89.7 / 55.1 (three reviews, PMC11645249) | `trotsky`: exon Sn 0.84, Sp 0.82 "on held-out test chromosomes" | Different quantity and different split; not comparable. Follow the sourced value |
| GENSCAN, Burset–Guigó | exon Sn/Sp 0.78/0.81, 243/570 exact genes | all five agree | No conflict |
| Helixer 2026 transcript F1 | fungi 0.5386, plants 0.4618, vertebrates 0.1977, invertebrates 0.3066 (`engels`, `marx`, `stalin`, PMC13076211 Table 2) | `trotsky`: "exon Sn 0.78, Sp 0.74 (out-of-clade plant benchmark)" | Different quantity; no traceable source. Follow the sourced values |

**One qualification that applies to the sourced values too.** `engels`'s annex
finds that the wording "exact primary CDS comparison" overstates what Helixer's
published transcript scores establish under the documented GffCompare
procedure, and `stalin`'s annex independently finds that Helixer's UTR-stripped
evaluation still uses matching rules that do not require exact coding starts and
stops. Both qualifications must travel with those Helixer numbers wherever they
are quoted. They are recorded in [`README.md`](README.md) §5.8 and §7.

---

## 4. The evidentiary status of the charter's own figures

The charter states that 1,409 EGAPx jobs on Galaxy consumed about 416,000
CPU-hours with a 45% failure rate and 24% CPU efficiency, sourced to
`nekrut/scalingPaper`.

- `trotsky` restates the figures as established fact.
- `marx` quotes them and cites `relay/TASK.md` as the source, explicitly, rather
  than the underlying repository.
- `lenin` records that `nekrut/scalingPaper` returns 404 to an unauthenticated
  client, marks the figures unverified, and argues that if they hold, the
  interesting quantity is the **45% failure rate** rather than the CPU-hours —
  a robustness failure being reported as a cost failure.
- `engels` and `stalin` do not use the figures.

**This is not resolvable by the agents.** `nekrut/axomeme` — the charter's
source for the HyphAeon design pattern that `T-human-011` is meant to transfer —
is unreadable for the same reason. A `question` has been open with the
coordinator since
[20260909T012944Z-lenin-0002](../../relay/messages/20260909T012944Z-lenin-0002.md).
**Open question 1**, below.

Partial mitigation exists for the cost half: EGAPx's own README supplies its
hardware floor and per-genome runtimes [ncbi2026egapx], so `T-human-009` is not
wholly dependent on the unreachable repository. It cannot supply the failure
rate, because documented runtimes are by construction the runs that succeeded.

---

## 5. Design disagreements

These are the substantive ones. All five reviews read overlapping evidence and
reached different conclusions; none of them is obviously wrong. They are stated
here as decisions `T-human-011` must make, with the strongest argument on each
side and, where one exists, the cheapest experiment that would settle it.

### 5.1 What decodes the structure: differentiable HMM, CRF, or splice graph?

- **Differentiable HMM, end-to-end** — `lenin` and `marx`. The argument is
  Tiberius's: the end-to-end integration, not Helixer's separate HMM
  post-processor, is what produced the gene-level jump [gabriel2024tiberius].
  `marx` does not keep the HMM as it stands: the variant proposed is
  candidate #5 of [`candidates.md`](candidates.md), which makes the duration
  and emission distributions **functions of per-genome covariates** instead of
  the fixed mammalian geometric durations. That is a direct answer to the
  defect `trotsky` cites as the reason to abandon the HMM, so this is a
  three-way disagreement about how to fix the length model, not two positions
  against one about whether it is broken.
- **A structured decoder, family unspecified, chosen by experiment** —
  `engels` and `stalin`. Both explicitly decline to pick, and both add the same
  requirement: whatever generates candidates must have its **recall ceiling
  measured before end-to-end accuracy**, because a cheap decoder looks excellent
  after discarding the difficult genes.
- **No HMM at all** — `trotsky`, in the most specific proposal of the five: a
  grammar-constrained DAG over candidate junctions, with edges weighted by an
  intronic span potential and traversal pruning in-frame stops, decoded with
  top-*K* paths to give native multi-isoform output. The stated motivation is
  that geometric self-transitions penalize long introns exponentially while
  explicit duration submodels cost O(D) or O(D²) and force hard length cutoffs.

**Assessment.** The motivation is real and independently confirmed: Tiberius's
HMM layer does fix geometric length distributions at empirical mammalian means
[gabriel2024tiberius §2.4], and that is the single most-cited concrete defect in
the strongest current tool ([`README.md`](README.md) §5.3). The published
alternative is not a DAG, though — it is GeneCAD's chromosome-scale CRF, which
reports ~9% transcript-F1 over Helixer and BRAKER3 on angiosperms
[liu2025genecad], and which is a real comparator the DAG proposal does not
mention. Note also that GeneCAD's grammar constraints (canonical start/stop,
donor/acceptor motifs, minimum intron and exon lengths) were **added in v0.4.0,
after the neural model** — evidence that the grammar is separable from the
network, which is the DAG proposal's core assumption.

*Cheapest decision:* hold the encoder fixed and swap the decoder. This is the
same experiment `engels` and `stalin` ask for from the other direction, so one
harness serves both. **Open question 3.**

### 5.2 Fix the parameter budget first, or prove the benefit first?

- `lenin`: target 2–5M parameters and treat exceeding 20M as evidence the
  geometry is wrong rather than that more capacity is needed — because a small
  model that generalizes across clades is a scientific claim about gene
  structure and a large one is not.
- `trotsky`: under 2M, as an architectural constraint.
- `marx`: the HyphAeon pattern, which is ~2M by construction.
- **`engels`, on the order of operations, not against the budget**: "Reducing
  the parameter count before establishing that benefit would optimize the wrong
  thing." The first question is whether comparative information helps
  reconstruct complete genes on genuinely held-out species; shrinking parameter
  count is not a substitute for establishing structural benefit. `engels`
  favours compact models and the charter's compute limits throughout, and is
  not proposing an unconstrained first experiment
  ([review §6](../../relay/artifacts/T-human-003/review.md#6-my-opinion-what-i-would-build),
  [20260910T042437Z-engels-0027](../../relay/messages/20260910T042437Z-engels-0027.md)).
- `stalin`: follow the charter's small-model ambition, but report capacity and
  **input preparation cost separately** — which reframes the disagreement
  usefully, because a 2M-parameter model that needs a Cactus alignment is not
  cheap.

**Assessment.** `engels` is right about the order of operations and `lenin` is
right that the budget is the part of the hypothesis nobody in the 2026 cohort
is testing: OrionGeno needs a compute-capability-7.0+ GPU with `mamba-ssm`,
GeneCAD sits on a foundation model, and Evo 2 is 7B–40B parameters. The
reconciliation `stalin` implies is to report **total cost** — parameters plus
alignment construction plus feature generation plus inference — as the budgeted
quantity, at which point the two positions stop conflicting. **Open question 3.**

### 5.3 Does the tree help at all?

This is the sharpest disagreement in the set, and it goes to the charter's
central hypothesis.

- **The charter's position**, endorsed by `marx`, `lenin` and `trotsky`:
  injecting the tree as a metric removes the need to learn the phylogeny, and
  that is what should buy clade independence.
- **`marx` states the counter-evidence most clearly, and against its own
  recommendation**: CONTRAST "threw the tree away and still gained from more
  informants, which N-SCAN could not". The numbers support the reading —
  CONTRAST, phylogeny-free, goes from gene Sn 50.8 with one informant to 58.6
  with eleven, while N-SCAN, with an explicit phylogenetic model, "performs as
  well using mouse as its only informant as it does with any combination"
  [gross2007contrast]. In 2007, discarding the tree scaled better than modelling
  it. (The prose of `marx`'s §6 gives this as "35.6 to 58.6"; `marx` has since
  corrected it — 35.6 is the N-SCAN-with-mouse row, and CONTRAST with mouse
  alone is 50.8. The table in `marx`'s §2 and the figures above are the right
  ones; do not restore the §6 numbers.
  [20260910T035517Z-marx-0027](../../relay/messages/20260910T035517Z-marx-0027.md).)
- **`engels` and `stalin` are agnostic and want it tested.** `engels`: "an
  attractive geometric embedding is not itself a validated evolutionary model",
  and proposes comparing an explicit tree likelihood or learned CTMC layer, a
  distance-aware small encoder, and an encoder with no tree information, with
  the decoder held fixed. `stalin` names ClaMSA as the comparator that makes the
  novelty question sharp: it already trains phylogenetic likelihoods
  discriminatively, so improvement must come from better *structural*
  prediction, robustness or total cost, not from adding a tree to a classifier.
- **`lenin` adds a revealed-preference argument.** Tiberius ships nine model
  configurations; exactly one uses ClaMSA and it is mammals-only, while all six
  clade models added in 2026 are sequence-only. Either comparative input buys
  less than assumed once a good sequence model exists, or generating it per
  clade is too expensive to scale. Both readings are bad news for a naively
  comparative design, and they point at different fixes.

**Assessment.** Three independent arguments — CONTRAST's 2007 result,
Tiberius's configuration list, and ClaMSA's existence — all say the tree's
contribution is unproven rather than assumed. This does not refute the charter;
it means **the charter's central claim is exactly the thing that has not been
measured**, which is also the strongest available argument for doing the
measurement. *Cheapest first look:* reproduce Tiberius ab initio vs ClaMSA mode
on mammals. Both the code and the weights are installed and MIT-licensed;
`lenin` estimates a day's work — but only if the published ClaMSA features for
the evaluated regions are downloadable and compatible with the released
checkpoints. If they must be regenerated, the comparison carries a
preprocessing bill that this estimate does not contain and that must be
reported as a separate line.

**What that reproduction can and cannot settle.** It measures the supplied
comparative pipeline *as a whole* on the evaluated checkpoints and regions, and
cannot separate tree geometry from alignment information or either from the
feature generator's supervised human exposure. Read its outcome as
prioritization evidence for candidate #4, not as a verdict: a small gain on
these mammalian regions does not show that other comparative encoders, other
informant sets, or other clades have no benefit, and is not grounds for
revising the charter's hypothesis. Attribution to the alignment or to the tree
requires the fixed-decoder ablation ([`candidates.md` §4](candidates.md), with
its matched MSA/no-tree, MSA/tree-token and MSA/tree-metric arms), and even
that attributes only within the encoders, informants, species and split it
evaluates; extrapolation to untested encoders or clades needs its own evidence.
**Open question 2.**

### 5.4 Is tree-as-metric mathematically well-posed?

`stalin` raises an objection nobody else does and it is a real one: **MDS
coordinates are not unique under rotation or sign change**. If the encoder
depends on those arbitrary coordinates, the representation needs a stable
convention or genuinely invariant operations, and the embedding's distortion
should be measured rather than assumed harmless. `engels` independently asks for
stability tests under changed taxon order, removed taxa and rescaled branch
lengths.

**Not a disagreement about direction, but an unanswered technical question**
that `T-human-011` must answer before proposing a tree-as-metric arm. No other
review addresses it. Cheap to settle: apply a random rotation to the MDS
embedding and check whether predictions change.

Keep the two halves apart. Rotation, sign change and a consistent permutation
of the taxon axis are **representation symmetries**: the evidence is
unchanged, so the prediction should be too, and equality is the right test.
Removing taxa and rescaling branch lengths **change the evidence**; predictions
may legitimately move, and requiring equality there would reward an encoder
that ignores the comparative channel. For those, measure robustness,
calibration and the degradation curve. `stalin` asks for stability tests, not
for invariance, throughout
([20260910T033840Z-stalin-0027](../../relay/messages/20260910T033840Z-stalin-0027.md)).

### 5.5 UTRs in version 1?

- `marx`, explicitly against: Helixer 2026 strips UTRs from both reference and
  prediction before scoring, Tiberius scores only the longest CDS per gene, and
  no paper reports transcript-level UTR accuracy — so **UTRs cannot yet be
  benchmarked against anything**.
- `trotsky`'s architecture includes a UTR class in its span head and cites
  OrionGeno and TOGA2, both of which do predict UTRs.
- The charter says "and eventually UTRs".

*Resolution proposed:* `marx` is right for version 1 and the reason is
measurement, not difficulty. Predicting UTRs is cheap; scoring them is not
currently possible against the panel `T-human-007` is building. Predict them if
convenient, exclude them from the headline metric, and revisit when a UTR
reference with an error bar exists.

### 5.6 What is the deliverable — an architecture or a measurement?

- `marx` and `trotsky` deliver an architecture.
- `engels` and `stalin` deliver an experimental protocol and decline to commit
  to an architecture before it runs.
- `lenin` revised position mid-review: the comparative-geometry design space is
  occupied (OrionGeno, ANNEVO, Tiberius-with-ClaMSA), the measurement is not, so
  the architecture should be demoted to one of three arms and the first
  effort should go into the comparison.

**Assessment.** This is the most consequential disagreement of the six, because
it decides what `T-human-011` is *for*. The reviews' own evidence favours the
measurement: the field's evaluation is fragmented enough that a new model can
be made to look excellent by choosing its benchmark
([`README.md`](README.md) §5.8), and at least two well-resourced groups have
already placed the phylogeny-as-input bet. See
[`candidates.md`](candidates.md), where the shortlist is ranked on exactly this
basis. **Open question 4.**

---

## 6. Non-disagreements worth recording

Three things looked like disagreements and are not:

- **"Tiberius is one model" vs "six models".** Both are true of different
  papers. The 2024 paper is a mammal model; the 2026 extension is six
  lineage-specific models. Reviews describing it either way are describing
  different artifacts.
- **Whether Helixer beats classical tools.** `marx` reports that in fungi
  Helixer's transcript F1 (0.54) is below GeneMark-ES's (0.60), which reads as
  contradicting other reviews' "deep learning wins". It does not: the same table
  has Helixer far ahead in plants (0.46 vs 0.09) and vertebrates (0.20 vs 0.02).
  The finding is that **clade decides the winner**, which is the project's whole
  premise.
- **"EGAPx is the industry standard" vs "EGAPx has no benchmark".** Both hold.
  It produces RefSeq annotations at scale and has no published methods paper and
  no controlled gene-level benchmark; `TITLE:"Gnomon"` and `TITLE:"EGAPx"` each
  return zero hits in Europe PMC. Any benchmark including it is benchmarking
  software, not a published method, and `T-human-007` should record that
  asymmetry.

---

## 7. Open questions for the coordinator

1. **`nekrut/axomeme` and `nekrut/scalingPaper` return 404.** Can they be made
   readable, or should `T-human-011` proceed with the HyphAeon pattern known
   only from the charter's paragraph, and `T-human-009` treat the Galaxy figures
   as uncited? Open since 2026-09-09 (§4).
2. **`T-human-012`'s status.** Nine of its 26 bibliography entries carry DOIs
   that resolve to unrelated papers or nothing, its Helixer repository
   conclusion is drawn from a personal fork, and two of its code links do not
   exist (§1, §2). Nothing in the synthesis rests on it as sole evidence.
   Should it be re-verified, annotated, or left in the record as-is? This is
   the coordinator's call, not `lenin`'s.
3. **Does `T-human-011` get a compute budget for the two cheap decisive
   experiments** — Tiberius ab initio vs ClaMSA on mammals (§5.3), and one
   fixed-decoder encoder ablation (§5.1, §5.2)? Both fit the charter's "laptop
   or one consumer GPU", with the caveat that Tiberius's own stack needs a
   compute-capability-appropriate GPU or a 30-minute PTX JIT.
4. **Is the Phase 3 deliverable an architecture or a measurement (§5.6)?** The
   charter says "two or three candidate model designs … with a recommendation",
   which is compatible with either reading. A `decision` here would let
   `T-human-011` start from the right premise rather than relitigating it.
