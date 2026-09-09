# Independent review of eukaryotic gene prediction — `lenin`

Task: T-human-002 (slot 1 of 4). Status: **draft, run 1 of an expected 2–3.**
Author: `lenin` (Claude Code / Opus 5). Date of this draft: 2026-09-09.

This draft was written blind. I have not opened
`relay/artifacts/T-human-003/`, `T-human-004/`, or `T-human-005/`, and will
not until this task is in `review`, per charter §Independence.

Everything below that carries a number carries a citation or names the
artifact that produced it (`repos.tsv`, or the search log in §1). Where I
have a claim I could not verify in this run, it is marked **[unverified]**.
Nothing paywalled was retrieved; only abstracts and metadata were pulled for
non-OA items.

---

## 1. Search log

All searches run 2026-09-09 (UTC) from a laptop, scripted, no browser.

| # | Source | Interface | Queries | Notes |
|---|--------|-----------|---------|-------|
| 1 | Europe PMC | REST `/search`, `resultType=core`, sort by citations | 20 broad topical queries (see below) | Ranking by citation count pulled in irrelevant high-citation papers for short/ambiguous names (`CONTRAST`, `Exonerate`, `Liftoff`, `TWINSCAN`); those needed targeted re-queries. |
| 2 | Europe PMC | REST `/search`, `TITLE:` field queries | 25 targeted queries on named tools | Recovered GENSCAN, SNAP, GlimmerHMM, GeneID, miniprot, Splign, TOGA, CAT, SpliceAI, minisplice, Evo 2, AlphaGenome, phyloP, Cactus, OMArk, ANNEVO. |
| 3 | Europe PMC | REST `/search` by DOI | 19 DOIs | Pulled verified abstract, journal, PMCID, OA flag and licence for each headline method. |
| 4 | GitHub | REST API v3 (`/repos/{full}`, `/commits?since=`) | 32 repositories | Last push, commits in trailing 12 months, open issues, stars, language, SPDX licence, archived flag. Output: `repos.tsv`. |

Broad queries in pass 1 (verbatim): GENSCAN; AUGUSTUS ab initio eukaryotic;
BRAKER pipeline; GeneMark-ES self-training; SNAP Korf; Gnomon NCBI; EGAPx;
MAKER; TWINSCAN/N-SCAN; CONTRAST CRF; Helixer; Tiberius; DNA language model
annotation; SpliceAI; Ka/Ks comparative human-mouse; Zoonomia; TOGA; GALBA;
StringTie; benchmark comparison of gene predictors.

Hit counts for the useful ones: `TITLE:"TWINSCAN"` 4, `TITLE:"N-SCAN"` 2,
`TITLE:"Gnomon"` **0**, `TITLE:"EGAPx"` **0**, `TITLE:"ANNEVO"` 2,
`TITLE:"miniprot"` 3, `Helixer` 44, `Tiberius` 6, `EGAPx` (free text) 41 —
all 41 being genome papers that *used* it, not method papers.

**Finding from the search log itself, worth stating plainly:** the two tools
the charter singles out as production reality — NCBI's **Gnomon** and
**EGAPx** — have no methods paper. `TITLE:"Gnomon"` and `TITLE:"EGAPx"` each
return zero hits in Europe PMC. They are describable only through NCBI
database papers [kuhn2013refseq; sayers2026ncbi] and through the repository
documentation. Any benchmark that includes them is benchmarking software, not
a published method, and that asymmetry should be recorded in T-human-007.

### 1.1 What this run did not cover

Listed so the next run and the synthesis (T-human-006) know the shape of the
hole, rather than assuming coverage:

- No bioRxiv, arXiv, OpenAlex or Semantic Scholar pass yet. Europe PMC
  indexes bioRxiv, so preprints leaked in (Tiberius-multiclade, ANNEVO), but
  not systematically.
- No GitHub *code search* pass — only a hand-seeded repository list. Snowball
  from Tiberius/egapx dependents is still to do.
- **No install was attempted for any repository.** Every `install_tested`
  cell in `repos.tsv` reads `no`. This is the single largest gap and is the
  first item of the next run.
- Runtime and memory figures below are as-reported by authors. Nothing was
  measured here. Measurement belongs to T-human-009, not to this task, but
  the reported figures are not comparable across papers (different hardware,
  different genomes, different definitions of "annotation") and I flag that
  rather than tabulate them as if they were.

---

## 2. Publications table

Approach classes: **GHMM** = generalized HMM ab initio; **COMP** =
comparative/alignment-informed; **EVID** = evidence-based pipeline;
**DL** = deep learning; **HYB** = hybrid. Accuracy figures are as reported by
the cited authors on their own benchmark; they are *not* mutually comparable.

| Method | Year | DOI | Class | Inputs required | Clades trained/evaluated | Reported accuracy (benchmark, metric) | Runtime / hardware | Cross-species evidence | Code |
|---|---|---|---|---|---|---|---|---|---|
| GENSCAN | 1997 | 10.1006/jmbi.1997.0951 | GHMM | genomic DNA only | human/vertebrate | 75–80% of exons identified exactly, on standardized human/vertebrate sets [burge1997genscan] | not reported | Author reports consistent accuracy across C+G content and vertebrate groups [burge1997genscan] | web only |
| GeneID | 2000 | 10.1101/gr.10.4.511 | GHMM | DNA | Drosophila, human | see paper | very fast, C | per-species parameter files | guigolab/geneid |
| GlimmerHMM / TigrScan | 2004 | 10.1093/bioinformatics/bth315 | GHMM | DNA | plant, human, misc. | see paper | fast | per-species training | see `repos.tsv` note |
| SNAP | 2004 | 10.1186/1471-2105-5-59 | GHMM | DNA + training set | multiple novel genomes | conclusion: "every genome needs a dedicated gene finder" [korf2004snap] | fast | **Explicitly negative**: foreign gene finders are "highly inaccurate"; the nearest phylogenetic neighbour is not necessarily the best donor [korf2004snap] | KorfLab/SNAP |
| GeneMark-ES | 2005 | 10.1093/nar/gki937 | GHMM (self-training) | DNA only | fungi + novel euk. | comparable to or better than supervised training [lomsadze2005genemarkes] | not reported | self-training removes the per-species labelled-data requirement, not the per-species *fit* | non-free (see §3) |
| TWINSCAN / N-SCAN | 2003 / 2007 | 10.1101/gr.830003 **[unverified]** / 10.1002/0471250953.bi0408s20 | COMP+GHMM | target DNA + informant genome(s) | human/mouse, rat, worm | rat and C. elegans predictions validated by RT-PCR and sequencing [10.1101/gr.1959604; 10.1101/gr.3329005] | not reported | The historical proof that one informant genome buys real accuracy | unmaintained |
| CONTRAST | 2007 | 10.1186/gb-2007-8-12-r269 **[unverified]** | COMP+DL(discriminative) | target + multiple informants | vertebrate | see paper | not reported | phylogeny-*free* multiple-informant design — directly relevant to our geometry question | unmaintained |
| MAKER / MAKER2 | 2008 / 2011 | 10.1101/gr.6743907 / 10.1186/1471-2105-12-491 | EVID | DNA + EST/protein/RNA-seq + ab initio predictors | any, community-driven | outperformed by BRAKER3 [gabriel2024braker3] | heavy, MPI cluster | pipeline, not a model | canonical repo **not found** (see `repos.tsv`) |
| EVidenceModeler | 2008 | 10.1186/gb-2008-9-1-r7 | EVID (combiner) | multiple predictions + evidence | any | see paper | light | combiner | active |
| AUGUSTUS (+hints) | 2006 | 10.1186/1471-2105-7-62 | GHMM (+EVID via hints) | DNA, optional hints | many species, per-species parameter sets | the long-standing reference ab initio; superseded ab initio by Tiberius [gabriel2024tiberius] | CPU, hours | per-species parameter sets shipped | Gaius-Augustus/Augustus |
| BRAKER1 / 2 / 3 | 2016 / 2021 / 2024 | 10.1093/bioinformatics/btv661 / 10.1093/nargab/lqaa108 / 10.1101/gr.278090.123 | EVID | RNA-seq (1), proteins (2), both (3) | 11-species benchmark | BRAKER3 raises average **transcript-level F1 by ~20 percentage points** over BRAKER1/2; largest gain on large complex genomes; beats MAKER2, Funannotate, FINDER [gabriel2024braker3] | CPU cluster, hours–days | benchmark assumes a stated level of proteome relatedness [gabriel2024braker3] | Gaius-Augustus/BRAKER |
| GeneMark-ETP | 2024 | 10.1101/gr.278373.123 | EVID | RNA-seq + proteins | large eukaryotic genomes | "significantly improves" over ETP predecessors [bruna2024genemarketp] | CPU | — | non-free |
| GALBA | 2023 | 10.1186/s12859-023-05449-z | EVID (protein-only) | DNA + protein db | many | see paper | CPU | miniprot+AUGUSTUS | Gaius-Augustus/GALBA |
| Gnomon (NCBI) | — | **no methods paper** | EVID | alignments (Splign/ProSplign/miniprot), RNA-seq | RefSeq organisms | not independently published | NCBI production | described only via [kuhn2013refseq; sayers2026ncbi] | inside egapx |
| EGAPx | — | **no methods paper** | EVID pipeline | assembly + RNA-seq + proteins | explicitly excludes fungi, protists, nematodes (charter; repo docs **[unverified — confirm next run]**) | not independently published | 32 CPU / 256 GB RAM per charter **[unverified]**; charter cites 416k CPU-hours over 1,409 Galaxy jobs on 290 genomes, 45% failure, 24% CPU efficiency, from `nekrut/scalingPaper` — **repo returns 404, see §3** | ncbi/egapx |
| Helixer | 2021 / 2026 | 10.1093/bioinformatics/btaa1044 / 10.1038/s41592-025-02939-1 | DL (CNN+bLSTM) + HMM post-processor | DNA only | 2021: one vertebrate model over 186 animal genomes, one land-plant model over 51 plant genomes. 2026: fungal, plant, vertebrate, invertebrate | 2021: predictions "much less sensitive to genome length" than the then state of the art; outputs base-wise probabilities, not complete gene models [stiehler2021helixer]. 2026: "on par with or exceeding current tools", pretrained models usable without retraining [holst2026helixer] | GPU | **The first serious cross-species claim.** One model, many genomes | weberlab-hhu/Helixer, GPL-3.0 |
| Tiberius | 2024 | 10.1093/bioinformatics/btae685 | DL end-to-end (CNN + LSTM + differentiable HMM) | DNA only | trained on mammals; evaluated on human + 2 | **gene-level F1 62% on human vs 21% for the next best ab initio**; exon-intron structure of 2 of 3 human genes exactly right in de novo mode; ab initio accuracy *matches BRAKER3*, which uses RNA-seq + a protein database [gabriel2024tiberius] | **human genome in under 2 hours**; "fastest state-of-the-art"; GPU [gabriel2024tiberius] | 2024 version: mammals only | Gaius-Augustus/Tiberius, MIT |
| Tiberius, multi-clade | 2026 | 10.64898/2026.04.24.720536 | DL | DNA only | lineage-specific models for Mesangiospermae, Fungi, Vertebrata, Insecta, Chlorophyta, Bacillariophyta → **92% of available eukaryotic assemblies** | across **33 species**: gene-level F1 **+12 to +37 points over Helixer**, **+10 to +22 over ANNEVO**; approaches BRAKER3 in Mesangiospermae, Fungi, Bacillariophyta, Chlorophyta while being **~80× faster on GPU**; backend rewrite cut runtime 31% [gabriel2026tiberiusclades] | GPU | Strong — but note it is *six lineage-specific models*, not one general model | same |
| ANNEVO | 2026 | 10.1038/s41592-026-03036-7 | DL | DNA | multiple | beaten by Tiberius-multiclade by 10–22 F1 points at gene level [gabriel2026tiberiusclades] — no independent read yet, not OA | — | — | **[to fetch next run]** |
| TOGA | 2023 | 10.1126/science.abn3107 | COMP (alignment projection + orthology) | whole-genome alignment chains + reference annotation | 488 placental mammals, 501 birds | improves ortholog detection and annotation of conserved genes vs. state of the art; handles fragmented assemblies; also yields a genome-quality measure [kirilenko2023toga] | scales to hundreds of genomes | **Only for what is conserved relative to a reference** — cannot find clade-specific or fast-evolving genes | hillerlab/TOGA, MIT |
| CAT | 2018 | 10.1101/gr.233460.117 | COMP (HAL projection) | Cactus/HAL alignment + reference annotation | clades, personal genomes | see paper | cluster | same limitation as TOGA | repo **unmaintained**, 0 commits/12 mo |
| Ka/Ks ratio test | 2002 | 10.1101/gr.200901 | COMP (single statistic) | one pairwise alignment (human/mouse) | human/mouse | **false-negative rate lower than most current gene prediction methods and false-positive rate lower than all of them**, at the time; especially good on long exons and single-exon genes, which were then the hard cases [nekrutenko2002kaks] | trivial | the charter's floor | — |
| SpliceAI | 2019 | 10.1016/j.cell.2018.12.015 | DL (CNN, 10 kb context) | pre-mRNA sequence | human | accurate splice-junction prediction; variant-level validation against RNA-seq in 21 of 28 patients [jaganathan2019spliceai] | GPU | human-trained; cross-species use is off-label | Illumina/SpliceAI — **archived** |
| Pangolin | 2022 | 10.1186/s13059-022-02664-4 | DL | sequence, multi-tissue/species | 4 species | see paper | GPU | multi-species splicing | active |
| minisplice | 2026 | 10.1186/s13015-025-00293-7 | DL (1D-CNN) | sequence | vertebrates + insects | **7,026 parameters**; captures splice signals conserved across phyla; reveals mammal/bird-specific GC-rich introns; improves junction accuracy in minimap2/miniprot for noisy long reads and distant-homology proteins [yang2026minisplice] | trivial | **trained across phyla with one tiny model** | lh3/minisplice |
| miniprot | 2023 | 10.1093/bioinformatics/btad014 | alignment | protein + genome | any | comparable accuracy to prior protein-to-genome aligners, **tens of times faster** [li2023miniprot] | CPU, light | — | lh3/miniprot, MIT |
| Evo 2 | 2026 | 10.1038/s41586-026-10176-5 | DL foundation model | DNA | all domains of life; 9 Tbp training; 1M-token context | interpretability analysis shows learned representations of **exon-intron boundaries** among other features; no gene-structure benchmark reported [brixi2026evo2] | very large | broad by construction | ArcInstitute/evo2 |
| AlphaGenome | 2026 | 10.1038/s41586-025-10014-0 | DL | DNA, 1 Mb context | human/mouse regulatory | regulatory variant effects; not gene structure [avsec2026alphagenome] | large | — | google-deepmind/alphagenome |
| AlphaFold-3 gene-model scoring | 2026 | 10.1093/nar/gkag369 | DL (structure as a QC signal) | predicted proteins from gene models | F. graminearum, T. gondii, A. fumigatus | AlphaFold-3 scores support **65–84% of manually curated changes**; combining AF3 + Foldseek is most discriminative; the far cheaper Protenix-Mini retains the same discriminatory power [davison2026alphafoldannotation] | GPU (or cheap, with Protenix-Mini) | 3 fungal/protist species | — |
| Review (independent) | 2025 | 10.1093/bioadv/vbaf222 | — | — | — | extends the G3PO benchmark over AUGUSTUS, GENSCAN, GeneID, GlimmerHMM, SNAP + a gene-model-free NN + Helixer [djossou2025overview] | — | proposes a gene-model-based / -free / hybrid taxonomy | benchmark repo linked in paper |
| OMArk | 2025 | 10.1038/s41587-024-02147-w | QC | annotated proteome | any | detects erroneous gene inference in existing annotations [nevers2025omark] | light | — | active |
| BUSCO | 2021 | 10.1093/molbev/msab199 | QC | assembly or proteome | any | completeness, not correctness [manni2021busco] | light | — | active |

Supporting infrastructure cited but not a gene predictor: Progressive Cactus
[armstrong2020cactus] (>600 amniote genomes aligned, reference-free), phyloP
/ phastCons [pollard2010phylop] (36-mammal power analysis; UCSC conservation
tracks), Splign [kapustin2008splign], StringTie [pertea2015stringtie],
Liftoff [shumate2021liftoff], Zoonomia [christmas2023zoonomia].

Full BibTeX: `refs.bib` (47 entries; two DOIs marked `TO VERIFY`).

---

## 3. Repository inventory

Full table in `repos.tsv` (32 rows, GitHub API, 2026-09-09). Columns: URL,
last commit, commits in trailing 12 months, open issues, stars, language,
SPDX licence, archived, install tested, notes.

What the numbers say:

**Alive and moving.** `ComparativeGenomicsToolkit/cactus` (485 commits/12 mo,
704 stars) and `EBI-Metagenomics/genomes-pipeline` (189) are the busiest.
Among gene predictors, **Tiberius is the only one under real development**:
174 commits in 12 months, last push 2026-09-07, MIT, 138 stars, 21 open
issues. `Helixer` is second (49 commits, GPL-3.0). `guigolab/geneid` is a
surprise at 105 commits — a 1990s-lineage GHMM in C still being worked on.

**Alive but slow.** `ncbi/egapx` 26 commits, 207 stars, licence
`NOASSERTION`. `Gaius-Augustus/BRAKER` 11 commits but **103 open issues** and
467 stars — heavy use, thin maintenance. `Gaius-Augustus/Augustus` 5 commits
and **175 open issues**. The BRAKER/AUGUSTUS stack is load-bearing for the
field and is being maintained at a fraction of the rate it is being used.

**Effectively dead, still widely used.** `Comparative-Annotation-Toolkit` 0
commits/12 mo with 116 open issues. `agshumate/Liftoff` last pushed
2023-08-01, 552 stars, 78 open issues. `KorfLab/SNAP` last pushed 2022 and
still a MAKER dependency. `Illumina/SpliceAI` is **archived** with 507 stars
— the field's de facto splice-site model is a read-only repository.
`hillerlab/TOGA` 1 commit/12 mo despite being the engine behind the largest
comparative gene resources published [kirilenko2023toga].

**Licensing.** MIT: Tiberius, miniprot, TOGA, StringTie, gffcompare.
GPL-3.0: Helixer, Liftoff. Apache-2.0: CAT, evo2, alphagenome, borzoi.
`NOASSERTION` (i.e. GitHub could not resolve a standard licence): egapx,
BRAKER, GALBA, SNAP, SpliceAI, miniprot's deps, nucleotide-transformer.
AUGUSTUS resolves to no licence at all. **The GeneMark components inside
BRAKER have historically been under a non-free academic licence
[unverified — confirm the current terms next run].** If true, that has direct
consequences for whether BRAKER3 can be a redistributable benchmark baseline,
and T-human-007 should know before it designs around it.

**Three repositories could not be reached and are blockers, not curiosities:**

- `nekrut/axomeme` → **404**. The charter names it as the HyphAeon design
  pattern this whole project is meant to transfer (2M parameters, 2D axial
  transformer, 4D MDS tree embeddings, Tree-RoPE). I cannot read it.
- `nekrut/scalingPaper` → **404**. The charter's EGAPx cost figures (416k
  CPU-hours, 1,409 jobs, 290 genomes, 45% failure, 24% CPU efficiency) come
  from it. I have restated them as charter claims, not as verified ones.
- `gmod/maker` → 404. MAKER is not at that path; the canonical location
  (Yandell-Lab / SourceForge) needs to be found next run.

Both `nekrut/*` repositories are private, renamed, or deleted from the
perspective of an unauthenticated client. I am raising this to the
coordinator as a `question` (see §7); it affects every agent's Phase 1 task,
not only mine, and it blocks Phase 3 more than Phase 1.

A fourth item: `salzberg-lab/GlimmerHMM` reports its primary language as
**Elixir** with 0 stars — almost certainly not the canonical GlimmerHMM.
Treat that row as unresolved. Two Ensembl annotation repositories were not
fetched because the GitHub API rate limit was reached at the end of the run.

---

## 4. Data sources noted in passing

Not an inventory — that is T-human-008 — just what these papers actually
consumed, so the inventory task starts from evidence rather than a list.

- **Whole-genome multiple alignments.** Progressive Cactus, reference-free,
  demonstrated at >600 amniote genomes [armstrong2020cactus]; HAL is the
  container CAT projects through [fiddes2018cat]. UCSC multiz/Cactus
  alignments and the chain/net files TOGA requires [kirilenko2023toga].
- **Conservation tracks.** phyloP and phastCons, with the 36-mammal power
  analysis and the UCSC 44-vertebrate Conservation tracks
  [pollard2010phylop]. Note the paper's own result: enough power for strong
  selection at single nucleotides, moderate selection in 3-bp elements, and
  weaker or clade-specific selection only in longer elements. That is
  directly the resolution question for splice-site scale features.
- **Zoonomia**: 240+ placental mammal alignments and constraint calls
  [christmas2023zoonomia]; the substrate TOGA's 488-mammal resource was
  built on.
- **Reference annotations.** NCBI RefSeq [kuhn2013refseq; sayers2026ncbi];
  Ensembl; the G3PO benchmark set used by the 2025 review
  [djossou2025overview], which is the only ready-made, published,
  gene-prediction-specific benchmark I found this run.
- **Expression.** StringTie-assembled transcriptomes [pertea2015stringtie]
  are what BRAKER3 and EGAPx actually consume, not raw reads.
- **Proteins.** OrthoDB partitions are BRAKER2/3's protein input; the
  benchmark in [gabriel2024braker3] is explicitly parameterized by *how
  related* the available proteome is — a leakage axis T-human-007 must
  control, because "protein evidence" silently smuggles in the answer when
  the informant is close.

---

## 5. Failure modes of current tools

**5.1 Per-species fitting is the field's founding assumption.** Korf states
it as the conclusion of the SNAP paper: "every genome needs a dedicated gene
finder", after showing that foreign gene finders are highly inaccurate and
that the nearest phylogenetic neighbour is not reliably the best parameter
donor [korf2004snap]. GeneMark-ES removed the *labelled training set*
requirement via self-training [lomsadze2005genemarkes] but not the per-genome
fit. Helixer was the first to claim a single cross-species model
[stiehler2021helixer], and Tiberius 2026 still ships **six lineage-specific
models** rather than one [gabriel2026tiberiusclades]. Twenty-two years after
Korf, the assumption has been dented, not overturned. That is the gap this
project exists to close, and it is the right gap.

**5.2 Clade exclusion is explicit, not accidental.** EGAPx declares fungi,
protists and nematodes out of scope (charter; **[unverified — confirm from
repo docs next run]**). Tiberius before 2026 was mammals-only
[gabriel2024tiberius]. Even the 2026 extension reaches "92% of currently
available eukaryotic assemblies" [gabriel2026tiberiusclades] — and the
remaining 8% is precisely the tail (protists, early-branching lineages) that
a species-independent method would be judged on. **A benchmark drawn from the
92% cannot measure species independence.** T-human-007 should deliberately
oversample the excluded tail.

**5.3 Genome size and intron length.** Helixer's headline generalization
claim is specifically that its predictions are "much less sensitive to the
length of the genome" than the then state of the art [stiehler2021helixer] —
which tells you that length sensitivity was the dominant failure mode of GHMM
methods. BRAKER3's gains are "most pronounced for species with large and
complex genomes" [gabriel2024braker3], same signal from the other direction:
the large-intron regime was where everything was breaking.

**5.4 GC composition.** GENSCAN's answer in 1997 was to fit *distinct
parameter sets per C+G compositional region* [burge1997genscan] — isochores
handled by partitioning the model, not by learning a representation. minisplice
finds GC-rich introns specific to mammals and birds [yang2026minisplice],
i.e. composition is still clade-entangled at the splice signal itself, in
2026. Any model claiming clade independence must show its splice-site
performance stratified by GC, or the claim is untested.

**5.5 Resource use.** The compute story is the clearest argument for this
project. BRAKER3 is a CPU-cluster pipeline; Tiberius annotates the human
genome in under two hours [gabriel2024tiberius] and is ~80× faster than
BRAKER3 on a GPU [gabriel2026tiberiusclades]. EGAPx's 32 CPU / 256 GB floor
and the Galaxy figures (416k CPU-hours, 45% failure rate, 24% CPU efficiency)
come from the charter and I could not verify them — **the source repository
is unreachable (§3)**. If those numbers hold, the interesting quantity is not
the CPU-hours but the *45% failure rate*: nearly half the compute bought
nothing. That is a robustness failure being reported as a cost failure, and
it should be framed that way in T-human-009.

**5.6 Comparative methods only find what is already known.** TOGA and CAT
project annotation from a reference through an alignment
[kirilenko2023toga; fiddes2018cat]. They are excellent on conserved genes
and structurally cannot discover clade-specific ones. Any comparative design
we propose inherits this unless it is explicitly built to predict *without* a
reference annotation on the informant side.

**5.7 Maintenance risk is a real failure mode.** §3: the field's splice-site
model is archived, its comparative annotation toolkit is at zero commits, its
most-used ab initio pipeline has 103 open issues and 11 commits a year. When
T-human-009 measures baselines, some of them will not build. Budget for that.

**5.8 Evaluation is not standardized.** No two rows in §2 share a benchmark.
GENSCAN reports exact exons on 1997 vertebrate sets, BRAKER3 reports
transcript-level F1 on 11 species under an assumed proteome-relatedness
level, Tiberius reports gene-level F1 on human, Tiberius-multiclade on 33
species over 6 clades, the 2025 review uses G3PO
[djossou2025overview]. Comparing published numbers across these is
meaningless, and I have deliberately not built a "which is best" column.
**Constructing a single benchmark all of them run on is the highest-value
thing this project can do first**, independent of any model design.

---

## 6. Opinion: what I would build

One page, and I am willing to be wrong on all of it.

**The bet.** The charter's hypothesis is that gene structure, given alignment
+ tree + codon geometry as coordinates, has a low-dimensional decision
surface. I think the evidence for this is stronger than the charter claims,
and it comes from an unexpected place: **minisplice learns splice signals
conserved across phyla with 7,026 parameters** [yang2026minisplice]. Not 2
million — seven thousand. And Nekrutenko's Ka/Ks test achieved a false
positive rate below every gene finder of its day using *one statistic on one
pairwise alignment* [nekrutenko2002kaks]. Both results say the same thing:
when the comparative signal is presented in the right coordinates, the
remaining decision is nearly trivial. What is large in current models is not
the decision surface. It is the machinery for *manufacturing the coordinates*
from raw DNA — which is exactly the work Tiberius's CNN+LSTM stack and
Helixer's bLSTM are doing, and exactly the work a comparative model gets
handed for free.

**What I would build.** A small axial model over a *codon-aware multiple
alignment window*, in the HyphAeon pattern: site axis × taxa axis, tree
injected as a metric (MDS embedding + Tree-RoPE) rather than learned. Output
per alignment column: a distribution over {intergenic, 5'UTR, CDS-phase-0/1/2,
intron, 3'UTR} plus explicit donor/acceptor/start/stop head. Decode with a
differentiable HMM, following Tiberius — the end-to-end integration, not
Helixer's separate post-processor, is the part of Tiberius I would keep
[gabriel2024tiberius].

**Inputs.** Cactus/HAL alignment window for the target locus with its induced
subtree [armstrong2020cactus]; the target sequence; the informant sequences;
the tree. Explicitly **not** RNA-seq and **not** a protein database, so that
the comparison against BRAKER3 is honest and the method stays usable on a
genome with no transcriptome — which is most genomes.

**The geometry I would encode.** Three inductive biases, and I would defend
each as a *constraint*, not a feature: (1) **phase is modular arithmetic** —
CDS length ≡ 0 mod 3 and phase composes across exons; make it a structural
constraint in the decoder, not something a network rediscovers per clade.
(2) **The tree is a metric, not a sequence** — inject it, as HyphAeon does;
this is what should buy clade independence, because a model that never
memorizes *which* taxa it saw cannot overfit to a clade's taxon set.
(3) **Intron length is a nuisance parameter, not a signal** — the single
clearest cross-clade failure mode in §5.3. I would make the model
length-equivariant over intron spans (predict boundaries, not spans) so that
a 60 bp fungal intron and a 100 kb mammalian intron are the same object.

**Parameter budget.** I would target 2–5M and treat exceeding 20M as evidence
the geometry is wrong rather than as evidence more capacity is needed.
Tiberius is already accurate at modest size; the argument for going smaller
is not efficiency, it is that a small model that generalizes across clades is
a *scientific claim* about gene structure, and a large one is not.

**The biggest risk, and I think it is underrated.** Not accuracy —
**alignment availability**. The whole design assumes a good multiple alignment
with a tree at every locus in a novel genome. For a newly sequenced genome in
a sparsely sampled clade, that alignment either does not exist or is the most
expensive part of the pipeline, at which point we have moved the cost rather
than removed it, and we have re-introduced exactly the clade dependence we
set out to kill: the method works where relatives have been sequenced. Ka/Ks
worked in 2002 because human *and mouse* existed [nekrutenko2002kaks]. Every
comparative method since inherits that precondition. **Before committing to a
comparative design, Phase 2 must measure how accuracy degrades as informant
density and divergence degrade** — with the honest possibility that the
answer sends us to a Tiberius-style single-sequence model with comparative
signal as an optional refinement, rather than the reverse. I would make that
degradation curve a required experiment in T-human-007, and I would rather
find out there than in Phase 4.

Second risk, smaller but real: the field's evaluation is fragmented enough
(§5.8) that a new model can look excellent by choosing its benchmark. We
should fix the benchmark before we have a model, precisely so we cannot.

---

## 7. State of this task and what run 2 does

Not part of the deliverable format; recorded so the work is resumable by
whoever holds the lease.

Done: search log, publications table (30 methods/resources), `refs.bib` (47
entries), `repos.tsv` (32 repositories with live GitHub metrics), failure
modes, opinion.

Run 2, in priority order:

1. **Attempt installs** for Tiberius, Helixer, egapx, BRAKER in a fresh
   environment and fill the `install_tested` column honestly. Largest gap.
2. Resolve the four unresolved repositories (§3) and the two `TO VERIFY`
   DOIs in `refs.bib` (TWINSCAN 2003, CONTRAST 2007).
3. Fetch ANNEVO (preprint, since the Nature Methods version is not OA) and
   complete its row.
4. Confirm EGAPx's documented hardware floor and clade exclusions from the
   repository docs rather than from the charter.
5. bioRxiv/arXiv and OpenAlex passes; GitHub code-search snowball.
6. Confirm the GeneMark licence terms.

Then move T-human-002 to `review` with a `note` to `all`.

**One `question` goes to `human` this run** (charter §Constraints puts
out-of-reach resources on the coordinator): `nekrut/axomeme` and
`nekrut/scalingPaper` both return 404 unauthenticated. Both are cited in the
charter as load-bearing evidence — the HyphAeon design pattern and the EGAPx
cost figures respectively. Every quantitative claim I make that traces to
them is marked as a charter claim, not a verified one, and will stay that way
until access exists.
