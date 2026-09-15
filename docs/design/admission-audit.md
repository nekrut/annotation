# Label admission audit on the train panel (T-human-013)

Phase 4 milestone 1, proposal section 3.6, run on the ten `train` species of
`benchmark/panel.tsv` with `model/labels/admission.py` (metadata audit on the
reference GFF3, then the FASTA checks on the checksummed genome). One admitted-label
manifest per species is committed under `model/labels/manifests/` with its summary
JSON; the checksums below are what later tasks train on. The audit does not repair the
GFF, does not change the benchmark denominator, and keeps every reason as an independent
flag; the manifest lists them all per representative.

Inputs are the panel's GFF3 files (md5 as in `panel.tsv`) and the NCBI genomic FASTA,
fetched with `benchmark/fetch.py` and verified against NCBI's `md5checksums.txt`.
Grammar parameters: m = 20 (minimum intron), translation table 1 for every train
species (`genetic_code` column). All ten ran on one laptop-class container (4 CPUs,
15 GB); wall clock and peak RSS per species are in the last table.

## 1. Admission order and what each stage means

1. **Benchmark filters** (`benchmark/score.py`): nuclear primary sequences
   (`select_seqids`: organelles, alt loci, sequences under 10 kb dropped) and
   protein-coding transcripts (`pseudo=true` and non-coding `gene_biotype` dropped).
2. **Representatives** (proposal 4.4): one transcript per gene, most CDS bases, then
   longest exonic span, then annotation order.
3. **Topology mask** (4.4): connected components of overlapping full CDS spans per
   sequence and strand; every representative in a component with more than one
   gene is masked, the whole component span becomes unconstrained.
4. **Metadata audit** (3.6) on the raw rows: `exception` and `transl_except` tags,
   partial declarations (`start_range`/`end_range`, swapped on the minus strand;
   `partial=true` without a range is *unlocated*; a range on an interior row is
   an interior partial), overlapping rows, positive gaps shorter than m, invalid or
   inconsistent phase (first row 0 unless 5'-partial, then `(-p) % 3` with p advanced
   by each row's length), `transl_table` conflicts, and a complete chain whose length
   is not a multiple of 3.
5. **FASTA audit** (3.6) on the spliced CDS in transcriptional orientation: any
   non-ACGT base, coordinates beyond the sequence, initiator (ATG; alternative
   initiators are not enabled), terminal stop, in-frame internal stop, length mod 3.
   A 5'-partial chain with first-row phase h starts with p = (-h) % 3 missing bases;
   the observed suffix of that codon is marginalized (every suffix has a non-stop
   completion), and no leading bases are trimmed from later exons.

A declared partial end is admissible only where it touches the true sequence edge
(the edge-initialized grammar of 3.1); `partial_5`/`partial_3` are recorded but do not
mask by themselves, `partial_away_from_edge` and `partial_unlocated` do. A
representative is **admitted** when it has no masking reason; every masked
representative contributes its full CDS span, introns included, to the unconstrained
mask, unioned with the topology components.

## 2. Per-species counts

| Train species | CDS transcripts in file | After benchmark filters | Loci = representatives | Non-representative omitted | Conflict components / masked | Metadata-masked | Sequence-masked | Masked (union) | **Admitted** | Admitted % of representatives |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| M. musculus | 98,005 | 97,324 | 22,183 | 75,141 | 137 / 370 | 108 | 72 | 489 | **21,694** | 97.8% |
| D. rerio | 94,046 | 93,990 | 28,402 | 65,588 | 315 / 892 | 8,061 | 207 | 8,828 | **19,574** | 68.9% |
| X. tropicalis | 45,171 | 45,054 | 21,788 | 23,266 | 96 / 217 | 3,926 | 139 | 4,128 | **17,660** | 81.1% |
| D. melanogaster | 30,811 | 30,746 | 13,965 | 16,781 | 234 / 530 | 489 | 443 | 991 | **12,974** | 92.9% |
| C. elegans | 30,560 | 28,590 | 19,971 | 8,619 | 65 / 132 | 30 | 29 | 187 | **19,784** | 99.1% |
| A. thaliana | 48,268 | 48,147 | 27,444 | 20,703 | 30 / 61 | 159 | 36 | 224 | **27,220** | 99.2% |
| Z. mays | 57,350 | 57,071 | 34,039 | 23,032 | 81 / 166 | 2,076 | 127 | 2,235 | **31,804** | 93.4% |
| S. cerevisiae | 6,027 | 6,002 | 6,002 | 0 | 69 / 138 | 47 | 0 | 144 | **5,858** | 97.6% |
| N. crassa | 10,812 | 10,784 | 9,729 | 1,055 | 0 / 0 | 7 | 0 | 7 | **9,722** | 99.9% |
| D. discoideum | 13,315 | 13,179 | 13,155 | 24 | 2 / 4 | 199 | 57 | 218 | **12,937** | 98.3% |

The four species the proposal's 4.4 table covered reproduce it exactly at this
stage: retained transcript IDs 6,002 / 28,590 / 30,746 / 97,324, loci 6,002 / 19,971 /
13,965 / 22,183, conflict components 69 / 65 / 234 / 137 and topology-masked
representatives 138 / 132 / 530 / 370 for *S. cerevisiae*, *C. elegans*,
*D. melanogaster* and *M. musculus* (the `topology` reason count for *D. melanogaster*
is 539 because nine FlyBase transcript IDs recur on more than one sequence; the
manifest keys rows by sequence, strand and ID). Metadata-masked counts a
representative with at least one masking metadata reason; sequence-masked one with at
least one FASTA reason; the union is the number of representatives masked for any
reason including topology.

### Benchmark filter reasons

| Species | Filtered before selection |
|---|---|
| M. musculus | biotype=C_region: 20; biotype=D_segment: 14; biotype=J_segment: 82; biotype=V_segment: 374; pseudogene: 176; sequence: organelle (genome=mitochondrion): 13; sequence: shorter than 10000 bp: 2 |
| D. rerio | biotype=C_region: 8; biotype=V_segment: 35; sequence: organelle (genome=mitochondrion): 13 |
| X. tropicalis | biotype=C_region: 7; biotype=V_segment: 65; sequence: organelle (genome=mitochondrion): 13; sequence: shorter than 10000 bp: 32 |
| D. melanogaster | biotype=segment: 50; sequence: organelle (genome=mitochondrion): 13; sequence: shorter than 10000 bp: 2 |
| C. elegans | pseudogene: 1,958; sequence: organelle (genome=mitochondrion): 12 |
| A. thaliana | sequence: organelle (genome=chloroplast): 86; sequence: organelle (genome=mitochondrion): 35 |
| Z. mays | sequence: organelle (genome=chloroplast): 111; sequence: organelle (genome=mitochondrion): 168 |
| S. cerevisiae | pseudogene: 6; sequence: organelle (genome=mitochondrion): 19 |
| N. crassa | sequence: organelle (genome=mitochondrion): 28 |
| D. discoideum | sequence: organelle (genome=mitochondrion): 42; sequence: shorter than 10000 bp: 94 |

### Reason flags per representative (independent; a representative may carry several)

| Reason | M. musculus | D. rerio | X. tropicalis | D. melanogaster | C. elegans | A. thaliana | Z. mays | S. cerevisiae | N. crassa | D. discoideum |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| `topology` | 370 | 892 | 217 | 539 | 132 | 61 | 166 | 138 | 0 | 4 |
| `exception` | 20 | 7,959 | 3,692 | 121 | 12 | 38 | 1,851 | 47 | 0 | 35 |
| `transl_except` | 50 | 175 | 108 | 368 | 1 | 0 | 122 | 0 | 0 | 5 |
| `partial_5` | 38 | 101 | 182 | 0 | 17 | 54 | 140 | 0 | 6 | 86 |
| `partial_3` | 20 | 63 | 158 | 3 | 0 | 38 | 67 | 0 | 1 | 74 |
| `partial_away_from_edge` | 49 | 107 | 259 | 3 | 17 | 86 | 190 | 0 | 7 | 153 |
| `partial_unlocated` | 2 | 26 | 27 | 0 | 0 | 0 | 3 | 0 | 0 | 0 |
| `overlapping_rows` | 1 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| `short_gap` | 11 | 204 | 141 | 1 | 12 | 39 | 225 | 47 | 0 | 7 |
| `phase_invalid` | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| `phase_inconsistent` | 0 | 0 | 0 | 18 | 0 | 0 | 0 | 0 | 0 | 0 |
| `table_conflict` | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| `frame_length` | 1 | 5 | 1 | 51 | 0 | 24 | 1 | 0 | 0 | 28 |
| `coords_out_of_range` | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| `ambiguous_base` | 0 | 0 | 0 | 0 | 0 | 4 | 0 | 0 | 0 | 18 |
| `no_initiator` | 27 | 26 | 5 | 41 | 28 | 0 | 4 | 0 | 0 | 0 |
| `no_stop` | 1 | 13 | 3 | 51 | 0 | 24 | 1 | 0 | 0 | 28 |
| `internal_stop` | 45 | 168 | 131 | 421 | 1 | 32 | 122 | 0 | 0 | 39 |
| `seq_frame_length` | 1 | 5 | 1 | 51 | 0 | 24 | 1 | 0 | 0 | 28 |

`partial_5` and `partial_3` are informational (declared ends); the masking partial
reasons are `partial_away_from_edge` and `partial_unlocated`. No representative had
`phase_invalid`, `table_conflict` or `coords_out_of_range` in any train species, and
no FASTA sequence referenced by a representative was missing.

### Masked oriented bases and the auxiliary boundary catalog

| Species | Masked oriented bases (union of masked spans, one strand each) | Starts retained / unknown | Stops retained / unknown | CDS donors | CDS acceptors |
|---|---:|---:|---:|---:|---:|
| M. musculus | 27,311,554 | 35,873 / 59 | 30,168 / 20 | 187,365 | 189,666 |
| D. rerio | 210,880,957 | 38,150 / 132 | 33,280 / 76 | 248,000 | 249,215 |
| X. tropicalis | 96,069,339 | 26,159 / 189 | 23,658 / 161 | 193,733 | 193,952 |
| D. melanogaster | 12,249,013 | 16,297 / 23 | 16,078 / 36 | 44,113 | 44,500 |
| C. elegans | 709,866 | 24,322 / 52 | 21,200 / 0 | 103,593 | 103,556 |
| A. thaliana | 352,898 | 31,777 / 67 | 32,328 / 65 | 117,134 | 118,322 |
| Z. mays | 10,524,942 | 37,429 / 144 | 36,722 / 68 | 140,738 | 142,881 |
| S. cerevisiae | 281,218 | 5,960 / 0 | 6,001 / 0 | 281 | 281 |
| N. crassa | 12,218 | 9,966 / 7 | 9,831 / 1 | 15,845 | 15,862 |
| D. discoideum | 366,013 | 13,082 / 86 | 13,057 / 102 | 16,729 | 16,730 |

Auxiliary sites are distinct (sequence, strand, position) targets over **all**
benchmark-accepted transcripts, not only representatives, so known alternative sites
stay positive even when their chains are masked (4.4). A start or stop is *unknown*,
suppressed from both positive and negative supervision, when its end is declared
partial, unlocated, or fails the initiator/stop check on the representative; it is
not counted as retained at any transcript. Donors and acceptors are the first and last
intron base between CDS rows with a gap of at least 20; the counts reproduce the
proposal's 4.4 *all* rows exactly for the four species it covered (281 / 281,
103,593 / 103,556, 44,113 / 44,500, 187,365 / 189,666). Start/stop totals differ from
4.4 by the unknown split only.

## 3. Findings the coordinator should see

- **Exception tags dominate the metadata mask in the RefSeq-annotated vertebrates and
  maize.** *D. rerio* masks 7,959 representatives (28.0%) on `exception`, *X. tropicalis*
  3,692 (16.9%), *Z. mays* 1,851 (5.4%). The tag is almost always
  `annotated by transcript or proteomic data` (97,129 of the zebrafish CDS rows carrying
  an exception, 24,667 in frog, 6,625 in maize), which says the RefSeq transcript
  disagrees with the assembly somewhere, not that the genomic ORF is broken. Of those
  exception-masked representatives only 156 (zebrafish), 90 (frog) and 72 (maize)
  also fail a FASTA check. The proposal's policy is to exclude them all from the first
  chain loss and expose the loss; this audit keeps that policy. If the coordinator
  wants the ORF-passing subset back, it is one flag in the manifest
  (`exception` without any sequence reason) and a `decision`, not a code change.
- **Sequence-level failures are rare and mostly tag-explained.** *D. melanogaster* has
  421 representatives with an in-frame stop, 368 of them carrying `transl_except`
  (stop-codon readthrough and selenocysteine are annotated this way in FlyBase);
  *D. rerio* 168 and *X. tropicalis* 131, mostly with `exception`. Chains that fail only
  a sequence check with no tag are the declaration/sequence disagreements the
  proposal asks to list: they are in the manifests with `internal_stop`, `no_stop`,
  `no_initiator` or `seq_frame_length` alone.
- **Partial declarations away from the sequence edge** are the other large
  reason class (259 in frog, 190 in maize, 153 in *D. discoideum*, 107 in zebrafish);
  these chains are masked, their reliable boundaries stay in the auxiliary catalog.
- **Short gaps under m.** Positive gaps of 1 to 19 bases between CDS rows occur in
  every species with exception tags (ribosomal slippage in yeast's 47 Ty ORFs, frame
  adjustments elsewhere); they are grammar exceptions and are masked, never filled.
- **Ambiguous CDS bases** are confined to *D. discoideum* (18) and *A. thaliana* (4).
- *N. crassa* has no exception tags and no sequence failure at all: 9,722 of 9,729
  representatives admitted.

## 4. Finite legal numerator on admitted chains

An admitted chain has, by construction, a legal path in candidate A's grammar: rows
in frame, every intron at least m, an initiator, a terminal stop, no internal stop, no
ambiguous base. `model/labels/numerator_check.py` checks that claim with the grammar
itself on a sample: it cuts each sampled admitted gene's window (span plus 10
unrewarded flanking bases) from the checksummed FASTA, orients it, rewards exactly the
annotated CDS and intron bases, and runs the delayed-entry decoder with seams every
1,024 bases (m = 20, three duration components). The numerator is finite when the log
partition and the Viterbi score are finite and the Viterbi chain reproduces the
annotated spliced CDS and intron intervals exactly. The pure-Python decoder runs at
about 0.7 ms per base, so the sample is bounded: ten spliced admitted genes per species
with span at most 5,000 bases, seed 20260915.

| Species | Eligible admitted (spliced, span <= 5 kb) | Sampled | Failures | Mean window (bp) | Mean decode (s) |
|---|---:|---:|---:|---:|---:|
| M. musculus | 3,583 | 10 | 0 | 2,528 | 1.80 |
| D. rerio | 4,338 | 10 | 0 | 3,367 | 2.48 |
| X. tropicalis | 3,070 | 10 | 0 | 2,894 | 2.21 |
| D. melanogaster | 8,489 | 10 | 0 | 1,710 | 1.23 |
| C. elegans | 16,527 | 10 | 0 | 1,418 | 1.00 |
| A. thaliana | 19,884 | 10 | 0 | 1,342 | 0.98 |
| Z. mays | 18,670 | 10 | 0 | 2,019 | 1.57 |
| S. cerevisiae | 259 | 10 | 0 | 819 | 0.56 |
| N. crassa | 7,401 | 10 | 0 | 1,933 | 1.35 |
| D. discoideum | 8,415 | 10 | 0 | 2,158 | 1.50 |

No sampled admitted chain failed. Every admitted chain that is not in the sample is
covered by the constructive argument above; a chain-scoring pass over all admitted
chains belongs to the tensor implementation of T-human-014, where it costs one
forward pass per training crop rather than a pure-Python decode.

## 5. Manifests and checksums

Each manifest is a gzipped TSV with one row per representative: `transcript`, `gene`,
`type`, `seqid`, `strand`, `cds_start`, `cds_end` (1-based inclusive genomic span),
`rows`, `cds_len`, `missing_prefix`, `status`, `reasons`. Later tasks train on the
`admitted` rows and mask the full span of every other row.

| Species | GFF md5 | FASTA md5 | Manifest | Manifest md5 | Wall (s) | Peak RSS (MB) |
|---|---|---|---|---|---:|---:|
| M. musculus | `f0bc4339da97f2301929d2028bbb35ce` | `c0b0c4c3f54d2b480efe68a18bf7e42b` | `Mus_musculus.manifest.tsv.gz` | `fa1f59516033f1263dec4268c2dbb0c0` | 51 | 3,395 |
| D. rerio | `4285795ae90599eec0163067a55ae965` | `0ede4b1f5c9d00b9287a90929858500f` | `Danio_rerio.manifest.tsv.gz` | `ae2dbd80d2506592d88493523bd3e0c5` | 43 | 3,154 |
| X. tropicalis | `620cd662848acfe02e262a858a4c1a19` | `0e5238f0b14d476f43f2f135f2640717` | `Xenopus_tropicalis.manifest.tsv.gz` | `fc6fc186fe20fa3681208957391fc1e0` | 27 | 1,943 |
| D. melanogaster | `1de22f17786c44ae98d7922116967b4e` | `869cf40e4f5c7ca27d04ca9ae7baee4f` | `Drosophila_melanogaster.manifest.tsv.gz` | `6e01677e2e4f2515fa1810eccc499730` | 8 | 617 |
| C. elegans | `1656bccd184bf57434724f35bea89738` | `177a91ec15063e305188306b7e709cb1` | `Caenorhabditis_elegans.manifest.tsv.gz` | `bf56d4b29ead0d28afd0f1dfe5b01b53` | 9 | 744 |
| A. thaliana | `b809ab966e03e2a2d438b78656d9b409` | `0fed2d7901bf1488f20e7ef6842cd84c` | `Arabidopsis_thaliana.manifest.tsv.gz` | `847661db65845733e3c1f67dc3859b81` | 12 | 1,145 |
| Z. mays | `d29888159bd6df944e622b825ac658fd` | `3fe51a7a675eb0ac410814a41d50c775` | `Zea_mays.manifest.tsv.gz` | `0228ad155220d48f75bd8183983fde33` | 32 | 1,829 |
| S. cerevisiae | `1fddbd976c4ce61e8a36a3908d46c25d` | `88c38b957b721dfc50e6c4df03b6242e` | `Saccharomyces_cerevisiae.manifest.tsv.gz` | `e283de282ded5fff0a284f6ea69438c9` | 1 | 74 |
| N. crassa | `511c0c3c3875ca39bc9ba298e05d97eb` | `dbd205b23f3e95461207f46c7cec6314` | `Neurospora_crassa.manifest.tsv.gz` | `28adfc914eb25b33268fbae85c83211a` | 3 | 143 |
| D. discoideum | `074ff2892412e4068538cf5df8dba350` | `45075a089e27703fdff5e574d3750525` | `Dictyostelium_discoideum.manifest.tsv.gz` | `0e50e8f1cd1de59b27213dfea8d8510b` | 3 | 160 |

Reproduce: `python3 benchmark/fetch.py --species <sp> --what gff,fasta --dest <dir>`,
then `python3 -m model.labels.admission --species <sp> --gff <gff> --fasta <fna> --out-dir
model/labels/manifests`. The summary JSON next to each manifest carries every count in
this document. Adapter fixtures (complete split codon, nonzero-phase edge partial,
interior partial, internal range tag, short-gap frame adjustment, translation
exception, an N, declaration/sequence disagreements, on both strands) are in
`tests/test_label_admission.py`.
