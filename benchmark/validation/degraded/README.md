# benchmark/validation/degraded/

The panel-wide **control** run: every one of the 20 reference annotations
scored against a synthetically degraded copy of itself, made by
[`benchmark/degrade.py`](../../degrade.py) at seed 20260909 with
`--drop-transcripts 0.10 --shift-cds 0.10`.

An identity run puts every metric at 1.0 by construction, so it cannot
distinguish a correct scorer from one that drops the same thing from both
sides. This run can: the perturbation is known exactly, so each metric's
direction and rough size is predictable before the run, and a metric that
does not move, or moves the wrong way, is a defect. The write-up and the
four findings are `docs/benchmark.md` §6.

| file | what it is |
|---|---|
| `degraded.yaml` | the §3.3 declaration all 20 runs were scored under |
| `<Species>.json` | `score.py` output |
| `<Species>.degrade.json` | what `degrade.py` actually did to that species: transcripts seen and deleted, CDS segments seen, shifted, split by strand, and refused at a sequence end |

Reproduce one species:

```
python3 benchmark/degrade.py \
    --reference /tmp/panel/Neurospora_crassa/*_genomic.gff.gz \
    --out /tmp/nc.degraded.gff3 --summary /tmp/nc.degrade.json
python3 benchmark/score.py \
    --reference /tmp/panel/Neurospora_crassa/*_genomic.gff.gz \
    --prediction /tmp/nc.degraded.gff3 --species Neurospora_crassa \
    --declaration benchmark/validation/degraded/degraded.yaml \
    --out /tmp/nc.json
```

The `reference` and `out` paths in the `.degrade.json` summaries are the
scratch paths of the machine that produced them; the files themselves were
deleted after scoring and are not in the repository.

These runs are a control, not an accuracy measurement: the input is the
reference with 10% of its transcripts removed, so the numbers say what the
scorer does to a known error, not what any predictor achieves.
