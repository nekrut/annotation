#!/usr/bin/env python3
"""Summarize a score-final directory: cost from measure_*.json, accuracy from score_*.json,
chain/intron geometry from the predicted GFF3 (plain or .gz). Usage: summarize.py DIR [DIR...]"""
import gzip, json, statistics, sys
from pathlib import Path

def med(xs): return statistics.median(xs) if xs else None
def gff_stats(path):
    op = gzip.open if str(path).endswith('.gz') else open
    cds = {}
    with op(path, 'rt') as f:
        for ln in f:
            if ln.startswith('#'): continue
            p = ln.rstrip('\n').split('\t')
            if len(p) < 9 or p[2] != 'CDS': continue
            attrs = dict(kv.split('=', 1) for kv in p[8].split(';') if '=' in kv)
            cds.setdefault(attrs.get('Parent'), []).append((int(p[3]), int(p[4])))
    spans, introns, nex = [], [], []
    for ex in cds.values():
        ex.sort(); spans.append(ex[-1][1] - ex[0][0] + 1); nex.append(len(ex))
        introns += [b[0] - a[1] - 1 for a, b in zip(ex, ex[1:])]
    return dict(chains=len(cds), span_min=min(spans), span_med=med(spans), span_max=max(spans),
                exons_med=med(nex), exons_max=max(nex), introns=len(introns),
                intron_min=min(introns) if introns else None, intron_med=med(introns), intron_max=max(introns) if introns else None)

for d in map(Path, sys.argv[1:]):
    print(f'== {d}')
    for m in sorted(d.glob('measure_*.json')):
        j = json.load(open(m)); n = m.name[len('measure_'):-5]
        print(f"{j['seqid']}: bases {j['sequence_bases']:,} tiles {j['tiles']} cpu_s {j['cpu_s']:.2f} ({j['cpu_s_per_mb']:.2f}/Mb) "
              f"enc {j['encoder_cpu_s']:.2f} dec {j['decode_cpu_s']:.2f} pre {j['preprocess_cpu_s']:.3f} io {j['io_cpu_s']:.3f} "
              f"out {j['output_cpu_s']:.4f} rss {j['peak_host_rss_gib']:.2f} GiB wall {j['wall_s']:.1f} chains {j['chains']} decoded {j['chains_decoded']} src {j['source']['source_sha256'][:8]} commit {j['commit'][:7]}")
        s = json.load(open(d / f'score_{n}.json'))
        nt, lo, ex, tr, sp, co = s['nucleotide'], s['locus'], s['exon']['all'], s['transcript'], s['splice'], s['codon']
        dn = sp['by_dinucleotide']
        print(f"  ref transcripts {s['reference_transcripts']} loci {lo['reference_loci']} | pred transcripts {s['predicted_transcripts']}")
        print(f"  nt sens/prec/F1/MCC {nt['sensitivity']:.3f}/{nt['precision']:.3f}/{nt['f1']:.3f}/{nt['mcc']:.3f}  pred/ref CDS bp {nt['predicted_cds_bp']:,}/{nt['reference_cds_bp']:,}")
        print(f"  locus TP/FP/FN {lo['tp']}/{lo['fp']}/{lo['fn']} fusion {lo['fusion']} split {lo['split']} | transcript exact {tr['tp']} | exon exact TP/FP/FN {ex['tp']}/{ex['fp']}/{ex['fn']}")
        print(f"  introns by dinuc: " + ', '.join(f"{k} tp {v['tp']} fp {v['fp']} fn {v['fn']}" for k, v in dn.items()) +
              f" | acceptor TP/FP/FN {sp['acceptor']['tp']}/{sp['acceptor']['fp']}/{sp['acceptor']['fn']} donor {sp['donor']['tp']}/{sp['donor']['fp']}/{sp['donor']['fn']}")
        print(f"  start TP/FP/FN {co['start']['tp']}/{co['start']['fp']}/{co['start']['fn']} stop {co['stop']['tp']}/{co['stop']['fp']}/{co['stop']['fn']}")
        g = next(p for p in [d / f'{n}.gff3', d / f'{n}.gff3.gz'] if p.exists())
        print('  gff: ' + ' '.join(f'{k}={v}' for k, v in gff_stats(g).items()))
