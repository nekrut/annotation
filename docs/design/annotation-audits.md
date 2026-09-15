# Annotation audit reproduction

Companion to [the design proposal](proposal.md), sections 4.3 and 4.4.
Run these read-only listings from the repository root with the checksummed
train GFFs at the paths named in the proposal. They are annotation audits,
not a model or structured-training implementation. The proposal records
the expected counts, exclusions and source hashes.

<a id="support-and-topology"></a>

## Support and topology (proposal section 4.3)

```python
from collections import defaultdict
from pathlib import Path
import csv, hashlib, json, runpy

score = runpy.run_path("benchmark/score.py")
panel = {r["species"]: r for r in csv.DictReader(open("benchmark/panel.tsv"), delimiter="\t")}
species = ["Saccharomyces_cerevisiae", "Caenorhabditis_elegans",
           "Drosophila_melanogaster", "Mus_musculus"]
print("score_sha256", hashlib.sha256(Path("benchmark/score.py").read_bytes()).hexdigest())
print("panel_sha256", hashlib.sha256(Path("benchmark/panel.tsv").read_bytes()).hexdigest())
for sp in species:
    row = panel[sp]
    assert row["split"] == "train"
    path = Path("/tmp/bench007/data") / sp / (row["ftp_dir"] + "_genomic.gff.gz")
    md5 = hashlib.md5(path.read_bytes()).hexdigest()
    assert md5 == row["md5_gff"], (sp, md5)
    ann = score["load_gff"](str(path))
    seqids = score["select_seqids"](ann)
    chains = score["select_transcripts"](
        {t: c for t, c in ann.chains().items() if c[0] in seqids}, ann)
    unique = set(chains.values())
    N = sum(ann.seq_len[s] for s in seqids)
    oriented, physical, spans = defaultdict(list), defaultdict(list), defaultdict(list)
    tiles = {"+": set(), "-": set()}
    for sid, strand, blocks in unique:
        assert strand in tiles
        L = ann.seq_len[sid]
        spans[sid, strand].append((min(a for a, b in blocks), max(b for a, b in blocks)))
        for lo, hi in blocks:
            assert 1 <= lo <= hi <= L
            oriented[sid, strand].append((lo, hi))
            physical[sid].append((lo, hi))
            a, b = (lo-1, hi) if strand == "+" else (L-hi, L-lo+1)
            tiles[strand].update((sid, t) for t in range(a//384, (b-1)//384+1))
    def union_bp(groups):
        return sum(e-s+1 for iv in groups.values()
                   for s, e in score["merge_intervals"](iv))
    max_pack = 0
    for intervals in spans.values():
        last = 0
        for lo, hi in sorted(intervals, key=lambda v: (v[1], v[0])):
            if lo > last:
                max_pack += 1
                last = hi
    f = 384 * sum(map(len, tiles.values())) / (2*N)
    result = dict(species=sp, md5_gff=md5, scored_sequences=len(seqids), scored_bp=N,
                  transcripts=len(chains), unique_CDS_chains=len(unique),
                  physical_CDS_union_bp=union_bp(physical),
                  oriented_CDS_union_bp=union_bp(oriented),
                  plus_tiles=len(tiles["+"]), minus_tiles=len(tiles["-"]),
                  oracle_CDS_allocated_f=f, B8_conditional_CPU_s_per_Mb=7.71535+1914.592*f/30,
                  max_nonoverlapping_unique_chains=max_pack)
    print(json.dumps(result, sort_keys=True), flush=True)
```

<a id="representatives-metadata-and-quota"></a>

## Representatives, metadata and quota (proposal section 4.4)

```python
from collections import Counter, defaultdict
from pathlib import Path
import bisect, csv, gzip, hashlib, json, runpy, sys
sys.dont_write_bytecode = True
score = runpy.run_path("benchmark/score.py")
cut = runpy.run_path("scripts/data/cut_windows.py")
panel_path = Path("benchmark/panel.tsv")
panel = {r["species"]: r for r in csv.DictReader(panel_path.open(), delimiter="\t")}
for source in ("benchmark/score.py", "scripts/data/cut_windows.py", "benchmark/panel.tsv"):
    print(json.dumps({"source": source, "sha256": hashlib.sha256(Path(source).read_bytes()).hexdigest()}), flush=True)

def components(txs):
    # Zero-based, half-open FULL CDS spans, grouped by sequence AND strand.
    by = defaultdict(list)
    for t in txs:
        by[t["seqid"], t["strand"]].append(t)
    masked, masks = set(), defaultdict(list)
    for key, group in by.items():
        members, lo, hi = [], None, None
        for t in sorted(group, key=lambda t: (t["cds"][0][0], t["cds"][-1][1], t["id"])):
            a, b = t["cds"][0][0], t["cds"][-1][1]
            if members and a >= hi:
                if len(members) > 1:
                    masked.update(members); masks[key].append((lo, hi))
                members, lo, hi = [], None, None
            if not members:
                lo, hi = a, b
            members.append(t["id"]); hi = max(hi, b)
        if len(members) > 1:
            masked.update(members); masks[key].append((lo, hi))
    return masked, masks

def sites(txs):
    # CDS-chain boundary convention, including short-gap exclusions.
    out = set()
    for t in txs:
        cds_only = dict(t, exons=t["cds"])
        out.update((t["seqid"], kind, strand, pos)
                   for kind, strand, pos in cut["distinct_sites"]([cds_only])
                   if kind in ("start", "stop", "donor", "acceptor"))
    return out

def count_sites(ss):
    c = Counter(kind for sid, kind, strand, pos in ss)
    return {kind: c[kind] for kind in ("start", "stop", "donor", "acceptor")}

def load_txs(path, ann, chains):
    # Keep scorer gene/sequence/strand namespacing; supplement full exon lengths
    # and first appearance in GFF order for the cutter's real tie breakers.
    extras = {}
    with gzip.open(path, "rt") as fh:
        for line_no, line in enumerate(fh):
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) != 9:
                continue
            raw = (score["_attr"](f[8], "Parent") or score["_attr"](f[8], "transcript_id")) if f[2] in ("exon", "CDS") else (score["_attr"](f[8], "ID") or score["_attr"](f[8], "transcript_id"))
            if not raw:
                continue
            ids = raw.split(",") if f[2] == "exon" else [raw.split(",")[0]]
            for raw_id in ids:
                tid = "\x00".join((f[0], f[6], raw_id))
                if tid not in chains:
                    continue
                meta = extras.setdefault(tid, {"order": line_no, "exons": [],
                                               "raw_cds": [], "tags": set()})
                if f[2] == "exon":
                    meta["exons"].append((int(f[3])-1, int(f[4])))
                if f[2] == "CDS":
                    meta["raw_cds"].append((int(f[3])-1, int(f[4]), f[7], f[8]))
                if f[2] == "CDS" or f[2] in score["TX_TYPES"]:
                    for key in ("exception", "transl_except"):
                        if score["_attr"](f[8], key) is not None:
                            meta["tags"].add(key)
    assert set(extras) == set(chains), "missing GFF transcript order"
    txs = []
    for tid in sorted(chains, key=lambda tid: extras[tid]["order"]):
        sid, strand, blocks = chains[tid]
        cds = [(a-1, b) for a, b in blocks]
        exons = sorted(extras[tid]["exons"])
        assert exons, ("missing full exons for longest-cds tie breaker", tid)
        txs.append(dict(id=tid, gene=ann.gene_of[tid], seqid=sid, strand=strand,
                        start=exons[0][0], end=exons[-1][1], cds=cds, exons=exons,
                        raw_cds=extras[tid]["raw_cds"], tags=extras[tid]["tags"],
                        _partial5=tid in ann.partial5, _partial3=tid in ann.partial3))
    return txs

def metadata_audit(txs, lengths):
    counts, flagged = Counter(), 0
    for t in txs:
        flags = set(t["tags"])
        rr = sorted(t["raw_cds"], reverse=t["strand"] == "-")
        assert rr, ("missing raw CDS rows", t["id"])
        L = lengths[t["seqid"]]
        iv = [(L-b, L-a) if t["strand"] == "-" else (a,b)
              for a,b,phase,attrs in rr]
        if t["_partial5"] or t["_partial3"]:
            flags.add("declared_partial")
        if ((t["_partial5"] and iv[0][0] != 0) or
                (t["_partial3"] and iv[-1][1] != L)):
            flags.add("nonedge_partial")
        for i, (a,b,phase,attrs) in enumerate(rr):
            sr = score["_attr"](attrs, "start_range") is not None
            er = score["_attr"](attrs, "end_range") is not None
            if t["strand"] == "-":
                sr, er = er, sr
            if (sr and i != 0) or (er and i != len(rr)-1):
                flags.add("internal_range")
        phases = [int(r[2]) if r[2] in ("0","1","2") else None for r in rr]
        if None in phases:
            flags.add("invalid_phase")
        else:
            p = (-phases[0]) % 3
            if p and not t["_partial5"]:
                flags.add("nonzero_initial_phase_without_partial5")
            for (a,b), phase in zip(iv, phases):
                if phase != (-p) % 3:
                    flags.add("phase_discontinuity")
                p = (p+b-a) % 3
            if not t["_partial3"] and p:
                flags.add("incomplete_terminal_frame")
        for (_,b), (a,_) in zip(iv, iv[1:]):
            if a < b:
                flags.add("overlapping_rows")
            if 0 < a-b < 20:
                flags.add("short_gap")
        flagged += bool(flags - {"declared_partial"})
        counts.update(flags)
    return dict(flag_counts=dict(sorted(counts.items())), flag_union=flagged,
                no_metadata_flag=len(txs)-flagged)

def tile_audit(txs, lengths, flank):
    tiles = defaultdict(set)
    for t in txs:
        sid, strand = t["seqid"], t["strand"]
        L = lengths[sid]
        for a, b in t["cds"]:
            a, b = max(0, a-flank), min(L, b+flank)
            if strand == "-":
                a, b = L-b, L-a
            tiles[sid, strand].update(range(a//384, (b-1)//384+1))
    served_local = served_carried = served_pooled = slots = groups_exceeding = 0
    for sid, L in lengths.items():
        for strand in ("+", "-"):
            counts = Counter(t//256 for t in tiles[sid, strand])
            credit = local = carried = qtotal = 0
            for g, lo in enumerate(range(0, L, 98304)):
                n = min(98304, L-lo)
                q = n//7680  # exact floor(0.05*n/384)
                demand = counts[g]
                qtotal += q; local += min(q, demand)
                credit += q
                take = min(credit, demand)
                carried += take; credit -= take
                assert carried <= qtotal
                groups_exceeding += int(demand > q)
            served_local += local; served_carried += carried
            served_pooled += min(qtotal, len(tiles[sid, strand]))
            slots += qtotal
    total = sum(map(len, tiles.values()))
    return dict(flank=flank, tiles=total, oracle_f=384*total/(2*sum(lengths.values())),
                slots=slots, served_local=served_local, served_carried=served_carried,
                served_pooled=served_pooled, groups_over_local_quota=groups_exceeding)

for sp in ("Saccharomyces_cerevisiae", "Caenorhabditis_elegans",
           "Drosophila_melanogaster", "Mus_musculus"):
    row = panel[sp]
    assert row["split"] == "train"
    path = Path("/tmp/bench007/data") / sp / (row["ftp_dir"] + "_genomic.gff.gz")
    md5 = hashlib.md5(path.read_bytes()).hexdigest()
    assert md5 == row["md5_gff"], (sp, md5)
    ann = score["load_gff"](str(path))
    seqids = score["select_seqids"](ann)
    chains = score["select_transcripts"]({t:c for t,c in ann.chains().items() if c[0] in seqids}, ann)
    txs = load_txs(path, ann, chains)
    selected, dropped, fallback = cut["select_isoforms"](txs, "longest-cds")
    assert not fallback
    assert len(selected) == len(score["loci_of"](ann, chains, seqids))
    masked, masks = components(selected)
    retained = [t for t in selected if t["id"] not in masked]
    assert not components(retained)[0]
    all_sites, selected_sites, retained_sites = sites(txs), sites(selected), sites(retained)
    assert retained_sites <= selected_sites <= all_sites
    lengths = {sid:ann.seq_len[sid] for sid in seqids}
    N = sum(lengths.values())
    cds_bp_masked = sum(sum(b-a for a,b in t["cds"]) for t in selected if t["id"] in masked)
    selected_cds_bp = sum(sum(b-a for a,b in t["cds"]) for t in selected)
    starts = {key:[a for a,b in iv] for key,iv in masks.items()}
    def in_mask(site):
        sid, kind, strand, pos = site
        key = sid, "+" if strand == 1 else "-"
        i = bisect.bisect_right(starts.get(key, []), pos)-1
        return i >= 0 and pos < masks[key][i][1]
    result = dict(species=sp, md5_gff=md5, transcripts=len(txs), loci=len(selected),
                  nonrepresentative_transcripts=len(dropped),
                  overlapping_components=sum(map(len,masks.values())),
                  masked_representatives=len(masked), retained_representatives=len(retained),
                  mask_oriented_span_bp=sum(b-a for iv in masks.values() for a,b in iv),
                  scored_oriented_bp=2*N, selected_cds_bp_sum=selected_cds_bp,
                  masked_selected_cds_bp_sum=cds_bp_masked,
                  all_sites=count_sites(all_sites), selected_sites=count_sites(selected_sites),
                  retained_sites=count_sites(retained_sites),
                  auxiliary_sites_in_mask=count_sites({s for s in all_sites if in_mask(s)}),
                  selected_partial_transcripts=sum(t["_partial5"] or t["_partial3"] for t in selected),
                  retained_partial_transcripts=sum(t["_partial5"] or t["_partial3"] for t in retained),
                  raw_cds_metadata=metadata_audit(retained, lengths),
                  quota=[tile_audit(txs,lengths,flank) for flank in (0,12)])
    print(json.dumps(result,sort_keys=True),flush=True)
```
