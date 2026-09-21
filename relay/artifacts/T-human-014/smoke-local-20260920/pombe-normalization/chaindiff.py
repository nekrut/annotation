import sys, collections
def chains(path):
    by = collections.defaultdict(list)
    for line in open(path):
        if line.startswith('#'): continue
        f = line.rstrip('\n').split('\t')
        if f[2] != 'CDS': continue
        attrs = dict(a.split('=',1) for a in f[8].split(';') if '=' in a)
        by[attrs['Parent']].append((f[0], f[6], int(f[3]), int(f[4]), f[7]))
    return collections.Counter(tuple(sorted(v)) for v in by.values())
a, b = chains(sys.argv[1]), chains(sys.argv[2])
print("chains", sum(a.values()), sum(b.values()), "identical", sum((a & b).values()), "only_a", sum((a - b).values()), "only_b", sum((b - a).values()))
