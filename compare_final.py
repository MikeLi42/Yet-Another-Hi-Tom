"""Validate pipeline output against reference example TSV files.

Usage:
    python compare_final.py                                    # compares all
    python compare_final.py 01-genotype.tsv 01-genotype_example.tsv
"""
import sys
import os

VALID_ROWS = {'A','B','C','D','E','F','G','H'}

def parse_genotype(path):
    rows = {r: {} for r in 'ABCDEFGH'}
    with open(path) as f:
        f.readline()
        for line in f:
            line = line.rstrip('\n')
            if not line:
                continue
            parts = line.split('\t')
            if parts[0] not in VALID_ROWS:
                continue
            row = parts[0]
            for i in range(12):
                if i + 1 < len(parts):
                    rows[row][i + 1] = parts[i + 1]
    return rows

def parse_sequence(path):
    samples = {}
    cur, haps = None, []
    with open(path) as f:
        f.readline()
        for line in f:
            line = line.rstrip('\n')
            if not line:
                continue
            parts = line.split('\t')
            if len(parts) == 1 and parts[0].strip():
                if cur and haps:
                    samples[cur] = haps
                cur, haps = parts[0], []
            elif len(parts) >= 9 and parts[0].isdigit():
                haps.append({
                    'rank': int(parts[0]), 'count': int(parts[1]),
                    'ratio': parts[2], 'l_var': parts[3], 'r_var': parts[4],
                })
        if cur and haps:
            samples[cur] = haps
    return samples

def compare_genotypes(ex_path, gen_path, label=''):
    gen_ex = parse_genotype(ex_path)
    gen_new = parse_genotype(gen_path)
    mm = 0
    for row in 'ABCDEFGH':
        for col in range(1, 13):
            if gen_ex[row].get(col) != gen_new[row].get(col):
                mm += 1
    pct = (96 - mm) / 96 * 100
    print(f"  [{label}] Genotype: {96-mm}/96 = {pct:.1f} %")
    if mm:
        for row in 'ABCDEFGH':
            for col in range(1, 13):
                va = gen_ex[row].get(col, '?'); vb = gen_new[row].get(col, '?')
                if va != vb:
                    print(f"    {row}{col:02d}  example={va:<10}  yours={vb}")
    return mm

def compare_sequence_variants(ex_path, gen_path, label=''):
    s_ex = parse_sequence(ex_path)
    s_new = parse_sequence(gen_path)
    mm = 0
    for sid in sorted(set(list(s_ex) + list(s_new))):
        te = sorted(set((h['l_var'], h['r_var']) for h in s_ex.get(sid, [])))
        tn = sorted(set((h['l_var'], h['r_var']) for h in s_new.get(sid, [])))
        if te != tn:
            mm += 1
    pct = (96 - mm) / 96 * 100
    print(f"  [{label}] Seq var types: {96-mm}/96 = {pct:.1f} %")
    return mm

if __name__ == '__main__':
    if len(sys.argv) >= 3:
        g = compare_genotypes(sys.argv[2], sys.argv[1])
        sys.exit(0 if g == 0 else 1)

    pairs = [
        ('01', 'example_results/01-genotype_example.tsv', 'demo/01-genotype.tsv',
              'example_results/01-sequence_example.tsv', 'demo/01-sequence.tsv'),
        ('27', 'example_results/27_genotype_example.tsv', 'demo/27-genotype.tsv',
              'example_results/27_sequence_example.tsv', 'demo/27-sequence.tsv'),
    ]
    for ds, gx, gn, sx, sn in pairs:
        if not os.path.exists(gn):
            print(f"[{ds}] SKIP — generated output not found. Run pipeline first.")
            continue
        compare_genotypes(gx, gn, ds)
        compare_sequence_variants(sx, sn, ds)
        print()
