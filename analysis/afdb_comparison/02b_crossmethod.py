#!/usr/bin/env python3
"""Predictor-vs-predictor concordance, computed with the same code path as the
AFDB comparison so that a single 5x5 agreement matrix (AFDB, AF3, Boltz-1,
Boltz-2, Chai-1) can be reported on one scale.

Writes crossmethod.tsv (one row per protein per unordered method pair).
"""
import os
import sys
import csv
import itertools
import numpy as np
from multiprocessing import Pool

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from structlib import (parse_ca, align_by_resnum, needleman_wunsch, kabsch_rt,
                       apply_rt, tm_score, lddt_ca, pearson)

MANIFEST = os.path.join(HERE, 'manifest.tsv')
OUT = os.path.join(HERE, 'crossmethod.tsv')
METHODS = ['afdb', 'af3', 'boltz', 'boltz2', 'chai']
FIELDS = ['ml_id', 'uniprot', 'gene', 'method_a', 'method_b', 'n_common',
          'tm_score', 'rmsd_global', 'lddt', 'plddt_r']


def path_for(row, m):
    return row['afdb_path'] if m == 'afdb' else row.get(f'{m}_path')


def corresp(a, b):
    if a['seq'] == b['seq']:
        return np.arange(len(a['seq'])), np.arange(len(b['seq']))
    ia, ib = align_by_resnum(a, b)
    if len(ia) >= 3:
        same = np.mean([a['seq'][x] == b['seq'][y] for x, y in zip(ia, ib)])
        if same >= 0.98 and len(ia) >= 0.9 * min(len(a['seq']), len(b['seq'])):
            return ia, ib
    return needleman_wunsch(a['seq'], b['seq'])


def one(row):
    st = {}
    for m in METHODS:
        p = path_for(row, m)
        if not p:
            continue
        try:
            st[m] = parse_ca(p)
        except Exception:
            pass
    out = []
    for a, b in itertools.combinations(METHODS, 2):
        if a not in st or b not in st:
            continue
        ia, ib = corresp(st[a], st[b])
        if len(ia) < 3:
            continue
        # normalise by the shorter chain, symmetric convention
        Q = st[a]['xyz'][ia]
        P = st[b]['xyz'][ib]
        Ln = min(len(st[a]['seq']), len(st[b]['seq']))
        tm, R, t = tm_score(P, Q, L_norm=Ln)
        R2, t2 = kabsch_rt(P, Q)
        rms = float(np.sqrt((np.linalg.norm(apply_rt(P, R2, t2) - Q, axis=1) ** 2).mean()))
        ld, _ = lddt_ca(P, Q)
        pr = pearson(st[a]['plddt'][ia], st[b]['plddt'][ib])
        out.append({'ml_id': row['ml_id'], 'uniprot': row['uniprot'],
                    'gene': row['gene'], 'method_a': a, 'method_b': b,
                    'n_common': len(ia), 'tm_score': round(tm, 4),
                    'rmsd_global': round(rms, 3), 'lddt': round(ld, 4),
                    'plddt_r': round(pr, 4) if pr == pr else ''})
    return out


def main():
    with open(MANIFEST) as fh:
        rows = list(csv.DictReader(fh, delimiter='\t'))
    print(f'{len(rows)} proteins')
    res = []
    with Pool(max(1, (os.cpu_count() or 4) - 1)) as pool:
        for i, r in enumerate(pool.imap_unordered(one, rows, chunksize=4), 1):
            res.extend(r)
            if i % 200 == 0:
                print(f'  {i}/{len(rows)}', flush=True)
    res.sort(key=lambda r: (r['ml_id'], r['method_a'], r['method_b']))
    with open(OUT, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDS, delimiter='\t')
        w.writeheader()
        w.writerows(res)
    print(f'wrote {OUT} ({len(res)} rows)')


if __name__ == '__main__':
    main()
