#!/usr/bin/env python3
"""Validate the self-contained TM-score implementation against the reference
TM-align code (tmtools wrapper of Zhang lab C++), on a random sample of pairs
spanning all four methods.

Note that TM-align performs a sequence-INDEPENDENT structural alignment and can
therefore only equal or exceed the sequence-dependent TM-score used throughout
this note; a non-positive mean difference is the expected result.
"""
import os
import sys
import csv
import json
import random
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from structlib import parse_ca, tm_score, align_by_resnum, needleman_wunsch
from tmtools import tm_align

METHODS = ['af3', 'boltz', 'boltz2', 'chai']
LAB = {'af3': 'AlphaFold 3', 'boltz': 'Boltz-1', 'boltz2': 'Boltz-2', 'chai': 'Chai-1'}
N_PER_METHOD = 100

rows = [r for r in csv.DictReader(open(os.path.join(HERE, 'manifest.tsv')), delimiter='\t')
        if r['afdb_path']]
random.seed(20260726)

out_rows = []
res = {}
for m in METHODS:
    pool = [r for r in rows if r[f'{m}_path']]
    samp = random.sample(pool, min(N_PER_METHOD, len(pool)))
    mine, offi = [], []
    for r in samp:
        a = parse_ca(r['afdb_path'])
        b = parse_ca(r[f'{m}_path'])
        if a['seq'] == b['seq']:
            ia = np.arange(len(a['seq'])); ib = ia
        else:
            ia, ib = align_by_resnum(a, b)
            if len(ia) < 3 or np.mean([a['seq'][x] == b['seq'][y]
                                       for x, y in zip(ia, ib)]) < 0.98:
                ia, ib = needleman_wunsch(a['seq'], b['seq'])
        tm, _, _ = tm_score(b['xyz'][ib], a['xyz'][ia], L_norm=len(a['seq']))
        ref = tm_align(b['xyz'], a['xyz'], b['seq'], a['seq']).tm_norm_chain2
        mine.append(tm); offi.append(ref)
    mine = np.array(mine); offi = np.array(offi)
    d = mine - offi
    res[m] = {'n': len(mine), 'pearson_r': float(np.corrcoef(mine, offi)[0, 1]),
              'mean_diff': float(d.mean()), 'median_diff': float(np.median(d)),
              'max_abs_diff': float(np.abs(d).max()),
              'frac_within_0.02': float(np.mean(np.abs(d) < 0.02)),
              'frac_within_0.05': float(np.mean(np.abs(d) < 0.05))}
    out_rows.append([LAB[m], len(mine), f"{res[m]['pearson_r']:.4f}",
                     f"{res[m]['mean_diff']:+.4f}", f"{res[m]['median_diff']:+.4f}",
                     f"{res[m]['max_abs_diff']:.3f}",
                     f"{100*res[m]['frac_within_0.02']:.0f}%",
                     f"{100*res[m]['frac_within_0.05']:.0f}%"])
    print(LAB[m], res[m])

with open(os.path.join(HERE, 'tables', 'S14_tm_validation.tsv'), 'w', newline='') as fh:
    w = csv.writer(fh, delimiter='\t')
    w.writerow(['Method', 'n pairs', 'Pearson r vs TM-align',
                'Mean difference', 'Median difference', 'Max |difference|',
                'Within 0.02', 'Within 0.05'])
    w.writerows(out_rows)

sm_path = os.path.join(HERE, 'summary.json')
S = json.load(open(sm_path))
S['tm_validation'] = res
json.dump(S, open(sm_path, 'w'), indent=2, default=float)
print('wrote tables/S14_tm_validation.tsv and updated summary.json')
