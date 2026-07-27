#!/usr/bin/env python3
"""Worked example: domain-resolved comparison of PknB (ML0016, P54744).

PknB is the highest-ranked prioritised target (rank 50, high-priority tier) that
falls below TM = 0.5 against AFDB for every method, and it illustrates why a
whole-chain score can be misleading for a modular protein. Domain boundaries are
taken from the UniProt feature annotation stored in the HANSEN database.
"""
import os
import sys
import csv
import json
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from structlib import parse_ca, tm_score, kabsch_rt, apply_rt, lddt_ca

METHODS = ['af3', 'boltz', 'boltz2', 'chai']
LAB = {'af3': 'AlphaFold 3', 'boltz': 'Boltz-1', 'boltz2': 'Boltz-2', 'chai': 'Chai-1'}

# UniProt P54744 features: DOMAIN 11..273 kinase; TRANSMEM 329..349;
# DOMAIN 352..418 / 419..486 / 487..553 / 554..622 PASTA 1-4;
# REGION 381..404 disordered
SEGMENTS = [
    ('Kinase domain (11-273)', 11, 273),
    ('Juxtamembrane linker (274-328)', 274, 328),
    ('Transmembrane helix (329-349)', 329, 349),
    ('PASTA 1-4 (352-622)', 352, 622),
    ('Whole chain (1-622)', 1, 10 ** 6),
]

man = {r['ml_id']: r for r in
       csv.DictReader(open(os.path.join(HERE, 'manifest.tsv')), delimiter='\t')}
r = man['ML0016']
a = parse_ca(r['afdb_path'])

rows = []
store = {}
for name, lo, hi in SEGMENTS:
    sel = (a['resnum'] >= lo) & (a['resnum'] <= hi)
    row = [name, int(sel.sum())]
    for m in METHODS:
        b = parse_ca(r[f'{m}_path'])
        P = b['xyz'][sel]
        Q = a['xyz'][sel]
        tm, _, _ = tm_score(P, Q, L_norm=int(sel.sum()))
        R, t = kabsch_rt(P, Q)
        rms = float(np.sqrt((np.linalg.norm(apply_rt(P, R, t) - Q, axis=1) ** 2).mean()))
        ld, _ = lddt_ca(P, Q)
        row += [f'{tm:.2f}', f'{rms:.1f}', f'{ld:.3f}']
        store.setdefault(name, {})[m] = {'tm': tm, 'rmsd': rms, 'lddt': ld,
                                         'n': int(sel.sum())}
    row.append(f"{np.nanmean(a['plddt'][sel]):.1f}")
    rows.append(row)

hdr = ['Segment', 'Residues']
for m in METHODS:
    hdr += [f'TM {LAB[m]}', f'RMSD {LAB[m]} (A)', f'lDDT {LAB[m]}']
hdr.append('AFDB pLDDT')

out = os.path.join(HERE, 'tables', 'S16_pknb_domains.tsv')
with open(out, 'w', newline='') as fh:
    w = csv.writer(fh, delimiter='\t')
    w.writerow(hdr)
    w.writerows(rows)
print('wrote', out)
for row in rows:
    print('  ', row)

sm = os.path.join(HERE, 'summary.json')
S = json.load(open(sm))
S['pknb_worked_example'] = store
json.dump(S, open(sm, 'w'), indent=2, default=float)
print('updated summary.json')
