#!/usr/bin/env python3
"""Figure S6: per-residue deviation profiles for the worked examples.

Shows, for each of four instructive proteins, how far each HANSEN model sits from
the AFDB reference at every residue after TM-superposition, with the AFDB
per-residue pLDDT underneath. These make the difference between a hinge and a
genuinely wrong fold visible at a glance.
"""
import os
import sys
import csv
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from structlib import parse_ca, tm_score, apply_rt, align_by_resnum

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402

METHODS = ['af3', 'boltz', 'boltz2', 'chai']
LABEL = {'af3': 'AlphaFold 3', 'boltz': 'Boltz-1', 'boltz2': 'Boltz-2', 'chai': 'Chai-1'}
COLOR = {'af3': '#1f77b4', 'boltz': '#ff7f0e', 'boltz2': '#2ca02c', 'chai': '#d62728'}

plt.rcParams.update({'figure.dpi': 150, 'savefig.dpi': 300, 'font.size': 8.5,
                     'axes.spines.top': False, 'axes.spines.right': False,
                     'axes.grid': True, 'grid.alpha': 0.25, 'legend.frameon': False})

CASES = [
    ('ML0016', 'PknB (ML0016), Kinase and PASTA Domains With a Disordered '
     'Juxtamembrane Linker',
     [(11, 273, 'Kinase'), (274, 328, 'Linker'), (329, 349, 'TM Helix'),
      (352, 622, 'PASTA')]),
    ('ML1629', 'Smc (ML1629), Antiparallel Coiled-Coil Hinge', []),
    ('ML2594', 'Mce1F (ML2594), MCE Domain and Helical Extension', []),
    ('ML2593', 'LprK (ML2593), mce1 Operon Lipoprotein', []),
]

man = {r['ml_id']: r for r in
       csv.DictReader(open(os.path.join(HERE, 'manifest.tsv')), delimiter='\t')}

fig, axes = plt.subplots(len(CASES), 1, figsize=(9, 2.55 * len(CASES)))
for ax, (ml, title, segs) in zip(axes, CASES):
    r = man[ml]
    a = parse_ca(r['afdb_path'])
    for m in METHODS:
        if not r.get(f'{m}_path'):
            continue
        b = parse_ca(r[f'{m}_path'])
        if a['seq'] == b['seq']:
            ia = np.arange(len(a['seq'])); ib = ia
        else:
            ia, ib = align_by_resnum(a, b)
        P, Q = b['xyz'][ib], a['xyz'][ia]
        tm, R, t = tm_score(P, Q, L_norm=len(a['seq']))
        d = np.linalg.norm(apply_rt(P, R, t) - Q, axis=1)
        ax.plot(a['resnum'][ia], np.maximum(d, 0.05), color=COLOR[m], lw=0.9,
                label=f'{LABEL[m]} (TM {tm:.2f})')
    ax.set_yscale('log')
    ax.set_ylim(0.1, 4000)
    ax.axhline(5, color='k', ls=':', lw=0.7)
    ax.set_ylabel('C\u03b1 Deviation (\u00c5)')
    ax.set_title(title, loc='left', fontweight='bold', fontsize=8.5)
    ax.legend(fontsize=6.8, ncol=4, loc='upper left', bbox_to_anchor=(0.0, 1.0))

    ax2 = ax.twinx()
    ax2.fill_between(a['resnum'], a['plddt'], color='#999999', alpha=0.18, lw=0)
    ax2.set_ylim(0, 100)
    ax2.set_ylabel('AFDB pLDDT', color='#666666', fontsize=7.5)
    ax2.tick_params(labelsize=7, colors='#666666')
    ax2.grid(False)
    for lo, hi, name in segs:
        ax.axvspan(lo, hi, color='#000000', alpha=0.04, lw=0)
        ax.text((lo + hi) / 2, 900, name, ha='center', fontsize=6.5, color='#555555')
axes[-1].set_xlabel('Residue Number (UniProt Numbering)')
fig.tight_layout()
for ext in ('png', 'pdf'):
    fig.savefig(os.path.join(HERE, 'figures', f'FigS6_deviation_profiles.{ext}'),
                bbox_inches='tight')
print('wrote figures/FigS6_deviation_profiles.{png,pdf}')
