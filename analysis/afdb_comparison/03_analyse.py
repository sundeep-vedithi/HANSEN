#!/usr/bin/env python3
"""Aggregate the AFDB-vs-HANSEN comparison into report tables + figures."""
import os
import sys
import csv
import json
import sqlite3
import collections
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from structlib import pearson, spearman  # noqa: E402

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.ticker import PercentFormatter  # noqa: E402

from paths import DB, PRIOR_TSV, FIGDIR, TABDIR  # noqa: E402
os.makedirs(FIGDIR, exist_ok=True)
os.makedirs(TABDIR, exist_ok=True)

METHODS = ['af3', 'boltz', 'boltz2', 'chai']
LABEL = {'af3': 'AlphaFold 3', 'boltz': 'Boltz-1', 'boltz2': 'Boltz-2',
         'chai': 'Chai-1', 'afdb': 'AlphaFold DB (v6)'}
COLOR = {'af3': '#1f77b4', 'boltz': '#ff7f0e', 'boltz2': '#2ca02c',
         'chai': '#d62728', 'afdb': '#7f7f7f'}

plt.rcParams.update({
    'figure.dpi': 150, 'savefig.dpi': 300, 'font.size': 9,
    'axes.spines.top': False, 'axes.spines.right': False,
    'axes.grid': True, 'grid.alpha': 0.25, 'grid.linewidth': 0.5,
    'legend.frameon': False, 'font.family': 'sans-serif',
})


def fnum(v):
    try:
        x = float(v)
        return x if np.isfinite(x) else np.nan
    except (TypeError, ValueError):
        return np.nan


def savefig(fig, name):
    for ext in ('png', 'pdf'):
        fig.savefig(os.path.join(FIGDIR, f'{name}.{ext}'), bbox_inches='tight')
    plt.close(fig)
    print('  figure:', name)


def write_tsv(path, header, rows):
    with open(path, 'w', newline='') as fh:
        w = csv.writer(fh, delimiter='\t')
        w.writerow(header)
        w.writerows(rows)
    print('  table:', os.path.basename(path))


def qstats(a):
    a = np.asarray(a, dtype=float)
    a = a[np.isfinite(a)]
    if a.size == 0:
        return dict(n=0)
    return dict(n=int(a.size), mean=float(a.mean()), sd=float(a.std(ddof=1)) if a.size > 1 else 0.0,
                median=float(np.median(a)), q1=float(np.percentile(a, 25)),
                q3=float(np.percentile(a, 75)), p5=float(np.percentile(a, 5)),
                p95=float(np.percentile(a, 95)), min=float(a.min()), max=float(a.max()))


# --------------------------------------------------------------- load inputs
pairs = list(csv.DictReader(open(os.path.join(HERE, 'pairs.tsv')), delimiter='\t'))
manifest = {r['ml_id']: r for r in
            csv.DictReader(open(os.path.join(HERE, 'manifest.tsv')), delimiter='\t')}
cross = []
cx_path = os.path.join(HERE, 'crossmethod.tsv')
if os.path.exists(cx_path):
    cross = list(csv.DictReader(open(cx_path), delimiter='\t'))

con = sqlite3.connect(f'file:{DB}?mode=ro', uri=True)
con.row_factory = sqlite3.Row
ess = {r['ml_id']: dict(r) for r in con.execute(
    'SELECT ml_id, proteomelm_ess_probability, proteomelm_ess_class, '
    'target_priority_score, priority_tier FROM hansen_essentiality_predictions')}
con.close()

prior = {}
if os.path.exists(PRIOR_TSV):
    for r in csv.DictReader(open(PRIOR_TSV), delimiter='\t'):
        prior[r['ml_id']] = {'rank': int(r['rank']),
                             'tier': r['priority_tier'],
                             'score': fnum(r['target_priority_score'])}

byml = collections.defaultdict(dict)
for r in pairs:
    byml[r['ml_id']][r['method']] = r

NUMCOLS = ['tm_score', 'rmsd_global', 'rmsd_core', 'gdt_ts', 'gdt_ha', 'lddt',
           'frac_lt2', 'frac_lt5', 'dev_p95', 'tm_core', 'lddt_core',
           'plddt_ref', 'plddt_mod', 'plddt_r', 'plddt_rho', 'rg_ratio',
           'ss_q3', 'seq_pid', 'helix_ref', 'strand_ref', 'helix_mod',
           'strand_mod']
V = {m: {c: np.array([fnum(r[c]) for r in pairs if r['method'] == m])
         for c in NUMCOLS} for m in METHODS}
NREF = {m: np.array([fnum(r['n_ref']) for r in pairs if r['method'] == m]) for m in METHODS}
MLS = {m: [r['ml_id'] for r in pairs if r['method'] == m] for m in METHODS}
# chain coverage: fraction of the AFDB reference that the model actually spans.
# TM-score is normalised by the AFDB length, so a truncated model is capped at
# roughly its coverage and must not be read as a fold difference.
COV = {m: np.array([fnum(r['n_common']) / fnum(r['n_ref'])
                    for r in pairs if r['method'] == m]) for m in METHODS}
FULLCOV = {m: COV[m] >= 0.95 for m in METHODS}

summary = {}

# ============================================================ S1 coverage
print('\n[1] coverage / correspondence')
n_afdb = sum(1 for r in manifest.values() if r['afdb_path'])
cov_rows = []
for m in METHODS:
    n = len(V[m]['tm_score'])
    cov_rows.append([LABEL[m], n, n_afdb, f'{100*n/n_afdb:.1f}%'])
write_tsv(os.path.join(TABDIR, 'S1_coverage.tsv'),
          ['Method', 'Comparisons', 'AFDB references', 'Coverage'], cov_rows)

seq_ident = collections.Counter(r['align_mode'] for r in pairs)
len_mismatch = []
for ml, d in byml.items():
    r = next(iter(d.values()))
    if int(r['n_mod']) != int(r['n_ref']):
        mrow = manifest[ml]
        len_mismatch.append([ml, r['uniprot'], r['gene'], int(r['n_ref']),
                             int(r['n_mod']), int(r['n_mod']) - int(r['n_ref']),
                             mrow['protein_name'][:80]])
len_mismatch.sort(key=lambda x: -abs(x[5]))
write_tsv(os.path.join(TABDIR, 'S2_construct_length_differences.tsv'),
          ['ML ID', 'UniProt', 'Gene', 'AFDB residues', 'HANSEN residues',
           'Difference', 'Protein'], len_mismatch)
summary['coverage'] = {'afdb_entries': 1602, 'afdb_paired': n_afdb,
                       'per_method': {m: len(V[m]['tm_score']) for m in METHODS},
                       'align_modes': dict(seq_ident),
                       'length_mismatch_proteins': len(len_mismatch)}

# ======================================================== S3 global agreement
print('\n[2] global agreement')
glob_rows = []
for m in METHODS:
    tm = V[m]['tm_score']
    row = [LABEL[m], len(tm)]
    for c, dec in (('tm_score', 3), ('gdt_ts', 3), ('gdt_ha', 3), ('lddt', 3),
                   ('rmsd_global', 2), ('rmsd_core', 2)):
        s = qstats(V[m][c])
        row += [f"{s['mean']:.{dec}f}", f"{s['median']:.{dec}f}",
                f"{s['q1']:.{dec}f}-{s['q3']:.{dec}f}"]
    for thr in (0.5, 0.7, 0.9):
        row.append(f'{100*np.nanmean(tm>=thr):.1f}%')
    glob_rows.append(row)
hdr = ['Method', 'n']
for c in ('TM-score', 'GDT-TS', 'GDT-HA', 'lDDT', 'RMSD (A)', 'Core RMSD (A)'):
    hdr += [f'{c} mean', f'{c} median', f'{c} IQR']
hdr += ['TM>=0.5', 'TM>=0.7', 'TM>=0.9']
write_tsv(os.path.join(TABDIR, 'S3_global_agreement.tsv'), hdr, glob_rows)
summary['global'] = {m: {c: qstats(V[m][c]) for c in
                         ('tm_score', 'gdt_ts', 'gdt_ha', 'lddt', 'rmsd_global',
                          'rmsd_core', 'frac_lt2', 'frac_lt5', 'dev_p95')}
                     for m in METHODS}
summary['global_thresholds'] = {
    m: {f'tm_ge_{t}': float(np.nanmean(V[m]['tm_score'] >= t))
        for t in (0.5, 0.7, 0.8, 0.9, 0.95)} for m in METHODS}

# Figure 1: TM-score distributions
fig, axes = plt.subplots(1, 3, figsize=(11, 3.2))
ax = axes[0]
bins = np.linspace(0, 1, 51)
for m in METHODS:
    ax.hist(V[m]['tm_score'], bins=bins, histtype='step', lw=1.6,
            color=COLOR[m], label=LABEL[m])
ax.set_xlabel('TM-Score vs AFDB'); ax.set_ylabel('Proteins')
ax.legend(fontsize=7.5); ax.set_title('a  TM-Score Distribution', loc='left', fontweight='bold')
ax = axes[1]
for m in METHODS:
    x = np.sort(V[m]['tm_score'][np.isfinite(V[m]['tm_score'])])
    ax.plot(x, np.linspace(0, 1, len(x)), color=COLOR[m], lw=1.6, label=LABEL[m])
ax.axvline(0.5, ls=':', c='k', lw=0.8)
ax.set_xlabel('TM-Score vs AFDB'); ax.set_ylabel('Cumulative Fraction')
ax.yaxis.set_major_formatter(PercentFormatter(1.0))
ax.set_title('b  Cumulative Distribution', loc='left', fontweight='bold')
ax = axes[2]
data = [V[m]['lddt'][np.isfinite(V[m]['lddt'])] for m in METHODS]
bp = ax.boxplot(data, tick_labels=[LABEL[m] for m in METHODS], showfliers=False,
                patch_artist=True, widths=0.6)
for p, m in zip(bp['boxes'], METHODS):
    p.set_facecolor(COLOR[m]); p.set_alpha(0.55)
for med in bp['medians']:
    med.set_color('k')
ax.set_ylabel('C\u03b1-lDDT vs AFDB')
ax.tick_params(axis='x', rotation=20)
ax.set_title('c  Local Agreement', loc='left', fontweight='bold')
fig.tight_layout()
savefig(fig, 'FigS1_global_agreement')

# ================================================= S4 confidence relationship
print('\n[3] confidence')
conf_rows = []
for m in METHODS:
    s_ref = qstats(V[m]['plddt_ref']); s_mod = qstats(V[m]['plddt_mod'])
    r_tm = pearson(V[m]['plddt_mod'], V[m]['tm_score'])
    r_reftm = pearson(V[m]['plddt_ref'], V[m]['tm_score'])
    conf_rows.append([LABEL[m], f"{s_ref['mean']:.1f}", f"{s_mod['mean']:.1f}",
                      f"{s_mod['mean']-s_ref['mean']:+.1f}",
                      f"{qstats(V[m]['plddt_r'])['median']:.3f}",
                      f"{r_tm:.3f}", f"{r_reftm:.3f}",
                      f"{spearman(np.minimum(V[m]['plddt_ref'],V[m]['plddt_mod']), V[m]['tm_score']):.3f}"])
write_tsv(os.path.join(TABDIR, 'S4_confidence.tsv'),
          ['Method', 'AFDB mean pLDDT', 'Model mean pLDDT', 'Delta',
           'Median per-residue pLDDT r', 'r(model pLDDT, TM)',
           'r(AFDB pLDDT, TM)', 'rho(min pLDDT, TM)'], conf_rows)
summary['confidence'] = {
    m: {'afdb_plddt_mean': qstats(V[m]['plddt_ref'])['mean'],
        'model_plddt_mean': qstats(V[m]['plddt_mod'])['mean'],
        'per_protein_plddt_r_median': qstats(V[m]['plddt_r'])['median'],
        'r_modelplddt_tm': pearson(V[m]['plddt_mod'], V[m]['tm_score']),
        'r_afdbplddt_tm': pearson(V[m]['plddt_ref'], V[m]['tm_score'])}
    for m in METHODS}

# per-residue pooled
res = np.load(os.path.join(HERE, 'per_residue.npz'))
pooled = {}
for m in METHODS:
    pr = res[f'{m}__plddt_ref']; pm = res[f'{m}__plddt_mod']; dv = res[f'{m}__dev']
    pooled[m] = (pr, pm, dv)
    mins = np.minimum(pr, pm)
    bins = {}
    for lo, hi, key in ((0, 50, 'lt50'), (50, 70, '50_70'),
                        (70, 90, '70_90'), (90, 101, 'ge90')):
        sel = (mins >= lo) & (mins < hi)
        bins[key] = {'n': int(sel.sum()),
                     'median_dev': float(np.median(dv[sel])) if sel.any() else None}
    summary.setdefault('per_residue', {})[m] = {
        'n_residues': int(pr.size),
        'r_plddt_afdb_vs_model': pearson(pr, pm),
        'median_dev_plddt_ge90': float(np.median(dv[(pr >= 90) & (pm >= 90)])) if ((pr >= 90) & (pm >= 90)).any() else None,
        'median_dev_plddt_lt50': float(np.median(dv[(pr < 50) | (pm < 50)])) if ((pr < 50) | (pm < 50)).any() else None,
        'dev_by_min_plddt_bin': bins,
    }

fig, axes = plt.subplots(1, 3, figsize=(11, 3.2))
ax = axes[0]
edges = np.array([0, 50, 70, 90, 100])
w = 0.2
for k, m in enumerate(METHODS):
    pr, pm, dv = pooled[m]
    mins = np.minimum(pr, pm)
    med = [np.median(dv[(mins >= edges[i]) & (mins < edges[i + 1])])
           if ((mins >= edges[i]) & (mins < edges[i + 1])).any() else np.nan
           for i in range(len(edges) - 1)]
    ax.bar(np.arange(4) + (k - 1.5) * w, med, width=w, color=COLOR[m],
           label=LABEL[m], alpha=0.85)
ax.set_xticks(range(4))
ax.set_xticklabels(['< 50', '50\u201370', '70\u201390', '\u2265 90'])
ax.set_xlabel('Minimum Per-Residue pLDDT of AFDB and Model')
ax.set_ylabel('Median C\u03b1 Deviation (\u00c5)')
ax.set_yscale('log'); ax.legend(fontsize=7.5)
ax.set_title('a  Deviation vs Confidence', loc='left', fontweight='bold')

ax = axes[1]
for m in METHODS:
    x = V[m]['plddt_mod']; y = V[m]['tm_score']
    ax.scatter(x, y, s=3, alpha=0.18, color=COLOR[m], edgecolors='none')
ax.set_xlabel('Model Mean pLDDT'); ax.set_ylabel('TM-Score vs AFDB')
ax.set_title('b  Per-Protein Confidence vs Agreement', loc='left', fontweight='bold')

ax = axes[2]
for m in METHODS:
    ax.scatter(V[m]['plddt_ref'], V[m]['plddt_mod'], s=3, alpha=0.18,
               color=COLOR[m], edgecolors='none', label=LABEL[m])
lims = [20, 100]
ax.plot(lims, lims, 'k--', lw=0.8)
ax.set_xlim(lims); ax.set_ylim(lims)
ax.set_xlabel('AFDB Mean pLDDT'); ax.set_ylabel('Model Mean pLDDT')
ax.legend(fontsize=7.5, markerscale=3)
ax.set_title('c  Confidence Calibration', loc='left', fontweight='bold')
fig.tight_layout()
savefig(fig, 'FigS2_confidence')

# ================================================ S5 disorder / core analysis
print('\n[4] confident-core restriction')
core_rows = []
for m in METHODS:
    full = V[m]['tm_score']; core = V[m]['tm_core']
    ld, ldc = V[m]['lddt'], V[m]['lddt_core']
    ok = np.isfinite(core)
    core_rows.append([LABEL[m], int(ok.sum()),
                      f'{np.nanmean(full[ok]):.3f}', f'{np.nanmean(core[ok]):.3f}',
                      f'{np.nanmean(core[ok])-np.nanmean(full[ok]):+.3f}',
                      f'{np.nanmean(ld[ok]):.3f}', f'{np.nanmean(ldc[ok]):.3f}',
                      f'{np.nanmean(ldc[ok])-np.nanmean(ld[ok]):+.3f}',
                      f'{100*np.nanmean(core[ok]>=0.5):.1f}%'])
write_tsv(os.path.join(TABDIR, 'S5_confident_core.tsv'),
          ['Method', 'n (>=20 core residues)', 'TM full chain', 'TM core',
           'Delta TM', 'lDDT full chain', 'lDDT core', 'Delta lDDT',
           'core TM>=0.5'], core_rows)
summary['core'] = {m: {'tm_full': float(np.nanmean(V[m]['tm_score'][np.isfinite(V[m]['tm_core'])])),
                       'tm_core': float(np.nanmean(V[m]['tm_core'])),
                       'lddt_full': float(np.nanmean(V[m]['lddt'][np.isfinite(V[m]['lddt_core'])])),
                       'lddt_core': float(np.nanmean(V[m]['lddt_core']))}
                   for m in METHODS}

# ==================================================== S6 length stratification
print('\n[5] chain-length stratification')
lbins = [(0, 100), (100, 200), (200, 400), (400, 700), (700, 5000)]
len_rows = []
for lo, hi in lbins:
    row = [f'{lo}-{hi if hi<5000 else "max"}']
    n_shown = None
    for m in METHODS:
        sel = (NREF[m] >= lo) & (NREF[m] < hi)
        n_shown = int(sel.sum()) if n_shown is None else n_shown
        row.append(f'{np.nanmean(V[m]["tm_score"][sel]):.3f}' if sel.any() else '-')
    row.insert(1, n_shown)
    for m in METHODS:
        sel = (NREF[m] >= lo) & (NREF[m] < hi)
        row.append(f'{100*np.nanmean(V[m]["tm_score"][sel]>=0.5):.0f}%' if sel.any() else '-')
    len_rows.append(row)
write_tsv(os.path.join(TABDIR, 'S6_length_strata.tsv'),
          ['Chain length', 'n'] + [f'{LABEL[m]} mean TM' for m in METHODS]
          + [f'{LABEL[m]} TM>=0.5' for m in METHODS], len_rows)

fig, axes = plt.subplots(1, 3, figsize=(11, 3.2))
ax = axes[0]
xs = np.arange(len(lbins))
for k, m in enumerate(METHODS):
    ys = [np.nanmean(V[m]['tm_score'][(NREF[m] >= lo) & (NREF[m] < hi)])
          for lo, hi in lbins]
    ax.plot(xs, ys, 'o-', color=COLOR[m], label=LABEL[m], lw=1.5, ms=4)
ax.set_xticks(xs); ax.set_xticklabels([f'{lo}-{hi if hi<5000 else "+"}' for lo, hi in lbins])
ax.set_xlabel('Chain Length (Residues)'); ax.set_ylabel('Mean TM-Score vs AFDB')
ax.legend(fontsize=7.5)
ax.set_title('a  Length Dependence', loc='left', fontweight='bold')

ax = axes[1]
xs = np.arange(len(METHODS))
full = [np.nanmean(V[m]['tm_score'][np.isfinite(V[m]['tm_core'])]) for m in METHODS]
core = [np.nanmean(V[m]['tm_core']) for m in METHODS]
ax.bar(xs - 0.2, full, 0.4, label='Full Chain', color='#999999')
ax.bar(xs + 0.2, core, 0.4, label='Confident Core (pLDDT \u2265 70)', color='#4c72b0')
ax.set_xticks(xs); ax.set_xticklabels([LABEL[m] for m in METHODS], rotation=20)
ax.set_ylabel('Mean TM-Score'); ax.set_ylim(0.6, 1.0); ax.legend(fontsize=7.5)
ax.set_title('b  Effect of Low-Confidence Regions', loc='left', fontweight='bold')

ax = axes[2]
for m in METHODS:
    ax.scatter(V[m]['ss_q3'], V[m]['tm_score'], s=3, alpha=0.15, color=COLOR[m],
               edgecolors='none')
ax.set_xlabel('Secondary-Structure Agreement (Q3)'); ax.set_ylabel('TM-Score vs AFDB')
ax.set_title('c  Fold vs Local Geometry', loc='left', fontweight='bold')
fig.tight_layout()
savefig(fig, 'FigS3_stratification')

# ===================================================== S7 secondary structure
print('\n[6] secondary structure and compactness')
ss_rows = []
for m in METHODS:
    ss_rows.append([LABEL[m],
                    f"{np.nanmean(V[m]['ss_q3']):.3f}",
                    f"{np.nanmedian(V[m]['ss_q3']):.3f}",
                    f"{100*np.nanmean(V[m]['helix_ref']):.1f}%",
                    f"{100*np.nanmean(V[m]['helix_mod']):.1f}%",
                    f"{100*np.nanmean(V[m]['strand_ref']):.1f}%",
                    f"{100*np.nanmean(V[m]['strand_mod']):.1f}%",
                    f"{np.nanmean(V[m]['rg_ratio']):.3f}",
                    f"{np.nanmedian(V[m]['rg_ratio']):.3f}"])
write_tsv(os.path.join(TABDIR, 'S7_secondary_structure.tsv'),
          ['Method', 'Q3 mean', 'Q3 median', 'AFDB helix', 'Model helix',
           'AFDB strand', 'Model strand', 'Rg ratio mean', 'Rg ratio median'],
          ss_rows)
summary['ss'] = {m: {'q3_mean': float(np.nanmean(V[m]['ss_q3'])),
                     'q3_median': float(np.nanmedian(V[m]['ss_q3'])),
                     'r_q3_tm': pearson(V[m]['ss_q3'], V[m]['tm_score']),
                     'q3_tm_ge_0.9': float(np.nanmean(V[m]['ss_q3'][V[m]['tm_score'] >= 0.9])),
                     'q3_tm_lt_0.5': float(np.nanmean(V[m]['ss_q3'][V[m]['tm_score'] < 0.5])),
                     'n_tm_lt_0.5': int((V[m]['tm_score'] < 0.5).sum()),
                     'helix_ref': float(np.nanmean(V[m]['helix_ref'])),
                     'helix_mod': float(np.nanmean(V[m]['helix_mod'])),
                     'strand_ref': float(np.nanmean(V[m]['strand_ref'])),
                     'strand_mod': float(np.nanmean(V[m]['strand_mod'])),
                     'rg_ratio_median': float(np.nanmedian(V[m]['rg_ratio']))}
                 for m in METHODS}

# ================================================== S8 discordance catalogue
print('\n[7] discordance catalogue')
disc_rows = []
consensus_fail = []
for ml, d in sorted(byml.items()):
    tms = {m: fnum(d[m]['tm_score']) for m in METHODS if m in d}
    if not tms:
        continue
    vals = np.array(list(tms.values()))
    n_fail = int(np.nansum(vals < 0.5))
    mrow = manifest[ml]
    e = ess.get(ml, {})
    rec = [ml, d[list(tms)[0]]['uniprot'], d[list(tms)[0]]['gene'],
           d[list(tms)[0]]['n_ref']]
    rec += [f'{tms.get(m, np.nan):.3f}' if m in tms else '' for m in METHODS]
    rec += [f'{np.nanmax(vals):.3f}', f'{np.nanmin(vals):.3f}',
            f'{np.nanmax(vals)-np.nanmin(vals):.3f}', n_fail,
            d[list(tms)[0]]['plddt_ref'],
            f'{np.nanmean([fnum(d[m]["plddt_mod"]) for m in tms]):.1f}',
            e.get('priority_tier') or '',
            mrow['protein_name'][:90]]
    if n_fail > 0 or (np.nanmax(vals) - np.nanmin(vals)) > 0.3:
        disc_rows.append(rec)
    if n_fail == len(tms):
        consensus_fail.append(rec)
disc_rows.sort(key=lambda r: (-int(r[11]), float(r[8])))
write_tsv(os.path.join(TABDIR, 'S8_discordant_proteins.tsv'),
          ['ML ID', 'UniProt', 'Gene', 'Residues'] + [f'TM {LABEL[m]}' for m in METHODS]
          + ['TM max', 'TM min', 'TM range', 'n methods TM<0.5', 'AFDB pLDDT',
             'mean model pLDDT', 'Priority tier', 'Protein'], disc_rows)
summary['discordance'] = {
    'proteins_flagged': len(disc_rows),
    'all_methods_below_0.5': len(consensus_fail),
    'method_specific_failures': {
        m: int(np.nansum((V[m]['tm_score'] < 0.5))) for m in METHODS},
}

# unique failures per method (only that method fails)
uniq = collections.Counter()
for ml, d in byml.items():
    tms = {m: fnum(d[m]['tm_score']) for m in METHODS if m in d}
    fails = [m for m, v in tms.items() if v < 0.5]
    if len(fails) == 1 and len(tms) == len(METHODS):
        uniq[fails[0]] += 1
summary['discordance']['method_unique_failures'] = dict(uniq)

# priority tier of the proteins that every method places differently
tier_cf = collections.Counter()
for ml, d in byml.items():
    tms = [fnum(d[m]['tm_score']) for m in METHODS if m in d]
    if len(tms) == len(METHODS) and max(tms) < 0.5:
        tier_cf[(prior.get(ml) or {}).get('tier', 'not ranked')] += 1
summary['tier_of_consensus_failures'] = dict(tier_cf)

# confidence of the shortest chains, quoted in the length-stratification section
short_pl, all_pl = [], []
for m in METHODS:
    sel = NREF[m] < 100
    short_pl.append(V[m]['plddt_ref'][sel])
    all_pl.append(V[m]['plddt_ref'])
summary['short_chain_confidence'] = {
    'n_under_100': int(sum(int((NREF[m] < 100).sum()) for m in METHODS) / len(METHODS)),
    'mean_afdb_plddt_under_100': float(np.nanmean(np.concatenate(short_pl))),
    'mean_afdb_plddt_all': float(np.nanmean(np.concatenate(all_pl))),
}

# ============================ S8b fold-conservation vs domain-packing split
print('\n[7b] local fold vs global packing decomposition')
split_rows = []
for m in METHODS:
    tm = V[m]['tm_score']; ld = V[m]['lddt']
    # partial-coverage pairs are excluded: their TM-score is capped by the
    # missing part of the chain, not by a structural difference
    ok = np.isfinite(tm) & np.isfinite(ld) & FULLCOV[m]
    rearr = ok & (ld >= 0.80) & (tm < 0.70)      # domains right, packing differs
    both_bad = ok & (ld < 0.80) & (tm < 0.70)     # genuinely different fold
    both_good = ok & (ld >= 0.80) & (tm >= 0.70)
    split_rows.append([LABEL[m], int(ok.sum()),
                       f'{int(both_good.sum())} ({100*both_good.sum()/ok.sum():.1f}%)',
                       f'{int(rearr.sum())} ({100*rearr.sum()/ok.sum():.1f}%)',
                       f'{int(both_bad.sum())} ({100*both_bad.sum()/ok.sum():.1f}%)',
                       f'{np.nanmean(V[m]["rmsd_global"][rearr]):.1f}' if rearr.any() else '-',
                       f'{np.nanmean(V[m]["rmsd_core"][rearr]):.2f}' if rearr.any() else '-',
                       f'{np.nanmean(NREF[m][rearr]):.0f}' if rearr.any() else '-'])
write_tsv(os.path.join(TABDIR, 'S8b_fold_vs_packing.tsv'),
          ['Method', 'n', 'Same fold and packing (lDDT>=0.8, TM>=0.7)',
           'Domain rearrangement (lDDT>=0.8, TM<0.7)',
           'Divergent local structure (lDDT<0.8, TM<0.7)',
           'Mean RMSD of rearranged (A)', 'Mean core RMSD of rearranged (A)',
           'Mean length of rearranged'], split_rows)
summary['fold_vs_packing'] = {
    m: {'n': int(FULLCOV[m].sum()),
        'domain_rearrangement': int(((V[m]['lddt'] >= 0.80) & (V[m]['tm_score'] < 0.70) & FULLCOV[m]).sum()),
        'divergent_local': int(((V[m]['lddt'] < 0.80) & (V[m]['tm_score'] < 0.70) & FULLCOV[m]).sum()),
        'concordant': int(((V[m]['lddt'] >= 0.80) & (V[m]['tm_score'] >= 0.70) & FULLCOV[m]).sum()),
        'partial_coverage_excluded': int((~FULLCOV[m]).sum())}
    for m in METHODS}
summary['coverage_sensitivity'] = {
    m: {'median_tm_all': float(np.nanmedian(V[m]['tm_score'])),
        'median_tm_full_coverage': float(np.nanmedian(V[m]['tm_score'][FULLCOV[m]])),
        'n_partial': int((~FULLCOV[m]).sum())} for m in METHODS}

# characterisation of the proteins that every method places differently
cf_ml = [r[0] for r in disc_rows if int(r[11]) == len([m for m in METHODS])]
cf_len = np.array([float(byml[ml][list(byml[ml])[0]]['n_ref']) for ml in cf_ml])
cf_pl = np.array([float(byml[ml][list(byml[ml])[0]]['plddt_ref']) for ml in cf_ml])
cf_ld = np.array([np.nanmean([fnum(byml[ml][m]['lddt']) for m in byml[ml]]) for ml in cf_ml])
all_len = np.concatenate([NREF[m] for m in METHODS])
summary['consensus_failures'] = {
    'n': len(cf_ml),
    'median_length': float(np.median(cf_len)) if len(cf_ml) else None,
    'median_length_all': float(np.median(all_len)),
    'median_afdb_plddt': float(np.median(cf_pl)) if len(cf_ml) else None,
    'median_afdb_plddt_all': float(np.median(np.concatenate([V[m]['plddt_ref'] for m in METHODS]))),
    'median_mean_lddt': float(np.median(cf_ld)) if len(cf_ml) else None,
    'frac_under_150aa': float(np.mean(cf_len < 150)) if len(cf_ml) else None,
    'frac_afdb_plddt_under_70': float(np.mean(cf_pl < 70)) if len(cf_ml) else None,
}

# ================================================= S15 named case studies
print('\n[7c] case studies')
case_rows = []


def prow(ml):
    d = byml[ml]
    r0 = next(iter(d.values()))
    cov = min(fnum(d[m]['n_common']) / fnum(d[m]['n_ref']) for m in d)
    return d, r0, cov


for ml, d in byml.items():
    if len(d) < len(METHODS):
        continue
    d, r0, cov = prow(ml)
    if cov < 0.95:
        continue
    tms = np.array([fnum(d[m]['tm_score']) for m in METHODS])
    lds = np.array([fnum(d[m]['lddt']) for m in METHODS])
    plr = fnum(r0['plddt_ref'])
    plm = np.nanmean([fnum(d[m]['plddt_mod']) for m in METHODS])
    cls = None
    if lds.min() >= 0.85 and tms.max() < 0.70 and plr >= 80:
        cls = 'Domain rearrangement (all four methods)'
    elif tms.max() < 0.50 and plr >= 80:
        cls = 'Whole-proteome discordance despite confident AFDB model'
    else:
        fails = [m for m in METHODS if fnum(d[m]['tm_score']) < 0.5]
        others = [fnum(d[m]['tm_score']) for m in METHODS if m not in fails]
        if len(fails) == 1 and others and min(others) >= 0.70:
            cls = f'Single-method outlier ({LABEL[fails[0]]})'
    if cls:
        case_rows.append([cls, ml, r0['uniprot'], r0['gene'], int(r0['n_ref'])]
                         + [f'{fnum(d[m]["tm_score"]):.2f}' for m in METHODS]
                         + [f'{np.nanmean(lds):.3f}', f'{plr:.1f}', f'{plm:.1f}',
                            manifest[ml]['protein_name'][:90]])
case_rows.sort(key=lambda r: (r[0], -int(r[4])))
write_tsv(os.path.join(TABDIR, 'S15_case_studies.tsv'),
          ['Class', 'ML ID', 'UniProt', 'Gene', 'Residues']
          + [f'TM {LABEL[m]}' for m in METHODS]
          + ['mean lDDT', 'AFDB pLDDT', 'mean model pLDDT', 'Protein'], case_rows)
summary['case_studies'] = collections.Counter(r[0] for r in case_rows)

# ==================================================== S9 5-way concordance
print('\n[8] five-way concordance matrix')
if cross:
    cm = collections.defaultdict(list)
    for r in cross:
        cm[(r['method_a'], r['method_b'])].append(fnum(r['tm_score']))
    order = ['afdb'] + METHODS
    mat = np.full((5, 5), np.nan)
    nmat = np.zeros((5, 5), dtype=int)
    for i, a in enumerate(order):
        for j, b in enumerate(order):
            if i == j:
                mat[i, j] = 1.0
                continue
            key = (a, b) if (a, b) in cm else (b, a)
            if key in cm:
                mat[i, j] = np.nanmedian(cm[key])
                nmat[i, j] = len(cm[key])
    rows = []
    for i, a in enumerate(order):
        rows.append([LABEL[a]] + [f'{mat[i,j]:.3f}' if np.isfinite(mat[i, j]) else '-'
                                  for j in range(5)]
                    + [f'{np.nanmean([mat[i,j] for j in range(5) if j!=i]):.3f}'])
    write_tsv(os.path.join(TABDIR, 'S9_five_way_concordance.tsv'),
              ['Method'] + [LABEL[o] for o in order] + ['mean vs others'], rows)
    summary['five_way'] = {
        'order': [LABEL[o] for o in order],
        'median_tm_matrix': [[None if not np.isfinite(v) else round(float(v), 4)
                              for v in mat[i]] for i in range(5)],
        'n_matrix': nmat.tolist(),
        'mean_vs_others': {LABEL[a]: float(np.nanmean([mat[i, j] for j in range(5) if j != i]))
                           for i, a in enumerate(order)},
    }

    fig, axes = plt.subplots(1, 2, figsize=(9.5, 4))
    ax = axes[0]
    im = ax.imshow(mat, vmin=0.6, vmax=1.0, cmap='viridis')
    ax.set_xticks(range(5)); ax.set_xticklabels([LABEL[o] for o in order], rotation=35, ha='right')
    ax.set_yticks(range(5)); ax.set_yticklabels([LABEL[o] for o in order])
    for i in range(5):
        for j in range(5):
            if np.isfinite(mat[i, j]):
                ax.text(j, i, f'{mat[i,j]:.2f}', ha='center', va='center',
                        color='w' if mat[i, j] < 0.9 else 'k', fontsize=8)
    ax.grid(False)
    fig.colorbar(im, ax=ax, shrink=0.8, label='Median Pairwise TM-Score')
    ax.set_title('a  Five-Way Concordance', loc='left', fontweight='bold')

    ax = axes[1]
    means = [np.nanmean([mat[i, j] for j in range(5) if j != i]) for i in range(5)]
    ax.barh([LABEL[o] for o in order], means,
            color=[COLOR[o] for o in order], alpha=0.85)
    ax.set_xlim(0.6, 1.0)
    ax.set_xlabel('Mean Median TM-Score to the Other Four Methods')
    ax.set_title('b  Centrality in the Ensemble', loc='left', fontweight='bold')
    fig.tight_layout()
    savefig(fig, 'FigS4_five_way')

# ============================================ S10 target-relevance stratum
print('\n[9] drug-target relevance')
tier_rows = []
tiers = ['High-priority', 'Strong candidate', 'Moderate candidate', 'Exploratory']
for t in tiers:
    mls = {ml for ml, v in prior.items() if v['tier'] == t}
    row = [t, len(mls & set(byml))]
    for m in METHODS:
        sel = np.array([ml in mls for ml in MLS[m]])
        row.append(f'{np.nanmean(V[m]["tm_score"][sel]):.3f}' if sel.any() else '-')
    for m in METHODS:
        sel = np.array([ml in mls for ml in MLS[m]])
        row.append(f'{np.nanmean(V[m]["plddt_mod"][sel]):.1f}' if sel.any() else '-')
    tier_rows.append(row)
write_tsv(os.path.join(TABDIR, 'S10_target_tiers.tsv'),
          ['Priority tier', 'n proteins'] + [f'{LABEL[m]} mean TM' for m in METHODS]
          + [f'{LABEL[m]} mean pLDDT' for m in METHODS], tier_rows)

# top-50 targets, per protein
top50 = sorted([ml for ml in prior if prior[ml]['rank'] <= 50], key=lambda x: prior[x]['rank'])
t50_rows = []
for ml in top50:
    d = byml.get(ml)
    if not d:
        continue
    r0 = next(iter(d.values()))
    t50_rows.append([prior[ml]['rank'], ml, r0['uniprot'], r0['gene'], r0['n_ref']]
                    + [d[m]['tm_score'] if m in d else '' for m in METHODS]
                    + [d[m]['lddt'] if m in d else '' for m in METHODS]
                    + [r0['plddt_ref'], manifest[ml]['protein_name'][:80]])
write_tsv(os.path.join(TABDIR, 'S11_top50_targets.tsv'),
          ['Priority rank', 'ML ID', 'UniProt', 'Gene', 'Residues']
          + [f'TM {LABEL[m]}' for m in METHODS]
          + [f'lDDT {LABEL[m]}' for m in METHODS]
          + ['AFDB pLDDT', 'Protein'], t50_rows)

# membrane / secreted strata
strat_rows = []
for name, key in (('Transmembrane annotated', 'transmembrane'),
                  ('Signal peptide annotated', 'signal_peptide')):
    mls = {ml for ml, r in manifest.items() if r.get(key) == '1'}
    row = [name, len(mls & set(byml))]
    for m in METHODS:
        sel = np.array([ml in mls for ml in MLS[m]])
        row.append(f'{np.nanmean(V[m]["tm_score"][sel]):.3f}' if sel.any() else '-')
    strat_rows.append(row)
mls_other = {ml for ml, r in manifest.items()
             if r.get('transmembrane') != '1' and r.get('signal_peptide') != '1'}
row = ['Soluble (neither annotation)', len(mls_other & set(byml))]
for m in METHODS:
    sel = np.array([ml in mls_other for ml in MLS[m]])
    row.append(f'{np.nanmean(V[m]["tm_score"][sel]):.3f}')
strat_rows.append(row)
write_tsv(os.path.join(TABDIR, 'S12_localisation_strata.tsv'),
          ['Stratum', 'n'] + [f'{LABEL[m]} mean TM' for m in METHODS], strat_rows)

# =============================================== S13 full per-protein table
print('\n[10] full per-protein table')
full_hdr = (['ML ID', 'UniProt', 'Gene', 'Protein', 'AFDB residues', 'AFDB pLDDT']
            + [f'{c} {LABEL[m]}' for m in METHODS
               for c in ('TM', 'GDT-TS', 'lDDT', 'RMSD', 'core TM', 'Q3', 'pLDDT')]
            + ['TM best method', 'TM best', 'TM worst', 'TM range',
               'Priority tier', 'Priority rank', 'ProteomeLM class'])
full_rows = []
for ml, d in sorted(byml.items()):
    r0 = next(iter(d.values()))
    e = ess.get(ml, {})
    tms = {m: fnum(d[m]['tm_score']) for m in METHODS if m in d}
    best = max(tms, key=lambda k: tms[k]) if tms else ''
    rec = [ml, r0['uniprot'], r0['gene'], manifest[ml]['protein_name'][:120],
           r0['n_ref'], r0['plddt_ref']]
    for m in METHODS:
        if m in d:
            rec += [d[m]['tm_score'], d[m]['gdt_ts'], d[m]['lddt'],
                    d[m]['rmsd_global'], d[m]['tm_core'], d[m]['ss_q3'],
                    d[m]['plddt_mod']]
        else:
            rec += [''] * 7
    rec += [LABEL.get(best, ''), f'{max(tms.values()):.3f}' if tms else '',
            f'{min(tms.values()):.3f}' if tms else '',
            f'{max(tms.values())-min(tms.values()):.3f}' if tms else '',
            (prior.get(ml) or {}).get('tier', ''),
            (prior.get(ml) or {}).get('rank', ''),
            e.get('proteomelm_ess_class', '')]
    full_rows.append(rec)
write_tsv(os.path.join(TABDIR, 'S13_per_protein_full.tsv'), full_hdr, full_rows)

# best-method counts
best_counts = collections.Counter()
for ml, d in byml.items():
    tms = {m: fnum(d[m]['tm_score']) for m in METHODS if m in d}
    if len(tms) == len(METHODS):
        best_counts[max(tms, key=lambda k: tms[k])] += 1
summary['best_method_counts'] = {LABEL[m]: best_counts[m] for m in METHODS}

# agreement-as-quality-proxy: cross-method agreement vs AFDB agreement
if cross:
    cmm = collections.defaultdict(list)
    for r in cross:
        if r['method_a'] != 'afdb' and r['method_b'] != 'afdb':
            cmm[r['ml_id']].append(fnum(r['tm_score']))
    x, y = [], []
    for ml, d in byml.items():
        if ml in cmm and cmm[ml]:
            tms = [fnum(d[m]['tm_score']) for m in METHODS if m in d]
            x.append(np.nanmean(cmm[ml])); y.append(np.nanmean(tms))
    summary['agreement_proxy'] = {
        'n': len(x),
        'pearson_crossmethod_vs_afdb': pearson(x, y),
        'spearman_crossmethod_vs_afdb': spearman(x, y),
    }
    fig, ax = plt.subplots(figsize=(4.2, 3.6))
    ax.scatter(x, y, s=4, alpha=0.25, color='#4c72b0', edgecolors='none')
    ax.plot([0, 1], [0, 1], 'k--', lw=0.8)
    ax.set_xlabel('Mean Pairwise TM-Score Among the Four Predictors')
    ax.set_ylabel('Mean TM-Score to AFDB')
    ax.set_title(f"r = {summary['agreement_proxy']['pearson_crossmethod_vs_afdb']:.3f}",
                 loc='left', fontsize=9)
    fig.tight_layout()
    savefig(fig, 'FigS5_agreement_proxy')

with open(os.path.join(HERE, 'summary.json'), 'w') as fh:
    json.dump(summary, fh, indent=2, default=float)
print('\nwrote summary.json')
print(json.dumps({k: summary[k] for k in ('coverage', 'best_method_counts',
                                          'discordance')}, indent=2, default=float))
