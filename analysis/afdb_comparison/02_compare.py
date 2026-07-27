#!/usr/bin/env python3
"""Pairwise structural comparison: AlphaFold DB (v6) vs HANSEN monomer models.

For every M. leprae protein with an AFDB entry, each HANSEN monomer model
(AF3, Boltz-1, Boltz-2, Chai-1) is compared against the AFDB reference.

Writes
  pairs.tsv        one row per (protein, method) comparison
  per_residue.npz  pooled per-residue arrays for confidence/deviation analysis
"""
import os
import sys
import csv
import json
import numpy as np
from multiprocessing import Pool

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from structlib import (parse_ca, align_by_resnum, needleman_wunsch, kabsch_rt,
                       apply_rt, tm_score, gdt, lddt_ca, psea_ss,
                       radius_of_gyration, pearson, spearman)

MANIFEST = os.path.join(HERE, 'manifest.tsv')
OUT_PAIRS = os.path.join(HERE, 'pairs.tsv')
OUT_RES = os.path.join(HERE, 'per_residue.npz')
METHODS = ['af3', 'boltz', 'boltz2', 'chai']

FIELDS = [
    'ml_id', 'uniprot', 'gene', 'method', 'n_ref', 'n_mod', 'n_common',
    'seq_identical', 'align_mode', 'seq_pid',
    'tm_score', 'rmsd_global', 'rmsd_core', 'gdt_ts', 'gdt_ha', 'lddt',
    'frac_lt2', 'frac_lt5', 'dev_p95',
    'tm_core', 'lddt_core', 'n_core',
    'plddt_ref', 'plddt_mod', 'plddt_r', 'plddt_rho',
    'rg_ref', 'rg_mod', 'rg_ratio',
    'ss_q3', 'helix_ref', 'strand_ref', 'coil_ref',
    'helix_mod', 'strand_mod', 'coil_mod',
]


def read_manifest():
    with open(MANIFEST) as fh:
        return list(csv.DictReader(fh, delimiter='\t'))


def compare_one(args):
    row = args
    ml = row['ml_id']
    afdb = row['afdb_path']
    if not afdb:
        return [], None
    try:
        ref = parse_ca(afdb)
    except Exception as e:
        return [], ('parse_fail_ref', ml, str(e))

    ref_ss = psea_ss(ref['xyz'])
    ref_rg = radius_of_gyration(ref['xyz'])

    out = []
    res_pack = {}
    for m in METHODS:
        p = row.get(f'{m}_path')
        if not p:
            continue
        try:
            mod = parse_ca(p)
        except Exception as e:
            out.append({'ml_id': ml, 'uniprot': row['uniprot'],
                        'gene': row['gene'], 'method': m,
                        'align_mode': f'parse_fail:{e}'})
            continue

        seq_ident = (ref['seq'] == mod['seq'])
        if seq_ident:
            ia = np.arange(len(ref['seq']))
            ib = np.arange(len(mod['seq']))
            mode = 'identity'
        else:
            ia, ib = align_by_resnum(ref, mod)
            same = sum(1 for x, y in zip(ia, ib)
                       if ref['seq'][x] == mod['seq'][y])
            if len(ia) >= 0.9 * min(len(ref['seq']), len(mod['seq'])) and \
               same >= 0.98 * max(1, len(ia)):
                mode = 'resnum'
            else:
                ia, ib = needleman_wunsch(ref['seq'], mod['seq'])
                mode = 'nw'
        if len(ia) < 3:
            out.append({'ml_id': ml, 'uniprot': row['uniprot'],
                        'gene': row['gene'], 'method': m,
                        'align_mode': 'too_short'})
            continue

        Q = ref['xyz'][ia]           # reference (AFDB)
        P = mod['xyz'][ib]           # model
        pid = float(np.mean([ref['seq'][x] == mod['seq'][y]
                             for x, y in zip(ia, ib)]))

        tm, R, t = tm_score(P, Q, L_norm=len(ref['seq']))
        Pt_tm = apply_rt(P, R, t)
        d_tm = np.linalg.norm(Pt_tm - Q, axis=1)

        Rg_, tg = kabsch_rt(P, Q)
        d_glob = np.linalg.norm(apply_rt(P, Rg_, tg) - Q, axis=1)
        rmsd_global = float(np.sqrt((d_glob ** 2).mean()))

        core_sel = d_tm < 5.0
        if core_sel.sum() >= 3:
            rmsd_core = float(np.sqrt((d_tm[core_sel] ** 2).mean()))
        else:
            rmsd_core = float('nan')

        gts, _ = gdt(P, Q, (1.0, 2.0, 4.0, 8.0))
        gha, _ = gdt(P, Q, (0.5, 1.0, 2.0, 4.0))
        ld, ld_res = lddt_ca(P, Q)

        pl_ref = ref['plddt'][ia]
        pl_mod = mod['plddt'][ib]

        # high-confidence core: both predictors confident
        cm = (pl_ref >= 70) & (pl_mod >= 70)
        if cm.sum() >= 20:
            tmc, _, _ = tm_score(P[cm], Q[cm], L_norm=int(cm.sum()))
            ldc, _ = lddt_ca(P[cm], Q[cm])
        else:
            tmc, ldc = float('nan'), float('nan')

        mod_ss = psea_ss(mod['xyz'])
        ss_r = ''.join(ref_ss[x] for x in ia)
        ss_m = ''.join(mod_ss[y] for y in ib)
        q3 = float(np.mean([a == b for a, b in zip(ss_r, ss_m)]))

        rec = {
            'ml_id': ml, 'uniprot': row['uniprot'], 'gene': row['gene'],
            'method': m,
            'n_ref': len(ref['seq']), 'n_mod': len(mod['seq']),
            'n_common': len(ia),
            'seq_identical': int(seq_ident), 'align_mode': mode,
            'seq_pid': round(pid, 4),
            'tm_score': round(tm, 4),
            'rmsd_global': round(rmsd_global, 3),
            'rmsd_core': round(rmsd_core, 3) if rmsd_core == rmsd_core else '',
            'gdt_ts': round(gts, 4), 'gdt_ha': round(gha, 4),
            'lddt': round(ld, 4),
            'frac_lt2': round(float((d_tm < 2).mean()), 4),
            'frac_lt5': round(float((d_tm < 5).mean()), 4),
            'dev_p95': round(float(np.percentile(d_tm, 95)), 3),
            'tm_core': round(tmc, 4) if tmc == tmc else '',
            'lddt_core': round(ldc, 4) if ldc == ldc else '',
            'n_core': int(cm.sum()),
            'plddt_ref': round(float(np.nanmean(pl_ref)), 2),
            'plddt_mod': round(float(np.nanmean(pl_mod)), 2),
            'plddt_r': round(pearson(pl_ref, pl_mod), 4),
            'plddt_rho': round(spearman(pl_ref, pl_mod), 4),
            # Rg computed over the corresponded residues only, so that
            # construct-length differences do not confound compactness
            'rg_ref': round(radius_of_gyration(Q), 2),
            'rg_mod': round(radius_of_gyration(P), 2),
            'rg_ratio': round(radius_of_gyration(P) / radius_of_gyration(Q), 4),
            'ss_q3': round(q3, 4),
            'helix_ref': round(ss_r.count('H') / len(ss_r), 4),
            'strand_ref': round(ss_r.count('E') / len(ss_r), 4),
            'coil_ref': round(ss_r.count('C') / len(ss_r), 4),
            'helix_mod': round(ss_m.count('H') / len(ss_m), 4),
            'strand_mod': round(ss_m.count('E') / len(ss_m), 4),
            'coil_mod': round(ss_m.count('C') / len(ss_m), 4),
        }
        out.append(rec)
        res_pack[m] = {
            'plddt_ref': pl_ref.astype(np.float32),
            'plddt_mod': pl_mod.astype(np.float32),
            'dev': d_tm.astype(np.float32),
            'lddt_res': ld_res.astype(np.float32),
        }
    return out, (ml, res_pack)


def main():
    rows = [r for r in read_manifest() if r['afdb_path']]
    print(f'comparing {len(rows)} proteins x up to {len(METHODS)} methods')
    nproc = max(1, (os.cpu_count() or 4) - 1)
    results, packs = [], {}
    with Pool(nproc) as pool:
        for i, (recs, pk) in enumerate(pool.imap_unordered(compare_one, rows, chunksize=4), 1):
            results.extend(recs)
            if pk and isinstance(pk, tuple) and isinstance(pk[1], dict):
                packs[pk[0]] = pk[1]
            if i % 100 == 0:
                print(f'  {i}/{len(rows)}', flush=True)

    results.sort(key=lambda r: (r['ml_id'], r['method']))
    with open(OUT_PAIRS, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDS, delimiter='\t',
                           extrasaction='ignore')
        w.writeheader()
        for r in results:
            w.writerow(r)
    print(f'wrote {OUT_PAIRS} ({len(results)} rows)')

    # pooled per-residue arrays (subsampled cap keeps the npz manageable)
    pooled = {}
    for m in METHODS:
        for key in ('plddt_ref', 'plddt_mod', 'dev', 'lddt_res'):
            pooled[f'{m}__{key}'] = np.concatenate(
                [packs[ml][m][key] for ml in packs if m in packs[ml]]
                or [np.zeros(0, dtype=np.float32)])
    np.savez_compressed(OUT_RES, **pooled)
    print(f'wrote {OUT_RES}')


if __name__ == '__main__':
    main()
