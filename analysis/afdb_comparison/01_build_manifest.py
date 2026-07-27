#!/usr/bin/env python3
"""Build the AFDB <-> HANSEN monomer manifest.

Outputs manifest.tsv with one row per M. leprae protein:
  ml_id, uniprot, entry_name, gene, protein_name, length, reviewed,
  afdb_path, af3_path, boltz_path, boltz2_path, chai_path, sequence
"""
import os
import sqlite3
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from paths import STATIC, DB, AFDB_DIR as AFDB, MANIFEST as OUT, check_inputs

check_inputs()

METHODS = ['af3', 'boltz', 'boltz2', 'chai']

con = sqlite3.connect(f'file:{DB}?mode=ro', uri=True)
con.row_factory = sqlite3.Row

cand = {}
for r in con.execute("""
        SELECT ml_id, entry, entry_name, gene_names, protein_names, length,
               reviewed, sequence, transmembrane, signal_peptide, protein_families
        FROM protein_characteristics
        WHERE ml_id IS NOT NULL AND TRIM(ml_id) <> ''"""):
    cand.setdefault(r['ml_id'], []).append(dict(r))

models = {}
for r in con.execute("""
        SELECT ml_id, method_key, model, path, ranking_score
        FROM protein_structure_models
        WHERE assembly_type = 'monomer'"""):
    models.setdefault(r['ml_id'], {})[r['method_key']] = dict(r)
con.close()

afdb_by_acc = {}
for fn in os.listdir(AFDB):
    if fn.endswith('.cif.gz'):
        acc = fn.split('-')[1]
        afdb_by_acc[acc] = os.path.join(AFDB, fn)

# --- disambiguate ML loci that carry more than one UniProt accession -------
# The HANSEN monomer models were built from one specific sequence; pick the
# accession whose UniProt sequence matches the modelled chain exactly.
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from structlib import parse_ca  # noqa: E402

chars = {}
ambiguous = []
for ml, lst in cand.items():
    if len(lst) == 1:
        chars[ml] = lst[0]
        continue
    modseq = None
    for m in METHODS:
        p = models.get(ml, {}).get(m, {}).get('path')
        if p and os.path.exists(os.path.join(STATIC, p)):
            try:
                modseq = parse_ca(os.path.join(STATIC, p))['seq']
                break
            except Exception:
                pass
    pick = None
    if modseq:
        for c in lst:
            if (c.get('sequence') or '') == modseq:
                pick = c
                break
        if pick is None:
            for c in lst:
                if c.get('length') == len(modseq):
                    pick = c
                    break
    if pick is None:
        pick = sorted(lst, key=lambda c: (c.get('reviewed') != 'reviewed',
                                          -(c.get('length') or 0)))[0]
    ambiguous.append((ml, [c['entry'] for c in lst], pick['entry'],
                      len(modseq) if modseq else None))
    chars[ml] = pick
    # the non-selected accessions of an ambiguous locus still have an AFDB
    # model but no HANSEN monomer of their own -> recorded, not compared
    for c in lst:
        if c['entry'] != pick['entry']:
            chars.setdefault('__unpaired__', []).append(
                (ml, c['entry'], c['entry_name'], c['length']))

unpaired = chars.pop('__unpaired__', [])

rows = []
stats = {'no_afdb': [], 'no_uniprot': [], 'missing_files': []}

all_ml = sorted(set(list(chars.keys()) + list(models.keys())))
for ml in all_ml:
    c = chars.get(ml, {})
    acc = c.get('entry')
    afdb = afdb_by_acc.get(acc) if acc else None
    mm = models.get(ml, {})
    row = {
        'ml_id': ml,
        'uniprot': acc or '',
        'entry_name': c.get('entry_name') or '',
        'gene': (c.get('gene_names') or '').split(' ')[0],
        'protein_name': (c.get('protein_names') or '').replace('\t', ' ').replace('\n', ' '),
        'length': c.get('length') or '',
        'reviewed': c.get('reviewed') or '',
        'transmembrane': 1 if (c.get('transmembrane') or '').strip() not in ('', 'None') else 0,
        'signal_peptide': 1 if (c.get('signal_peptide') or '').strip() not in ('', 'None') else 0,
        'protein_families': (c.get('protein_families') or '').replace('\t', ' '),
        'afdb_path': afdb or '',
        'sequence': c.get('sequence') or '',
    }
    for m in METHODS:
        p = mm.get(m, {}).get('path')
        full = os.path.join(STATIC, p) if p else ''
        if full and not os.path.exists(full):
            stats['missing_files'].append((ml, m, p))
            full = ''
        row[f'{m}_path'] = full
        row[f'{m}_model'] = mm.get(m, {}).get('model') or ''
        row[f'{m}_rank'] = mm.get(m, {}).get('ranking_score') or ''
    if not acc:
        stats['no_uniprot'].append(ml)
    elif not afdb:
        stats['no_afdb'].append((ml, acc))
    rows.append(row)

cols = (['ml_id', 'uniprot', 'entry_name', 'gene', 'protein_name', 'length',
         'reviewed', 'transmembrane', 'signal_peptide', 'protein_families',
         'afdb_path']
        + [f'{m}_{k}' for m in METHODS for k in ('path', 'model', 'rank')]
        + ['sequence'])

with open(OUT, 'w') as fh:
    fh.write('\t'.join(cols) + '\n')
    for r in rows:
        fh.write('\t'.join(str(r.get(c, '')) for c in cols) + '\n')

n_full = sum(1 for r in rows if r['afdb_path'] and all(r[f'{m}_path'] for m in METHODS))
print(f'proteins in manifest      : {len(rows)}')
print(f'with AFDB model           : {sum(1 for r in rows if r["afdb_path"])}')
for m in METHODS:
    print(f'with {m:7s} monomer      : {sum(1 for r in rows if r[f"{m}_path"])}')
print(f'complete 5-way (AFDB+4)   : {n_full}')
print(f'no UniProt accession      : {stats["no_uniprot"]}')
print(f'ambiguous ML loci resolved: {len(ambiguous)}')
for a in ambiguous:
    print('   ', a)
print(f'AFDB entries with no paired HANSEN monomer (2nd accession of a shared locus): {len(unpaired)}')
for u in unpaired:
    print('   ', u)
print(f'UniProt but no AFDB model : {stats["no_afdb"]}')
print(f'DB path missing on disk   : {len(stats["missing_files"])}')
for x in stats['missing_files'][:20]:
    print('   ', x)
print(f'\nwrote {OUT}')
