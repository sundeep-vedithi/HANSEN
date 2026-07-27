"""Input and output locations for the AFDB comparison pipeline.

Every path is resolved from an environment variable so that the pipeline can be
run against any local copy of the HANSEN data. Defaults are relative to this
directory, so a fresh checkout works once the two inputs are placed alongside
the scripts or the variables are exported.

    HANSEN_BASE   directory holding protein_characteristics.db and static/
                  (default ./new-hansen)
    AFDB_DIR      directory of extracted AlphaFold DB .cif.gz files
                  (default ./afdb_cif)
    WORK_DIR      where manifest.tsv, pairs.tsv, tables/ and figures/ are
                  written (default: this directory)
"""
import os

HERE = os.path.dirname(os.path.abspath(__file__))

HANSEN_BASE = os.environ.get('HANSEN_BASE', os.path.join(HERE, 'new-hansen'))
AFDB_DIR = os.environ.get('AFDB_DIR', os.path.join(HERE, 'afdb_cif'))
WORK_DIR = os.environ.get('WORK_DIR', HERE)

DB = os.path.join(HANSEN_BASE, 'protein_characteristics.db')
STATIC = os.path.join(HANSEN_BASE, 'static')
PRIOR_TSV = os.path.join(STATIC, 'proteomelm_ess',
                         'hansen_ai_target_prioritisation.tsv')

MANIFEST = os.path.join(WORK_DIR, 'manifest.tsv')
PAIRS = os.path.join(WORK_DIR, 'pairs.tsv')
CROSSMETHOD = os.path.join(WORK_DIR, 'crossmethod.tsv')
PER_RESIDUE = os.path.join(WORK_DIR, 'per_residue.npz')
SUMMARY = os.path.join(WORK_DIR, 'summary.json')
TABDIR = os.path.join(WORK_DIR, 'tables')
FIGDIR = os.path.join(WORK_DIR, 'figures')


def check_inputs():
    """Fail early with a useful message rather than deep inside a parser."""
    missing = []
    if not os.path.isfile(DB):
        missing.append(f'HANSEN database not found: {DB}')
    if not os.path.isdir(AFDB_DIR):
        missing.append(f'AlphaFold DB directory not found: {AFDB_DIR}')
    if missing:
        raise SystemExit('\n'.join(missing) + '\n\nSet HANSEN_BASE and AFDB_DIR, '
                         'for example:\n'
                         '  export HANSEN_BASE=/path/to/new-hansen\n'
                         '  export AFDB_DIR=/path/to/afdb_cif')
