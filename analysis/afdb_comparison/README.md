# AFDB comparison

Code for Supplementary Note S1 of the HANSEN manuscript, which compares the
AlphaFold Protein Structure Database (AFDB) proteome release for
*Mycobacterium leprae* TN with the monomeric models held in HANSEN.

The comparison covers 1,598 proteins and 6,379 model pairs, one for each
combination of protein and prediction method (AlphaFold 3, Boltz-1, Boltz-2 and
Chai-1). It also computes all ten pairwise combinations of the five model
sources so that the AFDB models can be placed within the predictor ensemble.

## Inputs

Two inputs are required and neither is held in this repository.

| Variable | Contents |
|---|---|
| `HANSEN_BASE` | Directory holding `protein_characteristics.db` and `static/` from the HANSEN deployment |
| `AFDB_DIR` | Directory of extracted `AF-*-model_v6.cif.gz` files from the AFDB proteome archive `UP000000806_272631_MYCLE_v6.tar` |

The AFDB archive is available from <https://alphafold.ebi.ac.uk/download>.
Outputs are written to `WORK_DIR`, which defaults to this directory.

## Running

```bash
python -m venv .venv && . .venv/bin/activate
pip install -r requirements.txt
export HANSEN_BASE=/path/to/new-hansen
export AFDB_DIR=/path/to/afdb_cif
./run_all.sh
```

`run_all.sh` enforces the dependency order. `03_analyse.py` rewrites
`summary.json`, so the scripts that append to it must run after it and before
the report is built; `04_report.py` aborts if any expected block is missing
rather than producing a silently shortened document.

| Script | Purpose |
|---|---|
| `paths.py` | Resolves input and output locations from the environment |
| `structlib.py` | mmCIF Cα parser and all structural metrics |
| `01_build_manifest.py` | Pairs AFDB entries with HANSEN monomers by UniProt accession |
| `02_compare.py` | Per-pair comparison, writes `pairs.tsv` and `per_residue.npz` |
| `02b_crossmethod.py` | All ten pairwise combinations of the five model sources |
| `03_analyse.py` | Aggregates into tables and figures |
| `03b_worked_example.py` | Domain-resolved comparison for PknB (ML0016) |
| `03c_profiles.py` | Per-residue deviation profiles |
| `04_report.py` | Builds the supplementary note and the tables workbook |
| `05_validate.py` | Checks the TM-score implementation against TM-align |

## Metrics

`structlib.py` implements every metric directly from its primary definition and
depends only on NumPy.

- **TM-score**, Zhang and Skolnick (2004), sequence-dependent and normalised by
  the reference chain length, using the standard iterative fragment-seed search
- **GDT-TS and GDT-HA**, Zemla (2003), cutoff sets (1, 2, 4, 8) Å and
  (0.5, 1, 2, 4) Å
- **lDDT**, Mariani *et al.* (2013), on Cα atoms with a 15 Å inclusion radius
  and the four standard tolerance thresholds
- **Cα-RMSD** after Kabsch superposition
- **Secondary structure**, P-SEA (Labesse *et al.*, 1997), assigned from Cα
  geometry alone and compared as a three-state agreement
- **Needleman-Wunsch** global alignment, used only where construct boundaries
  differ between the two model sets

The TM-score implementation is checked against the reference TM-align program on
400 seeded random pairs by `05_validate.py`. Correlation is r = 0.996 to 0.998
with a median difference of zero. Residual differences run in one direction,
with TM-align scoring higher, which follows from TM-align optimising a
sequence-independent alignment.

## Notes on the comparison

Three properties of the data affect interpretation and are handled explicitly.

1. Four *M. leprae* loci carry two UniProt accessions each. The accession whose
   sequence matches the modelled chain exactly is selected, otherwise four AFDB
   entries are silently dropped.
2. Around 680 proteins differ from the AFDB sequence at position 1 only.
   UniProt records the initiator as methionine while the genomic translation
   retains valine or leucine at GTG and TTG start codons. This does not affect
   the residue correspondence.
3. 84 proteins have different modelled chain lengths in the two sets. TM-score
   is normalised by the reference length, so a truncated model is capped near
   its coverage fraction. Filter on `n_common / n_ref >= 0.95` before reading a
   low score as a fold difference.

## Reproducibility

The pipeline is deterministic apart from the seeded sample drawn for the
TM-align validation. Running it end to end regenerates every number, table and
figure in Supplementary Note S1 from the deposited coordinates.
