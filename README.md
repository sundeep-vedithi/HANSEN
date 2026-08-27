# HANSEN

## An Integrated Structural and Functional Resource for the Proteome of *Mycobacterium leprae*

HANSEN is named after Gerhard Henrik Armauer Hansen, who first identified the
causative organism of leprosy. The database integrates genomic identifiers,
curated functional annotations, protein structure models, ligand information,
binding pocket predictions, B cell epitope predictions, gene essentiality
estimates and target prioritisation results for the *Mycobacterium leprae*
proteome. It supports biomarker research, molecular epidemiology, pathobiology
and the selection of diagnostic and therapeutic targets using structural
evidence.

The HANSEN Database is available at
<https://hansen-leprosy.medschl.cam.ac.uk/>.

This repository contains an archived HANSEN data extract and the Python
programs required to generate eight tab separated data files. It does not
contain the web application, protein sequences, structure coordinate files or
model parameters.

## Installation and data export

Python 3.10 or a later version is required. The export program uses only the
Python standard library.

```bash
git clone https://github.com/sundeep-vedithi/HANSEN.git
cd HANSEN
python -m venv .venv
source .venv/bin/activate
python -m pip install -e .
hansen-export --output outputs
```

The command creates eight TSV files and `outputs/manifest.json`. The manifest
records the number of rows and the SHA 256 value for each file.

## HANSEN Database exports

### S1 All Models

`S1_All_Models.tsv` contains 7132 structure model records. Models are ordered
by the recorded model score.

### S2 Oligomer Models

`S2_Oligomer_Models.tsv` contains 726 oligomer model records.

### S3 Oligomer Ligands

`S3_Oligomer_Ligands.tsv` contains 2439 curated ligand and cofactor records for
oligomer models.

### S4 Pockets

`S4_Pockets.tsv` contains AF2Bind, P2Rank and fpocket summaries for 1603
proteins.

### S5 B Cell Epitopes

`S5_BCell_Epitopes.tsv` contains DiscoTope 3.0 summaries for 1603 proteins.

### S6 Target Priority

`S6_Target_Priority.tsv` contains the HANSEN Target Priority Score and the
associated component measurements for 1603 proteins.

### S7 Essentiality

`S7_Essentiality.tsv` contains ProteomeLM essentiality estimates for 1603
proteins.

### S8 Combined

`S8_Combined.tsv` contains an integrated summary of the available measurements
for 1603 proteins.

The expected row counts and SHA 256 values are recorded in
`reference/expected_manifest.json`.

## Programs

The `hansen-export` command reads the archived HANSEN data and creates the
eight TSV files together with the output manifest.

The `src/hansen_database/export_hansen_data.py` program defines the command
options and accepts either the archived data supplied with this repository or
a compatible file provided by the user.

The `src/hansen_database/generate_exports.py` program calculates the stated
summary values, orders the HANSEN records, writes the TSV files and calculates
their SHA 256 values.

The `scripts/prepare_hansen_data.py` program prepares the archived data from
the HANSEN server database, statistics file and validation workbook. Protein
sequences, structure coordinates, residue measurements, server logs, access
credentials and web application code are excluded.

The `tests/test_hansen_exports.py` program confirms the number of rows and SHA
256 value for every export and checks that excluded fields are absent.

## Verification

```bash
python -m unittest discover -s tests -v
```

The tests generate all eight exports and compare their row counts and SHA 256
values with the HANSEN reference manifest.

## Archived HANSEN data

The archived data are stored in
`src/hansen_database/data/hansen_snapshot.json.gz.b64.part-*`. The four text files
contain a Base64 representation of one compressed JSON file. The export
program joins and decodes these files without changing the recorded bytes.

The archive contains protein display information, 7122 model records, 2439
curated oligomer ligand records and the recorded measurements for 1603
proteins. It excludes amino acid sequences, structure coordinate files,
residue measurements, access credentials, logs, backups and web application
code. The archive also records the SHA 256 values for the source database, the
statistics cache and the validation workbook.

The construction of each export is described in
[`docs/DATABASE_EXPORTS.md`](docs/DATABASE_EXPORTS.md). Definitions of the
columns are provided in
[`docs/DATA_DICTIONARY.md`](docs/DATA_DICTIONARY.md).

## Use of HANSEN data

The TSV files can be opened in Excel or read with R, Python and other
statistical software. The ML locus identifier and UniProt accession should be
retained when HANSEN records are joined to external resources. Missing values
are represented by empty fields and must not be interpreted as zero.

HANSEN values are computational predictions intended for hypothesis
generation and research prioritisation. They are not experimental findings
and they must not be used for clinical, diagnostic or therapeutic decisions.
Biological conclusions require independent experimental assessment. Users
must cite the relevant source databases and computational methods when HANSEN
data are analysed or redistributed.

## GitHub and Zenodo

GitHub should be used for the current programs, documentation and change
history. Zenodo should be connected to the repository when a fixed version is
ready. A tagged GitHub release can then be preserved by Zenodo and assigned a
DOI. The two services therefore have complementary purposes.

The GitHub guidance for preserving a repository with Zenodo is available at
<https://docs.github.com/en/repositories/archiving-a-github-repository/referencing-and-citing-content>.

## Citation

Citation information is provided in `CITATION.cff`. A Zenodo DOI can be added
to that file after the first archived release has been created.

## Licences

The programs are released under the MIT License. The archived HANSEN data and
the generated tables are provided under the conditions stated in
[`DATA_LICENSE.md`](DATA_LICENSE.md). Conditions imposed by the source
databases and computational methods remain applicable.
