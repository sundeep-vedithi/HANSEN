# HANSEN Database exports and data provenance

## Purpose

This repository generates eight HANSEN Database exports from an archived and
versioned HANSEN data extract. The recorded values, column order, row order and
number of rows are preserved.

The repository does not regenerate the original protein structure predictions.
A complete calculation requires the relevant software licences, model
parameters, computing environment and structure coordinate archive.

## Source records

The archived data extract was prepared from the HANSEN server directory
`/data/new-hansen`. The source database is
`protein_characteristics.db` with SHA 256 value
`0c3f786daa193530f2005f4f4009ed88f1624d2fe7231d1988687d0ed23db741`.
The statistics file is
`static/statistics_cache/hansen_statistics_advanced_v4.json` with SHA 256 value
`64f70f38c12741854efec0c0537699d41c7f35628d529188d13434b7852bd1cc`.
The validation workbook has SHA 256 value
`3bd3d5eb859ee8e36f00e63423f0756e601d30a38d409cf144108864ed1310d7`.

These values are stored within the archived data extract.
`scripts/prepare_hansen_data.py` reads the three sources and prepares the data used
for exports S1 to S8.

## Construction of each export

### S1 All Models

The 7132 model records are linked to their protein annotations. The selection
value uses pLDDT when it is available, followed by confidence and then the
model ranking score. Records are ordered from the highest available value to
the lowest. Records without a score are retained at the end.

### S2 Oligomer Models

Oligomer records from S1 are retained and ordered by method, chain count, ML
locus and database row. The export includes 22 AlphaFold 3 homomer records
that are oligomers in the source database.

### S3 Oligomer Ligands

The export contains 2439 curated ligand and cofactor records. Duplicate ligand
records are removed. The retained information includes template and assembly
identifiers, Chemical Component Dictionary codes, ligand names, SMILES,
assembly evidence and the source record.

### S4 Pockets

The mean pocket value for each of the 1603 proteins is the arithmetic mean of
the available AF2Bind, P2Rank and fpocket averages. A method without a value is
not included in the denominator. Proteins are ordered by the mean pocket value.

### S5 B Cell Epitopes

All 1603 proteins are ordered by the maximum DiscoTope 3.0 value and then by
the mean value. The number of residues examined and the number above the
threshold are retained.

### S6 Target Priority

All proteins are ordered by the recorded HANSEN Target Priority Score. The
export retains the component measurements, priority category, ProteomeLM
essentiality estimate, model quality values, pocket measurements, DiscoTope
measurements, structure methods and assembly classes.

### S7 Essentiality

All proteins are ordered by `proteomelm_ess_probability` from the HANSEN
essentiality table. The classification recorded in HANSEN is retained.

### S8 Combined

Records are ordered by ML locus. `Data Types Present` reports the number of
available values among the Target Priority Score, Essentiality Score, mean
pLDDT, AF2Bind average, P2Rank average, fpocket average and DiscoTope 3.0
average.

The HANSEN example score is calculated from the available values using the
following expression.

```text
Target Priority Score
+ 100 × Essentiality Score
+ mean pLDDT
+ 10 × (AF2Bind + P2Rank + fpocket + DiscoTope 3.0)
```

The example score is provided only to assist exploration of the database. It
is not a validated biological or clinical measure.

## Complete computational generation

A complete calculation begins with the *Mycobacterium leprae* strain TN
UniProt reference proteome UP000000806, which contains 1603 protein entries.
Monomer and oligomer models are generated with the recorded AlphaFold 3,
Boltz, Boltz 2 and Chai configurations. Model ranking, pLDDT and PAE values are
then calculated and entered into the database. AF2Bind, P2Rank and fpocket
provide the binding pocket measurements. DiscoTope 3.0 provides the B cell
epitope measurements. ProteomeLM L uses the recorded ESMC 600M representation
to estimate gene essentiality. HANSEN then calculates the Target Priority
Score and updates the statistics data before the eight exports are generated.

The database records P2Rank 2.5.1, fpocket 4.1, DiscoTope 3.0, ProteomeLM L and
ESMC 600M for versioned records. Some earlier binding pocket records do not
contain a software version. The relevant HANSEN computing environments should
therefore be retained when a complete calculation is archived.

Structure coordinate archives and residue measurements are large scientific
data objects. They should be preserved in a suitable data repository and not
stored directly in GitHub.
