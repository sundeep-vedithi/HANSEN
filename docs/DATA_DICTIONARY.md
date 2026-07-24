# HANSEN data dictionary

## Identifiers and annotations

### `ML ID`

This is the HANSEN *Mycobacterium leprae* locus identifier.

### `Gene`

This is the gene name associated with the locus.

### `Protein`

This is the protein description recorded in HANSEN.

### `UniProt`

This is the UniProt accession.

### `Pathway`

This is the pathway annotation recorded in the HANSEN source table.

### `GO Terms`

These are the Gene Ontology identifiers recorded for the protein. Multiple
identifiers are separated by semicolons within the data field.

## Structure model fields

### `Method`

This is the structure prediction method.

### `Assembly`

This is the monomer or oligomer classification.

### `Chain Count`

This is the number of chains resolved from the model.

### `pLDDT`

This is the model confidence summary when a model value is available.

### `PAE`

This is the predicted aligned error summary when a model value is available.

### `Confidence`

This is the confidence value supplied by the structure prediction method.

### `Ranking`

This is the ranking value supplied by the structure prediction method.

### `Model Path`

This is the path relative to the HANSEN model store.

### `Database Row`

This is the source model row identifier.

## Ligand fields

### `Oligomer Kind`

This identifies a homooligomer or heterooligomer. The number of subunits is
included when it is available.

### `Template PDB`

This is the Protein Data Bank template identifier.

### `Assembly ID`

This is the biological assembly identifier for the template.

### `Ligand CCD`

This is the Protein Data Bank Chemical Component Dictionary code.

### `Ligand Name`

This is the ligand or cofactor name.

### `Ligand SMILES`

This is the SMILES representation when it is available.

### `Source Assembly (YAML)`

This records the assembly evidence used during curation.

### `Record Source`

This identifies the HANSEN source table or associated source file.

## Protein measurements

### `Target Priority Score`

This is the integrated HANSEN prioritisation value between zero and one
hundred.

### `Target Priority Tier`

This is the HANSEN priority category.

### `Essentiality Score`

This is the ProteomeLM essentiality probability.

### `Essentiality Call`

This is the essentiality classification recorded in HANSEN.

### `Mean pLDDT`

This is the mean model pLDDT for the protein.

### `Mean PAE`

This is the mean predicted aligned error for the protein.

### `AF2Bind Average`

This is the mean AF2Bind value for the protein.

### `P2Rank Average`

This is the mean P2Rank value for the protein.

### `F-Pocket Average`

This is the mean fpocket value for the protein.

### `Mean Pocket Score`

This is the arithmetic mean of the available AF2Bind, P2Rank and fpocket
averages.

### `DiscoTope3 Average`

This is the mean DiscoTope 3.0 residue value.

### `DiscoTope3 Max`

This is the maximum DiscoTope 3.0 residue value.

### `Residues Scored`

This is the number of residues examined by DiscoTope 3.0.

### `Above-Threshold Residues`

This is the number of DiscoTope residues above the stated threshold.

### `Methods`

These are the structure methods represented for the protein.

### `Assemblies`

These are the assembly classes represented for the protein.

An empty field means that the corresponding value was not available in the
archived HANSEN data. An empty field must not be interpreted as zero.
