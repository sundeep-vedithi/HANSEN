# Archived HANSEN data

The archived HANSEN data are stored in
`src/hansen_repro/data/hansen_snapshot.json.gz.b64.part-*`.

These four text files contain a Base64 representation of one compressed JSON
file. The Python package joins and decodes them without changing the recorded
bytes. The SHA 256 value is checked before the exports are generated.

The compressed file is approximately 0.6 MB and contains only the fields
required for the eight HANSEN TSV exports. The source records and the
construction of each export are described in
`docs/DATABASE_EXPORTS.md`.
