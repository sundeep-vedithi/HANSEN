from __future__ import annotations

import base64
import csv
import gzip
import hashlib
import json
import math
import tempfile
from contextlib import contextmanager
from importlib.resources import files
from pathlib import Path
from typing import Any, Iterable, Iterator


TABLE_COLUMNS = {
    "S1_All_Models": [
        "ML ID",
        "Gene",
        "Protein",
        "UniProt",
        "Method",
        "Assembly",
        "Chain Count",
        "pLDDT",
        "PAE",
        "Confidence",
        "Ranking",
        "Ligand Metadata",
        "Model Path",
        "Pathway",
        "GO Terms",
        "Database Row",
        "Source Table",
        "Example Selection Score",
    ],
    "S2_Oligomer_Models": [
        "ML ID",
        "Gene",
        "Protein",
        "UniProt",
        "Method",
        "Assembly",
        "Chain Count",
        "pLDDT",
        "PAE",
        "Confidence",
        "Ranking",
        "Ligand Metadata",
        "Model Path",
        "Pathway",
        "GO Terms",
        "Database Row",
        "Source Table",
    ],
    "S3_Oligomer_Ligands": [
        "ML ID",
        "Gene",
        "Protein",
        "UniProt",
        "Oligomer Kind",
        "Template PDB",
        "Assembly ID",
        "Ligand CCD",
        "Ligand Name",
        "Ligand SMILES",
        "Source Assembly (YAML)",
        "Record Source",
        "Pathway",
        "GO Terms",
    ],
    "S4_Pockets": [
        "ML ID",
        "Gene",
        "Protein",
        "UniProt",
        "Target Priority Score",
        "Mean pLDDT",
        "AF2Bind Average",
        "P2Rank Average",
        "F-Pocket Average",
        "Mean Pocket Score",
        "Pathway",
        "GO Terms",
        "Methods",
    ],
    "S5_BCell_Epitopes": [
        "ML ID",
        "Gene",
        "Protein",
        "UniProt",
        "Target Priority Score",
        "DiscoTope3 Average",
        "DiscoTope3 Max",
        "Residues Scored",
        "Above-Threshold Residues",
        "Pathway",
        "GO Terms",
        "Methods",
    ],
    "S6_Target_Priority": [
        "ML ID",
        "Gene",
        "Protein",
        "UniProt",
        "Target Priority Score",
        "Target Priority Tier",
        "Essentiality Score",
        "Essentiality Call",
        "Mean pLDDT",
        "AF2Bind Average",
        "P2Rank Average",
        "F-Pocket Average",
        "DiscoTope3 Average",
        "Pathway",
        "GO Terms",
        "Mean PAE",
        "Methods",
        "Assemblies",
    ],
    "S7_Essentiality": [
        "ML ID",
        "Gene",
        "Protein",
        "UniProt",
        "Target Priority Score",
        "Essentiality Score",
        "Essentiality Call",
        "Pathway",
        "GO Terms",
        "Methods",
    ],
    "S8_Combined": [
        "ML ID",
        "Gene",
        "Protein",
        "UniProt",
        "Target Priority Score",
        "Essentiality Score",
        "Mean pLDDT",
        "AF2Bind Average",
        "P2Rank Average",
        "F-Pocket Average",
        "DiscoTope3 Average",
        "Pathway",
        "GO Terms",
        "Example Score",
        "Data Types Present",
        "Methods",
        "Assemblies",
    ],
}


def clean(value: Any) -> str:
    if value is None:
        return ""
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, float):
        if math.isnan(value) or math.isinf(value):
            return ""
        return format(value, ".15g")
    return str(value).replace("\r", " ").replace("\n", " ").strip()


def number(value: Any) -> float | None:
    text = clean(value)
    if not text or text.lower() in {"none", "null", "nan", "na", "n/a", "[]", "{}"}:
        return None
    try:
        result = float(text)
    except ValueError:
        return None
    return result if math.isfinite(result) else None


def rounded_mean(values: Iterable[Any]) -> float | str:
    numeric = [value for item in values if (value := number(item)) is not None]
    return round(sum(numeric) / len(numeric), 4) if numeric else ""


def load_snapshot(path: Path) -> dict[str, Any]:
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        snapshot = json.load(handle)
    if snapshot.get("schema_version") != "1.0":
        raise ValueError(
            f"The data schema {snapshot.get('schema_version')!r} is not supported"
        )
    return snapshot


@contextmanager
def bundled_snapshot() -> Iterator[Path]:
    data_directory = files("hansen_repro").joinpath("data")
    parts = sorted(
        item
        for item in data_directory.iterdir()
        if item.name.startswith("hansen_snapshot.json.gz.b64.part-")
    )
    if not parts:
        raise FileNotFoundError("Bundled HANSEN snapshot parts were not found")

    encoded = "".join(part.read_text(encoding="ascii").strip() for part in parts)
    payload = base64.b64decode(encoded, validate=True)
    with tempfile.NamedTemporaryFile(suffix=".json.gz") as handle:
        handle.write(payload)
        handle.flush()
        yield Path(handle.name)


def model_tables(snapshot: dict[str, Any]) -> tuple[list[dict], list[dict]]:
    models = []
    for source in snapshot["models"]:
        row = dict(source)
        score = (
            number(row.get("pLDDT"))
            or number(row.get("Confidence"))
            or number(row.get("Ranking"))
        )
        row["Example Selection Score"] = score if score is not None else ""
        models.append(row)

    models.sort(
        key=lambda row: (
            number(row.get("Example Selection Score")) is not None,
            number(row.get("Example Selection Score")) or float("-inf"),
        ),
        reverse=True,
    )

    oligomers = [
        dict(row)
        for row in snapshot["models"]
        if "oligomer" in clean(row.get("Assembly")).lower()
    ]
    oligomers.sort(
        key=lambda row: (
            clean(row.get("Method")),
            number(row.get("Chain Count")) or 0,
            clean(row.get("ML ID")),
            number(row.get("Database Row")) or 0,
        )
    )
    return models, oligomers


def metric_tables(snapshot: dict[str, Any]) -> dict[str, list[dict]]:
    metadata = {row["ML ID"]: row for row in snapshot["metadata"]}
    metrics = snapshot["protein_metrics"]

    pockets: list[dict] = []
    epitopes: list[dict] = []
    priorities: list[dict] = []
    essentiality: list[dict] = []
    combined: list[dict] = []

    for raw in metrics:
        ml_id = clean(raw.get("ml_id")).upper()
        meta = {
            "ML ID": ml_id,
            "Gene": clean(raw.get("gene")),
            "Protein": clean(raw.get("protein")),
            "UniProt": clean(raw.get("entry")),
            "Pathway": "",
            "GO Terms": "",
        }
        meta.update(metadata.get(ml_id, {}))

        priority = number(raw.get("target_priority_score"))
        essentiality_score = number(raw.get("essentiality_score"))
        plddt = number(raw.get("mean_plddt"))
        af2bind = number(raw.get("af2bind_avg"))
        p2rank = number(raw.get("p2rank_avg"))
        fpocket = number(raw.get("fpocket_avg"))
        discotope = number(raw.get("discotope3_avg"))

        pockets.append(
            {
                **meta,
                "Target Priority Score": priority if priority is not None else "",
                "Mean pLDDT": plddt if plddt is not None else "",
                "AF2Bind Average": af2bind if af2bind is not None else "",
                "P2Rank Average": p2rank if p2rank is not None else "",
                "F-Pocket Average": fpocket if fpocket is not None else "",
                "Mean Pocket Score": rounded_mean([af2bind, p2rank, fpocket]),
                "Methods": clean(raw.get("methods")),
            }
        )

        epitopes.append(
            {
                **meta,
                "Target Priority Score": priority if priority is not None else "",
                "DiscoTope3 Average": (
                    discotope if discotope is not None else ""
                ),
                "DiscoTope3 Max": clean(raw.get("discotope3_max")),
                "Residues Scored": clean(raw.get("discotope3_residue_count")),
                "Above-Threshold Residues": clean(
                    raw.get("discotope3_above_threshold_count")
                ),
                "Methods": clean(raw.get("methods")),
            }
        )

        priorities.append(
            {
                **meta,
                "Target Priority Score": priority if priority is not None else "",
                "Target Priority Tier": clean(raw.get("target_priority_tier")),
                "Essentiality Score": (
                    essentiality_score if essentiality_score is not None else ""
                ),
                "Essentiality Call": clean(raw.get("essentiality_call")),
                "Mean pLDDT": plddt if plddt is not None else "",
                "AF2Bind Average": af2bind if af2bind is not None else "",
                "P2Rank Average": p2rank if p2rank is not None else "",
                "F-Pocket Average": fpocket if fpocket is not None else "",
                "DiscoTope3 Average": (
                    discotope if discotope is not None else ""
                ),
                "Mean PAE": clean(raw.get("mean_pae")),
                "Methods": clean(raw.get("methods")),
                "Assemblies": clean(raw.get("assemblies")),
            }
        )

        essentiality.append(
            {
                **meta,
                "Target Priority Score": priority if priority is not None else "",
                "Essentiality Score": (
                    essentiality_score if essentiality_score is not None else ""
                ),
                "Essentiality Call": clean(raw.get("essentiality_call")),
                "Methods": clean(raw.get("methods")),
            }
        )

        evidence = [priority, essentiality_score, plddt, af2bind, p2rank, fpocket, discotope]
        evidence_count = sum(value is not None for value in evidence)
        example_score = 0.0
        if priority is not None:
            example_score += priority
        if essentiality_score is not None:
            example_score += 100 * essentiality_score
        if plddt is not None:
            example_score += plddt
        for value in (af2bind, p2rank, fpocket):
            if value is not None:
                example_score += 10 * value
        if discotope is not None:
            example_score += 10 * discotope

        combined.append(
            {
                **meta,
                "Target Priority Score": priority if priority is not None else "",
                "Essentiality Score": (
                    essentiality_score if essentiality_score is not None else ""
                ),
                "Mean pLDDT": plddt if plddt is not None else "",
                "AF2Bind Average": af2bind if af2bind is not None else "",
                "P2Rank Average": p2rank if p2rank is not None else "",
                "F-Pocket Average": fpocket if fpocket is not None else "",
                "DiscoTope3 Average": (
                    discotope if discotope is not None else ""
                ),
                "Example Score": round(example_score, 3),
                "Data Types Present": evidence_count,
                "Methods": clean(raw.get("methods")),
                "Assemblies": clean(raw.get("assemblies")),
            }
        )

    pockets.sort(
        key=lambda row: (
            number(row.get("Mean Pocket Score")) is not None,
            number(row.get("Mean Pocket Score")) or float("-inf"),
        ),
        reverse=True,
    )
    epitopes.sort(
        key=lambda row: (
            number(row.get("DiscoTope3 Max")) or float("-inf"),
            number(row.get("DiscoTope3 Average")) or float("-inf"),
        ),
        reverse=True,
    )
    priorities.sort(
        key=lambda row: number(row.get("Target Priority Score")) or float("-inf"),
        reverse=True,
    )
    essentiality.sort(
        key=lambda row: number(row.get("Essentiality Score")) or float("-inf"),
        reverse=True,
    )
    combined.sort(key=lambda row: clean(row.get("ML ID")))

    return {
        "S4_Pockets": pockets,
        "S5_BCell_Epitopes": epitopes,
        "S6_Target_Priority": priorities,
        "S7_Essentiality": essentiality,
        "S8_Combined": combined,
    }


def write_tsv(path: Path, columns: list[str], rows: Iterable[dict[str, Any]]) -> int:
    count = 0
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=columns,
            delimiter="\t",
            extrasaction="ignore",
            lineterminator="\n",
        )
        writer.writeheader()
        for row in rows:
            writer.writerow({column: clean(row.get(column, "")) for column in columns})
            count += 1
    return count


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def generate_exports(snapshot_path: Path, output_dir: Path) -> dict[str, Any]:
    snapshot = load_snapshot(snapshot_path)
    output_dir.mkdir(parents=True, exist_ok=True)

    all_models, oligomers = model_tables(snapshot)
    tables = {
        "S1_All_Models": all_models,
        "S2_Oligomer_Models": oligomers,
        "S3_Oligomer_Ligands": snapshot["oligomer_ligands"],
        **metric_tables(snapshot),
    }

    manifest = {
        "snapshot": {
            "path": snapshot_path.name,
            "sha256": sha256(snapshot_path),
            "source": snapshot["source"],
        },
        "tables": {},
    }
    for name, columns in TABLE_COLUMNS.items():
        path = output_dir / f"{name}.tsv"
        rows = write_tsv(path, columns, tables[name])
        manifest["tables"][name] = {
            "file": path.name,
            "rows": rows,
            "sha256": sha256(path),
        }

    manifest_path = output_dir / "manifest.json"
    manifest_path.write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    return manifest
