from __future__ import annotations

import json
import sys
import tempfile
import unittest
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from hansen_repro.export import (  # noqa: E402
    bundled_snapshot,
    generate_exports,
    load_snapshot,
)


class ExportTest(unittest.TestCase):
    def test_all_tables_match_reference_manifest(self) -> None:
        expected = json.loads(
            (PROJECT_ROOT / "reference" / "expected_manifest.json").read_text()
        )
        with (
            bundled_snapshot() as snapshot,
            tempfile.TemporaryDirectory() as temporary_directory,
        ):
            actual = generate_exports(snapshot, Path(temporary_directory))

        self.assertEqual(expected["tables"], actual["tables"])
        self.assertEqual(expected["snapshot"]["sha256"], actual["snapshot"]["sha256"])

    def test_archived_data_are_complete_and_restricted(self) -> None:
        with bundled_snapshot() as snapshot:
            data = load_snapshot(snapshot)

        self.assertEqual(7132, len(data["models"]))
        self.assertEqual(2439, len(data["oligomer_ligands"]))
        self.assertEqual(1603, len(data["protein_metrics"]))

        serialised_keys = json.dumps(
            sorted(
                {
                    key
                    for section in (
                        data["metadata"],
                        data["models"],
                        data["oligomer_ligands"],
                        data["protein_metrics"],
                    )
                    for row in section
                    for key in row
                }
            )
        ).lower()
        for excluded in ("sequence", "password", "api_key", "private_key"):
            self.assertNotIn(excluded, serialised_keys)


if __name__ == "__main__":
    unittest.main()
