from __future__ import annotations

import argparse
import json
from pathlib import Path

from .export import bundled_snapshot, generate_exports


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Generate eight HANSEN Database exports as TSV files."
    )
    parser.add_argument(
        "--snapshot",
        type=Path,
        default=None,
        help="Input compressed JSON data file. The archived HANSEN data are used by default.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("outputs"),
        help="Output directory. The default directory is outputs.",
    )
    args = parser.parse_args()

    if args.snapshot is None:
        with bundled_snapshot() as snapshot:
            manifest = generate_exports(snapshot, args.output)
    else:
        if not args.snapshot.exists():
            raise SystemExit(f"The HANSEN data file was not found at {args.snapshot}")
        manifest = generate_exports(args.snapshot, args.output)
    print(json.dumps(manifest["tables"], indent=2))


if __name__ == "__main__":
    main()
