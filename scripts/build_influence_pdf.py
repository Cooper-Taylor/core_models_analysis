#!/usr/bin/env python3
"""Deprecated shim: the report is now built by build_influence_report.py.

The PDF and the Markdown are two renderings of one content definition, so they
cannot disagree. This file keeps the older command working.

    python3 scripts/build_influence_report.py              # both
    python3 scripts/build_influence_report.py --format pdf
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from build_influence_report import main as _main  # noqa: E402

if __name__ == "__main__":
    print("note: build_influence_pdf.py is a shim for build_influence_report.py --format pdf")
    sys.argv = [sys.argv[0], "--format", "pdf", *sys.argv[1:]]
    sys.exit(_main())
