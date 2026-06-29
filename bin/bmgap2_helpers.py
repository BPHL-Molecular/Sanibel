#!/usr/bin/env python3
"""
bmgap2_helpers.py — shared boilerplate for the run_bmgap2_*.py host scripts.

The three BMGAP2 steps (AMR, LocusExtractor, BMScan) each derive the sample name
from the output path and read the MLST scheme to decide whether to run, skipping
(exit 0) any sample that is not Neisseria meningitidis or H. influenzae. This
per-script re-check is the intended safety net even though the pipeline already
gates BMGAP2 on genus upstream — keep it.
"""

import sys

MENINGITIS_SCHEMES = ('neisseria', 'hinfluenzae')


def read_scheme_and_sample(mypath, mlst_file, tool_tag):
    """Return (scheme, sample_name) for a BMGAP2 step.

    Exits 1 if the MLST file cannot be read; exits 0 with a skip message if the
    sample's MLST scheme is not a meningitis species. Behavior and messages match
    the original per-script logic exactly (tool_tag is the per-step prefix, e.g.
    'BMGAP2-AMR').
    """
    sample_name = mypath.strip().split('/')[-1]

    scheme = ''
    try:
        with open(mlst_file) as f:
            fields = f.readline().strip().split()
            if len(fields) > 1:
                scheme = fields[1]
    except Exception as e:
        print(f"{tool_tag}: Error reading MLST file - {e}", file=sys.stderr)
        sys.exit(1)

    if scheme not in MENINGITIS_SCHEMES:
        print(f"{tool_tag}: Skipping {sample_name} "
              f"- not a meningitis species (scheme={scheme})")
        sys.exit(0)

    return scheme, sample_name
