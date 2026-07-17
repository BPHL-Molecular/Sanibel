#!/usr/bin/env python3
"""
bmgap2_helpers.py - shared boilerplate for the run_bmgap2_*.py host scripts.

The three BMGAP2 steps (AMR, LocusExtractor, BMScan) each derive the sample name
from the output path and read the MLST scheme to decide whether to run, skipping
(exit 0) any sample that is not Neisseria meningitidis or H. influenzae.
"""

import subprocess
import sys

MENINGITIS_SCHEMES = ('neisseria', 'hinfluenzae')


def run_tool(cmd, tool_tag, cwd=None, on_fail='fail'):
    try:
        result = subprocess.run(cmd, cwd=cwd, capture_output=True, text=True, check=False)
    except Exception as e:
        print(f"{tool_tag}: Unexpected error - {e}", file=sys.stderr)
        sys.exit(1)

    if result.returncode != 0:
        print(f"{tool_tag}: subtool exited with code {result.returncode}", file=sys.stderr)
        if result.stderr:
            print(result.stderr, file=sys.stderr)
        if on_fail == 'fail':
            sys.exit(1)
        if on_fail == 'continue':
            print(f"{tool_tag}: continuing pipeline despite failure")
            sys.exit(0)

    if result.stdout:
        print(result.stdout)
    return result


def read_scheme_and_sample(mypath, mlst_file, tool_tag):
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
