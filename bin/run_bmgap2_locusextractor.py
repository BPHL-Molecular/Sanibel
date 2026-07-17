#!/usr/bin/env python3
"""
run_bmgap2_locusextractor.py - vaccine antigen typing for Neisseria / H. influenzae samples.

Usage:
    run_bmgap2_locusextractor.py <mypath> <mlst_file> <bmgap2_db>

Reads MLST scheme from mlst_file to determine whether to run.
Exits cleanly (code 0) if the sample is not a meningitis species.
"""

import glob
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.realpath(__file__)))
from bmgap2_helpers import read_scheme_and_sample, run_tool


def main():
    if len(sys.argv) != 4:
        print(
            f"Usage: {sys.argv[0]} <mypath> <mlst_file> <bmgap2_db>",
            file=sys.stderr
        )
        sys.exit(1)

    mypath    = sys.argv[1]
    mlst_file = sys.argv[2]
    bmgap2_db = sys.argv[3]

    _scheme, sample_name = read_scheme_and_sample(
        mypath, mlst_file, "BMGAP2-LocusExtractor"
    )

    assembly_dir  = f"{mypath}/assembly"
    output_dir    = f"{mypath}/bmgap2_locusextractor"
    le_dir        = os.path.join(bmgap2_db, "locusextractor")
    le_script     = os.path.join(le_dir, "LocusExtractor_RHEL8.py")
    settings_dir  = os.path.join(le_dir, "settings")

    if not os.path.isdir(assembly_dir):
        print(
            f"BMGAP2-LocusExtractor: Error - Assembly directory not found: {assembly_dir}",
            file=sys.stderr
        )
        sys.exit(1)

    if not os.path.isfile(le_script):
        print(
            f"BMGAP2-LocusExtractor: Error - LocusExtractor script not found: {le_script}",
            file=sys.stderr
        )
        sys.exit(1)

    os.makedirs(output_dir, exist_ok=True)

    print(f"BMGAP2-LocusExtractor: Running vaccine antigen typing for {sample_name}")

    cmd = [
        sys.executable,
        le_script,
        "--no_update",
        "-s", settings_dir,
        "-p", sample_name,
        assembly_dir,
    ]

    run_tool(cmd, "BMGAP2-LocusExtractor", cwd=output_dir, on_fail='warn')

    le_dirs = glob.glob(os.path.join(output_dir, f"LE_*_{sample_name}_*"))
    if le_dirs:
        expected = os.path.join(
            le_dirs[0], "Results_text", f"molecular_data_{sample_name}.csv"
        )
        if os.path.isfile(expected):
            print(f"BMGAP2-LocusExtractor: Successfully generated {expected}")
        else:
            print(
                f"BMGAP2-LocusExtractor: Warning - Expected output not found: {expected}",
                file=sys.stderr
            )
    else:
        print(
            "BMGAP2-LocusExtractor: Warning - LocusExtractor output directory not found",
            file=sys.stderr
        )


if __name__ == "__main__":
    main()
