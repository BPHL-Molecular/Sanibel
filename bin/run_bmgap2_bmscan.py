#!/usr/bin/env python3
"""
run_bmgap2_bmscan.py - species identification for Neisseria / H. influenzae samples.

Usage:
    run_bmgap2_bmscan.py <mypath> <mlst_file> <bmgap2_db> <cpus>

Reads MLST scheme from mlst_file to determine whether to run.
Exits cleanly (code 0) if the sample is not a meningitis species.
"""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.realpath(__file__)))
from bmgap2_helpers import read_scheme_and_sample, run_tool


def main():
    if len(sys.argv) != 5:
        print(
            f"Usage: {sys.argv[0]} <mypath> <mlst_file> <bmgap2_db> <cpus>",
            file=sys.stderr
        )
        sys.exit(1)

    mypath    = sys.argv[1]
    mlst_file = sys.argv[2]
    bmgap2_db = sys.argv[3]
    cpus      = sys.argv[4]

    _scheme, sample_name = read_scheme_and_sample(mypath, mlst_file, "BMGAP2-BMScan")

    assembly_dir = f"{mypath}/assembly"
    output_dir   = f"{mypath}/bmgap2_bmscan"
    bmscan       = os.path.join(bmgap2_db, "SpeciesDB", "bin", "identify_species.py")

    if not os.path.isdir(assembly_dir):
        print(
            f"BMGAP2-BMScan: Error - Assembly directory not found: {assembly_dir}",
            file=sys.stderr
        )
        sys.exit(1)

    if not os.path.isfile(bmscan):
        print(f"BMGAP2-BMScan: Error - BMScan script not found: {bmscan}", file=sys.stderr)
        sys.exit(1)

    os.makedirs(output_dir, exist_ok=True)

    print(f"BMGAP2-BMScan: Running species identification for {sample_name}")

    cmd = [
        sys.executable,
        bmscan,
        "-d", assembly_dir,
        "-o", output_dir,
        "-t", cpus,
        "-j",
    ]

    run_tool(cmd, "BMGAP2-BMScan", on_fail='fail')
    print(f"BMGAP2-BMScan: Successfully completed species identification for {sample_name}")


if __name__ == "__main__":
    main()
