#!/usr/bin/env python3
"""
run_bmgap2_bmscan.py — species identification for Neisseria / H. influenzae samples.

Usage:
    run_bmgap2_bmscan.py <mypath> <mlst_file> <bmgap2_db> <cpus>

Reads MLST scheme from mlst_file to determine whether to run.
Exits cleanly (code 0) if the sample is not a meningitis species.
"""

import os
import subprocess
import sys


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

    items       = mypath.strip().split("/")
    sample_name = items[-1]

    scheme = ""
    try:
        with open(mlst_file) as f:
            fields = f.readline().strip().split()
            if len(fields) > 1:
                scheme = fields[1]
    except Exception as e:
        print(f"BMGAP2-BMScan: Error reading MLST file - {e}", file=sys.stderr)
        sys.exit(1)

    if scheme not in ['neisseria', 'hinfluenzae']:
        print(
            f"BMGAP2-BMScan: Skipping {sample_name} "
            f"- not a meningitis species (scheme={scheme})"
        )
        sys.exit(0)

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

    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=False)

        if result.returncode != 0:
            print("BMGAP2-BMScan: Error running BMScan", file=sys.stderr)
            print(f"STDOUT: {result.stdout}", file=sys.stderr)
            print(f"STDERR: {result.stderr}", file=sys.stderr)
            sys.exit(1)

        print(
            f"BMGAP2-BMScan: Successfully completed species identification for {sample_name}"
        )
        print(result.stdout)

    except Exception as e:
        print(f"BMGAP2-BMScan: Exception - {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
