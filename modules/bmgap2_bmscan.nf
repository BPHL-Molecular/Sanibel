process bmgap2_bmscan {
    input:
        val mypath
        path pyoutputs
    output:
        val mypath
        path pyoutputs

    when:
        params.meningitis

    """
    #!/usr/bin/env python3

    import subprocess
    import os
    import sys

    items = "${mypath}".strip().split("/")
    sample_name = items[-1]

    pyoutputs_file = "${pyoutputs}"
    scheme = ""

    try:
        with open(pyoutputs_file, 'r') as f:
            content = f.read().strip()
            fields = content.split(',')
            if len(fields) > 16:
                scheme = fields[16]
    except Exception as e:
        print(f"BMGAP2-BMScan: Error reading pyoutputs - {str(e)}", file=sys.stderr)
        sys.exit(0)
    if scheme not in ['neisseria', 'hinfluenzae']:
        print(f"BMGAP2-BMScan: Skipping {sample_name} - not a meningitis species (scheme={scheme})")
        sys.exit(0)

    assembly_dir = f"${mypath}/{sample_name}_assembly"
    output_dir = f"${mypath}/bmgap2_bmscan"
    meningitis_scripts_dir = "${params.bmgap2_db}"
    bmscan_script = os.path.join(meningitis_scripts_dir, "SpeciesDB", "bin", "identify_species.py")

    if not os.path.isdir(assembly_dir):
        print(f"BMGAP2-BMScan: Error - Assembly directory not found: {assembly_dir}", file=sys.stderr)
        sys.exit(0)

    if not os.path.isfile(bmscan_script):
        print(f"BMGAP2-BMScan: Error - BMScan script not found: {bmscan_script}", file=sys.stderr)
        sys.exit(0)

    os.makedirs(output_dir, exist_ok=True)

    print(f"BMGAP2-BMScan: Running species identification for {sample_name}")

    cmd = [
        sys.executable,
        bmscan_script,
        "-d", assembly_dir,      # Input directory with FASTA files
        "-o", output_dir,         # Output directory
        "-t", "${task.cpus}",     # Number of threads
        "-j"                      # JSON output only
    ]

    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            check=False
        )

        if result.returncode != 0:
            print(f"BMGAP2-BMScan: Error running BMScan", file=sys.stderr)
            print(f"STDOUT: {result.stdout}", file=sys.stderr)
            print(f"STDERR: {result.stderr}", file=sys.stderr)
            sys.exit(0)

        print(f"BMGAP2-BMScan: Successfully completed species identification for {sample_name}")
        print(result.stdout)

    except Exception as e:
        print(f"BMGAP2-BMScan: Exception - {str(e)}", file=sys.stderr)
        sys.exit(0)
    """
}
