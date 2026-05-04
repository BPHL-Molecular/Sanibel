#!/usr/bin/env python3
"""
Parse mash distances and QUAST report into a single CSV line.

Usage: parse_assembly.py <distances_file> <quast_report>

Output (stdout):
  GENUS,SPECIES,DIST,ACCESSION,ASM_NAME,CONTIGS,LARGEST,N50,L50,TOTAL,GC
"""

import re
import sys


def parse_mash(distances_file):
    with open(distances_file) as fh:
        line = fh.readline().strip()

    fields = line.split('\t')
    ref_id = fields[0]
    dist   = fields[2] if len(fields) > 2 else 'NA'

    if '-.-' in ref_id:
        segs = ref_id.split('-.-')
        if len(segs) >= 3:
            # Format: pre -.- accession -.- organism.fna (e.g. some Listeria refs)
            acc_str  = segs[-2]
            org_part = segs[-1]
            m = re.search(r'((?:GC[FA]|N[CZ])_[A-Za-z0-9]+(?:\.[0-9]+)?)', acc_str)
            accession = m.group(1) if m else (acc_str.rsplit('-', 1)[-1] or 'Unknown')
        else:
            # Standard format: acc_part -.- organism.fna
            acc_part, org_part = segs
            m = re.search(r'((?:GC[FA]|N[CZ])_[A-Za-z0-9]+(?:\.[0-9]+)?)', acc_part)
            if m:
                accession = m.group(1)
            else:
                accession = acc_part.rsplit('-', 1)[-1] or 'Unknown'
        org_seg  = re.sub(r'\.fna.*', '', org_part)
        parts    = [p for p in org_seg.split('_') if p]
        genus    = parts[0] if parts else 'Unknown'
        species  = parts[1] if len(parts) > 1 else 'unknown'
        asm_name = '_'.join(parts[2:]) if len(parts) > 2 else '.'
    else:
        accession = 'Unknown'
        genus     = 'Unknown'
        species   = 'unknown'
        asm_name  = '.'

    return genus, species, dist, accession, asm_name


def parse_quast(quast_file):
    want = {
        '# contigs':      'NA',
        'Largest contig': 'NA',
        'N50':            'NA',
        'L50':            'NA',
        'Total length':   'NA',
        'GC (%)':         'NA',
    }
    with open(quast_file) as fh:
        for row in fh:
            parts = row.rstrip('\n').split('\t')
            if len(parts) >= 2 and parts[0] in want:
                want[parts[0]] = parts[1]
    return (
        want['# contigs'],
        want['Largest contig'],
        want['N50'],
        want['L50'],
        want['Total length'],
        want['GC (%)'],
    )


if __name__ == '__main__':
    if len(sys.argv) != 3:
        sys.exit(f"Usage: {sys.argv[0]} <distances_file> <quast_report>")

    genus, species, dist, accession, asm_name = parse_mash(sys.argv[1])
    contigs, largest, n50, l50, total, gc     = parse_quast(sys.argv[2])

    print(f"{genus},{species},{dist},{accession},{asm_name},"
          f"{contigs},{largest},{n50},{l50},{total},{gc}")
