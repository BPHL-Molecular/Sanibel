#!/usr/bin/env python3
"""
Parse mash distances and QUAST report into a single CSV line.

Usage: parse_assembly.py <distances_file> <quast_report>

Output (stdout):
  GENUS,SPECIES,DIST,ACCESSION,CONTIGS,LARGEST,N50,L50,TOTAL,GC
"""

import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.realpath(__file__)))
from sanibel_taxonomy import find_accession


def parse_mash(distances_file):
    with open(distances_file) as fh:
        line = fh.readline().strip()

    fields = line.split('\t')
    ref_id = fields[0]
    dist   = fields[2] if len(fields) > 2 else 'NA'

    if '-.-' in ref_id:
        segs     = ref_id.split('-.-')
        org_part = segs[-1]
        org_seg  = re.sub(r'\.fna.*', '', org_part)
        parts    = [p for p in org_seg.split('_') if p]
        genus    = parts[0] if parts else 'Unknown'
        species  = parts[1] if len(parts) > 1 else 'unknown'
        accession = 'Unknown'
        for seg in segs[:-1]:
            acc = find_accession(seg)
            if acc:
                accession = acc
                break
    else:
        org_seg   = re.sub(r'\.fna.*', '', ref_id)
        parts     = [p for p in org_seg.split('_') if p]
        genus     = parts[0] if parts else 'Unknown'
        species   = parts[1] if len(parts) > 1 else 'unknown'
        accession = find_accession(ref_id) or 'Unknown'

    return genus, species, dist, accession


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

    genus, species, dist, accession        = parse_mash(sys.argv[1])
    contigs, largest, n50, l50, total, gc  = parse_quast(sys.argv[2])

    print(f"{genus},{species},{dist},{accession},"
          f"{contigs},{largest},{n50},{l50},{total},{gc}")
