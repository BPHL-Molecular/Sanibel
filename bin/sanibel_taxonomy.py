#!/usr/bin/env python3
"""
sanibel_taxonomy.py — shared taxonomy / contamination primitives.

"""

import re
from collections import namedtuple

# Minimum identity/length for a 16S hit to count as a contamination candidate.
CONTAM_PIDENT = 99.0
CONTAM_LENGTH = 1400

# Minimum identity/length for a 16S hit to count as a species-ID candidate
BLAST16S_MIN_LENGTH = 400
BLAST16S_MIN_PIDENT = 97.0

# RefSeq/GenBank accession pattern (GCF_*, GCA_*, NC_*, NZ_*)
ACCESSION_RE = re.compile(
    r'GC[FA]_[A-Za-z0-9]+(?:\.[0-9]+)?|N[CZ]_[A-Za-z0-9]+(?:\.[0-9]+)?'
)


def find_accession(text):
    """First RefSeq/GenBank accession in `text`, or None."""
    if not text:
        return None
    m = ACCESSION_RE.search(text)
    return m.group(0) if m else None


# Genera indistinguishable by 16S; never flag each other as contamination.
SYNONYMOUS_PAIRS = {
    frozenset({'escherichia',   'shigella'}),
    frozenset({'klebsiella',    'enterobacter'}),
    frozenset({'klebsiella',    'raoultella'}),
    frozenset({'salmonella',    'citrobacter'}),
    frozenset({'yersinia',      'serratia'}),
    frozenset({'hafnia',        'escherichia'}),
    frozenset({'haemophilus',   'aggregatibacter'}),
    frozenset({'haemophilus',   'pasteurella'}),
    frozenset({'neisseria',     'kingella'}),
    frozenset({'neisseria',     'eikenella'}),
    frozenset({'streptococcus', 'lactococcus'}),
    frozenset({'streptococcus', 'enterococcus'}),
    frozenset({'staphylococcus', 'macrococcus'}),
    frozenset({'staphylococcus', 'mammaliicoccus'}),
    frozenset({'mycobacterium', 'mycobacteroides'}),
    frozenset({'mycobacterium', 'mycolicibacterium'}),
    frozenset({'campylobacter', 'arcobacter'}),
    frozenset({'campylobacter', 'aliarcobacter'}),
    frozenset({'campylobacter', 'helicobacter'}),
    frozenset({'listeria',      'brochothrix'}),
    frozenset({'bacillus',      'paenibacillus'}),
}


def genus_of(label):
    """First token of a 'Genus species' or 'Genus_species' string."""
    if not label:
        return None
    tokens = label.replace('_', ' ').split()
    return tokens[0] if tokens else None


def are_synonymous(genus_a, genus_b):
    if not genus_a or not genus_b:
        return False
    return frozenset({genus_a.lower(), genus_b.lower()}) in SYNONYMOUS_PAIRS


def ranges_overlap(s1, e1, s2, e2):
    lo1, hi1 = min(s1, e1), max(s1, e1)
    lo2, hi2 = min(s2, e2), max(s2, e2)
    return lo1 <= hi2 and lo2 <= hi1


# Parsed 16S BLAST hit
BlastHit = namedtuple(
    'BlastHit', 'qseqid pident length sstart send bitscore stitle genus species'
)


def iter_blast16s_rows(path):
    """Yield BlastHit for each valid 16S BLAST row (header/blank/malformed skipped)."""
    try:
        fh = open(path)
    except OSError:
        return
    with fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            parts = line.split('\t')
            if len(parts) < 13 or parts[0] == 'qseqid':
                continue
            try:
                pident   = float(parts[2])
                length   = int(parts[3])
                sstart   = int(parts[8])
                send     = int(parts[9])
                bitscore = float(parts[11])
            except (ValueError, IndexError):
                continue
            stitle  = parts[12].strip()
            words   = stitle.split()
            genus   = words[0] if words else None
            species = words[1] if len(words) > 1 else None
            yield BlastHit(parts[0], pident, length, sstart, send, bitscore,
                           stitle, genus, species)


def iter_kraken_species_rows(path):
    try:
        fh = open(path)
    except OSError:
        return
    with fh:
        for line in fh:
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 6 or parts[3].strip() != 'S':
                continue
            try:
                pct   = float(parts[0].strip())
                reads = int(parts[1].strip())
            except ValueError:
                continue
            words   = parts[5].strip().split()
            genus   = words[0] if words else None
            species = words[1] if len(words) > 1 else None
            yield genus, species, pct, reads


def extract_contam_candidates(blast16s_tsv_path):
    rows = []
    for hit in iter_blast16s_rows(blast16s_tsv_path):
        if hit.pident < CONTAM_PIDENT or hit.length < CONTAM_LENGTH:
            continue
        if not hit.genus:
            continue
        rows.append((hit.pident, hit.length, hit.bitscore, hit.genus,
                     hit.qseqid, hit.sstart, hit.send))

    rows.sort(key=lambda r: (-r[0], abs(r[1] - 1500), -r[2]))
    candidates = {}
    for _pi, _ln, _bs, genus, qseqid, sstart, send in rows:
        if qseqid not in candidates:
            candidates[qseqid] = (genus, sstart, send)
    return candidates


def detect_contamination(contam_candidates, self_genus, suppress_synonymous=True):
    if len(contam_candidates) < 2 or not self_genus:
        return 'None'

    def is_foreign(genus):
        if genus.lower() == self_genus.lower():
            return False
        if suppress_synonymous and are_synonymous(genus, self_genus):
            return False
        return True

    contigs = list(contam_candidates.items())
    for i in range(len(contigs)):
        qid_a, (genus_a, s_a, e_a) = contigs[i]
        for j in range(i + 1, len(contigs)):
            qid_b, (genus_b, s_b, e_b) = contigs[j]
            if ranges_overlap(s_a, e_a, s_b, e_b):
                if is_foreign(genus_a):
                    return f"Possible contamination: {genus_a} detected on contig {qid_a}"
                if is_foreign(genus_b):
                    return f"Possible contamination: {genus_b} detected on contig {qid_b}"
    return 'None'
