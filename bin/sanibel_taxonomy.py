#!/usr/bin/env python3
"""
sanibel_taxonomy.py — shared taxonomy / contamination primitives.

Single source of truth for the 16S-synonymous genus table and the contig-overlap
contamination logic, imported by aggregate_species_id.py (the candidate-pool vote)
and summary_report.py (the skani-anchored report). Keeping these in one place stops
the synonym table and the contamination thresholds from drifting between the two.
"""

# Minimum identity/length for a 16S hit to count as a contamination candidate.
CONTAM_PIDENT = 99.0
CONTAM_LENGTH = 1400


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


def extract_contam_candidates(blast16s_tsv_path):
    rows = []
    try:
        with open(blast16s_tsv_path) as fh:
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
                if pident < CONTAM_PIDENT or length < CONTAM_LENGTH:
                    continue
                stitle = parts[12].strip()
                genus  = stitle.split()[0] if stitle else None
                if not genus:
                    continue
                rows.append((pident, length, bitscore, genus, parts[0], sstart, send))
    except OSError:
        return {}

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
