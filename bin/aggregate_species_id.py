#!/usr/bin/env python3
"""
2-of-3 majority vote for species identification.

Usage: aggregate_species_id.py <assembly_stats_csv> <kraken_report> <blast_16s_tsv>

Output (stdout, TSV):
  genus  species  confidence  evidence  contamination_flag
  where evidence = mash:G/S|kraken:G/S|blast16s:G/S/Strain name

Confidence values: high (all 3 agree), medium (2 of 3 agree), low (no majority)
"""

import sys
from collections import Counter

# Thresholds for 16S BLAST filtering and contamination detection
MIN_LENGTH    = 400    # minimum alignment length (bp)
MIN_PIDENT    = 97.0   # genus-level identity threshold (%)
SP_PIDENT     = 98.7   # species-level identity threshold (%)
CONTAM_PIDENT = 99.0   # minimum pident for contamination candidates
CONTAM_LENGTH = 1400   # minimum length for contamination candidates


def parse_mash(assembly_stats):
    with open(assembly_stats) as fh:
        line = fh.readline().strip()
    if not line:
        return None, None
    fields = line.split(',')
    genus   = fields[0].strip() if len(fields) > 0 else None
    species = fields[1].strip() if len(fields) > 1 else None
    if not genus or genus == 'Unknown':
        return None, None
    return genus, species


def parse_kraken(kraken_report):
    with open(kraken_report) as fh:
        for line in fh:
            parts = line.strip().split('\t')
            if len(parts) < 6:
                continue
            rank = parts[3].strip()
            if rank == 'S':
                name = parts[5].strip()
                words = name.split()
                genus   = words[0] if len(words) > 0 else None
                species = words[1] if len(words) > 1 else None
                if genus:
                    return genus, species
    return None, None


def parse_blast16s(blast_tsv):
    rows = []
    with open(blast_tsv) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            parts = line.split('\t')
            if len(parts) < 13:
                continue
            if parts[0] == 'qseqid':
                continue
            qseqid = parts[0]
            try:
                pident   = float(parts[2])
                length   = int(parts[3])
                sstart   = int(parts[8])
                send     = int(parts[9])
                bitscore = float(parts[11])
            except ValueError:
                continue
            if length < MIN_LENGTH or pident < MIN_PIDENT:
                continue
            words   = parts[12].strip().split()
            genus   = words[0] if len(words) > 0 else None
            species = words[1] if len(words) > 1 else None
            stitle  = parts[12].strip()
            if genus:
                rows.append((genus, species, pident, length, stitle, qseqid, sstart, send, bitscore))

    if not rows:
        return None, None, None, {}, []

    rows.sort(key=lambda r: (-r[2], abs(r[3] - 1500), -r[8]))
    winner_genus, winner_species, winner_strain = _select_winner(rows)

    contam_candidates = {}
    for g, s, pi, ln, title, qseqid, ss, se, _bs in rows:
        if pi >= CONTAM_PIDENT and ln >= CONTAM_LENGTH:
            if qseqid not in contam_candidates:
                contam_candidates[qseqid] = (g, ss, se)

    return winner_genus, winner_species, winner_strain, contam_candidates, rows


# Genera indistinguishable by 16S; never flag each other as contamination
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


def _select_winner(rows, sp_pident=SP_PIDENT):
    if not rows:
        return None, None, None
    sorted_rows = sorted(rows, key=lambda r: (-r[2], abs(r[3] - 1500), -r[8]))
    best_g, best_s, best_pi, _ln, best_stitle, _qid, _ss, _se, _bs = sorted_rows[0]
    winner_genus   = best_g
    winner_species = best_s if best_pi >= sp_pident else None
    winner_strain  = best_stitle.split(' 16S ')[0].strip() if ' 16S ' in best_stitle else best_stitle

    top_pi_row   = max(sorted_rows, key=lambda r: r[2])
    top_pi_genus = top_pi_row[0]
    if (top_pi_genus.lower() != winner_genus.lower() and
            frozenset({winner_genus.lower(), top_pi_genus.lower()}) in SYNONYMOUS_PAIRS):
        winner_strain = f"{winner_genus}/{top_pi_genus} synonymous pair"

    return winner_genus, winner_species, winner_strain


def _ranges_overlap(s1, e1, s2, e2):
    lo1, hi1 = min(s1, e1), max(s1, e1)
    lo2, hi2 = min(s2, e2), max(s2, e2)
    return lo1 <= hi2 and lo2 <= hi1


def detect_contamination(contam_candidates, final_genus):
    if len(contam_candidates) < 2:
        return 'None'
    contigs = list(contam_candidates.items())
    for i in range(len(contigs)):
        qid_a, (genus_a, s_a, e_a) = contigs[i]
        for j in range(i + 1, len(contigs)):
            qid_b, (genus_b, s_b, e_b) = contigs[j]
            if _ranges_overlap(s_a, e_a, s_b, e_b):
                if genus_a.lower() != final_genus.lower():
                    return f"Possible contamination: {genus_a} detected on contig {qid_a}"
                if genus_b.lower() != final_genus.lower():
                    return f"Possible contamination: {genus_b} detected on contig {qid_b}"
    return 'None'


def vote(mash, kraken, blast16s):
    sources = {
        'mash':     mash,
        'kraken':   kraken,
        'blast16s': blast16s,
    }

    evidence = '|'.join(
        f"{name}:{g} {s}"
        for name, (g, s) in sources.items()
        if g is not None
    )

    genus_counts = Counter(
        g.lower() for g, _s in sources.values() if g is not None
    )
    if not genus_counts:
        return 'Unknown', 'unknown', 'low', evidence

    top_genus_lower, top_genus_count = genus_counts.most_common(1)[0]

    winner_genus = next(
        g for g, _s in sources.values()
        if g is not None and g.lower() == top_genus_lower
    )

    agreeing = {
        name: s
        for name, (g, s) in sources.items()
        if g is not None and g.lower() == top_genus_lower and s is not None
    }

    if not agreeing:
        winner_species = 'unknown'
        species_count  = 0
    else:
        species_counts = Counter(s.lower() for s in agreeing.values())
        top_sp_lower, species_count = species_counts.most_common(1)[0]
        winner_species = next(
            s for s in agreeing.values()
            if s is not None and s.lower() == top_sp_lower
        )

    if top_genus_count >= 2 and species_count >= 2:
        confidence = 'high' if top_genus_count == 3 and species_count == 3 else 'medium'
    elif top_genus_count >= 2:
        confidence = 'medium'
    else:
        confidence = 'low'

    return winner_genus, winner_species, confidence, evidence


def main():
    if len(sys.argv) != 4:
        sys.exit(f"Usage: {sys.argv[0]} <assembly_stats_csv> <kraken_report> <blast_16s_tsv>")

    assembly_stats = sys.argv[1]
    kraken_report  = sys.argv[2]
    blast_tsv      = sys.argv[3]

    mash     = parse_mash(assembly_stats)
    kraken   = parse_kraken(kraken_report)
    blast16s_genus, blast16s_species, blast16s_strain, contam_candidates, blast16s_rows = parse_blast16s(blast_tsv)
    blast16s = (blast16s_genus, blast16s_species)

    genus, species, confidence, evidence = vote(mash, kraken, blast16s)

    contam_flag = detect_contamination(contam_candidates, genus)

    if contam_flag != 'None' and blast16s_genus is not None and blast16s_genus.lower() != genus.lower():
        contig_best_genus = {}
        for r in blast16s_rows:
            if r[5] not in contig_best_genus:
                contig_best_genus[r[5]] = r[0]
        matching_rows = [
            r for r in blast16s_rows
            if contig_best_genus.get(r[5], '').lower() == genus.lower()
        ]
        if matching_rows:
            blast16s_genus, blast16s_species, blast16s_strain = _select_winner(matching_rows)

    orig_genus, orig_species = blast16s
    if orig_genus is not None:
        old_token  = f"blast16s:{orig_genus} {orig_species}"
        strain_val = blast16s_strain if blast16s_strain is not None else f"{blast16s_genus} {blast16s_species}"
        evidence   = evidence.replace(old_token, f"blast16s:{strain_val}")
    else:
        evidence = (evidence + '|' if evidence else '') + 'blast16s:No data'

    print('genus\tspecies\tconfidence\tevidence\tcontamination_flag')
    print(f"{genus}\t{species}\t{confidence}\t{evidence}\t{contam_flag}")


if __name__ == '__main__':
    main()
