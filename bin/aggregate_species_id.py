#!/usr/bin/env python3
"""

Usage: aggregate_species_id.py <assembly_stats_csv> <kraken_report> <blast_16s_tsv>

Output (stdout, TSV):
  genus  species  confidence  evidence  contamination_flag
  where evidence = mash:G/S|kraken:G/S|blast16s:G/S

Confidence values: high (all 3 agree on genus+species), medium (2 agree), low (otherwise)
"""

import os
import sys
from collections import Counter

sys.path.insert(0, os.path.dirname(os.path.realpath(__file__)))
from sanibel_taxonomy import (
    detect_contamination, extract_contam_candidates,
    iter_blast16s_rows, iter_kraken_species_rows,
    BLAST16S_MIN_LENGTH, BLAST16S_MIN_PIDENT,
)

# Species-level identity threshold (%); genus-level thresholds are shared via
# BLAST16S_MIN_LENGTH / BLAST16S_MIN_PIDENT in sanibel_taxonomy.
SP_PIDENT = 98.7


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
    for genus, species, _pct, _reads in iter_kraken_species_rows(kraken_report):
        if genus:
            return genus, species
    return None, None


def parse_blast16s(blast_tsv):
    rows = []
    for hit in iter_blast16s_rows(blast_tsv):
        if hit.length < BLAST16S_MIN_LENGTH or hit.pident < BLAST16S_MIN_PIDENT:
            continue
        if hit.genus:
            rows.append((hit.genus, hit.species, hit.pident, hit.length,
                         hit.stitle, hit.qseqid, hit.sstart, hit.send, hit.bitscore))

    if not rows:
        return None, None, {}, []

    rows.sort(key=lambda r: (-r[2], abs(r[3] - 1500), -r[8]))
    winner_genus, winner_species = _select_winner(rows)

    contam_candidates = extract_contam_candidates(blast_tsv)

    return winner_genus, winner_species, contam_candidates, rows


def _select_winner(rows, sp_pident=SP_PIDENT):
    if not rows:
        return None, None
    sorted_rows = sorted(rows, key=lambda r: (-r[2], abs(r[3] - 1500), -r[8]))
    best_g, best_s, best_pi, _ln, _stitle, _qid, _ss, _se, _bs = sorted_rows[0]
    winner_genus   = best_g
    winner_species = best_s if best_pi >= sp_pident else None
    return winner_genus, winner_species


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

    top_genus_lower = genus_counts.most_common(1)[0][0]

    winner_genus = next(
        g for g, _s in sources.values()
        if g is not None and g.lower() == top_genus_lower
    )

    agreeing = [
        s for g, s in sources.values()
        if g is not None and g.lower() == top_genus_lower and s is not None
    ]

    if not agreeing:
        winner_species = 'unknown'
    else:
        top_sp_lower   = Counter(s.lower() for s in agreeing).most_common(1)[0][0]
        winner_species = next(s for s in agreeing if s.lower() == top_sp_lower)

    # Confidence reflects genus + species agreement on the winner, not genus alone.
    if winner_species in (None, 'unknown'):
        confidence = 'low'
    else:
        agree = sum(
            1 for g, s in sources.values()
            if g is not None and s is not None
            and g.lower() == winner_genus.lower() and s.lower() == winner_species.lower()
        )
        confidence = 'high' if agree >= 3 else 'medium' if agree == 2 else 'low'

    return winner_genus, winner_species, confidence, evidence


def main():
    if len(sys.argv) != 4:
        sys.exit(f"Usage: {sys.argv[0]} <assembly_stats_csv> <kraken_report> <blast_16s_tsv>")

    assembly_stats = sys.argv[1]
    kraken_report  = sys.argv[2]
    blast_tsv      = sys.argv[3]

    mash     = parse_mash(assembly_stats)
    kraken   = parse_kraken(kraken_report)
    blast16s_genus, blast16s_species, contam_candidates, blast16s_rows = parse_blast16s(blast_tsv)
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
            blast16s_genus, blast16s_species = _select_winner(matching_rows)

    orig_genus, orig_species = blast16s
    if orig_genus is not None:
        old_token = f"blast16s:{orig_genus} {orig_species}"
        gs = f"{blast16s_genus} {blast16s_species}" if blast16s_species else f"{blast16s_genus} sp."
        evidence  = evidence.replace(old_token, f"blast16s:{gs}")
    else:
        evidence = (evidence + '|' if evidence else '') + 'blast16s:No data'

    disp_genus, disp_species = ('Inconclusive', 'unknown') if confidence == 'low' else (genus, species)

    print('genus\tspecies\tconfidence\tevidence\tcontamination_flag')
    print(f"{disp_genus}\t{disp_species}\t{confidence}\t{evidence}\t{contam_flag}")


if __name__ == '__main__':
    main()
