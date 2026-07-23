#!/usr/bin/env python3
"""
Build a ranked candidate species pool for Stage 2 multi-reference ANI.

Usage:
    build_candidate_pool.py <mash_distances_tab> <kraken_report> <blast16s_tsv>

Output (stdout, TSV):
    species  tools  accession  mash_distance  kraken_pct  blast16s_pident

Each row is one candidate species. Each tool's top-3 (mash by distance, kraken by
reads, 16S by pident) is seeded into the pool regardless of corroboration; the
remaining slots fill by number of supporting tools (desc), then mash_distance (asc),
then blast16s_pident (desc), capped at POOL_CAP.
"""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.realpath(__file__)))
from sanibel_taxonomy import (
    iter_blast16s_rows, iter_kraken_species_rows, find_accession, parse_mash_ref,
    BLAST16S_MIN_LENGTH, BLAST16S_MIN_PIDENT,
)

KRAKEN_ADAPT_FACTOR = 0.15
KRAKEN_MIN_PCT      = 5.0
KRAKEN_MIN_READS    = 10

SEED_TOP            = 3
POOL_CAP            = 15

NA = 'NA'

HEADER = ['species', 'tools', 'accession', 'mash_distance', 'kraken_pct', 'blast16s_pident']


# Mash parsing

def _parse_mash_ref(ref_name):
    try:
        genus, species, accession = parse_mash_ref(ref_name)
        if not genus:
            return None, None, None
        if accession is None:
            accession = find_accession(ref_name)
        return genus, species or NA, accession or NA
    except Exception:
        return None, None, None


def parse_mash_distances(filepath):
    best = {}

    try:
        with open(filepath) as fh:
            for line in fh:
                line = line.strip()
                if not line:
                    continue
                parts = line.split('\t')
                if len(parts) < 3:
                    continue
                ref_name = parts[0]
                try:
                    distance = float(parts[2])
                except ValueError:
                    continue
                genus, species, accession = _parse_mash_ref(ref_name)
                if genus is None:
                    continue
                key = f"{genus.lower()} {species.lower()}"
                if key not in best or distance < best[key]['mash_distance']:
                    best[key] = {
                        'genus':        genus,
                        'species':      species,
                        'species_key':  key,
                        'accession':    accession,
                        'mash_distance': distance,
                    }
    except OSError:
        pass

    sorted_all = sorted(best.values(), key=lambda d: d['mash_distance'])
    return sorted_all, best


# Kraken parsing

def parse_kraken_candidates(kraken_rows):
    species_rows = [(genus, species, pct) for genus, species, pct, _reads in kraken_rows]

    if not species_rows:
        return []

    species_rows.sort(key=lambda x: x[2], reverse=True)
    top_pct   = species_rows[0][2]
    threshold = max(top_pct * KRAKEN_ADAPT_FACTOR, KRAKEN_MIN_PCT)

    above   = [r for r in species_rows if r[2] >= threshold]

    included_genera = {r[0].lower() for r in above}
    seen_foreign    = set()
    foreign_top     = []
    for genus, species, pct in species_rows:
        gkey = genus.lower()
        if gkey not in included_genera and gkey not in seen_foreign:
            foreign_top.append((genus, species, pct))
            seen_foreign.add(gkey)

    return above + foreign_top


def kraken_seed_rows(kraken_rows, n):
    rows = sorted(kraken_rows, key=lambda x: x[3], reverse=True)
    return [(g, s, pct) for g, s, pct, _r in rows[:n]]


# 16S BLAST parsing

def parse_16s_candidates(filepath):
    qualifying = []
    for hit in iter_blast16s_rows(filepath):
        if hit.length < BLAST16S_MIN_LENGTH or hit.pident < BLAST16S_MIN_PIDENT:
            continue
        if hit.species is None:
            continue
        qualifying.append((hit.genus, hit.species, hit.pident))

    if not qualifying:
        return []

    best = {}
    for genus, species, pident in qualifying:
        key = f"{genus.lower()} {species.lower()}"
        if key not in best or pident > best[key][2]:
            best[key] = (genus, species, pident)

    sorted_cands = sorted(best.values(), key=lambda r: r[2], reverse=True)
    return sorted_cands


# Merge / rank

def merge_candidates(mash_cands, mash_all, kraken_cands, blast16s_cands, kraken_seeds, seed_keys):
    pool = {}  # species_key -> dict

    def _upsert(genus, species, tool, accession=NA,
                mash_distance=None, kraken_pct=None, blast16s_pident=None):
        key = f"{genus.lower()} {species.lower()}"
        if key not in pool:
            pool[key] = {
                'key':              key,
                'display_species':  f"{genus} {species}",
                'tools':            [],
                'accession':        accession,
                'mash_distance':    None,
                'kraken_pct':       None,
                'blast16s_pident':  None,
            }
        entry = pool[key]
        if tool not in entry['tools']:
            entry['tools'].append(tool)
        if accession != NA and entry['accession'] == NA:
            entry['accession'] = accession
        if mash_distance is not None:
            entry['mash_distance'] = mash_distance
        if kraken_pct is not None:
            entry['kraken_pct'] = kraken_pct
        if blast16s_pident is not None:
            if entry['blast16s_pident'] is None or blast16s_pident > entry['blast16s_pident']:
                entry['blast16s_pident'] = blast16s_pident

    for c in mash_cands:
        _upsert(c['genus'], c['species'], 'mash',
                accession=c['accession'],
                mash_distance=c['mash_distance'])

    for genus, species, pct in kraken_cands:
        key = f"{genus.lower()} {species.lower()}"
        acc = mash_all[key]['accession'] if key in mash_all else NA
        _upsert(genus, species, 'kraken2', accession=acc, kraken_pct=pct)

    for genus, species, pident in blast16s_cands:
        key = f"{genus.lower()} {species.lower()}"
        acc = mash_all[key]['accession'] if key in mash_all else NA
        _upsert(genus, species, '16S', accession=acc, blast16s_pident=pident)

    for genus, species, pct in kraken_seeds:
        key = f"{genus.lower()} {species.lower()}"
        if key not in pool:
            acc = mash_all[key]['accession'] if key in mash_all else NA
            _upsert(genus, species, 'kraken2', accession=acc, kraken_pct=pct)

    def _rank(e):
        return (
            -len(e['tools']),
            e['mash_distance'] if e['mash_distance'] is not None else 1.0,
            -(e['blast16s_pident'] if e['blast16s_pident'] is not None else 0.0),
        )

    seeds = sorted((e for e in pool.values() if e['key'] in seed_keys), key=_rank)
    fill  = sorted((e for e in pool.values() if e['key'] not in seed_keys), key=_rank)

    return (seeds + fill)[:POOL_CAP]


# Formatting helpers

def _fmt(value, fmt_spec):
    if value is None:
        return NA
    return format(value, fmt_spec)


def _row(entry):
    tools_str = ' + '.join(entry['tools'])
    return [
        entry['display_species'],
        tools_str,
        entry['accession'],
        _fmt(entry['mash_distance'],   '.7f'),
        _fmt(entry['kraken_pct'],      '.2f'),
        _fmt(entry['blast16s_pident'], '.3f'),
    ]


# Main

def main():
    if len(sys.argv) != 4:
        sys.exit(
            f"Usage: {sys.argv[0]} <mash_distances_tab> "
            "<kraken_report> <blast16s_tsv>"
        )

    mash_tab        = sys.argv[1]
    kraken_report   = sys.argv[2]
    blast16s_tsv    = sys.argv[3]

    mash_cands, mash_all = parse_mash_distances(mash_tab)
    kraken_rows = [
        (genus, species, pct, reads)
        for genus, species, pct, reads in iter_kraken_species_rows(kraken_report)
        if reads >= KRAKEN_MIN_READS and species is not None
    ]
    kraken_cands         = parse_kraken_candidates(kraken_rows)
    kraken_seeds         = kraken_seed_rows(kraken_rows, SEED_TOP)
    blast16s_cands       = parse_16s_candidates(blast16s_tsv)

    seed_keys = set()
    seed_keys.update(c['species_key'] for c in mash_cands[:SEED_TOP])
    seed_keys.update(f"{g.lower()} {s.lower()}" for g, s, _p in kraken_seeds)
    seed_keys.update(f"{g.lower()} {s.lower()}" for g, s, _pi in blast16s_cands[:SEED_TOP])

    ranked = merge_candidates(mash_cands, mash_all, kraken_cands, blast16s_cands,
                              kraken_seeds, seed_keys)

    print('\t'.join(HEADER))
    for entry in ranked:
        print('\t'.join(_row(entry)))


if __name__ == '__main__':
    main()
