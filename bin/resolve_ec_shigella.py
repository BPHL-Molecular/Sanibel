#!/usr/bin/env python3
"""resolve_ec_shigella.py - Pick Escherichia or Shigella from the ShigaTyper call."""

import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.realpath(__file__)))
from sanibel_taxonomy import genus_of

MIN_ANI = 95.0
MIN_AF  = 50.0


def read_prediction(path):
    if not os.path.isfile(path):
        return ''
    found_header = False
    with open(path) as f:
        for line in f:
            if found_header:
                fields = line.rstrip('\n').split('\t')
                return fields[1].strip() if len(fields) >= 2 else ''
            if line.startswith('sample\t'):
                found_header = True
    return ''


def resolve_target(prediction):
    p = prediction.lower()
    if 'not shigella' in p:
        return 'Escherichia', None
    if 'shigella' in p:
        tokens = prediction.split()
        named = f'{tokens[0]}_{tokens[1]}' if len(tokens) > 1 and tokens[1].islower() else None
        return 'Shigella', named
    if 'eiec' in p:
        return 'Escherichia', None
    return None, None


def read_skani(path):
    rows = []
    with open(path) as f:
        f.readline()
        for line in f:
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 5:
                continue
            try:
                ani = float(parts[2])
                af  = float(parts[4])
            except ValueError:
                continue
            if ani < MIN_ANI or af < MIN_AF:
                continue
            name = os.path.basename(parts[0])
            if name.endswith('.fna'):
                name = name[:-4]
            species = name.split('__', 1)[0]
            if species:
                rows.append((ani, species))
    rows.sort(key=lambda r: r[0], reverse=True)
    return rows


def main():
    ap = argparse.ArgumentParser(description='Resolve Escherichia vs Shigella from ShigaTyper.')
    ap.add_argument('--skani',      required=True)
    ap.add_argument('--shigatyper', required=True)
    ap.add_argument('--prefix',     required=True)
    args = ap.parse_args()

    rows       = read_skani(args.skani)
    prediction = read_prediction(args.shigatyper)
    genus, named_species = resolve_target(prediction)

    resolved = ''
    if not rows:
        print(f'resolve_ec_shigella: no skani row clears {MIN_ANI}/{MIN_AF}', file=sys.stderr)
    elif genus is None:
        print(f"resolve_ec_shigella: ShigaTyper prediction '{prediction}' is not usable, "
              f'keeping {rows[0][1]}', file=sys.stderr)
    else:
        pick = next((sp for _ani, sp in rows if sp == named_species), None) \
            or next((sp for _ani, sp in rows if genus_of(sp) == genus), None)
        if pick is None:
            print(f'resolve_ec_shigella: ShigaTyper says {genus} but no {genus} reference '
                  f'clears {MIN_ANI}/{MIN_AF}, keeping {rows[0][1]}', file=sys.stderr)
        elif pick != rows[0][1]:
            resolved = pick
            print(f"resolve_ec_shigella: {rows[0][1]} -> {pick} on ShigaTyper '{prediction}'",
                  file=sys.stderr)

    with open(f'{args.prefix}_species_resolved.txt', 'w') as fh:
        fh.write(f'{resolved}\n' if resolved else '')


if __name__ == '__main__':
    main()
