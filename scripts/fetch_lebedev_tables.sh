#!/usr/bin/env bash
# SPDX-License-Identifier: MIT
# Download Lebedev-Laikov 1999 public-domain quadrature tables for orders
# 9..41. The library does not bundle this data (it's ~50 KB of pure
# numerical tables with no source-code authorship question — but it
# lives upstream at John Burkardt's site in a clean ASCII format, and
# shipping a cached copy would just duplicate the source of truth).
#
# Usage:
#
#     scripts/fetch_lebedev_tables.sh [dest_dir]
#
# Default destination is `data/lebedev/`. Each file is
# `lebedev_<order>.txt`, one point per line, whitespace-separated
# `(theta_deg, phi_deg, weight)` triples. `examples/register_lebedev.c`
# parses this format and feeds it through `irrep_lebedev_register_rule`.
#
# Primary source: Burkardt, J. (2010), sphere_lebedev_rule dataset,
#   https://people.sc.fsu.edu/~jburkardt/datasets/sphere_lebedev_rule/
# Upstream: Lebedev, V.I. & Laikov, D.N. (1999), Doklady Mathematics 59,
#   477–481.
set -euo pipefail

DEST="${1:-data/lebedev}"
mkdir -p "$DEST"

BASE="https://people.sc.fsu.edu/~jburkardt/datasets/sphere_lebedev_rule"

# Lebedev-Laikov 1999 orders (odd, ≥ 9, up to 41). 31 is skipped in the
# paper (no rule published at that order) but Burkardt ships it as an
# additional table; include it since it's a valid positive-weight rule.
ORDERS=(9 11 13 15 17 19 21 23 25 27 29 31 35 41)

for o in "${ORDERS[@]}"; do
    # Burkardt's file naming is `lebedev_NNN.txt` with 3-digit zero-padded
    # order. e.g. order 9 → `lebedev_009.txt`, order 41 → `lebedev_041.txt`.
    filename=$(printf "lebedev_%03d.txt" "$o")
    dest="${DEST}/lebedev_${o}.txt"
    if [ -s "$dest" ]; then
        echo "skip ${filename} (already present at ${dest})"
        continue
    fi
    echo "fetch ${filename}"
    if command -v curl >/dev/null 2>&1; then
        curl -fsSL "${BASE}/${filename}" -o "$dest"
    elif command -v wget >/dev/null 2>&1; then
        wget -q "${BASE}/${filename}" -O "$dest"
    else
        echo "error: need curl or wget to fetch Lebedev tables" >&2
        exit 1
    fi
done

echo ""
echo "# fetched $(ls -1 "$DEST"/lebedev_*.txt 2>/dev/null | wc -l | tr -d ' ') Lebedev tables into ${DEST}"
echo "# load via examples/register_lebedev (or your own parser)"
