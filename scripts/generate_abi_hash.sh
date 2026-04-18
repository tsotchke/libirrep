#!/usr/bin/env bash
# SPDX-License-Identifier: MIT
# Emit SHA-256 over sorted exported symbol names of the static library under
# the given directory. Public struct layouts will fold in in M13 once we emit
# a canonical layout-dump.
set -euo pipefail

DIR="${1:-build/lib}"
STATIC_LIB="$(ls "${DIR}"/liblibirrep.a 2>/dev/null | head -n1 || true)"

if [ -z "${STATIC_LIB}" ]; then
    echo "no-static-lib"
    exit 0
fi

# Portable symbol extraction (works on BSD nm and GNU nm).
if nm -gU "${STATIC_LIB}" >/dev/null 2>&1; then
    SYMS="$(nm -gU "${STATIC_LIB}" 2>/dev/null | awk '{print $NF}' | sort -u)"
else
    SYMS="$(nm -g --defined-only "${STATIC_LIB}" 2>/dev/null | awk '{print $NF}' | sort -u)"
fi

# Filter to irrep-prefixed symbols only.
EXPORTED="$(printf "%s\n" "${SYMS}" | grep -E '^_?irrep_' || true)"

if command -v shasum >/dev/null 2>&1; then
    printf "%s" "${EXPORTED}" | shasum -a 256 | awk '{print $1}'
else
    printf "%s" "${EXPORTED}" | sha256sum | awk '{print $1}'
fi
