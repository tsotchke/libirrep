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

# Filter to public irrep-prefixed symbols only. Normalisations for
# portability across Mach-O / ELF / PE:
#   - strip the Mach-O leading underscore
#   - exclude internal SIMD dispatch kernels (suffix _neon / _avx2 /
#     _avx512 / _scalar): these are platform-specific implementation
#     details reached only through the internal function-pointer table,
#     so they must not perturb the public-API hash.
EXPORTED="$(printf "%s\n" "${SYMS}" \
    | grep -E '^_?irrep_' \
    | sed 's/^_//' \
    | grep -vE '_(neon|avx2|avx512|scalar)(_[a-z_0-9]*)?$' \
    | sort -u || true)"

if command -v shasum >/dev/null 2>&1; then
    printf "%s" "${EXPORTED}" | shasum -a 256 | awk '{print $1}'
else
    printf "%s" "${EXPORTED}" | sha256sum | awk '{print $1}'
fi
