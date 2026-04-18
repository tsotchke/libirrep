#!/usr/bin/env bash
# SPDX-License-Identifier: MIT
# Compare current build's ABI hash against the last tagged baseline.
# Full implementation lands in M13.
set -euo pipefail

CURRENT="$(bash scripts/generate_abi_hash.sh build/lib)"
BASELINE_FILE="release/BASELINE_ABI_HASH"

if [ ! -f "${BASELINE_FILE}" ]; then
    echo "# no baseline ABI hash yet; writing ${BASELINE_FILE}"
    echo "${CURRENT}" > "${BASELINE_FILE}"
    exit 0
fi

BASELINE="$(cat "${BASELINE_FILE}")"
if [ "${CURRENT}" = "${BASELINE}" ]; then
    echo "# ABI unchanged: ${CURRENT}"
    exit 0
fi

echo "# ABI drift detected"
echo "#   baseline: ${BASELINE}"
echo "#   current:  ${CURRENT}"
echo "# If intentional, bump MAJOR in VERSION and refresh ${BASELINE_FILE}."
exit 1
