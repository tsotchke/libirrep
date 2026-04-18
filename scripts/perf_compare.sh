#!/usr/bin/env bash
# SPDX-License-Identifier: MIT
# Compare two benchmark JSON files and flag regressions > 5%. M12 replaces the
# stub with structured diff; today it just echoes both files.
set -euo pipefail

BASELINE="${1:?usage: scripts/perf_compare.sh baseline.json candidate.json}"
CANDIDATE="${2:?}"

echo "# baseline: ${BASELINE}"
cat "${BASELINE}"
echo "# candidate: ${CANDIDATE}"
cat "${CANDIDATE}"
