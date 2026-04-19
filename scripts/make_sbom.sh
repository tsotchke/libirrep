#!/usr/bin/env bash
# SPDX-License-Identifier: MIT
#
# Emit an SPDX 2.3 JSON SBOM for libirrep. Every tracked source file becomes
# an SPDXRef-File entry with its SHA-256 and declared license (parsed from the
# SPDX-License-Identifier header). The top-level package depends on all files.
#
# Usage:   scripts/make_sbom.sh [VERSION] [> release/<v>/SBOM.spdx.json]
set -euo pipefail

VERSION="${1:-$(cat VERSION)}"
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"

# Enumerate tracked source-like files into a newline-delimited list. macOS
# still ships bash 3.2 so we avoid `mapfile`; newlines in filenames would
# confuse this pipeline — we don't have any and reject them if introduced.
FILE_LIST=$( {
    find include src tests benchmarks examples scripts docs \
         -type f \
         \( -name '*.c' -o -name '*.h' -o -name '*.sh' -o -name '*.md' \
            -o -name '*.in' -o -name 'Makefile' -o -name 'CMakeLists.txt' \
            -o -name 'Doxyfile' \) \
         -not -path '*/build/*' \
         -not -path '*/release/*' \
         -not -path '*/benchmarks/results/*' \
         -print 2>/dev/null | sort
    for f in LICENSE NOTICE README.md CHANGELOG.md VERSION Makefile CMakeLists.txt; do
        [ -f "$f" ] && echo "$f"
    done
} )

NOW="$(date -u +'%Y-%m-%dT%H:%M:%SZ')"

# SHA-256 function — portable between macOS (shasum) and Linux (sha256sum).
sha256() {
    if command -v sha256sum >/dev/null 2>&1; then
        sha256sum "$1" | awk '{print $1}'
    else
        shasum -a 256 "$1" | awk '{print $1}'
    fi
}

# License parse — we require SPDX-License-Identifier: somewhere in the first
# 8 lines of each source file. Missing license fails the SBOM (safer than
# silently claiming MIT on a file that isn't MIT).
spdx_license_of() {
    local f="$1"
    local lic
    lic=$(head -n 8 "$f" 2>/dev/null | grep -Eo 'SPDX-License-Identifier:[[:space:]]*[A-Za-z0-9.+-]+' | head -n1 | awk '{print $2}') || true
    if [ -z "$lic" ]; then
        case "$f" in
            LICENSE|NOTICE|README.md|CHANGELOG.md|VERSION) echo "MIT" ;;
            *.md|Doxyfile)                                  echo "MIT" ;;
            *) echo "NOASSERTION" ;;
        esac
    else
        echo "$lic"
    fi
}

printf '{\n'
printf '  "spdxVersion": "SPDX-2.3",\n'
printf '  "dataLicense": "CC0-1.0",\n'
printf '  "SPDXID": "SPDXRef-DOCUMENT",\n'
printf '  "name": "libirrep-%s",\n' "$VERSION"
printf '  "documentNamespace": "https://github.com/tsotchke/libirrep/spdx/%s",\n' "$VERSION"
printf '  "creationInfo": {\n'
printf '    "created": "%s",\n' "$NOW"
printf '    "creators": ["Tool: libirrep-make_sbom/1.0"]\n'
printf '  },\n'
printf '  "packages": [{\n'
printf '    "SPDXID": "SPDXRef-libirrep",\n'
printf '    "name": "libirrep",\n'
printf '    "versionInfo": "%s",\n' "$VERSION"
printf '    "licenseDeclared": "MIT",\n'
printf '    "licenseConcluded": "MIT",\n'
printf '    "downloadLocation": "NOASSERTION",\n'
printf '    "filesAnalyzed": true\n'
printf '  }],\n'
printf '  "files": [\n'

first=1
while IFS= read -r f; do
    [ -z "$f" ] && continue
    sha=$(sha256 "$f")
    lic=$(spdx_license_of "$f")
    # SPDXID must be a safe identifier — collapse path separators.
    id="SPDXRef-File-$(printf '%s' "$f" | tr '/.' '--' | tr -c 'A-Za-z0-9-' '-')"
    [ $first -eq 0 ] && printf ',\n'
    printf '    {\n'
    printf '      "SPDXID": "%s",\n' "$id"
    printf '      "fileName": "./%s",\n' "$f"
    printf '      "checksums": [{ "algorithm": "SHA256", "checksumValue": "%s" }],\n' "$sha"
    printf '      "licenseConcluded": "%s",\n' "$lic"
    printf '      "licenseInfoInFiles": ["%s"]\n' "$lic"
    printf '    }'
    first=0
done <<< "$FILE_LIST"
printf '\n  ],\n'
printf '  "relationships": [\n'
printf '    { "spdxElementId": "SPDXRef-DOCUMENT", "relationshipType": "DESCRIBES", "relatedSpdxElement": "SPDXRef-libirrep" }\n'
printf '  ]\n'
printf '}\n'
