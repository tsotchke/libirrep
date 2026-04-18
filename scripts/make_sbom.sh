#!/usr/bin/env bash
# SPDX-License-Identifier: MIT
# Emit a minimal SPDX 2.3 JSON SBOM for the current release tree.
# M13 replaces this stub with full source-file enumeration and dependency graph.
set -euo pipefail

VERSION="${1:-$(cat VERSION)}"
cat <<EOF
{
  "spdxVersion": "SPDX-2.3",
  "dataLicense": "CC0-1.0",
  "SPDXID": "SPDXRef-DOCUMENT",
  "name": "libirrep-${VERSION}",
  "documentNamespace": "https://example.com/libirrep/${VERSION}",
  "creationInfo": { "created": "$(date -u +'%Y-%m-%dT%H:%M:%SZ')", "creators": ["Tool: libirrep-make_sbom"] },
  "packages": [{
    "SPDXID": "SPDXRef-libirrep",
    "name": "libirrep",
    "versionInfo": "${VERSION}",
    "licenseDeclared": "MIT",
    "licenseConcluded": "MIT",
    "downloadLocation": "NOASSERTION"
  }]
}
EOF
