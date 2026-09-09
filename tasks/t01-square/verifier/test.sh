#!/usr/bin/env bash
# Task entry point.  Exits 0 only if every case passes every check.
set -euo pipefail
cd "$(dirname "$0")/.."
exec python3 verifier/checks.py "$@"
