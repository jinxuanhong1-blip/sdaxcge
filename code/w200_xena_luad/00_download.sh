#!/usr/bin/env bash
# Thin wrapper. Canonical downloader: scripts/xena_luad_download.sh
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
exec bash "$ROOT/scripts/xena_luad_download.sh" "$@"
