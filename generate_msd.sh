#!/usr/bin/env bash
# MSD source coordinates are millimetres; batch outputs are metres.
set -euo pipefail

repo_root="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
output="/shared/users/cjm6121/hybrid/MSD-RIT/Task03_Liver_NeoHookean/batch_$(date -u +%Y%m%dT%H%M%SZ)"

exec pixi run --manifest-path "$repo_root/pixi.toml" python \
  "$repo_root/batch.py" \
  /shared/users/cjm6121/hybrid/MSD-RIT/Task03_Liver \
  "$output" --preset all131_20 "$@"
