#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"

cmd=$'\
echo "===== $(date) ====="; \
cmake --build build -j; \
echo "===== RUN ====="; \
./build/main \
'

watchexec --restart \
  --watch src \
  --watch CMakeLists.txt \
  --ignore 'build/**' \
  --ignore '**/*.swp' \
  --ignore '**/*.swo' \
  --ignore '**/*~' \
  --ignore '**/.DS_Store' \
  -- /bin/bash -lc "$cmd"
