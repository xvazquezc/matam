#!/bin/bash
# Run upstream MATAM and this fork on the bundled 16-species example.
#
# Usage: run_compare.sh OUT_DIR
#
# Environment:
#   UPSTREAM_PREFIX  conda prefix with upstream MATAM, e.g. created with
#                    mamba create -p "$UPSTREAM_PREFIX" --channel-priority flexible \
#                        -c conda-forge -c bioconda matam=1.6.1
#                    (bioconda 1.6.1 ships the same scripts as upstream e7cd9b2)
#   FORK_PREFIX      conda prefix with this fork's dependencies (default: the
#                    active environment, i.e. the directory of `python`)
#   FORK_REPO        checkout of the fork to run (default: this repository);
#                    e.g. a `git worktree` of 264cf57 for the old pair rule
#   FORK_LABEL       run-name prefix for the fork (default "fork")
#   SKIP_UPSTREAM    set to 1 to run only the fork
#   SEEDS            seeds to run (default "1 2 3")
#   CPU              threads (default 8)
#
# Upstream gets its documented input (one interleaved file with /1, /2 names);
# the fork gets the same reads as separate R1/R2 files. Each version builds its
# own database from the same reference FASTA. Runs: OUT_DIR/up_<seed>,
# OUT_DIR/<FORK_LABEL>_<seed>, with logs and /usr/bin/time output next to them.
set -euo pipefail

OUT=${1:?usage: run_compare.sh OUT_DIR}
SKIP_UPSTREAM=${SKIP_UPSTREAM:-0}
if [ "$SKIP_UPSTREAM" != 1 ]; then
    UPSTREAM_PREFIX=${UPSTREAM_PREFIX:?set UPSTREAM_PREFIX to a conda prefix with matam=1.6.1}
fi
UPSTREAM_PREFIX=${UPSTREAM_PREFIX:-}
FORK_PREFIX=${FORK_PREFIX:-$(dirname "$(dirname "$(command -v python)")")}
SEEDS=${SEEDS:-"1 2 3"}
CPU=${CPU:-8}

REPO=$(cd "$(dirname "$0")/../.." && pwd)
FORK_REPO=$(cd "${FORK_REPO:-$REPO}" && pwd)
FORK_LABEL=${FORK_LABEL:-fork}
EX=$REPO/examples/16sp_simulated_dataset
UP_SCRIPTS=$UPSTREAM_PREFIX/opt/matam-1.6.1/scripts
UP_PATH=$UPSTREAM_PREFIX/opt/matam-1.6.1/ovgraphbuild/bin:$UPSTREAM_PREFIX/opt/matam-1.6.1/componentsearch:$UPSTREAM_PREFIX/bin:/usr/bin:/bin
FORK_PATH=$FORK_REPO/scripts:$FORK_PREFIX/bin:/usr/bin:/bin

mkdir -p "$OUT"
cd "$OUT"

# Shared inputs: SILVA-style reference headers (">ID taxonomy") and R1/R2 files.
"$FORK_PREFIX/bin/python" - "$REPO" <<'EOF'
import sys
sys.path.insert(0, sys.argv[1] + '/tests')
from pe_fixtures import write_pe_reads, write_reference_fasta
write_reference_fasta('refs.fasta')
write_pe_reads('reads', 'legacy')
EOF

if [ "$SKIP_UPSTREAM" != 1 ]; then
    rm -rf db_up
    PATH=$UP_PATH "$UPSTREAM_PREFIX/bin/python3" "$UP_SCRIPTS/matam_db_preprocessing.py" \
        -i refs.fasta -d db_up --cpu "$CPU" --max_memory 4000 --out_db_name refs > db_up.log 2>&1
fi
rm -rf "db_$FORK_LABEL"
PATH=$FORK_PATH "$FORK_PREFIX/bin/python" "$FORK_REPO/scripts/matam_db_preprocessing.py" \
    -i refs.fasta -d "db_$FORK_LABEL" --cpu "$CPU" --max_memory 4000 --out_db_name refs \
    > "db_$FORK_LABEL.log" 2>&1

for seed in $SEEDS; do
    if [ "$SKIP_UPSTREAM" != 1 ]; then
        rm -rf "up_$seed"
        /usr/bin/time -f "%e s wall, %M KB" -o "up_$seed.time" env PATH="$UP_PATH" \
            "$UPSTREAM_PREFIX/bin/python3" "$UP_SCRIPTS/matam_assembly.py" \
            -i "$EX/16sp.art_HS25_pe_100bp_50x.fq" -d db_up/refs -o "up_$seed" \
            --cpu "$CPU" --max_memory 4000 -v --seed "$seed" > "up_$seed.log" 2>&1
        echo "upstream seed $seed: done"
    fi
    run=${FORK_LABEL}_$seed
    rm -rf "$run"
    /usr/bin/time -f "%e s wall, %M KB" -o "$run.time" env PATH="$FORK_PATH" \
        "$FORK_PREFIX/bin/python" "$FORK_REPO/scripts/matam_assembly.py" \
        --forward reads/legacy_R1.fq --reverse reads/legacy_R2.fq -d "db_$FORK_LABEL/refs" -o "$run" \
        --cpu "$CPU" --max_memory 4000 -v --seed "$seed" > "$run.log" 2>&1
    echo "$FORK_LABEL seed $seed: done"
done
