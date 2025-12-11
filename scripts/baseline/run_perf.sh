#!/bin/bash
set -eu

# ============================
# CONFIG
# ============================

SRC="sw_baseline.c"
BIN="sw_baseline"
PATH_TO_CODE="/home/ritik/Downloads/hpc5b/codeBase/code/baseline"

# Sequence sizes you want to test
Ns=(1000 5000 10000 50000 100000)

# Thread counts (if your code uses OMP)
Threads=(1 2 4 8)

# Output file
OUTFILE="perf_results.txt"

# Redirect all output (stdout + stderr)
exec >> "$OUTFILE" 2>&1

echo "==============================================="
echo " Building $SRC  →  $BIN"
echo "==============================================="
echo

cd "$PATH_TO_CODE"

gcc -O3 -mavx2 -fopenmp -march=native -o "$BIN" "$SRC"

echo
echo "========== Running PERf for smith-waterman =========="
echo

for N in "${Ns[@]}"; do
    for T in "${Threads[@]}"; do

        echo "==============================================="
        echo "N = $N    Threads = $T"
        echo "==============================================="
        echo

        export OMP_NUM_THREADS=$T

        echo "---- Perf Basic Stats ----"
        sudo perf stat \
            -e cycles \
            -e instructions \
            -e branches \
            -e branch-misses \
            -e cache-references \
            -e cache-misses \
            ./"$BIN" "$N"

        echo
        echo "---- Perf L1 Cache ----"
        sudo perf stat \
            -e L1-dcache-loads \
            -e L1-dcache-load-misses \
            ./"$BIN" "$N"

        echo
        echo "---- Perf LLC / Memory ----"
        sudo perf stat \
            -e LLC-loads \
            -e LLC-load-misses \
            -e cache-references \
            -e cache-misses \
            ./"$BIN" "$N" || true

        echo -e "\n\n"

    done
done

