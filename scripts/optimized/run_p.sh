#!/bin/bash
set -eu

# ============================
# CONFIG
# ============================

SRC="sw_opt.c"
BIN="sw_opt"
N=100000
THREADS=(1 2 4 8)

OUTFILE="perf_N100000.txt"

echo "==============================================="
echo "  Compiling $SRC → $BIN"
echo "==============================================="

gcc -O3 -march=native -mavx2 -fopenmp -funroll-loops -o "$BIN" "$SRC"

echo
echo "==============================================="
echo " Running PERf for N = $N   (Threads: 1,2,4,8)"
echo " Results saved to $OUTFILE"
echo "==============================================="
echo

# Redirect stdout+stderr → OUTFILE
exec > "$OUTFILE" 2>&1

for T in "${THREADS[@]}"; do

    echo "==============================================="
    echo "      N = $N   |   Threads = $T"
    echo "==============================================="
    echo

    # FIX: Pass variable specifically to the command running under sudo
    
    echo "---- Basic Perf Counters ----"
    sudo OMP_NUM_THREADS=$T perf stat \
        -e task-clock,cycles,instructions,branches,branch-misses \
        -e cache-references,cache-misses \
        "./$BIN" "$N"

    echo
    echo "---- L1 Cache ----"
    sudo OMP_NUM_THREADS=$T perf stat \
        -e L1-dcache-loads,L1-dcache-load-misses \
        "./$BIN" "$N"

    echo
    echo "---- LLC / Memory ----"
    sudo OMP_NUM_THREADS=$T perf stat \
        -e LLC-loads,LLC-load-misses \
        -e cache-references,cache-misses \
        "./$BIN" "$N" || true

    echo -e "\n\n"

done
