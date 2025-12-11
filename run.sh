#!/usr/bin/env bash
set -euo pipefail

# ============================================================
#  run.sh — Build + run + perf for sw_opt.c
#
#  Directory structure:
#    code/run.sh
#    code/optimized/sw_opt.c
#    code/optimized/sw_opt         (binary after compile)
#
#  Usage:
#    ./run.sh build
#    ./run.sh run <N> <threads>
#    ./run.sh perf <N> <threads>
#
#  Example:
#    ./run.sh build
#    ./run.sh run 2048 8
#    ./run.sh perf 4096 8
# ============================================================

# Paths
SRC="optimized/sw_opt.c"
BIN="optimized/sw_opt"

# Compiler flags
CFLAGS="-O3 -march=native -funroll-loops -mavx2 -fopenmp -Wall -Wextra"

# ============================================================
# BUILD
# ============================================================
build_sw_opt() {
    echo "====================================================="
    echo " Compiling $SRC → $BIN"
    echo "====================================================="

    gcc $CFLAGS "$SRC" -o "$BIN"
    echo "Build complete."
}

# ============================================================
# RUN (no perf)
# ============================================================
run_sw_opt() {
    if [ $# -lt 2 ]; then
        echo "Usage: ./run.sh run <N> <threads>"
        exit 1
    fi

    N=$1
    T=$2

    # Auto-build if binary missing
    if [ ! -x "$BIN" ]; then
        build_sw_opt
    fi

    echo "====================================================="
    echo " Running optimized Smith–Waterman"
    echo "   N = $N"
    echo "   Threads = $T"
    echo "====================================================="

    "$BIN" "$N" "$T"
}

# ============================================================
# PERF MODE
# ============================================================
run_perf() {
    if [ $# -lt 2 ]; then
        echo "Usage: ./run.sh perf <N> <threads>"
        exit 1
    fi

    N=$1
    T=$2

    if [ ! -x "$BIN" ]; then
        build_sw_opt
    fi

    echo "====================================================="
    echo " Running perf for optimized SW"
    echo "====================================================="

    sudo perf stat \
        -e task-clock,cycles,instructions,branches,branch-misses,page-faults,context-switches,cpu-migrations \
        "$BIN" "$N" "$T"
}

# ============================================================
# MODE SELECTOR
# ============================================================

MODE=${1:-}

case "$MODE" in
    build)
        build_sw_opt
        ;;
    run)
        shift
        run_sw_opt "$@"
        ;;
    perf)
        shift
        run_perf "$@"
        ;;
    *)
        echo "Unknown mode: $MODE"
        echo "Usage:"
        echo "  ./run.sh build"
        echo "  ./run.sh run <N> <threads>"
        echo "  ./run.sh perf <N> <threads>"
        exit 1
        ;;
esac

