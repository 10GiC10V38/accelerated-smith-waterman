#!/usr/bin/env bash
set -euo pipefail

# ============================================================
#  run.sh — Wrapper for Smith-Waterman Optimization
# ============================================================

# Paths
BIN="optimized/sw_opt"

# ============================================================
# BUILD (Delegates to Makefile)
# ============================================================
build_sw_opt() {
    echo "====================================================="
    echo " Building project via Makefile..."
    echo "====================================================="
    
    # This ensures we use the exact same flags defined in Makefile
    make optimized_bin
    
    echo "Build complete."
}

# ============================================================
# RUN (Standard Execution)
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
# PERF MODE (Linux Only)
# ============================================================
run_perf() {
    if [ $# -lt 2 ]; then
        echo "Usage: ./run.sh perf <N> <threads>"
        exit 1
    fi

    # Check if user is on Linux before running perf
    if [[ "$OSTYPE" != "linux-gnu"* ]]; then
        echo "Error: 'perf' is a Linux-specific tool."
        echo "On macOS, try using 'Instruments' or 'time' instead."
        exit 1
    fi

    N=$1
    T=$2

    if [ ! -x "$BIN" ]; then
        build_sw_opt
    fi

    echo "====================================================="
    echo " Running perf (requires sudo)"
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
