#!/usr/bin/env bash
set -euo pipefail

# ============================================================
#  run.sh — Smith-Waterman Benchmark Runner
# ============================================================

# Paths
OPT_BIN="optimized/sw_opt"
BASE_BIN="baseline/sw_baseline"

# ============================================================
# BUILD (Delegates to Makefile)
# ============================================================
build_all() {
    echo "====================================================="
    echo " Building all targets via Makefile..."
    echo "====================================================="
    make all
    echo "Build complete."
}

# ============================================================
# RUN FUNCTION
# ============================================================
run_program() {
    TYPE=$1    # "opt" or "base"
    N=$2       # Sequence Length
    T=$3       # Threads (only used for opt)

    if [ "$TYPE" == "opt" ]; then
        BINARY="$OPT_BIN"
        LABEL="Optimized (AVX2 + OpenMP)"
        ARGS="$N $T"
    else
        BINARY="$BASE_BIN"
        LABEL="Baseline (Scalar)"
        # Baseline might not take thread arg, but we pass N. 
        # Adjust if sw_baseline.c arguments differ!
        ARGS="$N" 
    fi

    # Auto-build if missing
    if [ ! -x "$BINARY" ]; then
        build_all
    fi

    echo "====================================================="
    echo " Running $LABEL"
    echo "   N = $N"
    if [ "$TYPE" == "opt" ]; then echo "   Threads = $T"; fi
    echo "====================================================="

    # Time the execution (simple user-friendly timing)
    time "$BINARY" $ARGS
}

# ============================================================
# PERF FUNCTION (Linux Only)
# ============================================================
run_perf() {
    TYPE=$1
    N=$2
    T=$3

    if [[ "$OSTYPE" != "linux-gnu"* ]]; then
        echo "Error: 'perf' is Linux-only."
        exit 1
    fi

    if [ "$TYPE" == "opt" ]; then
        BINARY="$OPT_BIN"
        ARGS="$N $T"
    else
        BINARY="$BASE_BIN"
        ARGS="$N"
    fi

    if [ ! -x "$BINARY" ]; then build_all; fi

    echo "====================================================="
    echo " Running PERF on $TYPE (sudo required)"
    echo "====================================================="
    
    sudo perf stat \
        -e task-clock,cycles,instructions,branches,branch-misses,L1-dcache-loads,L1-dcache-load-misses \
        "$BINARY" $ARGS
}

# ============================================================
# INPUT PARSING
# ============================================================

MODE=${1:-help}

case "$MODE" in
    build)
        build_all
        ;;
    run)
        # Usage: ./run.sh run opt 2048 8
        if [ $# -lt 3 ]; then
            echo "Usage: ./run.sh run <opt|base> <N> [threads]"
            exit 1
        fi
        TYPE=$2
        N=$3
        THREADS=${4:-1} # Default to 1 if not provided
        run_program "$TYPE" "$N" "$THREADS"
        ;;
    perf)
        if [ $# -lt 3 ]; then
            echo "Usage: ./run.sh perf <opt|base> <N> [threads]"
            exit 1
        fi
        TYPE=$2
        N=$3
        THREADS=${4:-1}
        run_perf "$TYPE" "$N" "$THREADS"
        ;;
    *)
        echo "Usage:"
        echo "  ./run.sh build"
        echo "  ./run.sh run opt <N> <threads>   (Run Optimized)"
        echo "  ./run.sh run base <N>            (Run Baseline)"
        echo "  ./run.sh perf opt <N> <threads>  (Profile Optimized)"
        exit 1
        ;;
esac
