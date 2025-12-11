# Accelerated Smith-Waterman Algorithm (AVX2 + OpenMP)

![Language](https://img.shields.io/badge/language-C-blue.svg) ![Standard](https://img.shields.io/badge/SIMD-AVX2-green.svg) ![Parallelism](https://img.shields.io/badge/OpenMP-Wavefront-orange.svg)

## 🚀 Project Overview
This project focuses on accelerating the **Smith-Waterman local sequence alignment algorithm**, a fundamental but computationally intensive ($O(N^2)$) bioinformatics tool. Starting from a naive scalar C implementation, this work applies advanced architectural optimizations to achieve a **~65x speedup** over the baseline.

## 📊 Performance Highlights
* **Throughput:** Increased from **0.11 GCUPS** (Baseline) to **7.10 GCUPS** (Optimized).
* **Speedup:** Execution time for $N=100,000$ dropped from timeout (>infinite) to **1.4 seconds**.
* **Efficiency:** Reduced L1 Cache misses by **99.4%** and branch mispredictions by **99.5%**.

## 🧠 The Challenge
The standard algorithm suffers from three main bottlenecks:
1.  **Data Dependencies:** Every cell $H[i][j]$ depends on neighbors, preventing simple row-parallelization.
2.  **Memory Bound:** Accessing large matrices causes frequent cache misses.
3.  **Control Hazards:** The `max()` function creates frequent branch mispredictions.

## ⚡ The Solution: Hybrid Architecture
The implementation (`sw_opt.c`) uses a hybrid strategy that adapts based on sequence length ($N$):

### 1. Cache Blocking (Tiling)
The matrix is divided into $256 \times 256$ tiles that fit entirely in the L1 Cache, minimizing slow RAM access.

### 2. Wavefront Parallelism
Instead of row-by-row processing, we process **diagonals (waves)** of blocks. Blocks in the same diagonal are independent and computed in parallel using **OpenMP** threads.

### 3. SIMD Vectorization (AVX2)
* **Strategy A (Fastest):** For $N \le 100,000$, uses **16-bit unsigned integers** (16 cells/cycle). Handles negative gap penalties via unsigned saturation arithmetic.
* **Strategy B (Safest):** For $N > 100,000$, falls back to **32-bit signed integers** (8 cells/cycle) to prevent score overflow.

## 📂 Repository Structure
```text
.
├── baseline/           # Naive scalar implementation and analysis
├── optimized/          # AVX2/OpenMP optimized implementation
├── report.pdf          # Detailed architectural analysis
└── run.sh              # Master execution script
