# 1. Compiler and Flags
CC = gcc
# Common flags: Optimize code (-O3) and show warnings (-Wall)
CFLAGS = -O3 -Wall 
# Optimization flags: Enable AVX2 SIMD and OpenMP threading
OPT_FLAGS = -mavx2 -fopenmp

# 2. Targets
# 'all' is the default target when you just type 'make'
all: baseline_bin optimized_bin

# 3. Rules
# To build 'baseline_bin', we need 'baseline/sw_baseline.c'
baseline_bin: baseline/sw_baseline.c
	$(CC) $(CFLAGS) baseline/sw_baseline.c -o baseline/sw_baseline

# To build 'optimized_bin', we need 'optimized/sw_opt.c'
optimized_bin: optimized/sw_opt.c
	$(CC) $(CFLAGS) $(OPT_FLAGS) optimized/sw_opt.c -o optimized/sw_opt

# 4. Clean
# Typing 'make clean' deletes the created binaries (keeps folder clean)
clean:
	rm -f baseline/sw_baseline optimized/sw_opt
