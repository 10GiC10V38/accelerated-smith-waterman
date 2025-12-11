#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <immintrin.h>
#include <omp.h>
#include <stdint.h>

// --- CONFIGURATION ---
#define MATCH 2
#define MISMATCH -1 // For 32-bit
#define GAP -2      // For 32-bit

// For 16-bit Unsigned Logic (Magnitudes to subtract)
#define MISMATCH_U 1 
#define GAP_U 2

#define BLOCK_SIZE 256 
#define ALIGN 32
#define MAX(a,b) ((a) > (b) ? (a) : (b))
#define CUTOFF_N 100000 // Threshold to switch from 16-bit to 32-bit

double get_time() {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + ts.tv_nsec * 1e-9;
}

void* aligned_malloc(size_t size) {
    void* ptr;
    if (posix_memalign(&ptr, ALIGN, size) != 0) return NULL;
    return ptr;
}

void generate_sequence(char* restrict seq, int n) {
    const char alphabet[] = "ACGT";
    for (int i = 0; i < n; i++) seq[i] = alphabet[rand() % 4];
    seq[n] = '\0';
}

// ============================================================================
// STRATEGY 1: 16-BIT UNSIGNED AVX2 (Fastest, Max Score 65535)
// ============================================================================

int process_block_avx2_16bit(
    const char* restrict seq1, const char* restrict seq2,
    int rows, int cols,
    unsigned short* restrict top_buf, unsigned short* restrict left_buf,
    unsigned short corner_val
) {
    unsigned short buf1[BLOCK_SIZE] __attribute__((aligned(ALIGN)));
    unsigned short buf2[BLOCK_SIZE] __attribute__((aligned(ALIGN)));
    unsigned short* prev = buf1;
    unsigned short* curr = buf2;
    
    memcpy(prev, top_buf, cols * sizeof(unsigned short));

    const __m256i vZERO = _mm256_setzero_si256();
    const __m256i vGap1 = _mm256_set1_epi16((short)GAP_U);
    const __m256i vGap2 = _mm256_set1_epi16((short)(GAP_U * 2));
    const __m256i vGap4 = _mm256_set1_epi16((short)(GAP_U * 4));
    
    // Left Offsets
    unsigned short left_offs[16];
    for(int k=0; k<16; k++) left_offs[k] = (unsigned short)(GAP_U * (k + 1));
    const __m256i vLeftOffsets = _mm256_loadu_si256((__m256i*)left_offs);
    
    // Cross Offsets
    unsigned short cross_offs[16];
    for(int k=0; k<8; k++) cross_offs[k] = 0;
    for(int k=8; k<16; k++) cross_offs[k] = (unsigned short)(GAP_U * (k - 7));
    const __m256i vCrossOffsets = _mm256_loadu_si256((__m256i*)cross_offs);

    __m256i vBlockMax = vZERO;
    unsigned short block_corner_prev = corner_val;

    for (int i = 0; i < rows; i++) {
        __m256i vBase = _mm256_set1_epi16((short)seq1[i]);
        unsigned short current_left = left_buf[i];
        unsigned short diag_val = block_corner_prev;
        block_corner_prev = left_buf[i];

        int j = 0;
        for (; j <= cols - 16; j += 16) {
            __m256i vPrev = _mm256_loadu_si256((__m256i*)&prev[j]);
            __m128i vSeq2Raw = _mm_loadu_si128((__m128i*)&seq2[j]);
            __m256i vSeq2 = _mm256_cvtepu8_epi16(vSeq2Raw);

            __m256i vDiag;
            if (j == 0) {
                unsigned short temp[16]; temp[0] = diag_val; memcpy(&temp[1], &prev[0], 15 * sizeof(unsigned short));
                vDiag = _mm256_loadu_si256((__m256i*)temp);
            } else {
                vDiag = _mm256_loadu_si256((__m256i*)&prev[j - 1]);
            }
            
            __m256i vMatchVal = _mm256_adds_epu16(vDiag, _mm256_set1_epi16(MATCH));
            __m256i vMismatchVal = _mm256_subs_epu16(vDiag, _mm256_set1_epi16(MISMATCH_U));
            __m256i vMatchScore = _mm256_blendv_epi8(vMismatchVal, vMatchVal, _mm256_cmpeq_epi16(vBase, vSeq2));

            __m256i vUpScore = _mm256_subs_epu16(vPrev, vGap1);
            __m256i vMax = _mm256_max_epu16(vMatchScore, vUpScore);

            __m256i vLeftProp = _mm256_subs_epu16(_mm256_set1_epi16(current_left), vLeftOffsets);
            vMax = _mm256_max_epu16(vMax, vLeftProp);

            __m256i vShift = _mm256_bslli_epi128(vMax, 2);
            vMax = _mm256_max_epu16(vMax, _mm256_subs_epu16(vShift, vGap1));
            vShift = _mm256_bslli_epi128(vMax, 4);
            vMax = _mm256_max_epu16(vMax, _mm256_subs_epu16(vShift, vGap2));
            vShift = _mm256_bslli_epi128(vMax, 8);
            vMax = _mm256_max_epu16(vMax, _mm256_subs_epu16(vShift, vGap4));

            unsigned short lane0_end = (unsigned short)_mm256_extract_epi16(vMax, 7);
            __m256i vCross = _mm256_subs_epu16(_mm256_set1_epi16(lane0_end), vCrossOffsets);
            __m256i vUpperUpdate = _mm256_max_epu16(vMax, vCross);
            vMax = _mm256_blend_epi32(vMax, vUpperUpdate, 0xF0);

            _mm256_storeu_si256((__m256i*)&curr[j], vMax);
            vBlockMax = _mm256_max_epu16(vBlockMax, vMax);
            
            current_left = (unsigned short)_mm256_extract_epi16(vMax, 15);
            diag_val = prev[j + 15]; 
        }

        for (; j < cols; j++) {
            unsigned short up_val = prev[j];
            int d = (int)diag_val + ((seq1[i] == seq2[j]) ? MATCH : -MISMATCH_U);
            int u = (int)up_val - GAP_U;
            int l = (int)current_left - GAP_U;
            int val = MAX(0, MAX(d, MAX(u, l)));
            if (val > 65535) val = 65535;
            unsigned short uval = (unsigned short)val;
            __m256i vVal = _mm256_set1_epi16(uval);
            vBlockMax = _mm256_max_epu16(vBlockMax, vVal);
            diag_val = up_val;
            current_left = uval;
            curr[j] = uval;
        }
        left_buf[i] = current_left;
        unsigned short* tmp = prev; prev = curr; curr = tmp;
    }
    memcpy(top_buf, prev, cols * sizeof(unsigned short));
    unsigned short temp_max[16];
    _mm256_storeu_si256((__m256i*)temp_max, vBlockMax);
    unsigned short local_max = 0;
    for(int k=0; k<16; k++) if(temp_max[k] > local_max) local_max = temp_max[k];
    return (int)local_max;
}

int run_sw_16bit(const char* restrict seq1, const char* restrict seq2, int n1, int n2) {
    int n_blocks_row = (n1 + BLOCK_SIZE - 1) / BLOCK_SIZE;
    int n_blocks_col = (n2 + BLOCK_SIZE - 1) / BLOCK_SIZE;
    int global_max = 0;

    unsigned short* H_horizontal = (unsigned short*)calloc(n2 + BLOCK_SIZE, sizeof(unsigned short)); 
    unsigned short* H_vertical   = (unsigned short*)calloc(n1 + BLOCK_SIZE, sizeof(unsigned short));
    unsigned short* corners = (unsigned short*)calloc((n_blocks_row + 1) * (n_blocks_col + 1), sizeof(unsigned short));

    int num_waves = n_blocks_row + n_blocks_col - 1;

    for (int w = 0; w < num_waves; ++w) {
        int r_start = (w - n_blocks_col + 1 > 0) ? (w - n_blocks_col + 1) : 0;
        int r_end   = (w + 1 < n_blocks_row) ? (w + 1) : n_blocks_row;

        #pragma omp parallel for reduction(max:global_max) schedule(dynamic)
        for (int r = r_start; r < r_end; ++r) {
            int c = w - r;
            int i_start = r * BLOCK_SIZE;
            int j_start = c * BLOCK_SIZE;
            int i_len = (i_start + BLOCK_SIZE > n1) ? (n1 - i_start) : BLOCK_SIZE;
            int j_len = (j_start + BLOCK_SIZE > n2) ? (n2 - j_start) : BLOCK_SIZE;

            unsigned short* top_ptr  = &H_horizontal[j_start];
            unsigned short* left_ptr = &H_vertical[i_start];
            unsigned short corner_val = (r > 0 && c > 0) ? corners[(r-1) * n_blocks_col + (c-1)] : 0;

            int block_max = process_block_avx2_16bit(
                &seq1[i_start], &seq2[j_start],
                i_len, j_len,
                top_ptr, left_ptr,
                corner_val
            );

            if (block_max > global_max) global_max = block_max;
            corners[r * n_blocks_col + c] = top_ptr[j_len - 1];
        }
    }
    free(H_horizontal); free(H_vertical); free(corners);
    return global_max;
}

// ============================================================================
// STRATEGY 2: 32-BIT SIGNED AVX2 (Safe, Max Score 2 Billion)
// ============================================================================

int process_block_avx2_32bit(
    const char* restrict seq1, const char* restrict seq2,
    int rows, int cols,
    int* restrict top_buf, int* restrict left_buf,
    int corner_val
) {
    int buf1[BLOCK_SIZE] __attribute__((aligned(ALIGN)));
    int buf2[BLOCK_SIZE] __attribute__((aligned(ALIGN)));
    int* prev = buf1;
    int* curr = buf2;
    memcpy(prev, top_buf, cols * sizeof(int));

    const __m256i vGAP      = _mm256_set1_epi32(GAP);
    const __m256i vZERO     = _mm256_setzero_si256();
    const __m256i vMATCH    = _mm256_set1_epi32(MATCH);
    const __m256i vMISMATCH = _mm256_set1_epi32(MISMATCH);
    const __m256i vGap1 = _mm256_set1_epi32(GAP);
    const __m256i vGap2 = _mm256_set1_epi32(GAP * 2);
    const __m256i vGap4 = _mm256_set1_epi32(GAP * 4);
    const __m256i vLeftOffsets = _mm256_setr_epi32(GAP, GAP*2, GAP*3, GAP*4, GAP*5, GAP*6, GAP*7, GAP*8);
    const __m256i vCrossOffsets = _mm256_setr_epi32(0, 0, 0, 0, GAP, GAP*2, GAP*3, GAP*4);

    __m256i vBlockMax = vZERO;
    int block_corner_prev = corner_val;

    for (int i = 0; i < rows; i++) {
        __m256i vBase = _mm256_set1_epi32(seq1[i]);
        int current_left = left_buf[i];
        int diag_val = block_corner_prev;
        block_corner_prev = left_buf[i];

        int j = 0;
        for (; j <= cols - 8; j += 8) {
            __m256i vPrev = _mm256_loadu_si256((__m256i*)&prev[j]);
            __m256i vSeq2 = _mm256_cvtepi8_epi32(_mm_loadl_epi64((__m128i*)&seq2[j]));
            __m256i vMatchScore = _mm256_blendv_epi8(vMISMATCH, vMATCH, _mm256_cmpeq_epi32(vBase, vSeq2));

            __m256i vDiag;
            if (j == 0) {
                int temp[8]; temp[0] = diag_val; memcpy(&temp[1], &prev[0], 7 * sizeof(int));
                vDiag = _mm256_loadu_si256((__m256i*)temp);
            } else {
                vDiag = _mm256_loadu_si256((__m256i*)&prev[j - 1]);
            }

            __m256i vMax = _mm256_max_epi32(_mm256_add_epi32(vDiag, vMatchScore), _mm256_add_epi32(vPrev, vGAP));
            vMax = _mm256_max_epi32(vMax, vZERO);

            __m256i vLeftProp = _mm256_add_epi32(_mm256_set1_epi32(current_left), vLeftOffsets);
            vMax = _mm256_max_epi32(vMax, vLeftProp);

            __m256i vShift = _mm256_bslli_epi128(vMax, 4);
            vMax = _mm256_max_epi32(vMax, _mm256_add_epi32(vShift, vGap1));
            vShift = _mm256_bslli_epi128(vMax, 8);
            vMax = _mm256_max_epi32(vMax, _mm256_add_epi32(vShift, vGap2));
            
            int lane0_end = _mm256_extract_epi32(vMax, 3);
            __m256i vCross = _mm256_add_epi32(_mm256_set1_epi32(lane0_end), vCrossOffsets);
            __m256i vUpperUpdate = _mm256_max_epi32(vMax, vCross);
            vMax = _mm256_blend_epi32(vMax, vUpperUpdate, 0xF0);

            _mm256_storeu_si256((__m256i*)&curr[j], vMax);
            vBlockMax = _mm256_max_epi32(vBlockMax, vMax);
            current_left = _mm256_extract_epi32(vMax, 7);
            diag_val = prev[j + 7]; 
        }

        for (; j < cols; j++) {
            int up_val = prev[j];
            int score = (seq1[i] == seq2[j]) ? MATCH : MISMATCH;
            int val = MAX(0, MAX(diag_val + score, MAX(up_val + GAP, current_left + GAP)));
            __m256i vVal = _mm256_set1_epi32(val);
            vBlockMax = _mm256_max_epi32(vBlockMax, vVal);
            diag_val = up_val;
            current_left = val;
            curr[j] = val;
        }
        left_buf[i] = current_left;
        int* tmp = prev; prev = curr; curr = tmp;
    }
    memcpy(top_buf, prev, cols * sizeof(int));
    int temp_max[8];
    _mm256_storeu_si256((__m256i*)temp_max, vBlockMax);
    int local_max = 0;
    for(int k=0; k<8; k++) if(temp_max[k] > local_max) local_max = temp_max[k];
    return local_max;
}

int run_sw_32bit(const char* restrict seq1, const char* restrict seq2, int n1, int n2) {
    int n_blocks_row = (n1 + BLOCK_SIZE - 1) / BLOCK_SIZE;
    int n_blocks_col = (n2 + BLOCK_SIZE - 1) / BLOCK_SIZE;
    int global_max = 0;

    int* H_horizontal = (int*)calloc(n2 + BLOCK_SIZE, sizeof(int)); 
    int* H_vertical   = (int*)calloc(n1 + BLOCK_SIZE, sizeof(int));
    int* corners = (int*)calloc((n_blocks_row + 1) * (n_blocks_col + 1), sizeof(int));

    int num_waves = n_blocks_row + n_blocks_col - 1;

    for (int w = 0; w < num_waves; ++w) {
        int r_start = (w - n_blocks_col + 1 > 0) ? (w - n_blocks_col + 1) : 0;
        int r_end   = (w + 1 < n_blocks_row) ? (w + 1) : n_blocks_row;

        #pragma omp parallel for reduction(max:global_max) schedule(dynamic)
        for (int r = r_start; r < r_end; ++r) {
            int c = w - r;
            int i_start = r * BLOCK_SIZE;
            int j_start = c * BLOCK_SIZE;
            int i_len = (i_start + BLOCK_SIZE > n1) ? (n1 - i_start) : BLOCK_SIZE;
            int j_len = (j_start + BLOCK_SIZE > n2) ? (n2 - j_start) : BLOCK_SIZE;

            int* top_ptr  = &H_horizontal[j_start];
            int* left_ptr = &H_vertical[i_start];
            int corner_val = (r > 0 && c > 0) ? corners[(r-1) * n_blocks_col + (c-1)] : 0;

            int block_max = process_block_avx2_32bit(
                &seq1[i_start], &seq2[j_start],
                i_len, j_len,
                top_ptr, left_ptr,
                corner_val
            );

            if (block_max > global_max) global_max = block_max;
            corners[r * n_blocks_col + c] = top_ptr[j_len - 1];
        }
    }
    free(H_horizontal); free(H_vertical); free(corners);
    return global_max;
}

// ============================================================================
// MAIN DRIVER
// ============================================================================

int main(int argc, char **argv) {
    int N = 0;
    if (argc > 1) N = atoi(argv[1]);
    if (N <= 0) { printf("Enter N: "); if (scanf("%d", &N) != 1) N = 2048; }
    
    srand(42);
    char *seq1 = (char*)aligned_malloc(N + 64);
    char *seq2 = (char*)aligned_malloc(N + 64);
    generate_sequence(seq1, N); generate_sequence(seq2, N);

    int score = 0;
    double start = get_time();

    if (N <= CUTOFF_N) {
        printf("Strategy: 16-bit Unsigned (High Speed, Max Score 65535)\n");
        score = run_sw_16bit(seq1, seq2, N, N);
    } else {
        printf("Strategy: 32-bit Signed (Robust, Unlimited Score)\n");
        score = run_sw_32bit(seq1, seq2, N, N);
    }
    
    double end = get_time();

    printf("N=%d\n", N);
    printf("Optimized Score: %d\n", score);
    printf("Optimized Time:  %.4fs\n", end - start);
    printf("Throughput:      %.2f GCUPS\n", (double)N*N/(end-start)/1e9);
    
    free(seq1); free(seq2);
    return 0;
}
