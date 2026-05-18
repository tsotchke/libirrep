/* SPDX-License-Identifier: MIT */
/* Surface-code syndrome distribution under depolarising noise.
 *
 * Builds the rotated `[[d², 1, d]]` surface code at fixed distance and
 * samples random single-qubit Pauli errors at a parametric physical
 * error rate `p`. For each sample, computes the X- and Z-syndrome
 * weights (number of triggered stabilizers) and aggregates the
 * distribution.
 *
 * The relationship between (p, d) and syndrome weight is the
 * organising parameter of fault-tolerance threshold studies:
 *
 *   - At low p: most rounds are syndrome-empty; non-empty syndromes
 *     have small weight (2 — the two ends of a single-edge error
 *     string).
 *   - At threshold: syndrome weight crosses a percolation-like
 *     transition where decoders begin to fail systematically.
 *
 * This example does NOT decode. It demonstrates the building blocks:
 * surface-code construction, abstract stabilizer-group materialisation,
 * Pauli operator manipulation, and `irrep_stabilizer_syndrome`. A
 * full fault-tolerance pipeline would feed these syndromes to a
 * minimum-weight perfect matching (MWPM) or BP+OSD decoder; both are
 * straightforward to bolt on outside the library.
 *
 * Build / run:
 *   make examples
 *   ./build/bin/qec_surface_syndrome_distribution
 */

#include <irrep/css_code.h>
#include <irrep/stabilizer_group.h>
#include <irrep/surface_code.h>

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Deterministic xorshift64 — same seed → same output across machines. */
typedef struct {
    uint64_t s;
} xorshift64_t;

static uint64_t xorshift64_next(xorshift64_t *r) {
    uint64_t x = r->s;
    x ^= x << 13;
    x ^= x >> 7;
    x ^= x << 17;
    r->s = x;
    return x;
}

static double xorshift64_uniform(xorshift64_t *r) {
    return (xorshift64_next(r) >> 11) * (1.0 / (1ULL << 53));
}

/* Apply a depolarising Pauli error at rate p to each qubit:
 *   with probability p/3 each: X, Y, Z (mutually exclusive)
 *   with probability 1 - p: I  */
static void
inject_depolarising_error(xorshift64_t *r, double p, irrep_pauli_t *err)
{
    int n = err->n;
    for (int q = 0; q < n; ++q) {
        double u = xorshift64_uniform(r);
        if (u < p / 3.0) {
            irrep_pauli_set(err, q, IRREP_PAULI_LETTER_X);
        } else if (u < 2.0 * p / 3.0) {
            irrep_pauli_set(err, q, IRREP_PAULI_LETTER_Y);
        } else if (u < p) {
            irrep_pauli_set(err, q, IRREP_PAULI_LETTER_Z);
        }
    }
}

static int popcount_u64(uint64_t x) {
    int c = 0;
    while (x) { x &= x - 1; ++c; }
    return c;
}

static int hamming_weight(const uint64_t *bits, int n_bits) {
    int n_words = (n_bits + 63) / 64;
    int w = 0;
    for (int i = 0; i < n_words; ++i) w += popcount_u64(bits[i]);
    return w;
}

static void
print_histogram(const char *label, const int *hist, int max_w, int n_samples)
{
    printf("    %-12s : ", label);
    for (int w = 0; w <= max_w; ++w) {
        if (hist[w] > 0) {
            printf("[%d]:%4.1f%% ", w, 100.0 * (double)hist[w] / n_samples);
        }
    }
    printf("\n");
}

int main(void) {
    int distances[] = { 3, 5, 7 };
    double rates[] = { 0.005, 0.01, 0.02, 0.05 };
    int n_samples = 50000;

    printf("Surface-code syndrome distribution under depolarising noise\n");
    printf("(deterministic xorshift64, seed = 0x12345678, n_samples = %d)\n\n",
           n_samples);
    printf("Each row: distribution of syndrome Hamming weight at distance d,\n");
    printf("physical error rate p.\n\n");

    for (size_t di = 0; di < sizeof(distances)/sizeof(distances[0]); ++di) {
        int d = distances[di];
        irrep_surface_params_t p;
        if (irrep_surface_init(&p, d) != IRREP_OK) {
            fprintf(stderr, "init failed for d=%d\n", d);
            return 1;
        }
        irrep_css_code_t css;
        if (irrep_surface_build(&p, &css) != IRREP_OK) {
            fprintf(stderr, "build failed for d=%d\n", d);
            return 1;
        }
        irrep_stabilizer_group_t gx, gz;
        /* Build separate X-only / Z-only stabilizer groups for syndrome
         * extraction. Z-syndrome from X-stabilizers (detects Z errors);
         * X-syndrome from Z-stabilizers (detects X errors). */
        if (irrep_stabilizer_group_new(&gx, p.n_qubits, p.n_X_stabs) != IRREP_OK) return 1;
        for (int i = 0; i < p.n_X_stabs; ++i) {
            for (int j = 0; j < p.n_qubits; ++j) {
                if (irrep_parity_matrix_get(&css.H_X, i, j))
                    irrep_pauli_set(&gx.gens[i], j, IRREP_PAULI_LETTER_X);
            }
        }
        if (irrep_stabilizer_group_new(&gz, p.n_qubits, p.n_Z_stabs) != IRREP_OK) return 1;
        for (int i = 0; i < p.n_Z_stabs; ++i) {
            for (int j = 0; j < p.n_qubits; ++j) {
                if (irrep_parity_matrix_get(&css.H_Z, i, j))
                    irrep_pauli_set(&gz.gens[i], j, IRREP_PAULI_LETTER_Z);
            }
        }

        printf("d=%d  [[n=%d, 1, %d]]  n_X=%d  n_Z=%d\n",
               d, p.n_qubits, d, p.n_X_stabs, p.n_Z_stabs);

        for (size_t pi = 0; pi < sizeof(rates)/sizeof(rates[0]); ++pi) {
            double pe = rates[pi];
            xorshift64_t r = { .s = 0x12345678ULL + (uint64_t)d * 1000 + (uint64_t)(pe * 1e6) };

            int max_w = p.n_X_stabs > p.n_Z_stabs ? p.n_X_stabs : p.n_Z_stabs;
            int *hist_x = calloc((size_t)max_w + 1, sizeof(int));
            int *hist_z = calloc((size_t)max_w + 1, sizeof(int));

            irrep_pauli_t err;
            irrep_pauli_new(&err, p.n_qubits);
            int n_syndrome_words = (max_w + 63) / 64;
            uint64_t *sx = calloc((size_t)n_syndrome_words, sizeof(uint64_t));
            uint64_t *sz = calloc((size_t)n_syndrome_words, sizeof(uint64_t));

            for (int sample = 0; sample < n_samples; ++sample) {
                /* Reset err to identity. */
                for (int w = 0; w < err.n_words; ++w) {
                    err.x[w] = 0;
                    err.z[w] = 0;
                }
                err.sign = IRREP_PAULI_SIGN_POS;
                inject_depolarising_error(&r, pe, &err);

                /* X-syndrome: which Z-stabs anti-commute with err. */
                irrep_stabilizer_syndrome(&gz, &err, sz);
                /* Z-syndrome: which X-stabs anti-commute with err. */
                irrep_stabilizer_syndrome(&gx, &err, sx);

                int wx = hamming_weight(sx, p.n_X_stabs);
                int wz = hamming_weight(sz, p.n_Z_stabs);
                hist_x[wx]++;
                hist_z[wz]++;
            }

            printf("  p=%5.3f  (mean qubit errors ≈ %.1f)\n",
                   pe, pe * p.n_qubits);
            print_histogram("X-syndrome", hist_x, max_w, n_samples);
            print_histogram("Z-syndrome", hist_z, max_w, n_samples);

            free(hist_x);
            free(hist_z);
            free(sx);
            free(sz);
            irrep_pauli_free(&err);
        }
        printf("\n");

        irrep_stabilizer_group_free(&gx);
        irrep_stabilizer_group_free(&gz);
        irrep_css_code_free(&css);
    }

    return 0;
}
