/* SPDX-License-Identifier: MIT */
/** @file floquet_code.c
 *  @brief Implementation of the Floquet-code measurement-schedule abstraction. */
#include <irrep/floquet_code.h>
#include <irrep/types.h>

#include <stdlib.h>
#include <string.h>

irrep_status_t
irrep_floquet_code_new(irrep_floquet_code_t *out,
                       int n_qubits, int n_rounds)
{
    if (out == NULL || n_qubits <= 0 || n_rounds <= 0) return IRREP_ERR_INVALID_ARG;
    out->n_qubits = n_qubits;
    out->n_rounds = n_rounds;
    out->rounds = (irrep_floquet_round_t *)calloc(n_rounds, sizeof(irrep_floquet_round_t));
    if (out->rounds == NULL) return IRREP_ERR_OUT_OF_MEMORY;
    for (int t = 0; t < n_rounds; ++t) {
        out->rounds[t].n_qubits = n_qubits;
        out->rounds[t].n_meas = 0;
        out->rounds[t].meas = NULL;
    }
    return IRREP_OK;
}

void
irrep_floquet_code_free(irrep_floquet_code_t *f)
{
    if (f == NULL) return;
    if (f->rounds) {
        for (int t = 0; t < f->n_rounds; ++t) {
            for (int i = 0; i < f->rounds[t].n_meas; ++i) {
                irrep_pauli_free(&f->rounds[t].meas[i]);
            }
            free(f->rounds[t].meas);
        }
        free(f->rounds);
        f->rounds = NULL;
    }
    f->n_qubits = 0;
    f->n_rounds = 0;
}

irrep_status_t
irrep_floquet_round_alloc(irrep_floquet_code_t *f,
                          int round_index, int n_meas)
{
    if (f == NULL) return IRREP_ERR_INVALID_ARG;
    if (round_index < 0 || round_index >= f->n_rounds) return IRREP_ERR_INVALID_ARG;
    if (n_meas < 0) return IRREP_ERR_INVALID_ARG;
    irrep_floquet_round_t *r = &f->rounds[round_index];
    /* Free anything already there. */
    for (int i = 0; i < r->n_meas; ++i) irrep_pauli_free(&r->meas[i]);
    free(r->meas);
    r->n_meas = n_meas;
    if (n_meas == 0) {
        r->meas = NULL;
        return IRREP_OK;
    }
    r->meas = (irrep_pauli_t *)calloc(n_meas, sizeof(irrep_pauli_t));
    if (r->meas == NULL) {
        r->n_meas = 0;
        return IRREP_ERR_OUT_OF_MEMORY;
    }
    for (int i = 0; i < n_meas; ++i) {
        irrep_status_t s = irrep_pauli_new(&r->meas[i], f->n_qubits);
        if (s != IRREP_OK) {
            for (int j = 0; j < i; ++j) irrep_pauli_free(&r->meas[j]);
            free(r->meas);
            r->meas = NULL;
            r->n_meas = 0;
            return s;
        }
    }
    return IRREP_OK;
}

irrep_status_t
irrep_floquet_round_check_commutativity(const irrep_floquet_code_t *f,
                                        int round_index)
{
    if (f == NULL) return IRREP_ERR_INVALID_ARG;
    if (round_index < 0 || round_index >= f->n_rounds) return IRREP_ERR_INVALID_ARG;
    const irrep_floquet_round_t *r = &f->rounds[round_index];
    for (int i = 0; i < r->n_meas; ++i) {
        for (int j = i + 1; j < r->n_meas; ++j) {
            if (!irrep_pauli_commute(&r->meas[i], &r->meas[j])) {
                return IRREP_ERR_PRECONDITION;
            }
        }
    }
    return IRREP_OK;
}

irrep_status_t
irrep_floquet_code_check(const irrep_floquet_code_t *f)
{
    if (f == NULL) return IRREP_ERR_INVALID_ARG;
    for (int t = 0; t < f->n_rounds; ++t) {
        irrep_status_t s = irrep_floquet_round_check_commutativity(f, t);
        if (s != IRREP_OK) return s;
    }
    return IRREP_OK;
}

irrep_status_t
irrep_floquet_round_to_stabilizer_group(const irrep_floquet_code_t *f,
                                        int round_index,
                                        irrep_stabilizer_group_t *out)
{
    if (f == NULL || out == NULL) return IRREP_ERR_INVALID_ARG;
    if (round_index < 0 || round_index >= f->n_rounds) return IRREP_ERR_INVALID_ARG;
    const irrep_floquet_round_t *r = &f->rounds[round_index];
    irrep_status_t s = irrep_stabilizer_group_new(out, f->n_qubits, r->n_meas);
    if (s != IRREP_OK) return s;
    /* Deep-copy the bit-vectors. */
    for (int i = 0; i < r->n_meas; ++i) {
        memcpy(out->gens[i].x, r->meas[i].x,
               (size_t)r->meas[i].n_words * sizeof(uint64_t));
        memcpy(out->gens[i].z, r->meas[i].z,
               (size_t)r->meas[i].n_words * sizeof(uint64_t));
        out->gens[i].sign = r->meas[i].sign;
    }
    return IRREP_OK;
}
