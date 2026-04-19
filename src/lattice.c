/* SPDX-License-Identifier: MIT */
/* 2D Bravais-lattice primitives with periodic boundary conditions.
 *
 * Bonds are enumerated from canonical direction sets per lattice (chosen so
 * every NN or NNN pair is visited once per unit cell), canonicalised to
 * (min, max), sorted lexicographically, then deduplicated. This leaves the
 * bond lists insensitive to the lattice anisotropy that PBC wrap introduces
 * on small (Lx or Ly) clusters, including the degenerate case where the two
 * endpoints of a candidate bond happen to coincide through the wrap — such
 * self-bonds are dropped in the deduplication pass.
 *
 * Conventions (cf. Ashcroft-Mermin Chapter 7; Elser 1989 for kagome):
 *   square      a1=(1,0),       a2=(0,1);         A=(0,0)
 *   triangular  a1=(1,0),       a2=(1/2, √3/2);   A=(0,0)
 *   honeycomb   a1=(3/2,√3/2),  a2=(3/2,-√3/2);   A=(0,0), B=(1,0)
 *   kagome      a1=(2,0),       a2=(1,√3);        A=(0,0), B=(1,0), C=(1/2,√3/2)
 *
 * Chosen so every NN bond has length 1 in every lattice, so Heisenberg J
 * and Hubbard t transfer without a rescale.
 */

#include <math.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

#include <irrep/lattice.h>

#ifndef M_PI
#  define M_PI 3.14159265358979323846
#endif

extern void irrep_set_error_(const char *fmt, ...);

/* -------------------------------------------------------------------------- *
 * Internal lattice representation                                            *
 * -------------------------------------------------------------------------- */

struct irrep_lattice {
    irrep_lattice_kind_t kind;
    int Lx, Ly;
    int sites_per_cell;
    int num_cells;
    int num_sites;

    double a1[2], a2[2];
    double b1[2], b2[2];

    /* sublattice_xy[s*2 + 0/1]: cartesian offset of sublattice s within a cell */
    double sublattice_xy[6];

    /* cached bond lists */
    int *nn_i;
    int *nn_j;
    int  nn_count;

    int *nnn_i;
    int *nnn_j;
    int  nnn_count;
};

/* -------------------------------------------------------------------------- *
 * Bond candidate description: one per canonical direction.                   *
 * Each candidate says: "from site (ix, iy, src_sub), walk (dx, dy) unit cells
 * and land on sublattice dst_sub".                                           *
 * -------------------------------------------------------------------------- */

typedef struct {
    int src_sub;
    int dst_sub;
    int dx;
    int dy;
} bond_dir_t;

/* Directions are listed so each undirected NN / NNN bond is visited exactly
 * once per origin cell in an infinite lattice. PBC makes this untrue on
 * small clusters, which is why the builder canonicalises + dedups. */

static const bond_dir_t NN_SQUARE[] = {
    { 0, 0, +1,  0 },
    { 0, 0,  0, +1 }
};

static const bond_dir_t NNN_SQUARE[] = {
    { 0, 0, +1, +1 },
    { 0, 0, +1, -1 }
};

static const bond_dir_t NN_TRIANGULAR[] = {
    { 0, 0, +1,  0 },
    { 0, 0,  0, +1 },
    { 0, 0, +1, -1 }   /* -a1 + a2 reversed sign — or equivalently */
};

static const bond_dir_t NNN_TRIANGULAR[] = {
    { 0, 0, +1, +1 },
    { 0, 0, +2, -1 },
    { 0, 0, -1, +2 }
};

static const bond_dir_t NN_HONEYCOMB[] = {
    { 0, 1,  0,  0 },  /* A(ix,iy) - B(ix,iy) */
    { 0, 1, -1,  0 },  /* A(ix,iy) - B(ix-1,iy) */
    { 0, 1,  0, -1 }   /* A(ix,iy) - B(ix,iy-1) */
};

static const bond_dir_t NNN_HONEYCOMB[] = {
    /* same sublattice; hex ring of 6 per site — emit 3 per unit cell (for each
     * of the 2 sublattices, making 6 per unit cell ≡ 3/site net). */
    { 0, 0, +1,  0 },
    { 0, 0,  0, +1 },
    { 0, 0, +1, -1 },
    { 1, 1, +1,  0 },
    { 1, 1,  0, +1 },
    { 1, 1, +1, -1 }
};

static const bond_dir_t NN_KAGOME[] = {
    /* within cell */
    { 0, 1,  0,  0 },   /* A - B */
    { 0, 2,  0,  0 },   /* A - C */
    { 1, 2,  0,  0 },   /* B - C */
    /* between cells */
    { 1, 0, +1,  0 },   /* B(ix,iy) - A(ix+1,iy) */
    { 2, 0,  0, +1 },   /* C(ix,iy) - A(ix,iy+1) */
    { 2, 1, -1, +1 }    /* C(ix,iy) - B(ix-1,iy+1) */
};

static const bond_dir_t NNN_KAGOME[] = {
    /* Kagome NNN are the across-hexagon, mixed-sublattice bonds at distance
     * √3 (not the same-sublattice distance-2 bonds — those would be further
     * neighbours). Each site has 4 NNN; per unit cell, 4·3 / 2 = 6 directed
     * candidates, chosen so no bond is visited twice in the infinite lattice. */
    { 0, 1,  0, -1 },   /* A(ix,iy) - B(ix,iy-1) */
    { 0, 1, -1, +1 },   /* A(ix,iy) - B(ix-1,iy+1) */
    { 0, 2, -1,  0 },   /* A(ix,iy) - C(ix-1,iy) */
    { 0, 2, +1, -1 },   /* A(ix,iy) - C(ix+1,iy-1) */
    { 1, 2, +1,  0 },   /* B(ix,iy) - C(ix+1,iy) */
    { 1, 2,  0, -1 }    /* B(ix,iy) - C(ix,iy-1) */
};

/* -------------------------------------------------------------------------- *
 * Lattice metadata lookup                                                    *
 * -------------------------------------------------------------------------- */

static int lattice_sites_per_cell_(irrep_lattice_kind_t k) {
    switch (k) {
        case IRREP_LATTICE_SQUARE:      return 1;
        case IRREP_LATTICE_TRIANGULAR:  return 1;
        case IRREP_LATTICE_HONEYCOMB:   return 2;
        case IRREP_LATTICE_KAGOME:      return 3;
    }
    return 0;
}

static void lattice_primitive_(irrep_lattice_kind_t k, double a1[2], double a2[2]) {
    switch (k) {
        case IRREP_LATTICE_SQUARE:
            a1[0] = 1.0; a1[1] = 0.0;
            a2[0] = 0.0; a2[1] = 1.0;
            return;
        case IRREP_LATTICE_TRIANGULAR:
            a1[0] = 1.0; a1[1] = 0.0;
            a2[0] = 0.5; a2[1] = 0.5 * sqrt(3.0);
            return;
        case IRREP_LATTICE_HONEYCOMB:
            a1[0] = 1.5; a1[1] =  0.5 * sqrt(3.0);
            a2[0] = 1.5; a2[1] = -0.5 * sqrt(3.0);
            return;
        case IRREP_LATTICE_KAGOME:
            a1[0] = 2.0; a1[1] = 0.0;
            a2[0] = 1.0; a2[1] = sqrt(3.0);
            return;
    }
    a1[0] = a1[1] = a2[0] = a2[1] = 0.0;
}

static void lattice_sublattice_(irrep_lattice_kind_t k, double sub[][2]) {
    switch (k) {
        case IRREP_LATTICE_SQUARE:
        case IRREP_LATTICE_TRIANGULAR:
            sub[0][0] = 0.0; sub[0][1] = 0.0;
            return;
        case IRREP_LATTICE_HONEYCOMB:
            sub[0][0] = 0.0; sub[0][1] = 0.0;
            sub[1][0] = 1.0; sub[1][1] = 0.0;
            return;
        case IRREP_LATTICE_KAGOME:
            sub[0][0] = 0.0; sub[0][1] = 0.0;
            sub[1][0] = 1.0; sub[1][1] = 0.0;
            sub[2][0] = 0.5; sub[2][1] = 0.5 * sqrt(3.0);
            return;
    }
}

/* Reciprocal of a 2D primitive basis: a1·b1 = 2π, a1·b2 = 0, etc. Closed form
 * via the inverse of [[a1x, a2x],[a1y, a2y]] scaled by 2π. */
static void reciprocal_(const double a1[2], const double a2[2],
                        double b1[2], double b2[2]) {
    double det = a1[0] * a2[1] - a1[1] * a2[0];
    double s = (2.0 * M_PI) / det;
    b1[0] =  s * a2[1];
    b1[1] = -s * a2[0];
    b2[0] = -s * a1[1];
    b2[1] =  s * a1[0];
}

/* Non-negative modulo */
static int mod_(int x, int m) {
    int r = x % m;
    if (r < 0) r += m;
    return r;
}

static int pack_site_(const struct irrep_lattice *L, int ix, int iy, int sub) {
    ix = mod_(ix, L->Lx);
    iy = mod_(iy, L->Ly);
    return (iy * L->Lx + ix) * L->sites_per_cell + sub;
}

/* -------------------------------------------------------------------------- *
 * Bond list construction (canonicalise, sort, dedup)                         *
 * -------------------------------------------------------------------------- */

typedef struct {
    int i;
    int j;
} pair_t;

static int pair_cmp_(const void *a, const void *b) {
    const pair_t *pa = a;
    const pair_t *pb = b;
    if (pa->i != pb->i) return (pa->i < pb->i) ? -1 : 1;
    if (pa->j != pb->j) return (pa->j < pb->j) ? -1 : 1;
    return 0;
}

static int build_bonds_(const struct irrep_lattice *L,
                        const bond_dir_t *dirs, int n_dirs,
                        int **out_i, int **out_j) {
    int capacity = L->num_cells * n_dirs;
    if (capacity == 0) {
        *out_i = NULL;
        *out_j = NULL;
        return 0;
    }
    pair_t *buf = malloc((size_t)capacity * sizeof(*buf));
    if (!buf) {
        *out_i = NULL;
        *out_j = NULL;
        return -1;
    }

    int k = 0;
    for (int iy = 0; iy < L->Ly; ++iy) {
        for (int ix = 0; ix < L->Lx; ++ix) {
            for (int d = 0; d < n_dirs; ++d) {
                int i = pack_site_(L, ix,             iy,             dirs[d].src_sub);
                int j = pack_site_(L, ix + dirs[d].dx, iy + dirs[d].dy, dirs[d].dst_sub);
                if (i == j) continue;            /* self-bond from PBC wrap */
                if (i > j) { int t = i; i = j; j = t; }
                buf[k].i = i;
                buf[k].j = j;
                ++k;
            }
        }
    }

    qsort(buf, (size_t)k, sizeof(*buf), pair_cmp_);

    /* uniq in place */
    int w = 0;
    for (int r = 0; r < k; ++r) {
        if (w == 0 || buf[r].i != buf[w-1].i || buf[r].j != buf[w-1].j) {
            buf[w++] = buf[r];
        }
    }

    int *ii = malloc((size_t)(w > 0 ? w : 1) * sizeof(*ii));
    int *jj = malloc((size_t)(w > 0 ? w : 1) * sizeof(*jj));
    if (!ii || !jj) {
        free(ii); free(jj); free(buf);
        *out_i = NULL; *out_j = NULL;
        return -1;
    }
    for (int r = 0; r < w; ++r) {
        ii[r] = buf[r].i;
        jj[r] = buf[r].j;
    }
    free(buf);
    *out_i = ii;
    *out_j = jj;
    return w;
}

/* -------------------------------------------------------------------------- *
 * Public API                                                                 *
 * -------------------------------------------------------------------------- */

irrep_lattice_t *irrep_lattice_build(irrep_lattice_kind_t kind, int Lx, int Ly) {
    if (Lx < 2 || Ly < 2) {
        irrep_set_error_("irrep_lattice_build: Lx and Ly must be >= 2 (got %d, %d)",
                         Lx, Ly);
        return NULL;
    }
    int spc = lattice_sites_per_cell_(kind);
    if (spc <= 0) {
        irrep_set_error_("irrep_lattice_build: unknown lattice kind %d", (int)kind);
        return NULL;
    }

    struct irrep_lattice *L = calloc(1, sizeof(*L));
    if (!L) {
        irrep_set_error_("irrep_lattice_build: out of memory");
        return NULL;
    }
    L->kind           = kind;
    L->Lx             = Lx;
    L->Ly             = Ly;
    L->sites_per_cell = spc;
    L->num_cells      = Lx * Ly;
    L->num_sites      = L->num_cells * spc;

    lattice_primitive_(kind, L->a1, L->a2);
    reciprocal_(L->a1, L->a2, L->b1, L->b2);

    double sub[6][2] = {{0}};
    lattice_sublattice_(kind, sub);
    for (int s = 0; s < spc; ++s) {
        L->sublattice_xy[2*s + 0] = sub[s][0];
        L->sublattice_xy[2*s + 1] = sub[s][1];
    }

    const bond_dir_t *nn = NULL, *nnn = NULL;
    int n_nn = 0, n_nnn = 0;
    switch (kind) {
        case IRREP_LATTICE_SQUARE:
            nn  = NN_SQUARE;      n_nn  = sizeof(NN_SQUARE) / sizeof(NN_SQUARE[0]);
            nnn = NNN_SQUARE;     n_nnn = sizeof(NNN_SQUARE) / sizeof(NNN_SQUARE[0]);
            break;
        case IRREP_LATTICE_TRIANGULAR:
            nn  = NN_TRIANGULAR;  n_nn  = sizeof(NN_TRIANGULAR) / sizeof(NN_TRIANGULAR[0]);
            nnn = NNN_TRIANGULAR; n_nnn = sizeof(NNN_TRIANGULAR) / sizeof(NNN_TRIANGULAR[0]);
            break;
        case IRREP_LATTICE_HONEYCOMB:
            nn  = NN_HONEYCOMB;   n_nn  = sizeof(NN_HONEYCOMB) / sizeof(NN_HONEYCOMB[0]);
            nnn = NNN_HONEYCOMB;  n_nnn = sizeof(NNN_HONEYCOMB) / sizeof(NNN_HONEYCOMB[0]);
            break;
        case IRREP_LATTICE_KAGOME:
            nn  = NN_KAGOME;      n_nn  = sizeof(NN_KAGOME) / sizeof(NN_KAGOME[0]);
            nnn = NNN_KAGOME;     n_nnn = sizeof(NNN_KAGOME) / sizeof(NNN_KAGOME[0]);
            break;
    }

    int rc1 = build_bonds_(L, nn,  n_nn,  &L->nn_i,  &L->nn_j);
    int rc2 = build_bonds_(L, nnn, n_nnn, &L->nnn_i, &L->nnn_j);
    if (rc1 < 0 || rc2 < 0) {
        irrep_lattice_free(L);
        irrep_set_error_("irrep_lattice_build: bond construction out of memory");
        return NULL;
    }
    L->nn_count  = rc1;
    L->nnn_count = rc2;

    return L;
}

void irrep_lattice_free(irrep_lattice_t *L) {
    if (!L) return;
    free(L->nn_i);
    free(L->nn_j);
    free(L->nnn_i);
    free(L->nnn_j);
    free(L);
}

int irrep_lattice_num_sites     (const irrep_lattice_t *L) { return L ? L->num_sites      : 0; }
int irrep_lattice_num_cells     (const irrep_lattice_t *L) { return L ? L->num_cells      : 0; }
int irrep_lattice_sites_per_cell(const irrep_lattice_t *L) { return L ? L->sites_per_cell : 0; }
int irrep_lattice_Lx            (const irrep_lattice_t *L) { return L ? L->Lx             : 0; }
int irrep_lattice_Ly            (const irrep_lattice_t *L) { return L ? L->Ly             : 0; }
irrep_lattice_kind_t irrep_lattice_kind(const irrep_lattice_t *L) {
    return L ? L->kind : IRREP_LATTICE_SQUARE;
}

void irrep_lattice_primitive_vectors(const irrep_lattice_t *L,
                                     double a1[2], double a2[2]) {
    if (!L) return;
    a1[0] = L->a1[0]; a1[1] = L->a1[1];
    a2[0] = L->a2[0]; a2[1] = L->a2[1];
}

void irrep_lattice_reciprocal_vectors(const irrep_lattice_t *L,
                                      double b1[2], double b2[2]) {
    if (!L) return;
    b1[0] = L->b1[0]; b1[1] = L->b1[1];
    b2[0] = L->b2[0]; b2[1] = L->b2[1];
}

irrep_status_t irrep_lattice_site_position(const irrep_lattice_t *L,
                                           int site, double xy[2]) {
    if (!L || site < 0 || site >= L->num_sites) return IRREP_ERR_INVALID_ARG;
    int sub   = site % L->sites_per_cell;
    int cell  = site / L->sites_per_cell;
    int ix    = cell % L->Lx;
    int iy    = cell / L->Lx;
    xy[0] = ix * L->a1[0] + iy * L->a2[0] + L->sublattice_xy[2*sub + 0];
    xy[1] = ix * L->a1[1] + iy * L->a2[1] + L->sublattice_xy[2*sub + 1];
    return IRREP_OK;
}

int irrep_lattice_sublattice_of(const irrep_lattice_t *L, int site) {
    if (!L || site < 0 || site >= L->num_sites) return -1;
    return site % L->sites_per_cell;
}

irrep_status_t irrep_lattice_cell_of(const irrep_lattice_t *L, int site,
                                     int *ix, int *iy) {
    if (!L || site < 0 || site >= L->num_sites || !ix || !iy)
        return IRREP_ERR_INVALID_ARG;
    int cell = site / L->sites_per_cell;
    *ix = cell % L->Lx;
    *iy = cell / L->Lx;
    return IRREP_OK;
}

int irrep_lattice_site_index(const irrep_lattice_t *L, int ix, int iy, int sub) {
    if (!L) return -1;
    if (sub < 0 || sub >= L->sites_per_cell) return -1;
    return pack_site_(L, ix, iy, sub);
}

int irrep_lattice_translate(const irrep_lattice_t *L, int site, int dx, int dy) {
    if (!L || site < 0 || site >= L->num_sites) return -1;
    int sub  = site % L->sites_per_cell;
    int cell = site / L->sites_per_cell;
    int ix   = cell % L->Lx;
    int iy   = cell / L->Lx;
    return pack_site_(L, ix + dx, iy + dy, sub);
}

int irrep_lattice_num_bonds_nn (const irrep_lattice_t *L) { return L ? L->nn_count  : 0; }
int irrep_lattice_num_bonds_nnn(const irrep_lattice_t *L) { return L ? L->nnn_count : 0; }

void irrep_lattice_fill_bonds_nn(const irrep_lattice_t *L, int *i_out, int *j_out) {
    if (!L) return;
    if (i_out) memcpy(i_out, L->nn_i, (size_t)L->nn_count * sizeof(int));
    if (j_out) memcpy(j_out, L->nn_j, (size_t)L->nn_count * sizeof(int));
}

void irrep_lattice_fill_bonds_nnn(const irrep_lattice_t *L, int *i_out, int *j_out) {
    if (!L) return;
    if (i_out) memcpy(i_out, L->nnn_i, (size_t)L->nnn_count * sizeof(int));
    if (j_out) memcpy(j_out, L->nnn_j, (size_t)L->nnn_count * sizeof(int));
}

void irrep_lattice_k_grid(const irrep_lattice_t *L, double *kx, double *ky) {
    if (!L) return;
    for (int n2 = 0; n2 < L->Ly; ++n2) {
        for (int n1 = 0; n1 < L->Lx; ++n1) {
            double f1 = (double)n1 / (double)L->Lx;
            double f2 = (double)n2 / (double)L->Ly;
            int idx = n2 * L->Lx + n1;
            if (kx) kx[idx] = f1 * L->b1[0] + f2 * L->b2[0];
            if (ky) ky[idx] = f1 * L->b1[1] + f2 * L->b2[1];
        }
    }
}
