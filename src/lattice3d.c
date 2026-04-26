/* SPDX-License-Identifier: MIT */
/* 3D Bravais-lattice primitives with periodic boundary conditions.
 *
 * Conventions (Ashcroft-Mermin Chapter 4 / Kittel Chapter 1):
 *   simple cubic    a₁=(1,0,0),  a₂=(0,1,0),  a₃=(0,0,1)
 *   sublattices, in units of `a = 1`:
 *     SC      A=(0,0,0)
 *     BCC     A=(0,0,0), B=(½,½,½)
 *     FCC     A=(0,0,0), B=(0,½,½), C=(½,0,½), D=(½,½,0)
 *     Diamond FCC + the same FCC shifted by (¼,¼,¼); 8 sublattices.
 *
 * Bond candidates are auto-discovered: scan every (src_sub, dst_sub) pair
 * with cell offsets in [-1, +1]³ (sufficient — NN distance is below the
 * conventional cell edge for every supported family) and keep those whose
 * displacement matches the lattice's NN or NNN distance to within tol.
 * Then canonicalise (i < j), sort, dedup — the same shape as the 2D
 * pipeline. */

#include <math.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

#include <irrep/lattice3d.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

extern void irrep_set_error_(const char *fmt, ...);

#define MAX_SUBLATTICES 16
#define MAX_BOND_DIRS   512

struct irrep_lattice3d {
    irrep_lattice3d_kind_t kind;
    int                    Lx, Ly, Lz;
    int                    sites_per_cell;
    int                    num_cells;
    int                    num_sites;

    double a1[3], a2[3], a3[3];
    double b1[3], b2[3], b3[3];

    /* sublattice_xyz[s*3 + 0/1/2]: cartesian offset of sublattice s within a cell */
    double sublattice_xyz[MAX_SUBLATTICES * 3];

    double nn_distance;
    double nnn_distance;

    int *nn_i;
    int *nn_j;
    int  nn_count;

    int *nnn_i;
    int *nnn_j;
    int  nnn_count;
};

typedef struct {
    int src_sub;
    int dst_sub;
    int dx, dy, dz;
} bond_dir3_t;

typedef struct {
    int i, j;
} pair_t;

static int pair_cmp_(const void *a, const void *b) {
    const pair_t *pa = a;
    const pair_t *pb = b;
    if (pa->i != pb->i)
        return (pa->i < pb->i) ? -1 : 1;
    if (pa->j != pb->j)
        return (pa->j < pb->j) ? -1 : 1;
    return 0;
}

/* -------------------------------------------------------------------------- *
 * Lattice metadata                                                           *
 * -------------------------------------------------------------------------- */

static int sites_per_cell_(irrep_lattice3d_kind_t k) {
    switch (k) {
    case IRREP_LATTICE3D_SC:
        return 1;
    case IRREP_LATTICE3D_BCC:
        return 2;
    case IRREP_LATTICE3D_FCC:
        return 4;
    case IRREP_LATTICE3D_DIAMOND:
        return 8;
    case IRREP_LATTICE3D_PYROCHLORE:
        return 16;
    }
    return 0;
}

static void sublattice_offsets_(irrep_lattice3d_kind_t k, double sub[][3], int spc) {
    /* All entries pre-zeroed by caller. */
    switch (k) {
    case IRREP_LATTICE3D_SC:
        return; /* sub[0] = (0,0,0) already */
    case IRREP_LATTICE3D_BCC:
        sub[1][0] = 0.5;
        sub[1][1] = 0.5;
        sub[1][2] = 0.5;
        return;
    case IRREP_LATTICE3D_FCC:
        sub[1][0] = 0.0;
        sub[1][1] = 0.5;
        sub[1][2] = 0.5;
        sub[2][0] = 0.5;
        sub[2][1] = 0.0;
        sub[2][2] = 0.5;
        sub[3][0] = 0.5;
        sub[3][1] = 0.5;
        sub[3][2] = 0.0;
        return;
    case IRREP_LATTICE3D_DIAMOND: {
        /* FCC subs 0..3 plus FCC + (¼,¼,¼) at subs 4..7. */
        double fcc[4][3] = {
            {0.0, 0.0, 0.0}, {0.0, 0.5, 0.5}, {0.5, 0.0, 0.5}, {0.5, 0.5, 0.0}};
        for (int s = 0; s < 4; ++s) {
            sub[s][0] = fcc[s][0];
            sub[s][1] = fcc[s][1];
            sub[s][2] = fcc[s][2];
            sub[s + 4][0] = fcc[s][0] + 0.25;
            sub[s + 4][1] = fcc[s][1] + 0.25;
            sub[s + 4][2] = fcc[s][2] + 0.25;
        }
        return;
    }
    case IRREP_LATTICE3D_PYROCHLORE: {
        /* Each FCC sublattice carries a 4-site tetrahedral basis at
         * (0,0,0), (¼,¼,0), (¼,0,¼), (0,¼,¼). The 16-site basis is
         * the union over the 4 FCC sites in the conventional cubic cell. */
        double fcc[4][3] = {
            {0.0, 0.0, 0.0}, {0.0, 0.5, 0.5}, {0.5, 0.0, 0.5}, {0.5, 0.5, 0.0}};
        double tet[4][3] = {
            {0.0, 0.0, 0.0}, {0.25, 0.25, 0.0}, {0.25, 0.0, 0.25}, {0.0, 0.25, 0.25}};
        for (int f = 0; f < 4; ++f)
            for (int t = 0; t < 4; ++t) {
                int s = f * 4 + t;
                sub[s][0] = fcc[f][0] + tet[t][0];
                sub[s][1] = fcc[f][1] + tet[t][1];
                sub[s][2] = fcc[f][2] + tet[t][2];
            }
        return;
    }
    }
    (void)spc;
}

/* NN / NNN reference distances per lattice family. The values below match
 * the geometric truth for `a = 1` and are used to identify the bond shells
 * with a small tolerance. */
static void neighbour_distances_(irrep_lattice3d_kind_t k, double *nn, double *nnn) {
    switch (k) {
    case IRREP_LATTICE3D_SC:
        *nn = 1.0;
        *nnn = sqrt(2.0);
        return;
    case IRREP_LATTICE3D_BCC:
        *nn = 0.5 * sqrt(3.0);
        *nnn = 1.0;
        return;
    case IRREP_LATTICE3D_FCC:
        *nn = 0.5 * sqrt(2.0);
        *nnn = 1.0;
        return;
    case IRREP_LATTICE3D_DIAMOND:
        *nn = 0.25 * sqrt(3.0);
        *nnn = 0.5 * sqrt(2.0);
        return;
    case IRREP_LATTICE3D_PYROCHLORE:
        *nn = 0.25 * sqrt(2.0);  /* edge of the tetrahedral basis */
        *nnn = 0.25 * sqrt(6.0); /* next shell — verified by symmetry */
        return;
    }
    *nn = 0;
    *nnn = 0;
}

/* -------------------------------------------------------------------------- *
 * Reciprocal lattice                                                         *
 * -------------------------------------------------------------------------- */

static void reciprocal_(const double a1[3], const double a2[3], const double a3[3], double b1[3],
                        double b2[3], double b3[3]) {
    /* Volume of the primitive cell. */
    double v = a1[0] * (a2[1] * a3[2] - a2[2] * a3[1]) -
               a1[1] * (a2[0] * a3[2] - a2[2] * a3[0]) +
               a1[2] * (a2[0] * a3[1] - a2[1] * a3[0]);
    double s = (2.0 * M_PI) / v;

    /* b_i = (2π / V) · a_j × a_k cyclic. */
    b1[0] = s * (a2[1] * a3[2] - a2[2] * a3[1]);
    b1[1] = s * (a2[2] * a3[0] - a2[0] * a3[2]);
    b1[2] = s * (a2[0] * a3[1] - a2[1] * a3[0]);

    b2[0] = s * (a3[1] * a1[2] - a3[2] * a1[1]);
    b2[1] = s * (a3[2] * a1[0] - a3[0] * a1[2]);
    b2[2] = s * (a3[0] * a1[1] - a3[1] * a1[0]);

    b3[0] = s * (a1[1] * a2[2] - a1[2] * a2[1]);
    b3[1] = s * (a1[2] * a2[0] - a1[0] * a2[2]);
    b3[2] = s * (a1[0] * a2[1] - a1[1] * a2[0]);
}

/* -------------------------------------------------------------------------- *
 * Bond candidate auto-discovery                                              *
 * -------------------------------------------------------------------------- */

/* Walk all (src, dst, dx, dy, dz) triples, retain those at the requested
 * distance, and keep only one orientation per undirected bond.
 *
 * Canonicalisation: an undirected bond {(s, c), (t, c+(dx,dy,dz))} is
 * emitted only when (s, dx, dy, dz) < (t, -dx, -dy, -dz) lexicographically,
 * with `s` as the leading key. Same-sublattice bonds use (dx, dy, dz)
 * positive in the leading nonzero coordinate. */
static int discover_dirs_(const double sub_xyz[][3], int spc, double target_dist, double tol,
                          bond_dir3_t *out, int max_out) {
    int n = 0;
    for (int s = 0; s < spc; ++s) {
        for (int t = 0; t < spc; ++t) {
            for (int dx = -1; dx <= 1; ++dx)
                for (int dy = -1; dy <= 1; ++dy)
                    for (int dz = -1; dz <= 1; ++dz) {
                        if (s == t && dx == 0 && dy == 0 && dz == 0)
                            continue;

                        double rx = sub_xyz[t][0] + dx - sub_xyz[s][0];
                        double ry = sub_xyz[t][1] + dy - sub_xyz[s][1];
                        double rz = sub_xyz[t][2] + dz - sub_xyz[s][2];
                        double dist = sqrt(rx * rx + ry * ry + rz * rz);
                        if (fabs(dist - target_dist) > tol)
                            continue;

                        /* canonical orientation */
                        bool keep;
                        if (s != t) {
                            keep = (s < t);
                        } else {
                            /* same sublattice — pick one of {(dx,dy,dz), (-dx,-dy,-dz)}. */
                            if (dx != 0)
                                keep = (dx > 0);
                            else if (dy != 0)
                                keep = (dy > 0);
                            else
                                keep = (dz > 0);
                        }
                        if (!keep)
                            continue;

                        if (n >= max_out)
                            return -1;
                        out[n].src_sub = s;
                        out[n].dst_sub = t;
                        out[n].dx = dx;
                        out[n].dy = dy;
                        out[n].dz = dz;
                        ++n;
                    }
        }
    }
    return n;
}

static int mod_(int x, int m) {
    int r = x % m;
    if (r < 0)
        r += m;
    return r;
}

static int pack_site_(const struct irrep_lattice3d *L, int ix, int iy, int iz, int sub) {
    ix = mod_(ix, L->Lx);
    iy = mod_(iy, L->Ly);
    iz = mod_(iz, L->Lz);
    int cell = (iz * L->Ly + iy) * L->Lx + ix;
    return cell * L->sites_per_cell + sub;
}

static int build_bonds_(const struct irrep_lattice3d *L, const bond_dir3_t *dirs, int n_dirs,
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
    for (int iz = 0; iz < L->Lz; ++iz)
        for (int iy = 0; iy < L->Ly; ++iy)
            for (int ix = 0; ix < L->Lx; ++ix)
                for (int d = 0; d < n_dirs; ++d) {
                    int i = pack_site_(L, ix, iy, iz, dirs[d].src_sub);
                    int j = pack_site_(L, ix + dirs[d].dx, iy + dirs[d].dy, iz + dirs[d].dz,
                                       dirs[d].dst_sub);
                    if (i == j)
                        continue;
                    if (i > j) {
                        int t = i;
                        i = j;
                        j = t;
                    }
                    buf[k].i = i;
                    buf[k].j = j;
                    ++k;
                }

    qsort(buf, (size_t)k, sizeof(*buf), pair_cmp_);

    int w = 0;
    for (int r = 0; r < k; ++r) {
        if (w == 0 || buf[r].i != buf[w - 1].i || buf[r].j != buf[w - 1].j) {
            buf[w++] = buf[r];
        }
    }

    int *ii = malloc((size_t)(w > 0 ? w : 1) * sizeof(*ii));
    int *jj = malloc((size_t)(w > 0 ? w : 1) * sizeof(*jj));
    if (!ii || !jj) {
        free(ii);
        free(jj);
        free(buf);
        *out_i = NULL;
        *out_j = NULL;
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

irrep_lattice3d_t *irrep_lattice3d_build(irrep_lattice3d_kind_t kind, int Lx, int Ly, int Lz) {
    if (Lx < 1 || Ly < 1 || Lz < 1) {
        irrep_set_error_("irrep_lattice3d_build: Lx, Ly, Lz must be >= 1 (got %d, %d, %d)", Lx,
                         Ly, Lz);
        return NULL;
    }
    /* L = 1 along any axis collapses every cell offset to the same image,
     * which the bond builder handles via self-bond dedup. For SC this
     * yields 0 bonds at 1³ (every site is its own NN under PBC); for
     * lattices with a non-trivial basis (BCC, FCC, Diamond, Pyrochlore)
     * the intra-cell NN structure survives and 1×1×1 is meaningful. */
    int spc = sites_per_cell_(kind);
    if (spc <= 0) {
        irrep_set_error_("irrep_lattice3d_build: unknown lattice kind %d", (int)kind);
        return NULL;
    }

    struct irrep_lattice3d *L = calloc(1, sizeof(*L));
    if (!L) {
        irrep_set_error_("irrep_lattice3d_build: out of memory");
        return NULL;
    }
    L->kind = kind;
    L->Lx = Lx;
    L->Ly = Ly;
    L->Lz = Lz;
    L->sites_per_cell = spc;
    L->num_cells = Lx * Ly * Lz;
    L->num_sites = L->num_cells * spc;

    L->a1[0] = 1.0;
    L->a1[1] = 0.0;
    L->a1[2] = 0.0;
    L->a2[0] = 0.0;
    L->a2[1] = 1.0;
    L->a2[2] = 0.0;
    L->a3[0] = 0.0;
    L->a3[1] = 0.0;
    L->a3[2] = 1.0;
    reciprocal_(L->a1, L->a2, L->a3, L->b1, L->b2, L->b3);

    double sub[MAX_SUBLATTICES][3] = {{0}};
    sublattice_offsets_(kind, sub, spc);
    for (int s = 0; s < spc; ++s) {
        L->sublattice_xyz[3 * s + 0] = sub[s][0];
        L->sublattice_xyz[3 * s + 1] = sub[s][1];
        L->sublattice_xyz[3 * s + 2] = sub[s][2];
    }

    neighbour_distances_(kind, &L->nn_distance, &L->nnn_distance);

    bond_dir3_t nn_dirs[MAX_BOND_DIRS], nnn_dirs[MAX_BOND_DIRS];
    int         n_nn = discover_dirs_(sub, spc, L->nn_distance, 1e-9, nn_dirs, MAX_BOND_DIRS);
    int         n_nnn = discover_dirs_(sub, spc, L->nnn_distance, 1e-9, nnn_dirs, MAX_BOND_DIRS);
    if (n_nn < 0 || n_nnn < 0) {
        free(L);
        irrep_set_error_("irrep_lattice3d_build: bond direction table overflow");
        return NULL;
    }

    int rc1 = build_bonds_(L, nn_dirs, n_nn, &L->nn_i, &L->nn_j);
    int rc2 = build_bonds_(L, nnn_dirs, n_nnn, &L->nnn_i, &L->nnn_j);
    if (rc1 < 0 || rc2 < 0) {
        irrep_lattice3d_free(L);
        irrep_set_error_("irrep_lattice3d_build: bond construction out of memory");
        return NULL;
    }
    L->nn_count = rc1;
    L->nnn_count = rc2;

    return L;
}

void irrep_lattice3d_free(irrep_lattice3d_t *L) {
    if (!L)
        return;
    free(L->nn_i);
    free(L->nn_j);
    free(L->nnn_i);
    free(L->nnn_j);
    free(L);
}

int irrep_lattice3d_num_sites(const irrep_lattice3d_t *L) {
    return L ? L->num_sites : 0;
}
int irrep_lattice3d_num_cells(const irrep_lattice3d_t *L) {
    return L ? L->num_cells : 0;
}
int irrep_lattice3d_sites_per_cell(const irrep_lattice3d_t *L) {
    return L ? L->sites_per_cell : 0;
}
int irrep_lattice3d_Lx(const irrep_lattice3d_t *L) {
    return L ? L->Lx : 0;
}
int irrep_lattice3d_Ly(const irrep_lattice3d_t *L) {
    return L ? L->Ly : 0;
}
int irrep_lattice3d_Lz(const irrep_lattice3d_t *L) {
    return L ? L->Lz : 0;
}
irrep_lattice3d_kind_t irrep_lattice3d_kind(const irrep_lattice3d_t *L) {
    return L ? L->kind : IRREP_LATTICE3D_SC;
}
double irrep_lattice3d_nn_distance(const irrep_lattice3d_t *L) {
    return L ? L->nn_distance : 0.0;
}
double irrep_lattice3d_nnn_distance(const irrep_lattice3d_t *L) {
    return L ? L->nnn_distance : 0.0;
}

void irrep_lattice3d_primitive_vectors(const irrep_lattice3d_t *L, double a1[3], double a2[3],
                                       double a3[3]) {
    if (!L)
        return;
    memcpy(a1, L->a1, sizeof L->a1);
    memcpy(a2, L->a2, sizeof L->a2);
    memcpy(a3, L->a3, sizeof L->a3);
}

void irrep_lattice3d_reciprocal_vectors(const irrep_lattice3d_t *L, double b1[3], double b2[3],
                                        double b3[3]) {
    if (!L)
        return;
    memcpy(b1, L->b1, sizeof L->b1);
    memcpy(b2, L->b2, sizeof L->b2);
    memcpy(b3, L->b3, sizeof L->b3);
}

irrep_status_t irrep_lattice3d_site_position(const irrep_lattice3d_t *L, int site,
                                             double xyz[3]) {
    if (!L || site < 0 || site >= L->num_sites)
        return IRREP_ERR_INVALID_ARG;
    int sub = site % L->sites_per_cell;
    int cell = site / L->sites_per_cell;
    int ix = cell % L->Lx;
    int iy = (cell / L->Lx) % L->Ly;
    int iz = cell / (L->Lx * L->Ly);
    xyz[0] = ix * L->a1[0] + iy * L->a2[0] + iz * L->a3[0] + L->sublattice_xyz[3 * sub + 0];
    xyz[1] = ix * L->a1[1] + iy * L->a2[1] + iz * L->a3[1] + L->sublattice_xyz[3 * sub + 1];
    xyz[2] = ix * L->a1[2] + iy * L->a2[2] + iz * L->a3[2] + L->sublattice_xyz[3 * sub + 2];
    return IRREP_OK;
}

int irrep_lattice3d_sublattice_of(const irrep_lattice3d_t *L, int site) {
    if (!L || site < 0 || site >= L->num_sites)
        return -1;
    return site % L->sites_per_cell;
}

irrep_status_t irrep_lattice3d_cell_of(const irrep_lattice3d_t *L, int site, int *ix, int *iy,
                                       int *iz) {
    if (!L || site < 0 || site >= L->num_sites || !ix || !iy || !iz)
        return IRREP_ERR_INVALID_ARG;
    int cell = site / L->sites_per_cell;
    *ix = cell % L->Lx;
    *iy = (cell / L->Lx) % L->Ly;
    *iz = cell / (L->Lx * L->Ly);
    return IRREP_OK;
}

int irrep_lattice3d_site_index(const irrep_lattice3d_t *L, int ix, int iy, int iz, int sub) {
    if (!L)
        return -1;
    if (sub < 0 || sub >= L->sites_per_cell)
        return -1;
    return pack_site_(L, ix, iy, iz, sub);
}

int irrep_lattice3d_translate(const irrep_lattice3d_t *L, int site, int dx, int dy, int dz) {
    if (!L || site < 0 || site >= L->num_sites)
        return -1;
    int sub = site % L->sites_per_cell;
    int cell = site / L->sites_per_cell;
    int ix = cell % L->Lx;
    int iy = (cell / L->Lx) % L->Ly;
    int iz = cell / (L->Lx * L->Ly);
    return pack_site_(L, ix + dx, iy + dy, iz + dz, sub);
}

int irrep_lattice3d_num_bonds_nn(const irrep_lattice3d_t *L) {
    return L ? L->nn_count : 0;
}
int irrep_lattice3d_num_bonds_nnn(const irrep_lattice3d_t *L) {
    return L ? L->nnn_count : 0;
}

void irrep_lattice3d_fill_bonds_nn(const irrep_lattice3d_t *L, int *i_out, int *j_out) {
    if (!L)
        return;
    if (i_out)
        memcpy(i_out, L->nn_i, (size_t)L->nn_count * sizeof(int));
    if (j_out)
        memcpy(j_out, L->nn_j, (size_t)L->nn_count * sizeof(int));
}

void irrep_lattice3d_fill_bonds_nnn(const irrep_lattice3d_t *L, int *i_out, int *j_out) {
    if (!L)
        return;
    if (i_out)
        memcpy(i_out, L->nnn_i, (size_t)L->nnn_count * sizeof(int));
    if (j_out)
        memcpy(j_out, L->nnn_j, (size_t)L->nnn_count * sizeof(int));
}

void irrep_lattice3d_k_grid(const irrep_lattice3d_t *L, double *kx, double *ky, double *kz) {
    if (!L)
        return;
    for (int n3 = 0; n3 < L->Lz; ++n3)
        for (int n2 = 0; n2 < L->Ly; ++n2)
            for (int n1 = 0; n1 < L->Lx; ++n1) {
                double f1 = (double)n1 / (double)L->Lx;
                double f2 = (double)n2 / (double)L->Ly;
                double f3 = (double)n3 / (double)L->Lz;
                int    idx = (n3 * L->Ly + n2) * L->Lx + n1;
                if (kx)
                    kx[idx] = f1 * L->b1[0] + f2 * L->b2[0] + f3 * L->b3[0];
                if (ky)
                    ky[idx] = f1 * L->b1[1] + f2 * L->b2[1] + f3 * L->b3[1];
                if (kz)
                    kz[idx] = f1 * L->b1[2] + f2 * L->b2[2] + f3 * L->b3[2];
            }
}
