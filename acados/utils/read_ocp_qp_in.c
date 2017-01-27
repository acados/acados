/*
 *    This file is part of acados.
 *
 *    acados is free software; you can redistribute it and/or
 *    modify it under the terms of the GNU Lesser General Public
 *    License as published by the Free Software Foundation; either
 *    version 3 of the License, or (at your option) any later version.
 *
 *    acados is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *    Lesser General Public License for more details.
 *
 *    You should have received a copy of the GNU Lesser General Public
 *    License along with acados; if not, write to the Free Software Foundation,
 *    Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include "blasfeo/include/blasfeo_target.h"
#include "blasfeo/include/blasfeo_common.h"
#include "blasfeo/include/blasfeo_d_aux.h"
#include "blasfeo/include/blasfeo_i_aux.h"

#include "acados/utils/allocate_ocp_qp.h"
#include "acados/utils/read_ocp_qp_in.h"

// TODO(dimitris): IMPORTANT! Add asserts for all the data that are read in

static void transpose_matrix(real_t *mat, int m, int n) {
    int_t jj, ii;
    real_t *tmp;
    d_zeros(&tmp, m, n);

    for (jj = 0; jj < n; jj++) {
        for (ii = 0; ii < m; ii++) {
             tmp[ii*n+jj] = mat[jj*m + ii];
        }
    }
    for (ii = 1; ii < m*n; ii++) mat[ii] = tmp[ii];
    d_free(tmp);
}


static int_t read_int_vector_from_txt(int_t *vec, int_t n, const char *filename) {
    int_t ii;
    FILE *myFile;
    myFile = fopen(filename, "r");

    if (myFile == NULL) {
        printf("Error Reading File ! ! ! ! ! ! ! ! !\n");
        return -1;
    }

    for (ii = 0; ii < n; ii++) {
        fscanf(myFile, "%d,", &vec[ii]);
    }

    fclose(myFile);

    return 0;
}


static int_t read_double_vector_from_txt(real_t *vec, int_t n, const char *filename) {
    int_t ii;
    FILE *myFile;
    myFile = fopen(filename, "r");

    if (myFile == NULL) {
        printf("Error Reading File ! ! ! ! ! ! ! ! !\n");
        return -1;
    }

    for (ii = 0; ii < n; ii++) {
         fscanf(myFile, "%lf,", &vec[ii]);
    }

    fclose(myFile);

    return 0;
}


static int_t read_double_matrix_from_txt(real_t *mat, int_t m, int_t n, const char *filename) {
    int_t status;
    status = read_double_vector_from_txt(mat, m*n, filename);
    transpose_matrix(mat, m, n);
    return status;
}


static int_t read_ocp_qp_in_N(int_t *N, const char *fpath) {
    int_t status;
    char fname[256];

    snprintf(fname, sizeof(fname), "%s%s", fpath, "N.txt");
    status = read_int_vector_from_txt(N, 1, fname);
    return status;
}


static int_t read_ocp_qp_in_nx(int_t **nx, int_t N, const char *fpath) {
    int_t status;
    char fname[256];

    snprintf(fname, sizeof(fname), "%s%s", fpath, "nx.txt");
    status = read_int_vector_from_txt(*nx, N+1, fname);
    return status;
}


static int_t read_ocp_qp_in_nu(int_t **nu, int_t N, const char *fpath) {
    int_t status;
    char fname[256];

    snprintf(fname, sizeof(fname), "%s%s", fpath, "nu.txt");
    status = read_int_vector_from_txt(*nu, N, fname);
    return status;
}


static int_t read_ocp_qp_in_nb(int_t **nb, int_t N, const char *fpath) {
    int_t status;
    char fname[256];

    snprintf(fname, sizeof(fname), "%s%s", fpath, "nb.txt");
    status = read_int_vector_from_txt(*nb, N+1, fname);
    return status;
}


static int_t read_ocp_qp_in_nc(int_t **nc, int_t N, const char *fpath) {
    int_t status;
    char fname[256];

    snprintf(fname, sizeof(fname), "%s%s", fpath, "nc.txt");
    status = read_int_vector_from_txt(*nc, N+1, fname);
    return status;
}


static int_t read_ocp_qp_in_basic(ocp_qp_in *const in, const char *fpath) {
    int_t kk, N, status;
    char fname[256];
    char stage[16];
    N = in->N;

    for (kk = 0; kk < N; kk++) {
        snprintf(stage, sizeof(stage), "%d", kk);

        snprintf(fname, sizeof(fname), "%s%s%d%s", fpath, "Q", kk, ".txt");
        status = read_double_matrix_from_txt((real_t*)in->Q[kk], in->nx[kk], in->nx[kk], fname);
        if (status != 0) return status;

        snprintf(fname, sizeof(fname), "%s%s%d%s", fpath, "q", kk, ".txt");
        status = read_double_vector_from_txt((real_t*)in->q[kk], in->nx[kk], fname);
        if (status != 0) return status;

        snprintf(fname, sizeof(fname), "%s%s%d%s", fpath, "S", kk, ".txt");
        status = read_double_matrix_from_txt((real_t*)in->S[kk], in->nu[kk], in->nx[kk], fname);
        if (status != 0) return status;

        snprintf(fname, sizeof(fname), "%s%s%d%s", fpath, "R", kk, ".txt");
        status = read_double_matrix_from_txt((real_t*)in->R[kk], in->nu[kk], in->nu[kk], fname);
        if (status != 0) return status;

        snprintf(fname, sizeof(fname), "%s%s%d%s", fpath, "r", kk, ".txt");
        status = read_double_vector_from_txt((real_t*)in->r[kk], in->nu[kk], fname);
        if (status != 0) return status;

        snprintf(fname, sizeof(fname), "%s%s%d%s", fpath, "A", kk, ".txt");
        status = read_double_matrix_from_txt((real_t*)in->A[kk], in->nx[kk+1], in->nx[kk], fname);
        if (status != 0) return status;

        snprintf(fname, sizeof(fname), "%s%s%d%s", fpath, "B", kk, ".txt");
        status = read_double_matrix_from_txt((real_t*)in->B[kk], in->nx[kk+1], in->nu[kk], fname);
        if (status != 0) return status;

        snprintf(fname, sizeof(fname), "%s%s%d%s", fpath, "b", kk, ".txt");
        status = read_double_vector_from_txt((real_t*)in->b[kk], in->nx[kk+1], fname);
        if (status != 0) return status;
    }
    snprintf(stage, sizeof(stage), "%d", N);
    snprintf(fname, sizeof(fname), "%s%s%d%s", fpath, "Q", kk, ".txt");
    status = read_double_matrix_from_txt((real_t*)in->Q[N], in->nx[N], in->nx[N], fname);
    if (status != 0) return status;

    snprintf(fname, sizeof(fname), "%s%s%d%s", fpath, "q", kk, ".txt");
    status = read_double_vector_from_txt((real_t*)in->q[N], in->nx[N], fname);
    return status;
}


static int_t read_ocp_qp_in_bounds(ocp_qp_in *const in, const char *fpath) {
    char fname[256];
    char stage[16];
    int_t kk, status;

    for (kk = 0; kk <= in->N; kk++) {
        snprintf(stage, sizeof(stage), "%d", kk);

        snprintf(fname, sizeof(fname), "%s%s%d%s", fpath, "lb", kk, ".txt");
        status = read_double_vector_from_txt((real_t*)in->lb[kk], in->nb[kk], fname);
        if (status != 0) return status;

        snprintf(fname, sizeof(fname), "%s%s%d%s", fpath, "ub", kk, ".txt");
        status = read_double_vector_from_txt((real_t*)in->ub[kk], in->nb[kk], fname);
        if (status != 0) return status;

        snprintf(fname, sizeof(fname), "%s%s%d%s", fpath, "idxb", kk, ".txt");
        status = read_int_vector_from_txt((int_t*)in->idxb[kk], in->nb[kk], fname);
        if (status != 0) return status;
    }
    return status;
}


static int_t read_ocp_qp_in_polyhedral(ocp_qp_in *const in, const char *fpath) {
    char fname[256];
    char stage[16];
    int_t kk, status;

    for (kk = 0; kk <= in->N; kk++) {
        snprintf(stage, sizeof(stage), "%d", kk);

        snprintf(fname, sizeof(fname), "%s%s%d%s", fpath, "lc", kk, ".txt");
        status = read_double_vector_from_txt((real_t*)in->lc[kk], in->nc[kk], fname);
        if (status != 0) return status;

        snprintf(fname, sizeof(fname), "%s%s%d%s", fpath, "uc", kk, ".txt");
        status = read_double_vector_from_txt((real_t*)in->uc[kk], in->nc[kk], fname);
        if (status != 0) return status;

        snprintf(fname, sizeof(fname), "%s%s%d%s", fpath, "Cx", kk, ".txt");
        status = read_double_matrix_from_txt((real_t*)in->Cx[kk], in->nc[kk], in->nx[kk], fname);
        if (status != 0) return status;

        if (kk < in->N) {
            snprintf(fname, sizeof(fname), "%s%s%d%s", fpath, "Cu", kk, ".txt");
            status = read_double_matrix_from_txt((real_t*)in->Cu[kk], in->nc[kk], in->nu[kk],
            fname);
            if (status != 0) return status;
        }
    }
    return status;
}


static int_t read_ocp_qp_in_x0(ocp_qp_in *const in, const char *fpath) {
    char fname[256];
    int ii, status, *ptr;

    snprintf(fname, sizeof(fname), "%s%s", fpath, "x0.txt");
    status = read_double_vector_from_txt((real_t*)in->lb[0], in->nx[0], fname);
    if (status != 0) return status;
    status = read_double_vector_from_txt((real_t*)in->ub[0], in->nx[0], fname);
    if (status != 0) return status;

    ptr = (int_t*) in->idxb[0];
    for (ii = 0; ii < in->nx[0]; ii++) ptr[ii] = ii;

    return status;
}


void print_ocp_qp_in(ocp_qp_in const in) {
    int_t kk, ii;
    int_t N = in.N;
    printf("\n----------------------------------\n");
    printf("LOADED DATA FROM TXT FILES:\n\n");
    printf("N:\n%d\n\n", N);
    printf("nx:\n");
    int_print_mat(1, N+1, (int_t *)in.nx, 1);
    printf("nu:\n");
    int_print_mat(1, N, (int_t *)in.nu, 1);
    printf("nb:\n");
    int_print_mat(1, N+1, (int_t *)in.nb, 1);
    printf("nc:\n");
    int_print_mat(1, N+1, (int_t *)in.nc, 1);
    printf("bounds:\n");
    for (kk = 0; kk <= N; kk++) {
        for (ii = 0; ii < in.nb[kk]; ii++) {
            printf("%2.2f  <=  z_%d[%d]  <=  %2.2f\n",
            in.lb[kk][ii], kk, in.idxb[kk][ii], in.ub[kk][ii]);
        }
    }
    printf("\ninequality bounds:\n");
    for (kk = 0; kk <= N; kk++) {
        for (ii = 0; ii < in.nc[kk]; ii++) {
            printf("%2.2f  <=  [Cx_%d[%d] Cu_%d[%d]]*z_%d  <=  %2.2f\n",
            in.lc[kk][ii], kk, ii, kk, ii, kk, in.uc[kk][ii]);
        }
    }
    printf("\nobjective:");
    for (kk = 0; kk < N; kk++) {
        printf("\nQ[%d] =\n", kk);
        d_print_mat(in.nx[kk], in.nx[kk], (real_t*)in.Q[kk], in.nx[kk]);
        printf("\nR[%d] =\n", kk);
        d_print_mat(in.nu[kk], in.nu[kk], (real_t*)in.R[kk], in.nu[kk]);
        printf("\nS[%d] =\n", kk);
        d_print_mat(in.nu[kk], in.nx[kk], (real_t*)in.S[kk], in.nu[kk]);
        printf("\nq[%d] =\n", kk);
        d_print_mat(in.nx[kk], 1, (real_t*)in.q[kk], in.nx[kk]);
        printf("\nr[%d] =\n", kk);
        d_print_mat(in.nu[kk], 1, (real_t*)in.r[kk], in.nu[kk]);
    }
    printf("\nQ[%d] =\n", kk);
    d_print_mat(in.nx[N], in.nx[N], (real_t*)in.Q[N], in.nx[N]);
    printf("\nq[%d] =\n", kk);
    d_print_mat(in.nx[N], 1, (real_t*)in.q[N], in.nx[N]);

    printf("\nequalities:");
    for (kk = 0; kk < N; kk++) {
        printf("\nA[%d] =\n", kk);
        d_print_mat(in.nx[kk+1], in.nx[kk], (real_t*)in.A[kk], in.nx[kk+1]);
        printf("\nB[%d] =\n", kk);
        d_print_mat(in.nx[kk+1], in.nu[kk], (real_t*)in.B[kk], in.nx[kk+1]);
        printf("\nb[%d] =\n", kk);
        d_print_mat(in.nx[kk+1], 1, (real_t*)in.B[kk], in.nx[kk+1]);
    }
    printf("\ninequalities:");
    for (kk = 0; kk < N; kk++) {
        printf("\nCx[%d] =\n", kk);
        d_print_mat(in.nc[kk], in.nx[kk], (real_t*)in.Cx[kk], in.nc[kk]);
        printf("\nCu[%d] =\n", kk);
        d_print_mat(in.nc[kk], in.nu[kk], (real_t*)in.Cu[kk], in.nc[kk]);
    }
    printf("\nCx[%d] =\n", kk);
    d_print_mat(in.nc[N], in.nx[N], (real_t*)in.Cx[N], in.nc[N]);
    printf("\n----------------------------------\n");
}


int_t read_ocp_qp_in(ocp_qp_in *const in, const char *fpath_,
    int_t BOUNDS, int_t INEQUALITIES, int_t MPC, int_t QUIET) {
    char fpath[256];
    int_t ii, N, pathLength, status;
    int_t *nx, *nu, *nb, *nc;

    // add slash at the end if missing
    pathLength = (int_t) strlen(fpath_);
    // TODO(dimitris): add windows support..
    if (fpath_[pathLength-1] != '/') {
        snprintf(fpath, sizeof(fpath), "%s%c", fpath_, '/');
    } else {
        snprintf(fpath, sizeof(fpath), "%s", fpath_);
    }

    status = read_ocp_qp_in_N(&N, fpath);
    if (status != 0) return status;

    int_zeros(&nu, N, 1);
    int_zeros(&nx, N+1, 1);
    int_zeros(&nb, N+1, 1);
    int_zeros(&nc, N+1, 1);

    status = read_ocp_qp_in_nx(&nx, N, fpath);
    if (status != 0) return status;
    status = read_ocp_qp_in_nu(&nu, N, fpath);
    if (status != 0) return status;

    if (BOUNDS) status = read_ocp_qp_in_nb(&nb, N, fpath);
    if (status != 0) return status;
    if (INEQUALITIES) status = read_ocp_qp_in_nc(&nc, N, fpath);
    if (status != 0) return status;
    if (MPC && !BOUNDS) {
        nb[0] = nx[0];
    }

    allocate_ocp_qp_in(N, nx, nu, nb, nc, in);
    status = read_ocp_qp_in_basic(in, fpath);
    if (status != 0) return status;
    if (BOUNDS) status = read_ocp_qp_in_bounds(in, fpath);
    if (status != 0) return status;
    if (INEQUALITIES) status = read_ocp_qp_in_polyhedral(in, fpath);
    if (status != 0) return status;
    if (MPC) status = read_ocp_qp_in_x0(in, fpath);  // NOTE: call it AFTER the bounds
    if (status != 0) return status;

    if (BOUNDS && MPC) {
        // if BOUNDS == 1, we assume that x0 is already bounded..
        for (ii = 0; ii < in->nx[0]; ii++) {
            if (in->idxb[0][ii] != ii) {
                printf("\nERROR: Not implemented yet!\n");
                return -1;
            }
        }
    }

    if (!QUIET) print_ocp_qp_in(*in);

    int_free(nx);
    int_free(nu);
    int_free(nb);
    int_free(nc);
    return status;
}
