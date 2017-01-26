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

#include "acados/utils/allocate_ocp_qp_in.h"
#include "acados/utils/read_ocp_qp_in.h"


static void transpose_matrix(real_t *mat, int m, int n) {
    // real_t tmp[m*n];
    real_t *tmp;
    d_zeros(&tmp, m, n);

    int r, c, i;

    for (c = 0; c < n; c++) {
        for (r = 0; r < m; r++) {
             tmp[r*n+c] = mat[c*m + r];
        }
    }
    for (i = 1; i < m*n; i++) {
        mat[i] = tmp[i];
    }
    d_free(tmp);
}



static int_t read_int_vector_from_txt(int_t *vec, int_t n, const char *filename) {
    int_t ii;
    FILE *myFile;
    myFile = fopen(filename, "r");

    if (myFile == NULL) {
        printf("Error Reading File ! ! ! ! ! ! ! ! !  \n");
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
        printf("Error Reading File ! ! ! ! ! ! ! ! !  \n");
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


static int_t read_ocp_qp_in_dimensions(int_t *N, int_t **nx, int_t **nu,  const char *fpath) {
    int_t ii;
    int_t status = 0;
    int_t quiet = 1;
    char fname[256];

    snprintf(fname, sizeof(fname), "%s%s", fpath, "N.txt");
    status = read_int_vector_from_txt(N, 1, fname);
    if (status != 0) return status;

    *nx = (int_t*)malloc(sizeof(int_t)*(*N+1));
    *nu = (int_t*)malloc(sizeof(int_t)*(*N));

    snprintf(fname, sizeof(fname), "%s%s", fpath, "nx.txt");
    status = read_int_vector_from_txt(*nx, *N+1, fname);
    if (status != 0) return status;

    snprintf(fname, sizeof(fname), "%s%s", fpath, "nu.txt");
    status = read_int_vector_from_txt(*nu, *N, fname);
    if (status != 0) return status;

    if (!quiet) {
        printf("\nDIMENSIONS:\n");
        printf("\nN = %d\n", *N);
        printf("\nnx = [ ");
        for (ii = 0; ii <= *N; ii++) printf("%d ", nx[0][ii]);
        printf("]\n");
        printf("\nnu = [ ");
        for (ii = 0; ii < *N; ii++) printf("%d ", nu[0][ii]);
        printf("]\n\n");
    }
    // free(fname);

    return status;
}


int_t read_ocp_qp_in_unconstrained(ocp_qp_in *const in, const char *fpath) {
    int_t ii, kk, N;
    int_t status = 0;
    int_t quiet = 1;
    int_t *nx, *nu, *ptr;
    char fname[256];
    char stage[16];

    read_ocp_qp_in_dimensions(&N, &nx, &nu, fpath);
    allocate_ocp_qp_in_unconstrained(N, nx, nu, in);

    // read x0
    snprintf(fname, sizeof(fname), "%s%s", fpath, "x0.txt");
    status = read_double_vector_from_txt((real_t*)in->lb[0], in->nx[0], fname);
    if (status != 0) return status;
    status = read_double_vector_from_txt((real_t*)in->ub[0], in->nx[0], fname);
    if (status != 0) return status;
    ptr = (int_t*) in->idxb[0];
    for (ii = 0; ii < in->nx[0]; ii++) ptr[ii] = ii;

    for (kk = 0; kk < in->N; kk++) {
        snprintf(stage, sizeof(stage), "%d", kk);

        snprintf(fname, sizeof(fname), "%s%s%d%s", fpath, "Q", kk, ".txt");
        status = read_double_matrix_from_txt((real_t*)in->Q[kk], in->nx[kk], in->nx[kk], fname);
        if (status != 0) return status;

        snprintf(fname, sizeof(fname), "%s%s%d%s", fpath, "q", kk, ".txt");
        status = read_double_vector_from_txt((real_t*)in->q[kk], in->nx[kk], fname);
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
    snprintf(fname, sizeof(fname), "%s%s%d%s", fpath, "Q", kk, ".txt");
    status = read_double_matrix_from_txt((real_t*)in->Q[N], in->nx[N], in->nx[N], fname);
    if (status != 0) return status;

    snprintf(fname, sizeof(fname), "%s%s%d%s", fpath, "q", kk, ".txt");
    status = read_double_vector_from_txt((real_t*)in->q[N], in->nx[N], fname);
    if (status != 0) return status;

    if (!quiet) {
        for (kk = 0; kk < in->N; kk++) {
            printf("\nQ[%d] =\n", kk);
            d_print_mat(in->nx[kk], in->nx[kk], (real_t*)in->Q[kk], in->nx[kk]);
            printf("\nR[%d] =\n", kk);
            d_print_mat(in->nu[kk], in->nu[kk], (real_t*)in->R[kk], in->nu[kk]);
            printf("\nq[%d] =\n", kk);
            d_print_mat(in->nx[kk], 1, (real_t*)in->q[kk], in->nx[kk]);
            printf("\nr[%d] =\n", kk);
            d_print_mat(in->nu[kk], 1, (real_t*)in->r[kk], in->nu[kk]);

            printf("\nA[%d] =\n", kk);
            d_print_mat(in->nx[kk+1], in->nx[kk], (real_t*)in->A[kk], in->nx[kk+1]);
            printf("\nB[%d] =\n", kk);
            d_print_mat(in->nx[kk+1], in->nu[kk], (real_t*)in->B[kk], in->nx[kk+1]);
            printf("\nb[%d] =\n", kk);
            d_print_mat(in->nx[kk+1], 1, (real_t*)in->B[kk], in->nx[kk+1]);
        }
        printf("\nQ[%d] =\n", kk);
        d_print_mat(in->nx[N], in->nx[N], (real_t*)in->Q[N], in->nx[N]);
        printf("\nq[%d] =\n", kk);
        d_print_mat(in->nx[N], 1, (real_t*)in->q[N], in->nx[N]);
    }

    free(nx);
    free(nu);
    return status;
}
