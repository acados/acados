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

// external
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
// blasfeo
#include "blasfeo/include/blasfeo_target.h"
#include "blasfeo/include/blasfeo_common.h"
#include "blasfeo/include/blasfeo_d_aux.h"
#include "blasfeo/include/blasfeo_d_aux_ext_dep.h"
#include "blasfeo/include/blasfeo_i_aux_ext_dep.h"
// hpipm
#include "hpipm/include/hpipm_d_ocp_qp.h"
#include "hpipm/include/hpipm_d_ocp_qp_sol.h"
// acados
#include "acados/dense_qp/dense_qp_common.h"
#include "acados/ocp_qp/ocp_qp_common.h"
#include "acados/ocp_qp/ocp_qp_common_frontend.h"
#include "acados/utils/print.h"



void print_matrix(char *file_name, const real_t *matrix, const int_t nrows,
                  const int_t ncols) {
    FILE *output;
    if (strcmp(file_name, "stdout") == 0) {
        output = stdout;
    } else {
        output = fopen(file_name, "w");
    }
    if (output == NULL) {
        fprintf(stderr, "Opening of file `%s' failed!\n", file_name);
    }
    // Assumes column major ordering
    for (int_t i = 0; i < nrows; i++) {
        for (int_t j = 0; j < ncols; j++) {
            fprintf(output, "%+.3e ", matrix[j * nrows + i]);
        }
        fprintf(output, "\n");
    }
    if (output != stdout) fclose(output);
}



void print_matrix_name(char *file_name, char *name, const real_t *matrix,
                       const int_t nrows, const int_t ncols) {
    FILE *output;
    if (strcmp(file_name, "stdout") == 0) {
        output = stdout;
    } else {
        output = fopen(file_name, "w");
    }
    if (output == NULL) {
        fprintf(stderr, "Opening of file `%s' failed!\n", file_name);
    }
    fprintf(output, "%s:\n", name);
    // Assumes column major ordering
    for (int_t i = 0; i < nrows; i++) {
        for (int_t j = 0; j < ncols; j++) {
            fprintf(output, "%+.3e ", matrix[j * nrows + i]);
        }
        fprintf(output, "\n");
    }
    if (output != stdout) fclose(output);
}



void print_int_matrix(char *file_name, const int_t *matrix, const int_t nrows,
                      const int_t ncols) {
    FILE *output;
    if (strcmp(file_name, "stdout") == 0) {
        output = stdout;
    } else {
        output = fopen(file_name, "w");
    }
    if (output == NULL) {
        fprintf(stderr, "Opening of file `%s' failed!\n", file_name);
    }
    // Assumes column major ordering
    for (int_t i = 0; i < nrows; i++) {
        for (int_t j = 0; j < ncols; j++) {
            fprintf(output, "%d ", matrix[j * nrows + i]);
        }
        fprintf(output, "\n");
    }
    if (output != stdout) fclose(output);
}



void print_array(char *file_name, real_t *array, int_t size) {
    print_matrix(file_name, array, size, 1);
}



void print_int_array(char *file_name, const int_t *array, int_t size) {
    print_int_matrix(file_name, array, size, 1);
}



// Read space delimited file into column-major matrix
void read_matrix(const char *file_name, real_t *array, const int_t nrows,
                 const int_t ncols) {
    FILE *file;
    file = fopen(file_name, "r");

    if (file == NULL) {
        printf("Error opening file %s ! ! ! ! ! ! ! ! !\n", file_name);
        exit(1);
    }

    // Read numbers from file into buffer.
    for (int i = 0; i < nrows; i++) {
        for (int j = 0; j < ncols; j++) {
            if (!fscanf(file, "%lf", &array[nrows * j + i])) {
                break;
            }
        }
    }

    fclose(file);
}


void print_ocp_qp_dims(ocp_qp_dims *dims)
{
    int N = dims->N;

    printf("k\tnx\tnu\tnb\tnbx\tnbu\tng\tns\n");

    for (int kk = 0; kk < N+1; kk++)
    {
        printf("%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t\n", kk, dims->nx[kk], dims->nu[kk], dims->nb[kk],
            dims->nbx[kk], dims->nbu[kk], dims->ng[kk], dims->ns[kk]);
    }

    printf("\nmemsize = %d\n", dims->memsize);
}



void print_dense_qp_dims(dense_qp_dims *dims)
{
    printf("nv = %d\n", dims->nv);
    printf("ne = %d\n", dims->ne);
    printf("nb = %d\n", dims->nb);
    printf("ng = %d\n", dims->ng);
    printf("ns = %d\n", dims->ns);
}



void print_ocp_qp_in(ocp_qp_in *qp_in)
{
    int N = qp_in->dim->N;
    int *nx = qp_in->dim->nx;
    int *nu = qp_in->dim->nu;
    int *nb = qp_in->dim->nb;
    int *ng = qp_in->dim->ng;

    for (int ii = 0; ii < N+1; ii++)
    {
        printf("k = %d\n\n", ii);

        printf("RSQrq =\n");
        d_print_strmat(nu[ii]+nx[ii]+1, nu[ii]+nx[ii], &qp_in->RSQrq[ii], 0 , 0);

        printf("rq =\n");
        d_print_tran_strvec(nu[ii]+nx[ii], &qp_in->rq[ii], 0);


        if (ii < N)
        {
            printf("BAbt =\n");
            d_print_strmat(nu[ii]+nx[ii]+1, nx[ii+1], &qp_in->BAbt[ii], 0 , 0);

            printf("b =\n");
            d_print_tran_strvec(nx[ii+1], &qp_in->b[ii], 0);
        }

        printf("idxb = (nb = %d = %d + %d)\n", qp_in->dim->nb[ii], qp_in->dim->nbu[ii], qp_in->dim->nbx[ii]);
        int_print_mat(1, nb[ii], qp_in->idxb[ii], 1);

        printf("d =\n");
        d_print_tran_strvec(2*nb[ii]+2*ng[ii], &qp_in->d[ii], 0);
    }
}



void print_ocp_qp_out(ocp_qp_out *qp_out)
{
    int N = qp_out->dim->N;
    int *nx = qp_out->dim->nx;
    int *nu = qp_out->dim->nu;
    int *nb = qp_out->dim->nb;
    int *ng = qp_out->dim->ng;

    for (int ii = 0; ii < N+1; ii++)
    {
        printf("k = %d\n\n", ii);

        printf("ux =\n");
        d_print_tran_strvec(nu[ii]+nx[ii], &qp_out->ux[ii], 0);

        if (ii < N)
        {
            printf("pi =\n");
            d_print_tran_strvec(nx[ii], &qp_out->pi[ii], 0);
        }

        printf("lam =\n");
        d_print_tran_strvec(2*nb[ii]+2*ng[ii], &qp_out->lam[ii], 0);
    }
}



void print_colmaj_ocp_qp_in(colmaj_ocp_qp_in *qp)
{
    int_t N = qp->N;
    printf("ocp_qp structure with contents:\n");
    printf("N: %d\n", qp->N);
    printf("nx:\n");
    print_int_matrix("stdout", qp->nx, 1, N + 1);
    printf("nu:\n");
    print_int_matrix("stdout", qp->nu, 1, N + 1);
    printf("nb:\n");
    print_int_matrix("stdout", qp->nb, 1, N + 1);
    printf("nc:\n");
    print_int_matrix("stdout", qp->nc, 1, N + 1);
    for (int_t stage = 0; stage < N+1; stage++) {
        if (stage < N) {
            printf("A[%d]:\n", stage);
            print_matrix("stdout", qp->A[stage], qp->nx[stage], qp->nx[stage]);
            printf("B[%d]:\n", stage);
            print_matrix("stdout", qp->B[stage], qp->nx[stage], qp->nu[stage]);
            printf("b[%d]:\n", stage);
            print_matrix("stdout", qp->b[stage], qp->nx[stage], 1);
        }
        printf("Q[%d]:\n", stage);
        print_matrix("stdout", qp->Q[stage], qp->nx[stage], qp->nx[stage]);
        printf("R[%d]:\n", stage);
        print_matrix("stdout", qp->R[stage], qp->nu[stage], qp->nu[stage]);
        printf("S[%d]:\n", stage);
        print_matrix("stdout", qp->S[stage], qp->nu[stage], qp->nx[stage]);
        printf("q[%d]:\n", stage);
        print_matrix("stdout", qp->q[stage], qp->nx[stage], 1);
        printf("r[%d]:\n", stage);
        print_matrix("stdout", qp->r[stage], qp->nu[stage], 1);
        printf("lb[%d]:\n", stage);
        print_matrix("stdout", qp->lb[stage], qp->nb[stage], 1);
        printf("ub[%d]:\n", stage);
        print_matrix("stdout", qp->ub[stage], qp->nb[stage], 1);
        printf("Cx[%d]:\n", stage);
        print_matrix("stdout", qp->Cx[stage], qp->nc[stage], qp->nx[stage]);
        printf("Cu[%d]:\n", stage);
        print_matrix("stdout", qp->Cu[stage], qp->nc[stage], qp->nu[stage]);
        printf("lc[%d]:\n", stage);
        print_matrix("stdout", qp->lc[stage], qp->nc[stage], 1);
        printf("uc[%d]:\n", stage);
        print_matrix("stdout", qp->uc[stage], qp->nc[stage], 1);
    }
    printf("\n");
}



void print_colmaj_ocp_qp_in_to_file(colmaj_ocp_qp_in *qp)
{
    char filename[MAX_STR_LEN];
    for (int_t i = 0; i <= qp->N; i++) {
        snprintf(filename, sizeof(filename), "Qm%d.txt", i);
        print_matrix(filename, qp->Q[i], qp->nx[i], qp->nx[i]);
        snprintf(filename, sizeof(filename), "Sm%d.txt", i);
        print_matrix(filename, qp->S[i], qp->nu[i], qp->nx[i]);
        snprintf(filename, sizeof(filename), "Rm%d.txt", i);
        print_matrix(filename, qp->R[i], qp->nu[i], qp->nu[i]);
        snprintf(filename, sizeof(filename), "qv%d.txt", i);
        print_matrix(filename, qp->q[i], qp->nx[i], 1);
        snprintf(filename, sizeof(filename), "rv%d.txt", i);
        print_matrix(filename, qp->r[i], qp->nu[i], 1);
        if (i < qp->N) {
            snprintf(filename, sizeof(filename), "Am%d.txt", i);
            print_matrix(filename, qp->A[i], qp->nx[i+1], qp->nx[i+1]);
            snprintf(filename, sizeof(filename), "Bm%d.txt", i);
            print_matrix(filename, qp->B[i], qp->nx[i+1], qp->nu[i]);
            snprintf(filename, sizeof(filename), "bv%d.txt", i);
            print_matrix(filename, qp->b[i], qp->nx[i+1], 1);
        }
        snprintf(filename, sizeof(filename), "idxb%d.txt", i);
        print_int_matrix(filename, qp->idxb[i], qp->nb[i], 1);
        snprintf(filename, sizeof(filename), "lb%d.txt", i);
        print_matrix(filename, qp->lb[i], qp->nb[i], 1);
        snprintf(filename, sizeof(filename), "ub%d.txt", i);
        print_matrix(filename, qp->ub[i], qp->nb[i], 1);
        snprintf(filename, sizeof(filename), "Cx%d.txt", i);
        print_matrix(filename, qp->Cx[i], qp->nc[i], qp->nx[i]);
        snprintf(filename, sizeof(filename), "Cu%d.txt", i);
        print_matrix(filename, qp->Cu[i], qp->nc[i], qp->nu[i]);
    }
}



void print_colmaj_ocp_qp_out(char *filename, colmaj_ocp_qp_in *qp, colmaj_ocp_qp_out *out)
{
    for (int_t i = 0; i <= qp->N; i++) {
        printf("x[%d]:\n", i);
        print_matrix(filename, out->x[i], qp->nx[i], 1);
        printf("u[%d]:\n", i);
        print_matrix(filename, out->u[i], qp->nu[i], 1);
        if (i < qp->N) {
            printf("pi[%d]:\n", i);
            print_matrix(filename, out->pi[i], qp->nx[i], 1);
        }
        printf("lam[%d]:\n", i);
        print_matrix(filename, out->x[i], 2*qp->nb[i] + 2*qp->nc[i], 1);
    }
}



void print_dense_qp_in(dense_qp_in *qp_in)
{
    int nv = qp_in->dim->nv;

    printf("H =\n");
    d_print_strmat(nv, nv, qp_in->Hv, 0, 0);
    // TODO(dimitris): print all data
}



void print_ocp_qp_info(ocp_qp_info *info)
{
    double misc = info->total_time - info->condensing_time - info->solve_QP_time - info->interface_time;
    assert((misc >= 0 || fabs(misc) <= ACADOS_EPS) && "sum of timings larger than total time!");

    printf("\n***************************************************************\n");
    printf("total time \t=\t%7.3f ms \t=\t %6.2f %%\n", 1000*info->total_time, 100.0);
    printf("condensing time =\t%7.3f ms \t=\t %6.2f %%\n", 1000*info->condensing_time, 100*info->condensing_time/info->total_time);
    printf("QP time \t=\t%7.3f ms \t=\t %6.2f %%\n", 1000*info->solve_QP_time, 100*info->solve_QP_time/info->total_time);
    printf("interface time \t=\t%7.3f ms \t=\t %6.2f %%\n", 1000*info->interface_time, 100*info->interface_time/info->total_time);
    printf("misc \t\t=\t%7.3f ms \t=\t %6.2f %%\n", 1000*misc, 100*misc/info->total_time);
    printf("***************************************************************\n\n");
}