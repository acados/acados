/*
 * Copyright 2019 Gianluca Frison, Dimitris Kouzoupis, Robin Verschueren,
 * Andrea Zanelli, Niels van Duijkeren, Jonathan Frey, Tommaso Sartor,
 * Branimir Novoselnik, Rien Quirynen, Rezart Qelibari, Dang Doan,
 * Jonas Koenemann, Yutao Chen, Tobias Schöls, Jonas Schlagenhauf, Moritz Diehl
 *
 * This file is part of acados.
 *
 * The 2-Clause BSD License
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 * this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 * this list of conditions and the following disclaimer in the documentation
 * and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.;
 */


// external
#include <assert.h>
#if defined(__MABX2__)
#include <brtenv.h>
#define printf(...)                               \
    msg_info_printf(MSG_SM_USER, 0, __VA_ARGS__); \
    ds1401_tic_delay(0.01);
#else
#include <stdio.h>
#endif
#include <math.h>
#include <stdlib.h>
#include <string.h>

// blasfeo
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

// void print_matrix(char *file_name, const real_t *matrix, const int_t nrows, const int_t ncols)
// {
//     FILE *output;
//     if (strcmp(file_name, "stdout") == 0)
//     {
//         output = stdout;
//     }
//     else
//     {
//         output = fopen(file_name, "w");
//     }
//     if (output == NULL)
//     {
//         fprintf(stderr, "Opening of file `%s' failed!\n", file_name);
//     }
//     // Assumes column major ordering
//     for (int_t i = 0; i < nrows; i++)
//     {
//         for (int_t j = 0; j < ncols; j++)
//         {
//             fprintf(output, "%+.3e ", matrix[j * nrows + i]);
//         }
//         fprintf(output, "\n");
//     }
//     if (output != stdout) fclose(output);
// }

// void print_matrix_name(char *file_name, char *name, const real_t *matrix, const int_t nrows,
//                        const int_t ncols)
// {
//     FILE *output;
//     if (strcmp(file_name, "stdout") == 0)
//     {
//         output = stdout;
//     }
//     else
//     {
//         output = fopen(file_name, "w");
//     }
//     if (output == NULL)
//     {
//         fprintf(stderr, "Opening of file `%s' failed!\n", file_name);
//     }
//     fprintf(output, "%s:\n", name);
//     // Assumes column major ordering
//     for (int_t i = 0; i < nrows; i++)
//     {
//         for (int_t j = 0; j < ncols; j++)
//         {
//             fprintf(output, "%+.3e ", matrix[j * nrows + i]);
//         }
//         fprintf(output, "\n");
//     }
//     if (output != stdout) fclose(output);
// }

// void print_int_matrix(char *file_name, const int_t *matrix, const int_t nrows, const int_t ncols)
// {
//     FILE *output;
//     if (strcmp(file_name, "stdout") == 0)
//     {
//         output = stdout;
//     }
//     else
//     {
//         output = fopen(file_name, "w");
//     }
//     if (output == NULL)
//     {
//         fprintf(stderr, "Opening of file `%s' failed!\n", file_name);
//     }
//     // Assumes column major ordering
//     for (int_t i = 0; i < nrows; i++)
//     {
//         for (int_t j = 0; j < ncols; j++)
//         {
//             fprintf(output, "%d ", matrix[j * nrows + i]);
//         }
//         fprintf(output, "\n");
//     }
//     if (output != stdout) fclose(output);
// }

// void print_array(char *file_name, real_t *array, int_t size)
// {
//     print_matrix(file_name, array, size, 1);
// }

// void print_int_array(char *file_name, const int_t *array, int_t size)
// {
//     print_int_matrix(file_name, array, size, 1);
// }

// Read space delimited file into column-major matrix
void read_matrix(const char *file_name, real_t *array, const int_t nrows, const int_t ncols)
{
    FILE *file;
    file = fopen(file_name, "r");

    if (file == NULL)
    {
        printf("Error opening file %s ! ! ! ! ! ! ! ! !\n", file_name);
        exit(1);
    }

    // Read numbers from file into buffer.
    for (int i = 0; i < nrows; i++)
    {
        for (int j = 0; j < ncols; j++)
        {
            if (!fscanf(file, "%lf", &array[nrows * j + i]))
            {
                break;
            }
        }
    }

    fclose(file);
}

void ocp_nlp_dims_print(ocp_nlp_dims *dims)
{
    int N = dims->N;

    //    printf("k\tnx\tnu\tnb\tnbx\tnbu\tng\tns\n");
    printf("k\tnx\tnu\tni\n");

    for (int kk = 0; kk < N + 1; kk++)
    {
        //        printf("%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t\n", kk, dims->nx[kk], dims->nu[kk],
        //        dims->nb[kk],
        //            dims->nbx[kk], dims->nbu[kk], dims->ng[kk], dims->ns[kk]);
        printf("%d\t%d\t%d\t%d\n", kk, dims->nx[kk], dims->nu[kk], dims->ni[kk]);
    }
}

void print_ocp_qp_dims(ocp_qp_dims *dims)
{
    int N = dims->N;

    printf("k\tnx\tnu\tnb\tnbx\tnbu\tng\tns\n");

    for (int kk = 0; kk < N + 1; kk++)
    {
        printf("%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t\n", kk, dims->nx[kk], dims->nu[kk], dims->nb[kk],
               dims->nbx[kk], dims->nbu[kk], dims->ng[kk], dims->ns[kk]);
    }
}

// void print_dense_qp_dims(dense_qp_dims *dims)
// {
//     printf("nv = %d\n", dims->nv);
//     printf("ne = %d\n", dims->ne);
//     printf("nb = %d\n", dims->nb);
//     printf("ng = %d\n", dims->ng);
//     printf("ns = %d\n", dims->ns);
// }

void print_ocp_qp_in(ocp_qp_in *qp_in)
{
#ifndef BLASFEO_EXT_DEP_OFF
    int N = qp_in->dim->N;
    int *nx = qp_in->dim->nx;
    int *nu = qp_in->dim->nu;
    int *nb = qp_in->dim->nb;
    int *ng = qp_in->dim->ng;
    int *ns = qp_in->dim->ns;

#if 1
    printf("BAbt =\n");
    for (int ii = 0; ii < N; ii++)
    {
        blasfeo_print_dmat(nu[ii] + nx[ii] + 1, nx[ii + 1], &qp_in->BAbt[ii], 0, 0);
    }

    printf("b =\n");
    for (int ii = 0; ii < N; ii++)
    {
        blasfeo_print_tran_dvec(nx[ii + 1], &qp_in->b[ii], 0);
    }

    printf("RSQrq =\n");
    for (int ii = 0; ii < N + 1; ii++)
    {
        blasfeo_print_dmat(nu[ii] + nx[ii] + 1, nu[ii] + nx[ii], &qp_in->RSQrq[ii], 0, 0);
    }

    printf("rq =\n");
    for (int ii = 0; ii < N + 1; ii++)
    {
        blasfeo_print_tran_dvec(nu[ii] + nx[ii], &qp_in->rqz[ii], 0);
    }

    printf("d =\n");
    for (int ii = 0; ii < N + 1; ii++)
    {
        blasfeo_print_tran_dvec(2 * nb[ii] + 2 * ng[ii], &qp_in->d[ii], 0);
    }

    printf("idxb =\n");
    for (int ii = 0; ii < N + 1; ii++)
    {
        int_print_mat(1, nb[ii], qp_in->idxb[ii], 1);
    }

    printf("DCt =\n");
    for (int ii = 0; ii < N + 1; ii++)
    {
        if (ng[ii] > 0)
            blasfeo_print_dmat(nu[ii] + nx[ii], ng[ii], &qp_in->DCt[ii], 0, 0);
    }

    printf("d_s =\n");
    for (int ii = 0; ii <= N; ii++)
    {
        if (ns[ii] > 0)
            blasfeo_print_tran_dvec(2 * ns[ii], &qp_in->d[ii], 2 * nb[ii] + 2 * ng[ii]);
    }

    printf("idxs_rev =\n");
    for (int ii = 0; ii <= N; ii++)
    {
        int_print_mat(1, nb[ii]+ng[ii], qp_in->idxs_rev[ii], 1);
    }

    printf("Z =\n");
    for (int ii = 0; ii <= N; ii++)
    {
        if (ns[ii] > 0)
            blasfeo_print_tran_dvec(2 * ns[ii], &qp_in->Z[ii], 0);
    }
    printf("z =\n");
    for (int ii = 0; ii <= N; ii++)
    {
        if (ns[ii] > 0)
            blasfeo_print_tran_dvec(2 * ns[ii], &qp_in->rqz[ii], nu[ii] + nx[ii]);
    }

#else

    for (int ii = 0; ii <= N; ii++)
    {
        printf("k = %d\n\n", ii);

        printf("RSQrq =\n");
        blasfeo_print_dmat(nu[ii] + nx[ii] + 1, nu[ii] + nx[ii], &qp_in->RSQrq[ii], 0, 0);

        printf("rq =\n");
        blasfeo_print_tran_dvec(nu[ii] + nx[ii], &qp_in->rqz[ii], 0);

        if (ii < N)
        {
            printf("BAbt =\n");
            blasfeo_print_dmat(nu[ii] + nx[ii] + 1, nx[ii + 1], &qp_in->BAbt[ii], 0, 0);

            printf("b =\n");
            blasfeo_print_tran_dvec(nx[ii + 1], &qp_in->b[ii], 0);
        }

        printf("d =\n");
        blasfeo_print_tran_dvec(2 * nb[ii] + 2 * ng[ii], &qp_in->d[ii], 0);

        printf("idxb = (nb = %d = %d + %d)\n", qp_in->dim->nb[ii], qp_in->dim->nbu[ii],
               qp_in->dim->nbx[ii]);
        int_print_mat(1, nb[ii], qp_in->idxb[ii], 1);

        printf("DCt =\n");
        blasfeo_print_dmat(nu[ii] + nx[ii], ng[ii], &qp_in->DCt[ii], 0, 0);

        printf("d_s =\n");
        blasfeo_print_tran_dvec(2 * ns[ii], &qp_in->d[ii], 2 * nb[ii] + 2 * ng[ii]);

        printf("idxs = (ns = %d)\n", qp_in->dim->ns[ii]);
        int_print_mat(1, ns[ii], qp_in->idxs[ii], 1);

        printf("Z =\n");
        blasfeo_print_tran_dvec(2 * ns[ii], &qp_in->Z[ii], 0);

        printf("z =\n");
        blasfeo_print_tran_dvec(2 * ns[ii], &qp_in->rqz[ii], nu[ii] + nx[ii]);
    }

#endif
#endif
    return;
}

void print_ocp_qp_out(ocp_qp_out *qp_out)
{
#ifndef BLASFEO_EXT_DEP_OFF

    int ii;

    int N = qp_out->dim->N;
    int *nx = qp_out->dim->nx;
    int *nu = qp_out->dim->nu;
    int *nb = qp_out->dim->nb;
    int *ng = qp_out->dim->ng;
    int *ns = qp_out->dim->ns;

#if 1

    printf("ux =\n");
    for (ii = 0; ii <= N; ii++)
        blasfeo_print_tran_dvec(nu[ii] + nx[ii] + 2 * ns[ii], &qp_out->ux[ii], 0);

    printf("pi =\n");
    for (ii = 0; ii < N; ii++) blasfeo_print_tran_dvec(nx[ii + 1], &qp_out->pi[ii], 0);

    printf("lam =\n");
    for (ii = 0; ii <= N; ii++)
        blasfeo_print_tran_dvec(2 * nb[ii] + 2 * ng[ii] + 2 * ns[ii], &qp_out->lam[ii], 0);

    printf("t =\n");
    for (ii = 0; ii <= N; ii++)
        blasfeo_print_exp_tran_dvec(2 * nb[ii] + 2 * ng[ii] + 2 * ns[ii], &qp_out->t[ii], 0);

#else

    for (ii = 0; ii < N + 1; ii++)
    {
        printf("k = %d\n\n", ii);

        printf("ux =\n");
        blasfeo_print_tran_dvec(nu[ii] + nx[ii], &qp_out->ux[ii], 0);

        if (ii < N)
        {
            printf("pi =\n");
            blasfeo_print_tran_dvec(nx[ii], &qp_out->pi[ii], 0);
        }

        printf("lam =\n");
        blasfeo_print_tran_dvec(2 * nb[ii] + 2 * ng[ii], &qp_out->lam[ii], 0);
    }

#endif
#endif
    return;
}


void print_ocp_qp_out_to_file(FILE *file, ocp_qp_out *qp_out)
{
#ifndef BLASFEO_EXT_DEP_OFF

    int ii;

    int N = qp_out->dim->N;
    int *nx = qp_out->dim->nx;
    int *nu = qp_out->dim->nu;
    int *nb = qp_out->dim->nb;
    int *ng = qp_out->dim->ng;
    int *ns = qp_out->dim->ns;

    fprintf(file, "ux =\n");
    for (ii = 0; ii <= N; ii++)
        blasfeo_print_to_file_tran_dvec(file, nu[ii] + nx[ii] + 2 * ns[ii], &qp_out->ux[ii], 0);

    fprintf(file, "pi =\n");
    for (ii = 0; ii < N; ii++) blasfeo_print_to_file_tran_dvec(file, nx[ii + 1], &qp_out->pi[ii], 0);

    fprintf(file, "lam =\n");
    for (ii = 0; ii <= N; ii++)
        blasfeo_print_to_file_tran_dvec(file, 2 * nb[ii] + 2 * ng[ii] + 2 * ns[ii], &qp_out->lam[ii], 0);

    fprintf(file, "t =\n");
    for (ii = 0; ii <= N; ii++)
        blasfeo_print_to_file_tran_dvec(file, 2 * nb[ii] + 2 * ng[ii] + 2 * ns[ii], &qp_out->t[ii], 0);

#endif
    return;
}

void ocp_nlp_out_print(ocp_nlp_dims *dims, ocp_nlp_out *nlp_out)
{
#ifndef BLASFEO_EXT_DEP_OFF
    int ii;

    int N = dims->N;
    int *nv = dims->nv;
    int *nx = dims->nx;
    // int *nu = dims->nu;
    int *ni = dims->ni;

#if 1

    printf("\nN = %d\n", N);
    printf("uxs =\n");
    for (ii = 0; ii <= N; ii++)
    {
        blasfeo_print_tran_dvec(nv[ii], &nlp_out->ux[ii], 0);
    }

    printf("pi =\n");
    for (ii = 0; ii < N; ii++)
    {
        blasfeo_print_tran_dvec(nx[ii + 1], &nlp_out->pi[ii], 0);
    }

    printf("lam =\n");
    for (ii = 0; ii <= N; ii++)
    {
        blasfeo_print_tran_dvec(2 * ni[ii], &nlp_out->lam[ii], 0);
    }

// printf("t =\n");
// for (ii=0; ii<=N; ii++)
//        blasfeo_print_exp_tran_dvec(2*nb[ii]+2*ng[ii], &nlp_out->t[ii], 0);

#else

    for (ii = 0; ii < N + 1; ii++)
    {
        printf("k = %d\n\n", ii);

        printf("ux =\n");
        blasfeo_print_tran_dvec(nu[ii] + nx[ii], &nlp_out->ux[ii], 0);

        if (ii < N)
        {
            printf("pi =\n");
            blasfeo_print_tran_dvec(nx[ii], &nlp_out->pi[ii], 0);
        }

        printf("lam =\n");
        blasfeo_print_tran_dvec(2 * nb[ii] + 2 * ng[ii], &nlp_out->lam[ii], 0);
    }

#endif
#endif
    return;
}

void ocp_nlp_res_print(ocp_nlp_dims *dims, ocp_nlp_res *nlp_res)
{
#ifndef BLASFEO_EXT_DEP_OFF
    int ii;

    int N = dims->N;
    int *nv = dims->nv;
    int *nx = dims->nx;
    // int *nu = dims->nu;
    int *ni = dims->ni;

    printf("res_stat =\n");
    for (ii = 0; ii <= N; ii++)
    {
        blasfeo_print_exp_tran_dvec(nv[ii], &nlp_res->res_stat[ii], 0);
    }

    printf("res_eq =\n");
    for (ii = 0; ii < N; ii++)
    {
        blasfeo_print_exp_tran_dvec(nx[ii + 1], &nlp_res->res_eq[ii], 0);
    }

    printf("res_ineq =\n");
    for (ii = 0; ii <= N; ii++)
    {
        blasfeo_print_exp_tran_dvec(2 * ni[ii], &nlp_res->res_ineq[ii], 0);
    }

    printf("res_comp =\n");
    for (ii = 0; ii <= N; ii++)
    {
        blasfeo_print_exp_tran_dvec(2 * ni[ii], &nlp_res->res_comp[ii], 0);
    }

    printf("inf norm res=\t%e\t%e\t%e\t%e\n", nlp_res->inf_norm_res_stat, nlp_res->inf_norm_res_eq,
        nlp_res->inf_norm_res_ineq, nlp_res->inf_norm_res_comp);

#endif
    return;
}

void print_ocp_qp_res(ocp_qp_res *qp_res)
{
#ifndef BLASFEO_EXT_DEP_OFF
    int ii;

    int N = qp_res->dim->N;
    int *nx = qp_res->dim->nx;
    int *nu = qp_res->dim->nu;
    int *nb = qp_res->dim->nb;
    int *ng = qp_res->dim->ng;

    printf("res_g =\n");
    for (ii = 0; ii <= N; ii++) blasfeo_print_exp_tran_dvec(nu[ii] + nx[ii], &qp_res->res_g[ii], 0);

    printf("res_b =\n");
    for (ii = 0; ii < N; ii++) blasfeo_print_exp_tran_dvec(nx[ii + 1], &qp_res->res_b[ii], 0);

    printf("res_d =\n");
    for (ii = 0; ii <= N; ii++)
        blasfeo_print_exp_tran_dvec(2 * nb[ii] + 2 * ng[ii], &qp_res->res_d[ii], 0);

    printf("res_m =\n");
    for (ii = 0; ii <= N; ii++)
        blasfeo_print_exp_tran_dvec(2 * nb[ii] + 2 * ng[ii], &qp_res->res_m[ii], 0);
#endif

    return;
}


static void int_print_mat_to_file(FILE *file, int row, int col, int *A, int lda)
{
    int i, j;
    for(i=0; i<row; i++)
    {
        for(j=0; j<col; j++)
        {
            fprintf(file, "%d ", *(A+i+j*lda));
        }
        fprintf(file, "\n");
    }
    fprintf(file, "\n");
}



void print_ocp_qp_in_to_file(FILE *file, ocp_qp_in *qp_in)
{
#ifndef BLASFEO_EXT_DEP_OFF
    int N = qp_in->dim->N;
    int *nx = qp_in->dim->nx;
    int *nu = qp_in->dim->nu;
    int *nb = qp_in->dim->nb;
    int *ng = qp_in->dim->ng;
    int *ns = qp_in->dim->ns;

    fprintf(file, "BAbt =\n");
    for (int ii = 0; ii < N; ii++)
    {
        blasfeo_print_to_file_dmat(file, nu[ii] + nx[ii] + 1, nx[ii + 1], &qp_in->BAbt[ii], 0, 0);
    }

    fprintf(file, "b =\n");
    for (int ii = 0; ii < N; ii++)
    {
        blasfeo_print_to_file_tran_dvec(file, nx[ii + 1], &qp_in->b[ii], 0);
    }

    fprintf(file, "RSQrq =\n");
    for (int ii = 0; ii < N + 1; ii++)
    {
        blasfeo_print_to_file_dmat(file, nu[ii] + nx[ii] + 1, nu[ii] + nx[ii], &qp_in->RSQrq[ii], 0, 0);
    }

    fprintf(file, "rq =\n");
    for (int ii = 0; ii < N + 1; ii++)
    {
        blasfeo_print_to_file_tran_dvec(file, nu[ii] + nx[ii], &qp_in->rqz[ii], 0);
    }

    fprintf(file, "d =\n");
    for (int ii = 0; ii < N + 1; ii++)
    {
        blasfeo_print_to_file_tran_dvec(file, 2 * nb[ii] + 2 * ng[ii], &qp_in->d[ii], 0);
    }

    fprintf(file, "idxb =\n");
    for (int ii = 0; ii < N + 1; ii++)
    {
        int_print_mat_to_file(file, 1, nb[ii], qp_in->idxb[ii], 1);
    }

    fprintf(file, "DCt =\n");
    for (int ii = 0; ii < N + 1; ii++)
    {
        if (ng[ii] > 0)
            blasfeo_print_to_file_dmat(file, nu[ii] + nx[ii], ng[ii], &qp_in->DCt[ii], 0, 0);
    }

    fprintf(file, "d_s =\n");
    for (int ii = 0; ii <= N; ii++)
    {
        if (ns[ii] > 0)
            blasfeo_print_to_file_tran_dvec(file, 2 * ns[ii], &qp_in->d[ii], 2 * nb[ii] + 2 * ng[ii]);
    }

    fprintf(file, "idxs_rev =\n");
    for (int ii = 0; ii <= N; ii++)
    {
        int_print_mat_to_file(file, 1, nb[ii]+ng[ii], qp_in->idxs_rev[ii], 1);
    }

    fprintf(file, "Z =\n");
    for (int ii = 0; ii <= N; ii++)
    {
        if (ns[ii] > 0)
            blasfeo_print_to_file_tran_dvec(file, 2 * ns[ii], &qp_in->Z[ii], 0);
    }
    fprintf(file, "z =\n");
    for (int ii = 0; ii <= N; ii++)
    {
        if (ns[ii] > 0)
            blasfeo_print_to_file_tran_dvec(file, 2 * ns[ii], &qp_in->rqz[ii], nu[ii] + nx[ii]);
    }

#endif
    return;
}



// void print_ocp_qp_in_to_string(char *string_out, ocp_qp_in *qp_in)
// {
//     int N = qp_in->dim->N;
//     int *nx = qp_in->dim->nx;
//     int *nu = qp_in->dim->nu;
//     int *nb = qp_in->dim->nb;
//     int *ng = qp_in->dim->ng;

//     for (int ii = 0; ii < N + 1; ii++)
//     {
//         string_out += sprintf(string_out, "k = %d\n\n", ii);

//         string_out += sprintf(string_out, "RSQrq =\n");
//         blasfeo_print_to_string_dmat(&string_out, nu[ii] + nx[ii] + 1, nu[ii] + nx[ii],
//                                      &qp_in->RSQrq[ii], 0, 0);

//         string_out += sprintf(string_out, "rq =\n");
//         blasfeo_print_to_string_tran_dvec(&string_out, nu[ii] + nx[ii], &qp_in->rqz[ii], 0);

//         if (ii < N)
//         {
//             string_out += sprintf(string_out, "BAbt =\n");
//             blasfeo_print_to_string_dmat(&string_out, nu[ii] + nx[ii] + 1, nx[ii + 1],
//                                          &qp_in->BAbt[ii], 0, 0);

//             string_out += sprintf(string_out, "b =\n");
//             blasfeo_print_to_string_tran_dvec(&string_out, nx[ii + 1], &qp_in->b[ii], 0);
//         }

//         string_out += sprintf(string_out, "idxb = (nb = %d = %d + %d)\n", qp_in->dim->nb[ii],
//                               qp_in->dim->nbu[ii], qp_in->dim->nbx[ii]);
//         int_print_to_string_mat(&string_out, 1, nb[ii], qp_in->idxb[ii], 1);

//         string_out += sprintf(string_out, "d =\n");
//         blasfeo_print_to_string_tran_dvec(&string_out, 2 * nb[ii] + 2 * ng[ii], &qp_in->d[ii], 0);
//     }
// }

// void print_ocp_qp_out_to_string(char *string_out, ocp_qp_out *qp_out)
// {
//     int N = qp_out->dim->N;
//     int *nx = qp_out->dim->nx;
//     int *nu = qp_out->dim->nu;
//     int *nb = qp_out->dim->nb;
//     int *ng = qp_out->dim->ng;

//     for (int ii = 0; ii < N + 1; ii++)
//     {
//         string_out += sprintf(string_out, "k = %d\n\n", ii);

//         string_out += sprintf(string_out, "ux =\n");
//         blasfeo_print_to_string_tran_dvec(&string_out, nu[ii] + nx[ii], &qp_out->ux[ii], 0);

//         if (ii < N)
//         {
//             string_out += sprintf(string_out, "pi =\n");
//             blasfeo_print_to_string_tran_dvec(&string_out, nx[ii], &qp_out->pi[ii], 0);
//         }

//         string_out += sprintf(string_out, "lam =\n");
//         blasfeo_print_to_string_tran_dvec(&string_out, 2 * nb[ii] + 2 * ng[ii], &qp_out->lam[ii],
//                                           0);
//     }
// }

// void print_colmaj_ocp_qp_in(colmaj_ocp_qp_in *qp)
// {
//     int_t N = qp->N;
//     printf("ocp_qp structure with contents:\n");
//     printf("N: %d\n", qp->N);
//     printf("nx:\n");
//     print_int_matrix("stdout", qp->nx, 1, N + 1);
//     printf("nu:\n");
//     print_int_matrix("stdout", qp->nu, 1, N + 1);
//     printf("nb:\n");
//     print_int_matrix("stdout", qp->nb, 1, N + 1);
//     printf("nc:\n");
//     print_int_matrix("stdout", qp->nc, 1, N + 1);
//     for (int_t stage = 0; stage < N + 1; stage++)
//     {
//         if (stage < N)
//         {
//             printf("A[%d]:\n", stage);
//             print_matrix("stdout", qp->A[stage], qp->nx[stage], qp->nx[stage]);
//             printf("B[%d]:\n", stage);
//             print_matrix("stdout", qp->B[stage], qp->nx[stage], qp->nu[stage]);
//             printf("b[%d]:\n", stage);
//             print_matrix("stdout", qp->b[stage], qp->nx[stage], 1);
//         }
//         printf("Q[%d]:\n", stage);
//         print_matrix("stdout", qp->Q[stage], qp->nx[stage], qp->nx[stage]);
//         printf("R[%d]:\n", stage);
//         print_matrix("stdout", qp->R[stage], qp->nu[stage], qp->nu[stage]);
//         printf("S[%d]:\n", stage);
//         print_matrix("stdout", qp->S[stage], qp->nu[stage], qp->nx[stage]);
//         printf("q[%d]:\n", stage);
//         print_matrix("stdout", qp->q[stage], qp->nx[stage], 1);
//         printf("r[%d]:\n", stage);
//         print_matrix("stdout", qp->r[stage], qp->nu[stage], 1);
//         printf("lb[%d]:\n", stage);
//         print_matrix("stdout", qp->lb[stage], qp->nb[stage], 1);
//         printf("ub[%d]:\n", stage);
//         print_matrix("stdout", qp->ub[stage], qp->nb[stage], 1);
//         printf("Cx[%d]:\n", stage);
//         print_matrix("stdout", qp->Cx[stage], qp->nc[stage], qp->nx[stage]);
//         printf("Cu[%d]:\n", stage);
//         print_matrix("stdout", qp->Cu[stage], qp->nc[stage], qp->nu[stage]);
//         printf("lc[%d]:\n", stage);
//         print_matrix("stdout", qp->lc[stage], qp->nc[stage], 1);
//         printf("uc[%d]:\n", stage);
//         print_matrix("stdout", qp->uc[stage], qp->nc[stage], 1);
//     }
//     printf("\n");
// }

// void print_colmaj_ocp_qp_in_to_file(colmaj_ocp_qp_in *qp)
// {
//     char filename[MAX_STR_LEN];
//     for (int_t i = 0; i <= qp->N; i++)
//     {
//         snprintf(filename, sizeof(filename), "Qm%d.txt", i);
//         print_matrix(filename, qp->Q[i], qp->nx[i], qp->nx[i]);
//         snprintf(filename, sizeof(filename), "Sm%d.txt", i);
//         print_matrix(filename, qp->S[i], qp->nu[i], qp->nx[i]);
//         snprintf(filename, sizeof(filename), "Rm%d.txt", i);
//         print_matrix(filename, qp->R[i], qp->nu[i], qp->nu[i]);
//         snprintf(filename, sizeof(filename), "qv%d.txt", i);
//         print_matrix(filename, qp->q[i], qp->nx[i], 1);
//         snprintf(filename, sizeof(filename), "rv%d.txt", i);
//         print_matrix(filename, qp->r[i], qp->nu[i], 1);
//         if (i < qp->N)
//         {
//             snprintf(filename, sizeof(filename), "Am%d.txt", i);
//             print_matrix(filename, qp->A[i], qp->nx[i + 1], qp->nx[i + 1]);
//             snprintf(filename, sizeof(filename), "Bm%d.txt", i);
//             print_matrix(filename, qp->B[i], qp->nx[i + 1], qp->nu[i]);
//             snprintf(filename, sizeof(filename), "bv%d.txt", i);
//             print_matrix(filename, qp->b[i], qp->nx[i + 1], 1);
//         }
//         snprintf(filename, sizeof(filename), "idxb%d.txt", i);
//         print_int_matrix(filename, qp->idxb[i], qp->nb[i], 1);
//         snprintf(filename, sizeof(filename), "lb%d.txt", i);
//         print_matrix(filename, qp->lb[i], qp->nb[i], 1);
//         snprintf(filename, sizeof(filename), "ub%d.txt", i);
//         print_matrix(filename, qp->ub[i], qp->nb[i], 1);
//         snprintf(filename, sizeof(filename), "Cx%d.txt", i);
//         print_matrix(filename, qp->Cx[i], qp->nc[i], qp->nx[i]);
//         snprintf(filename, sizeof(filename), "Cu%d.txt", i);
//         print_matrix(filename, qp->Cu[i], qp->nc[i], qp->nu[i]);
//     }
// }

// void print_colmaj_ocp_qp_out(char *filename, colmaj_ocp_qp_in *qp, colmaj_ocp_qp_out *out)
// {
//     for (int_t i = 0; i <= qp->N; i++)
//     {
//         printf("x[%d]:\n", i);
//         print_matrix(filename, out->x[i], qp->nx[i], 1);
//         printf("u[%d]:\n", i);
//         print_matrix(filename, out->u[i], qp->nu[i], 1);
//         if (i < qp->N)
//         {
//             printf("pi[%d]:\n", i);
//             print_matrix(filename, out->pi[i], qp->nx[i], 1);
//         }
//         printf("lam[%d]:\n", i);
//         print_matrix(filename, out->x[i], 2 * qp->nb[i] + 2 * qp->nc[i], 1);
//     }
// }

void print_dense_qp_in(dense_qp_in *qp_in)
{
#ifndef BLASFEO_EXT_DEP_OFF
    int nv = qp_in->dim->nv;
    int ne = qp_in->dim->ne;
    int nb = qp_in->dim->nb;
    int ng = qp_in->dim->ng;
    int ns = qp_in->dim->ns;

    printf("H =\n");
    blasfeo_print_dmat(nv, nv, qp_in->Hv, 0, 0);
    printf("g =\n");
    blasfeo_print_dvec(nv, qp_in->gz, 0);
    printf("A =\n");
    blasfeo_print_dmat(ne, nv, qp_in->A, 0, 0);
    printf("b =\n");
    blasfeo_print_dvec(ne, qp_in->b, 0);
    printf("Ct =\n");
    blasfeo_print_dmat(nv, ng, qp_in->Ct, 0, 0);
    printf("d =\n");
    blasfeo_print_dvec(2*nb+2*ng+2*ns, qp_in->d, 0);
#endif
}

void print_qp_info(qp_info *info)
{
    double misc =
        info->total_time - info->condensing_time - info->solve_QP_time - info->interface_time;
    //    assert((misc >= 0 || fabs(misc) <= ACADOS_EPS) && "sum of timings larger than total
    //    time!");

    printf("\n***************************************************************\n");
    printf("total time \t=\t%7.3f ms \t=\t %6.2f %%\n", 1000 * info->total_time, 100.0);
    printf("condensing time =\t%7.3f ms \t=\t %6.2f %%\n", 1000 * info->condensing_time,
           100 * info->condensing_time / info->total_time);
    printf("QP time \t=\t%7.3f ms \t=\t %6.2f %%\n", 1000 * info->solve_QP_time,
           100 * info->solve_QP_time / info->total_time);
    printf("interface time \t=\t%7.3f ms \t=\t %6.2f %%\n", 1000 * info->interface_time,
           100 * info->interface_time / info->total_time);
    printf("misc \t\t=\t%7.3f ms \t=\t %6.2f %%\n", 1000 * misc, 100 * misc / info->total_time);
    printf("QP iter \t=\t%7d\n", info->num_iter);
    printf("***************************************************************\n\n");
}

// void acados_warning(char warning_string [])
// {
//     printf("\n-> acados warning: %s\n", warning_string);
// }

// void acados_error(char error_string [])
// {
//     printf("\n-> acados error: %s. Exiting. \n", error_string);
//     assert(0);
// }

// void acados_not_implemented(char feature_string [])
// {
//     printf("\n-> acados feature %s not yet implemented. Exiting.\n", feature_string);
//     assert(0);
// }

// void print_blasfeo_target()
// {
//     printf("\n");
// #if defined(LA_HIGH_PERFORMANCE)
//     printf("blasfeo compiled with LA = HIGH_PERFORMANCE\n");
// #elif defined(LA_REFERENCE)
//     printf("blasfeo compiled with LA = REFERENCE\n");
// #elif defined(LA_BLAS)
//     printf("blasfeo compiled with LA = BLAS\n");
// #endif
//     printf("\n");
// }
