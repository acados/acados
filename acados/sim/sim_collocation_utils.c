/*
 * Copyright (c) The acados authors.
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


#include "acados/sim/sim_collocation_utils.h"

#include <math.h>
#include <assert.h>
#include <stdlib.h>
#include <string.h>

#include "blasfeo/include/blasfeo_d_aux_ext_dep.h"
#include "blasfeo/include/blasfeo_d_blas.h"
#include "blasfeo/include/blasfeo_d_kernel.h"
#include "blasfeo/include/blasfeo_i_aux_ext_dep.h"

#include "acados/utils/print.h"
#include "acados/utils/mem.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif



// TODO(rien): replace these LU codes with blasfeo
static double lu_system_solve(double *A, double *b, int *perm, int dim, int dim_rhs, void *work)
{
    char *c_ptr = work;

    double *b_perm = (double *) c_ptr;
    c_ptr += dim * dim_rhs * sizeof(double);

    double det;
    double swap;
    double valueMax;

    int i, j, k;
    int index1;
    int indexMax;
    int intSwap;
    double tmp_var;

    for (i = 0; i < dim; ++i)
    {
        perm[i] = i;
    }
    det = 1.0000000000000000e+00;
    for (i = 0; i < (dim - 1); i++)
    {
        indexMax = i;
        valueMax = fabs(A[i * dim + i]);
        for (j = (i + 1); j < dim; j++)
        {
            swap = fabs(A[i * dim + j]);
            if (swap > valueMax)
            {
                indexMax = j;
                valueMax = swap;
            }
        }
        if (indexMax > i)
        {
            for (k = 0; k < dim; ++k)
            {
                swap = A[k * dim + i];
                A[k * dim + i] = A[k * dim + indexMax];
                A[k * dim + indexMax] = swap;
            }
            intSwap = perm[i];
            perm[i] = perm[indexMax];
            perm[indexMax] = intSwap;
        }
        for (j = i + 1; j < dim; j++)
        {
            A[i * dim + j] = -A[i * dim + j] / A[i * dim + i];
            for (k = i + 1; k < dim; k++)
            {
                A[k * dim + j] += A[i * dim + j] * A[k * dim + i];
            }
        }
    }

    for (i = 0; i < dim; ++i)
    {
        index1 = perm[i];
        for (j = 0; j < dim_rhs; ++j)
        {
            b_perm[j * dim + i] = b[j * dim + index1];
        }
    }
    for (j = 1; j < dim; ++j)
    {
        for (i = 0; i < j; ++i)
        {
            tmp_var = A[i * dim + j];
            for (k = 0; k < dim_rhs; ++k)
            {
                b_perm[k * dim + j] += tmp_var * b_perm[k * dim + i];
            }
        }
    }
    for (i = dim - 1; - 1 < i; --i)
    {
        for (j = dim - 1; i < j; --j)
        {
            tmp_var = A[j * dim + i];
            for (k = 0; k < dim_rhs; ++k)
            {
                b_perm[k * dim + i] -= tmp_var * b_perm[k * dim + j];
            }
        }
        tmp_var = 1.0 / A[i * (dim + 1)];
        for (k = 0; k < dim_rhs; ++k)
        {
            b_perm[k * dim + i] = tmp_var * b_perm[k * dim + i];
        }
    }
    for (k = 0; k < dim * dim_rhs; ++k)
    {
        b[k] = b_perm[k];
    }

    return det;
}



static acados_size_t gauss_legendre_nodes_work_calculate_size(int ns)
{
    int N = ns - 1;
    int N1 = N + 1;
    int N2 = N + 2;

    acados_size_t size = 0;

    size += 4 * N1 * sizeof(double);       // x_init, y, y_prev, der_lgvm
    size += 1 * N1 * N2 * sizeof(double);  // lgvm

    return size;
}


// calculates Gauss-Legendre nodes (c) for Butcher tableau of size ns
static void gauss_legendre_nodes(int ns, double *nodes, void *work)
{
    int N = ns - 1;
    int N1 = N + 1;
    int N2 = N + 2;
    double err = 1;
    double eps = 2e-16;

    char *c_ptr = work;

    // x_init
    double *x_init = (double *) c_ptr;
    c_ptr += N1 * sizeof(double);
    // y
    double *y = (double *) c_ptr;
    c_ptr += N1 * sizeof(double);
    // y_prev
    double *y_prev = (double *) c_ptr;
    c_ptr += N1 * sizeof(double);
    // lgvm // Legendre-Gauss Vandermonde Matrix
    double *lgvm = (double *) c_ptr;
    c_ptr += N1 * N2 * sizeof(double);
    // der_lgvm // derivative of LGVM
    double *der_lgvm = (double *) c_ptr;
    c_ptr += N1 * sizeof(double);

    assert((char *) work + gauss_legendre_nodes_work_calculate_size(ns) >= c_ptr);

    double a = 0.0;
    double b = 1.0;  // code for collocation interval [a,b]

    for (int i = 0; i < N1; i++)
    {
        if (N > 0)
        {
            x_init[i] = -1 + i * 2.0 / N;
        }
        else
        {
            x_init[i] = -1;
        }
        y[i] = cos((2 * i + 1) * M_PI / (2 * N + 2)) + (0.27 / N1) * sin(M_PI * x_init[i] * N / N2);
        y_prev[i] = 2.0;
    }

    while (err > eps)
    {  // iterate until step sufficiently small
        for (int i = 0; i < N1; i++) lgvm[i] = 1.0;
        for (int i = 0; i < N1; i++) lgvm[N1 + i] = y[i];
        for (int k = 2; k < N2; k++)
        {
            for (int i = 0; i < N1; i++)
                lgvm[k * N1 + i] = ((2 * k - 1) * y[i] * lgvm[(k - 1) * N1 + i] -
                                    (k - 1) * lgvm[(k - 2) * N1 + i]) /
                                   k;
        }
        for (int i = 0; i < N1; i++)
            der_lgvm[i] = N2 * (lgvm[N * N1 + i] - y[i] * lgvm[N1 * N1 + i]) / (1 - pow(y[i], 2));
        for (int i = 0; i < N1; i++) y_prev[i] = y[i];

        for (int i = 0; i < N1; i++)
            y[i] = y_prev[i] - lgvm[N1 * N1 + i] / der_lgvm[i];  // Newton step

        err = 0;
        for (int i = 0; i < N1; i++)
        {
            if (err < fabs(y[i] - y_prev[i])) err = fabs(y[i] - y_prev[i]);
        }
    }
    for (int i = 0; i < N1; i++)
        nodes[i] = (a * (1 - y[i]) + 0.5 * b * (1 + y[i]));

    return;
}

// looks up Gauss-Radau IIA nodes (c) for Butcher tableau of size ns
static void gauss_radau_iia_nodes(int ns, double *nodes, void *work)
{
    // Radau nodes can be calculated as zeros of polynomial, see:
    // Encyclopedia of Applied and Computational Mathematics
    // Full entry, version of August 2, 2011
    // Radau methods
    // Ernst Hairer 1 , Gerhard Wanner 1
    // case 1:
    //     nodes[0] = 1.0;
    //     break;
    // case 2:
    //     nodes[0] = 1.0 / 3;
    //     nodes[1] = 1.0;
    //     break;
    // case 3:
    //     nodes[0] = (4 - sqrt(6.0)) / 10;
    //     nodes[1] = (4 + sqrt(6.0)) / 10;
    //     nodes[2] = 1.0;
    //     break;

    /* hard code: Radau nodes, as done in CasADi */
    // https://github.com/casadi/casadi/issues/673

    if (ns < 1 || ns > 9)
    {
        printf("\n\nerror: gauss_radau_iia_nodes only available for num_stages = 1,...,9; got %d.\n", ns);
        exit(1);
    }
    // Radau collocation points
    nodes[ns-1] = 1.00000000000000000000;
    switch (ns)
    {
        case 2:
            nodes[0] = 0.33333333333333337034;
            break;
        case 3:
            nodes[0] = 0.15505102572168222297;
            nodes[1] = 0.64494897427831787695;
            break;
        case 4:
            nodes[0] = 0.08858795951270420632;
            nodes[1] = 0.40946686444073465694;
            nodes[2] = 0.78765946176084700170;
            break;
        case 5:
            nodes[0] = 0.05710419611451822419;
            nodes[1] = 0.27684301363812369168;
            nodes[2] = 0.58359043236891683382;
            nodes[3] = 0.86024013565621926247;
            break;
        case 6:
            nodes[0] = 0.03980985705146905529;
            nodes[1] = 0.19801341787360787761;
            nodes[2] = 0.43797481024738621480;
            nodes[3] = 0.69546427335363603106;
            nodes[4] = 0.90146491420117347282;
            break;
        case 7:
            nodes[0] = 0.02931642715978521885;
            nodes[1] = 0.14807859966848435640;
            nodes[2] = 0.33698469028115418666;
            nodes[3] = 0.55867151877155019069;
            nodes[4] = 0.76923386203005450490;
            nodes[5] = 0.92694567131974103802;
            break;
        case 8:
            nodes[0] = 0.02247938643871305597;
            nodes[1] = 0.11467905316090415413;
            nodes[2] = 0.26578982278458951338;
            nodes[3] = 0.45284637366944457959;
            nodes[4] = 0.64737528288683043876;
            nodes[5] = 0.81975930826310761113;
            nodes[6] = 0.94373743946307731001;
            break;
        case 9:
            nodes[0] = 0.01777991514736393386;
            nodes[1] = 0.09132360789979432347;
            nodes[2] = 0.21430847939563035798;
            nodes[3] = 0.37193216458327238438;
            nodes[4] = 0.54518668480342658000;
            nodes[5] = 0.71317524285556954666;
            nodes[6] = 0.85563374295785443735;
            nodes[7] = 0.95536604471003006012;
            break;
    }
}

// acados_size_t gauss_simplified_work_calculate_size(int ns)
// {
//     acados_size_t size = 0;

//     size += 1 * 2 * ns * sizeof(double);   // D
//     size += 3 * ns * ns * sizeof(double);  // T, T_inv, lu_work

//     size += 1 * ns * sizeof(int);  // perm

//     return size;
// }


// // TODO(all): understand how this works and leave a comment!
// void gauss_simplified(int ns, Newton_scheme *scheme, void *work)
// {
//     char *c_ptr = work;

//     // D
//     double *D = (double *) c_ptr;
//     c_ptr += 2 * ns * sizeof(double);
//     // T
//     double *T = (double *) c_ptr;
//     c_ptr += ns * ns * sizeof(double);
//     // T_inv
//     double *T_inv = (double *) c_ptr;
//     c_ptr += ns * ns * sizeof(double);
//     // lu_work
//     double *lu_work = (double *) c_ptr;
//     c_ptr += ns * ns * sizeof(double);
//     // perm
//     int *perm = (int *) c_ptr;
//     c_ptr += ns * sizeof(int);

//     assert((char *) work + gauss_simplified_work_calculate_size(ns) >= c_ptr);

//     char simplified[MAX_STR_LEN];

//     snprintf(simplified, sizeof(simplified), "simplified/GL%d_simpl_%s.txt", 2 * ns, "D");
//     read_matrix(simplified, D, ns, 2);

//     snprintf(simplified, sizeof(simplified), "simplified/GL%d_simpl_%s.txt", 2 * ns, "T");
//     read_matrix(simplified, T, ns, ns);

//     scheme->single = false;
//     scheme->low_tria = 0;
//     for (int i = 0; i < ns; i++)
//     {
//         scheme->eig[i] = D[i];
//     }
//     for (int i = 0; i < ns * ns; i++)
//     {
//         scheme->transf2[i] = T[i];
//     }
//     // transf1_T:
//     for (int i = 0; i < ns; i++)
//     {
//         if ((i + 1) < ns)
//         {  // complex conjugate pair of eigenvalues
//             for (int i1 = i; i1 < i + 2; i1++)
//             {
//                 for (int i2 = 0; i2 < ns; i2++)
//                 {
//                     scheme->transf1_T[i2 * ns + i1] = 0.0;
//                     for (int i3 = 0; i3 < 2; i3++)
//                     {
//                         scheme->transf1_T[i2 * ns + i1] +=
//                             D[(i1 - i) * ns + (i + i3)] * T[(i + i3) * ns + i2];
//                     }
//                 }
//             }
//             i++;
//         }
//         else
//         {  // real eigenvalue
//             for (int i2 = 0; i2 < ns; i2++)
//             {
//                 scheme->transf1_T[i2 * ns + i] = D[i] * T[i * ns + i2];
//             }
//         }
//     }

//     for (int i = 0; i < ns; i++)
//     {
//         T_inv[i * (ns + 1)] = 1.0;
//     }

//     lu_system_solve(T, T_inv, perm, ns, ns, lu_work);

//     // transf1:
//     for (int i = 0; i < ns; i++)
//     {
//         if ((i + 1) < ns)
//         {  // complex conjugate pair of eigenvalues
//             for (int i1 = i; i1 < i + 2; i1++)
//             {
//                 for (int i2 = 0; i2 < ns; i2++)
//                 {
//                     scheme->transf1[i2 * ns + i1] = 0.0;
//                     for (int i3 = 0; i3 < 2; i3++)
//                     {
//                         scheme->transf1[i2 * ns + i1] += D[i3 * ns + i1] * T_inv[i2 * ns + i + i3];
//                     }
//                 }
//             }
//             i++;
//         }
//         else
//         {  // real eigenvalue
//             for (int i2 = 0; i2 < ns; i2++)
//             {
//                 scheme->transf1[i2 * ns + i] = D[i] * T_inv[i2 * ns + i];
//             }
//         }
//     }
//     // transf2_T:
//     for (int i = 0; i < ns; i++)
//     {
//         for (int i2 = 0; i2 < ns; i2++)
//         {
//             scheme->transf2_T[i2 * ns + i] = T_inv[i * ns + i2];
//         }
//     }

//     return;
// }



acados_size_t butcher_tableau_work_calculate_size(int ns)
{
    acados_size_t size = 0;

    size += 3 * ns * ns * sizeof(double);  // can_vm, rhs, lu_work
    size += 1 * ns * sizeof(int);  // perm

    acados_size_t size_legendre = gauss_legendre_nodes_work_calculate_size(ns);

    size = size > size_legendre ? size : size_legendre;
    make_int_multiple_of(8, &size);

    return size;
}



void calculate_butcher_tableau_from_nodes(int ns, double *nodes, double *b, double *A, void *work)
{
    int i, j, k;

    char *c_ptr = work;

    // can_vm
    double *can_vm = (double *) c_ptr;
    c_ptr += ns * ns * sizeof(double);
    // rhs
    double *rhs = (double *) c_ptr;
    c_ptr += ns * ns * sizeof(double);
    // lu_work
    double *lu_work = (double *) c_ptr;
    c_ptr += ns * ns * sizeof(double);
    // perm
    int *perm = (int *) c_ptr;
    c_ptr += ns * sizeof(int);

    assert((char *) work + butcher_tableau_work_calculate_size(ns) >= c_ptr);

    for (j = 0; j < ns; j++)
    {
        for (i = 0; i < ns; i++) can_vm[i + j * ns] = pow(nodes[i], j);
    }

    for (i = 0; i < ns * ns; i++) rhs[i] = 0.0;
    for (i = 0; i < ns; i++) rhs[i * (ns + 1)] = 1.0;

    lu_system_solve(can_vm, rhs, perm, ns, ns, lu_work);

    for (k = 0; k < ns; k++)
    {
        for (i = 0; i < ns; i++)
        {
            A[i * ns + k] = 0.0;
            for (j = 0; j < ns; j++)
            {
                A[i * ns + k] = A[i * ns + k] + pow(nodes[k], j + 1) / (j + 1) * rhs[i * ns + j];
            }
        }
    }

    for (i = 0; i < ns; i++)
    {
        b[i] = 0.0;
        for (j = 0; j < ns; j++)
        {
            b[i] = b[i] + 1.0 / (j + 1) * rhs[i * ns + j];
        }
    }

    return;
}


void calculate_butcher_tableau(int ns, sim_collocation_type collocation_type, double *c_vec, double *b_vec, double *A_mat, void *work)
{
    // compute collocation nodes
    switch (collocation_type)
    {
        case GAUSS_LEGENDRE:
            gauss_legendre_nodes(ns, c_vec, work);
            calculate_butcher_tableau_from_nodes(ns, c_vec, b_vec, A_mat, work);
            break;
        case GAUSS_RADAU_IIA:
            gauss_radau_iia_nodes(ns, c_vec, work);
            calculate_butcher_tableau_from_nodes(ns, c_vec, b_vec, A_mat, work);
            break;
        case EXPLICIT_RUNGE_KUTTA:
            get_explicit_butcher_tableau(ns, A_mat, b_vec, c_vec);
            break;
        default:
            printf("\nerror: calculate_butcher_tableau: unsupported collocation_type\n");
            exit(1);
    }
}


void get_explicit_butcher_tableau(int ns, double *A, double *b, double *c)
{
    switch (ns)
    {
        case 1:
        {
            // A
            A[0 + ns * 0] = 0.0;
            // b
            b[0] = 1.0;
            // c
            c[0] = 0.0;
            break;
        }
        case 2:
        {
            // A
            A[0 + ns * 0] = 0.0;
            A[0 + ns * 1] = 0.0;
            A[1 + ns * 0] = 0.5;
            A[1 + ns * 1] = 0.0;
            // b
            b[0] = 0.0;
            b[1] = 1.0;
            // c
            c[0] = 0.0;
            c[1] = 0.5;
            break;
        }
        case 3:
        {
            //A
            A[0 + ns * 0] = 0.0;
            A[0 + ns * 1] = 0.0;
            A[0 + ns * 2] = 0.0;
            A[1 + ns * 0] = 0.5;
            A[1 + ns * 1] = 0.0;
            A[1 + ns * 2] = 0.0;
            A[2 + ns * 0] = -1.0;
            A[2 + ns * 1] = 2.0;
            A[2 + ns * 2] = 0.0;
            //b
            b[0] = 1.0 / 6.0;
            b[1] = 2.0 / 3.0;
            b[2] = 1.0 / 6.0;
            // c
            c[0] = 0.0;
            c[1] = 0.5;
            c[2] = 1.0;
            break;
        }
        case 4:
        {
            // A
            A[0 + ns * 0] = 0.0;
            A[0 + ns * 1] = 0.0;
            A[0 + ns * 2] = 0.0;
            A[0 + ns * 3] = 0.0;
            A[1 + ns * 0] = 0.5;
            A[1 + ns * 1] = 0.0;
            A[1 + ns * 2] = 0.0;
            A[1 + ns * 3] = 0.0;
            A[2 + ns * 0] = 0.0;
            A[2 + ns * 1] = 0.5;
            A[2 + ns * 2] = 0.0;
            A[2 + ns * 3] = 0.0;
            A[3 + ns * 0] = 0.0;
            A[3 + ns * 1] = 0.0;
            A[3 + ns * 2] = 1.0;
            A[3 + ns * 3] = 0.0;
            // b
            b[0] = 1.0 / 6.0;
            b[1] = 1.0 / 3.0;
            b[2] = 1.0 / 3.0;
            b[3] = 1.0 / 6.0;
            // c
            c[0] = 0.0;
            c[1] = 0.5;
            c[2] = 0.5;
            c[3] = 1.0;
            break;
        }
        default:
        {
            printf("\n error: ERK: num_stages = %d not available. Only number of stages = {1,2,3,4} implemented!\n",ns);
            exit(1);
        }
    }
}
