/*
 * Copyright (c) The acados authors.
 *
 * This file is part of acados.
 *
 * Licensed under the 2-Clause BSD License.
 */


#include "acados/sim/sim_collocation_utils.h"

#include <math.h>
#include <assert.h>
#include <stdlib.h>
#include <string.h>

#include "blasfeo_d_aux_ext_dep.h"
#include "blasfeo_d_blas.h"
#include "blasfeo_d_kernel.h"
#include "blasfeo_i_aux_ext_dep.h"

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
        case ESDIRK:
            get_esdirk_butcher_tableau(ns, A_mat, b_vec, c_vec);
            break;
        case SDIRK:
            get_sdirk_butcher_tableau(ns, A_mat, b_vec, c_vec);
            break;
        case DIRK:
            get_dirk_butcher_tableau(ns, A_mat, b_vec, c_vec);
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

void get_sdirk_butcher_tableau(int ns, double *A, double *b, double *c)
{
    switch (ns)
    {
        case 1:
        {
            // Implicit euler as sdirk method
            // A integral
            A[0 + ns * 0] = 1.0;
            // b integral
            b[0] = 1.0;
            // c integral
            c[0] = 1.0;

            break;
        }
        case 2:
        {
            // A stable SDIRK method of order 3, stage order 1, with 2 stages.
            // https://ntrs.nasa.gov/api/citations/20010075154/downloads/20010075154.pdf
            // equation (50)
            // A integral
            A[0 + ns * 0] = 1.2679491924311228;
            A[0 + ns * 1] = 0.0;
            A[1 + ns * 0] = 0.9282032302755092;
            A[1 + ns * 1] = 1.2679491924311228;
            // b integral
            b[0] = 1.098076211353316;
            b[1] = 0.6339745962155614;
            // c integral
            c[0] = 0.7886751345948129;
            c[1] = 0.2113248654051871;

            break;
        }
        case 3:
        {
            // best stability
            // A integral
            A[0 + ns * 0] = 0.9358222275240878;
            A[0 + ns * 1] = 0.0;
            A[0 + ns * 2] = 0.0;
            A[1 + ns * 0] = 0.497940606760015;
            A[1 + ns * 1] = 0.9358222275240878;
            A[1 + ns * 2] = 0.0;
            A[2 + ns * 0] = -0.34586591580096876;
            A[2 + ns * 1] = 2.867525668568206;
            A[2 + ns * 2] = 0.9358222275240878;
            // b integral
            b[0] = 0.4456224072877139;
            b[1] = 1.064177772475912;
            b[2] = 0.12061475842818323;
            // c integral
            c[0] = 1.0685790213016289;
            c[1] = 0.5;
            c[2] = -0.06857902130162881;
            /* // A integral */
            /* A[0 + ns * 0] = 3.3054072893322783; */
            /* A[0 + ns * 1] = 0.0; */
            /* A[0 + ns * 2] = 0.0; */
            /* A[1 + ns * 0] = -2.157451384853202; */
            /* A[1 + ns * 1] = 3.3054072893322783; */
            /* A[1 + ns * 2] = 0.0; */
            /* A[2 + ns * 0] = -8.109364584203572; */
            /* A[2 + ns * 1] = 2.295911808958153; */
            /* A[2 + ns * 2] = 3.3054072893322783; */
            /* // b integral */
            /* b[0] = -2.680044790717033; */
            /* b[1] = -1.3054072893322786; */
            /* b[2] = 3.5320888862379562; */
            /* // c integral */
            /* c[0] = 0.3025345781826508; */
            /* c[1] = 0.5; */
            /* c[2] = 0.6974654218173493; */
            /* // A integral */
            /* A[0 + ns * 0] = 7.7587704831436355; */
            /* A[0 + ns * 1] = 0.0; */
            /* A[0 + ns * 2] = 0.0; */
            /* A[1 + ns * 0] = -22.340489221906825; */
            /* A[1 + ns * 1] = 7.7587704831436355; */
            /* A[1 + ns * 2] = 0.0; */
            /* A[2 + ns * 0] = 68.45523050000462; */
            /* A[2 + ns * 1] = -29.163437477526376; */
            /* A[2 + ns * 2] = 7.7587704831436355; */
            /* // b integral */
            /* b[0] = 14.234422383429331; */
            /* b[1] = -5.758770483143636; */
            /* b[2] = 2.347296355333861; */
            /* // c integral */
            /* c[0] = 0.1288864005157204; */
            /* c[1] = 0.5; */
            /* c[2] = 0.8711135994842796; */
            break;
        }
        case 5:
        {
            // Stiffly Accurate, L-Stable, SDIRK method of order 4, stage order 1, with 5 stages.
            // Harrier and Wanner, Solving Ordinary Differential Equations II, page 100, (6.16)
            // A
            A[0 + ns * 0] = 1.0 / 4.0;
            A[0 + ns * 1] = 0.0;
            A[0 + ns * 2] = 0.0;
            A[0 + ns * 3] = 0.0;
            A[0 + ns * 4] = 0.0;
            A[1 + ns * 0] = 1.0 / 2.0;
            A[1 + ns * 1] = 1.0 / 4.0;
            A[1 + ns * 2] = 0.0;
            A[1 + ns * 3] = 0.0;
            A[1 + ns * 4] = 0.0;
            A[2 + ns * 0] = 17.0 / 50.0;
            A[2 + ns * 1] = -1.0 / 25.0;
            A[2 + ns * 2] = 1.0 / 4.0;
            A[2 + ns * 3] = 0.0;
            A[2 + ns * 4] = 0.0;
            A[3 + ns * 0] = 371.0 / 1360.0;
            A[3 + ns * 1] = -137.0 / 2720.0;
            A[3 + ns * 2] = 15.0 / 544.0;
            A[3 + ns * 3] = 1.0 / 4.0;
            A[3 + ns * 4] = 0.0;
            A[4 + ns * 0] = 25.0 / 24.0;
            A[4 + ns * 1] = -49.0 / 48.0;
            A[4 + ns * 2] = 125.0 / 16.0;
            A[4 + ns * 3] = -85.0 / 12.0;
            A[4 + ns * 4] = 1.0 / 4.0;
            // b
            b[0] = 25.0 / 24.0;
            b[1] = -49.0 / 48.0;
            b[2] = 125.0 / 16.0;
            b[3] = -85.0 / 12.0;
            b[4] = 1.0 / 4.0;
            // c
            c[0] = 1.0 / 4.0;
            c[1] = 3.0 / 4.0;
            c[2] = 11.0 / 20.0;
            c[3] = 1.0 / 2.0;
            c[4] = 1.0;

            break;
        }
        case 9:
        {
            // A integral
            A[0 + ns * 0] = 4.584468750764483;
            A[0 + ns * 1] = 0.0;
            A[0 + ns * 2] = 0.0;
            A[0 + ns * 3] = 0.0;
            A[0 + ns * 4] = 0.0;
            A[0 + ns * 5] = 0.0;
            A[0 + ns * 6] = 0.0;
            A[0 + ns * 7] = 0.0;
            A[0 + ns * 8] = 0.0;
            A[1 + ns * 0] = 1.8989491328422865;
            A[1 + ns * 1] = 4.584468750764483;
            A[1 + ns * 2] = 0.0;
            A[1 + ns * 3] = 0.0;
            A[1 + ns * 4] = 0.0;
            A[1 + ns * 5] = 0.0;
            A[1 + ns * 6] = 0.0;
            A[1 + ns * 7] = 0.0;
            A[1 + ns * 8] = 0.0;
            A[2 + ns * 0] = -0.5561893238812065;
            A[2 + ns * 1] = 7.432892450027895;
            A[2 + ns * 2] = 4.584468750764483;
            A[2 + ns * 3] = 0.0;
            A[2 + ns * 4] = 0.0;
            A[2 + ns * 5] = 0.0;
            A[2 + ns * 6] = 0.0;
            A[2 + ns * 7] = 0.0;
            A[2 + ns * 8] = 0.0;
            A[3 + ns * 0] = -11.049759380723843;
            A[3 + ns * 1] = -0.11794945196882703;
            A[3 + ns * 2] = 0.3027149332471444;
            A[3 + ns * 3] = 4.584468750764483;
            A[3 + ns * 4] = 0.0;
            A[3 + ns * 5] = 0.0;
            A[3 + ns * 6] = 0.0;
            A[3 + ns * 7] = 0.0;
            A[3 + ns * 8] = 0.0;
            A[4 + ns * 0] = 0.06850489003494811;
            A[4 + ns * 1] = -5.8252685493318515;
            A[4 + ns * 2] = -4.537843642344846;
            A[4 + ns * 3] = 0.4340491884966701;
            A[4 + ns * 4] = 4.584468750764483;
            A[4 + ns * 5] = 0.0;
            A[4 + ns * 6] = 0.0;
            A[4 + ns * 7] = 0.0;
            A[4 + ns * 8] = 0.0;
            A[5 + ns * 0] = -17.00435152896912;
            A[5 + ns * 1] = -5.8919841565920485;
            A[5 + ns * 2] = -3.4646265857898535;
            A[5 + ns * 3] = -0.7161723475289629;
            A[5 + ns * 4] = 7.374820515533151;
            A[5 + ns * 5] = 4.584468750764483;
            A[5 + ns * 6] = 0.0;
            A[5 + ns * 7] = 0.0;
            A[5 + ns * 8] = 0.0;
            A[6 + ns * 0] = -6.950755145922249;
            A[6 + ns * 1] = -5.6357037806741594;
            A[6 + ns * 2] = -9.787894730249818;
            A[6 + ns * 3] = 1.0859473919889138;
            A[6 + ns * 4] = 2.470596480516108;
            A[6 + ns * 5] = -0.44940495536564223;
            A[6 + ns * 6] = 4.584468750764483;
            A[6 + ns * 7] = 0.0;
            A[6 + ns * 8] = 0.0;
            A[7 + ns * 0] = -15.806046349411968;
            A[7 + ns * 1] = -22.415289767618436;
            A[7 + ns * 2] = 2.7926792928271835;
            A[7 + ns * 3] = -14.64310339699807;
            A[7 + ns * 4] = 26.9807455026971;
            A[7 + ns * 5] = 10.951713067000027;
            A[7 + ns * 6] = -7.952789324789637;
            A[7 + ns * 7] = 4.584468750764483;
            A[7 + ns * 8] = 0.0;
            A[8 + ns * 0] = 3.433922955161458;
            A[8 + ns * 1] = -24.029408301849863;
            A[8 + ns * 2] = -4.4537310739969005;
            A[8 + ns * 3] = -33.611082189404776;
            A[8 + ns * 4] = 24.376513447879983;
            A[8 + ns * 5] = 12.714242667842194;
            A[8 + ns * 6] = -1.2537906153730913;
            A[8 + ns * 7] = 6.7237108325106565;
            A[8 + ns * 8] = 4.584468750764483;
            // b integral
            b[0] = 0.0;
            b[1] = 0.0;
            b[2] = 0.0;
            b[3] = 0.0;
            b[4] = 0.0;
            b[5] = 0.0;
            b[6] = 0.0;
            b[7] = 0.0;
            b[8] = 1.0;
            // c integral
            c[0] = 0.218127781944908;
            c[1] = 0.127776296332966;
            c[2] = 0.0374248107210239;
            c[3] = 0.744688596518583;
            c[4] = 0.343766086047543;
            c[5] = 0.783024769450069;
            c[6] = 0.500923374474345;
            c[7] = 0.925980212553918;
            c[8] = 1.0;

            break;
        }
        default:
        {
            printf("\n error: SDIRK: num_stages = %d not available. Only number of stages = {2,3,5} implemented!\n",ns);
            exit(1);
        }
    }
}

void get_dirk_butcher_tableau(int ns, double *A, double *b, double *c)
{
    switch (ns)
    {
        case 6:
        {
            // A integral
            A[0 + ns * 0] = 3.295024883010639;
            A[0 + ns * 1] = 0.0;
            A[0 + ns * 2] = 0.0;
            A[0 + ns * 3] = 0.0;
            A[0 + ns * 4] = 0.0;
            A[0 + ns * 5] = 0.0;
            A[1 + ns * 0] = 1.8434903557385343;
            A[1 + ns * 1] = 1.9998710642297644;
            A[1 + ns * 2] = 0.0;
            A[1 + ns * 3] = 0.0;
            A[1 + ns * 4] = 0.0;
            A[1 + ns * 5] = 0.0;
            A[2 + ns * 0] = -0.53416359646691;
            A[2 + ns * 1] = 4.037102931086472;
            A[2 + ns * 2] = 4.602976382238963;
            A[2 + ns * 3] = 0.0;
            A[2 + ns * 4] = 0.0;
            A[2 + ns * 5] = 0.0;
            A[3 + ns * 0] = -7.000282320364896;
            A[3 + ns * 1] = -3.600471105615286;
            A[3 + ns * 2] = 4.957064833763692;
            A[3 + ns * 3] = 4.037661962977095;
            A[3 + ns * 4] = 0.0;
            A[3 + ns * 5] = 0.0;
            A[4 + ns * 0] = 0.4102327195731729;
            A[4 + ns * 1] = -2.1229097148735967;
            A[4 + ns * 2] = -2.840537042266239;
            A[4 + ns * 3] = 0.9795825031060229;
            A[4 + ns * 4] = 1.1537435869683463;
            A[4 + ns * 5] = 0.0;
            A[5 + ns * 0] = 2.6780621774808924;
            A[5 + ns * 1] = 0.7422115909150838;
            A[5 + ns * 2] = 0.21873314405788713;
            A[5 + ns * 3] = 0.14867974931749606;
            A[5 + ns * 4] = -1.0186787116908351;
            A[5 + ns * 5] = 1.4356472990415838;
            // b integral
            b[0] = -1.1075220742405483;
            b[1] = -0.9695102361041349;
            b[2] = 0.6486190903767947;
            b[3] = 1.3677771503132274;
            b[4] = 0.7302551199901542;
            b[5] = -0.3841304062808289;
            // c integral
            c[0] = 0.303487844706747;
            c[1] = 0.220275743310934;
            c[2] = 0.0592738886728798;
            c[3] = 0.89749279829076;
            c[4] = 0.548064850129365;
            c[5] = 0.303450087867971;


            break;
        }
        case 10:
        {
            // A integral
            A[0 + ns * 0] = 4.278905346916733;
            A[0 + ns * 1] = 0.0;
            A[0 + ns * 2] = 0.0;
            A[0 + ns * 3] = 0.0;
            A[0 + ns * 4] = 0.0;
            A[0 + ns * 5] = 0.0;
            A[0 + ns * 6] = 0.0;
            A[0 + ns * 7] = 0.0;
            A[0 + ns * 8] = 0.0;
            A[0 + ns * 9] = 0.0;
            A[1 + ns * 0] = 1.5813010569440609;
            A[1 + ns * 1] = 4.998579450054526;
            A[1 + ns * 2] = 0.0;
            A[1 + ns * 3] = 0.0;
            A[1 + ns * 4] = 0.0;
            A[1 + ns * 5] = 0.0;
            A[1 + ns * 6] = 0.0;
            A[1 + ns * 7] = 0.0;
            A[1 + ns * 8] = 0.0;
            A[1 + ns * 9] = 0.0;
            A[2 + ns * 0] = -6.164507830348227;
            A[2 + ns * 1] = -9.905836828209257;
            A[2 + ns * 2] = 7.504960640715964;
            A[2 + ns * 3] = 0.0;
            A[2 + ns * 4] = 0.0;
            A[2 + ns * 5] = 0.0;
            A[2 + ns * 6] = 0.0;
            A[2 + ns * 7] = 0.0;
            A[2 + ns * 8] = 0.0;
            A[2 + ns * 9] = 0.0;
            A[3 + ns * 0] = -1.5316612201417596;
            A[3 + ns * 1] = 9.814237732928433;
            A[3 + ns * 2] = 0.14523523745884595;
            A[3 + ns * 3] = 3.975847265182771;
            A[3 + ns * 4] = 0.0;
            A[3 + ns * 5] = 0.0;
            A[3 + ns * 6] = 0.0;
            A[3 + ns * 7] = 0.0;
            A[3 + ns * 8] = 0.0;
            A[3 + ns * 9] = 0.0;
            A[4 + ns * 0] = 8.115586383156664;
            A[4 + ns * 1] = 16.438781754518068;
            A[4 + ns * 2] = -16.052688382979568;
            A[4 + ns * 3] = 0.2529562402025452;
            A[4 + ns * 4] = 5.256417991472141;
            A[4 + ns * 5] = 0.0;
            A[4 + ns * 6] = 0.0;
            A[4 + ns * 7] = 0.0;
            A[4 + ns * 8] = 0.0;
            A[4 + ns * 9] = 0.0;
            A[5 + ns * 0] = -25.355197261299487;
            A[5 + ns * 1] = -12.800204401367228;
            A[5 + ns * 2] = 14.269560295130885;
            A[5 + ns * 3] = 5.741001773452404;
            A[5 + ns * 4] = -0.3063973123780201;
            A[5 + ns * 5] = 2.1638273472763068;
            A[5 + ns * 6] = 0.0;
            A[5 + ns * 7] = 0.0;
            A[5 + ns * 8] = 0.0;
            A[5 + ns * 9] = 0.0;
            A[6 + ns * 0] = -12.738301449222966;
            A[6 + ns * 1] = -9.41674517639761;
            A[6 + ns * 2] = 0.335534707546381;
            A[6 + ns * 3] = 1.2120119969408916;
            A[6 + ns * 4] = -0.9136261987143;
            A[6 + ns * 5] = 0.4795370048765861;
            A[6 + ns * 6] = 9.661323592475965;
            A[6 + ns * 7] = 0.0;
            A[6 + ns * 8] = 0.0;
            A[6 + ns * 9] = 0.0;
            A[7 + ns * 0] = 1.822736596460327;
            A[7 + ns * 1] = 20.013406903261796;
            A[7 + ns * 2] = -38.99600951197574;
            A[7 + ns * 3] = 0.8349487481118018;
            A[7 + ns * 4] = -2.257500501637924;
            A[7 + ns * 5] = 2.2588179888003688;
            A[7 + ns * 6] = 22.196465185327206;
            A[7 + ns * 7] = 6.304994521225663;
            A[7 + ns * 8] = 0.0;
            A[7 + ns * 9] = 0.0;
            A[8 + ns * 0] = 4.960016894639058;
            A[8 + ns * 1] = 18.422697972257026;
            A[8 + ns * 2] = -22.243752684294428;
            A[8 + ns * 3] = 1.1535942712254306;
            A[8 + ns * 4] = -1.358706661158572;
            A[8 + ns * 5] = 0.7509217045361036;
            A[8 + ns * 6] = 5.111417536175991;
            A[8 + ns * 7] = 2.2690279138428915;
            A[8 + ns * 8] = 4.887525081108087;
            A[8 + ns * 9] = 0.0;
            A[9 + ns * 0] = 7.882060894650045;
            A[9 + ns * 1] = -20.676096095585407;
            A[9 + ns * 2] = 44.78265443908609;
            A[9 + ns * 3] = -2.9956858673310607;
            A[9 + ns * 4] = 30.714288137116842;
            A[9 + ns * 5] = -4.607081553766025;
            A[9 + ns * 6] = -29.816690195466876;
            A[9 + ns * 7] = -16.47234088821905;
            A[9 + ns * 8] = -22.5076978125548;
            A[9 + ns * 9] = 4.984266587289817;
            // b integral
            b[0] = 0.0;
            b[1] = 0.0;
            b[2] = 0.0;
            b[3] = 0.0;
            b[4] = 0.0;
            b[5] = 0.0;
            b[6] = 0.0;
            b[7] = 0.0;
            b[8] = 0.0;
            b[9] = 1.0;
            // c integral
            c[0] = 0.233704632125264;
            c[1] = 0.126124356831163;
            c[2] = 0.491680304562551;
            c[3] = 0.0122571595694357;
            c[4] = 0.935943315417379;
            c[5] = 0.804308619035087;
            c[6] = 0.564545415166426;
            c[7] = 0.789591288449452;
            c[8] = 0.906473723240191;
            c[9] = 1.0;
            break;
        }
        default:
        {
            printf("\n error: DIRK: num_stages = %d not available. Only number of stages = {6} implemented!\n",ns);
            exit(1);
        }
    }
}


void get_esdirk_butcher_tableau(int ns, double *A, double *b, double *c)
{
    switch (ns)
    {
        case 4:
        {
            // Stiffly Accurate, L-Stable, ESDIRK method of order 3, stage order 2, with 4 stages.
            // https://ntrs.nasa.gov/api/citations/20010075154/downloads/20010075154.pdf
            // A
            A[0 + ns * 0] = 0.0;
            A[0 + ns * 1] = 0.0;
            A[0 + ns * 2] = 0.0;
            A[0 + ns * 3] = 0.0;
            A[1 + ns * 0] = 0.435866521508459;
            A[1 + ns * 1] = 0.435866521508459;
            A[1 + ns * 2] = 0.0;
            A[1 + ns * 3] = 0.0;
            A[2 + ns * 0] = 0.2576482460664272;
            A[2 + ns * 1] = -0.09351476757488625;
            A[2 + ns * 2] = 0.435866521508459;
            A[2 + ns * 3] = 0.0;
            A[3 + ns * 0] = 0.18764102434672383;
            A[3 + ns * 1] = -0.595297473576955;
            A[3 + ns * 2] = 0.9717899277217721;
            A[3 + ns * 3] = 0.435866521508459;
            // b
            b[0] = 0.18764102434672383;
            b[1] = -0.595297473576955;
            b[2] = 0.9717899277217721;
            b[3] = 0.435866521508459;
            // c
            c[0] = 0.0;
            c[1] = 0.871733043016918;
            c[2] = 0.6;
            c[3] = 1.0;
            break;
        }
        default:
        {
            printf("\n error: ESDIRK: num_stages = %d not available. Only number of stages = {4} implemented!\n",ns);
            exit(1);
        }
    }
}


inline bool is_dirk_tableau(sim_collocation_type collocation_type)
{
    return collocation_type == SDIRK || collocation_type == ESDIRK || collocation_type == DIRK;
}

inline bool is_stiffly_accurate(sim_collocation_type collocation_type, int ns)
{
    return (collocation_type == SDIRK && (ns == 9)) || (collocation_type == DIRK && (ns == 10));
}
