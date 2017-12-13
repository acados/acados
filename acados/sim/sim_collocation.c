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

#include "acados/sim/sim_collocation.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "blasfeo/include/blasfeo_target.h"
#include "blasfeo/include/blasfeo_common.h"
#include "blasfeo/include/blasfeo_d_aux_ext_dep.h"
#include "blasfeo/include/blasfeo_d_blas.h"
#include "blasfeo/include/blasfeo_d_kernel.h"
#include "blasfeo/include/blasfeo_i_aux_ext_dep.h"

#include "acados/utils/print.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

void get_Gauss_nodes(const int_t num_stages, real_t *nodes) {
    //    if ( num_stages == 1 ) {         // GL2
    //        nodes[0] = 1.0/2.0;
    //    } else if ( num_stages == 2 ) {  // GL4
    //        memcpy(nodes,
    //                ((real_t[]) {1.0/2.0+sqrt(3.0)/6.0,
    //                1.0/2.0-sqrt(3.0)/6.0}), sizeof(*nodes) * (num_stages));
    //    } else if ( num_stages == 3 ) {  // GL6
    //        memcpy(nodes,
    //                ((real_t[]) {1.0/2.0-sqrt(15.0)/10.0, 1.0/2.0,
    //                1.0/2.0+sqrt(15.0)/10.0}), sizeof(*nodes) * (num_stages));
    //    } else {
    //        // throw error somehow?
    //    }
    uint N = num_stages - 1;
    uint N1 = N + 1;
    uint N2 = N + 2;
    real_t *x_init;
    real_t *y;
    real_t *y_prev;
    real_t *lgvm;      // Legendre-Gauss Vandermonde Matrix
    real_t *der_lgvm;  // derivative of LGVM
    real_t err = 1;
    real_t eps = 2e-16;

    d_zeros(&x_init, N1, 1);
    d_zeros(&y, N1, 1);
    d_zeros(&y_prev, N1, 1);
    d_zeros(&lgvm, N1, N2);
    d_zeros(&der_lgvm, N1, 1);

    real_t a = 0.0;
    real_t b = 1.0;  // code for collocation interval [a,b]

    for (uint i = 0; i < N1; i++) {
        if (N > 0) {
            x_init[i] = -1 + i * 2.0 / N;
        } else {
            x_init[i] = -1;
        }
        y[i] = cos((2 * i + 1) * M_PI / (2 * N + 2)) +
               (0.27 / N1) * sin(M_PI * x_init[i] * N / N2);
        y_prev[i] = 2.0;
    }

    while (err > eps) {  // iterate until step sufficiently small
        for (uint i = 0; i < N1; i++) lgvm[i] = 1.0;
        for (uint i = 0; i < N1; i++) lgvm[N1 + i] = y[i];
        for (uint k = 2; k < N2; k++) {
            for (uint i = 0; i < N1; i++)
                lgvm[k * N1 + i] =
                    ((2 * k - 1) * y[i] * lgvm[(k - 1) * N1 + i] -
                     (k - 1) * lgvm[(k - 2) * N1 + i]) /
                    k;
        }
        for (uint i = 0; i < N1; i++)
            der_lgvm[i] = N2 * (lgvm[N * N1 + i] - y[i] * lgvm[N1 * N1 + i]) /
                          (1 - pow(y[i], 2));
        for (uint i = 0; i < N1; i++) y_prev[i] = y[i];

        for (uint i = 0; i < N1; i++)
            y[i] = y_prev[i] - lgvm[N1 * N1 + i] / der_lgvm[i];  // Newton step
        //        print_matrix("stdout", y, 1, num_stages);

        err = 0;
        for (uint i = 0; i < N1; i++) {
            if (err < fabs(y[i] - y_prev[i])) err = fabs(y[i] - y_prev[i]);
        }
    }
    for (uint i = 0; i < N1; i++)
        nodes[i] = (a * (1 - y[i]) + b * (1 + y[i])) / 2;

    free(x_init);
    free(y);
    free(y_prev);
    free(lgvm);
    free(der_lgvm);
}

void read_Gauss_simplified(const int_t num_stages, Newton_scheme *scheme) {
    real_t *D;
    D = (real_t *) calloc(2 * num_stages, sizeof(real_t));
    real_t *T;
    T = (real_t *) calloc(num_stages*num_stages, sizeof(real_t));
    char simplified[MAX_STR_LEN];
    int_t *perm;
    real_t *T_inv;

    snprintf(simplified, sizeof(simplified), "simplified/GL%d_simpl_%s.txt",
             2 * num_stages, "D");
    read_matrix(simplified, D, num_stages, 2);

    snprintf(simplified, sizeof(simplified), "simplified/GL%d_simpl_%s.txt",
             2 * num_stages, "T");
    read_matrix(simplified, T, num_stages, num_stages);

    //    print_matrix("stdout", D, num_stages, 2);
    //    print_matrix("stdout", T, num_stages, num_stages);

    scheme->single = false;
    scheme->low_tria = 0;
    for (int_t i = 0; i < num_stages; i++) {
        scheme->eig[i] = D[i];
    }
    for (int_t i = 0; i < num_stages * num_stages; i++) {
        scheme->transf2[i] = T[i];
    }
    // transf1_T:
    for (int_t i = 0; i < num_stages; i++) {
        if ((i + 1) < num_stages) {  // complex conjugate pair of eigenvalues
            for (int_t i1 = i; i1 < i + 2; i1++) {
                for (int_t i2 = 0; i2 < num_stages; i2++) {
                    scheme->transf1_T[i2 * num_stages + i1] = 0.0;
                    for (int_t i3 = 0; i3 < 2; i3++) {
                        scheme->transf1_T[i2 * num_stages + i1] +=
                            D[(i1 - i) * num_stages + (i + i3)] *
                            T[(i + i3) * num_stages + i2];
                    }
                }
            }
            i++;
        } else {  // real eigenvalue
            for (int_t i2 = 0; i2 < num_stages; i2++) {
                scheme->transf1_T[i2 * num_stages + i] =
                    D[i] * T[i * num_stages + i2];
            }
        }
    }

    int_zeros(&perm, num_stages, 1);
    d_zeros(&T_inv, num_stages, num_stages);
    for (int_t i = 0; i < num_stages; i++) {
        T_inv[i * (num_stages + 1)] = 1.0;
    }
    LU_system_solve(T, T_inv, perm, num_stages, num_stages);

    // transf1:
    for (int_t i = 0; i < num_stages; i++) {
        if ((i + 1) < num_stages) {  // complex conjugate pair of eigenvalues
            for (int_t i1 = i; i1 < i + 2; i1++) {
                for (int_t i2 = 0; i2 < num_stages; i2++) {
                    scheme->transf1[i2 * num_stages + i1] = 0.0;
                    for (int_t i3 = 0; i3 < 2; i3++) {
                        scheme->transf1[i2 * num_stages + i1] +=
                            D[i3 * num_stages + i1] *
                            T_inv[i2 * num_stages + i + i3];
                    }
                }
            }
            i++;
        } else {  // real eigenvalue
            for (int_t i2 = 0; i2 < num_stages; i2++) {
                scheme->transf1[i2 * num_stages + i] =
                    D[i] * T_inv[i2 * num_stages + i];
            }
        }
    }
    // transf2_T:
    for (int_t i = 0; i < num_stages; i++) {
        for (int_t i2 = 0; i2 < num_stages; i2++) {
            scheme->transf2_T[i2 * num_stages + i] = T_inv[i * num_stages + i2];
        }
    }

    //    print_matrix("stdout", scheme->transf1, num_stages, num_stages);
    //    print_matrix("stdout", T_inv, 1, 1);
    //    print_matrix("stdout", scheme->transf2, num_stages, num_stages);
    //    print_matrix("stdout", T_inv, 1, 1);
    //        print_matrix("stdout", scheme->transf1_T, num_stages, num_stages);
    //        print_matrix("stdout", T_inv, 1, 1);
    //        print_matrix("stdout", scheme->transf2_T, num_stages, num_stages);
    //        print_matrix("stdout", T_inv, 1, 1);

    free(perm);
    free(T_inv);
}

void create_Butcher_table(const int_t num_stages, const real_t *nodes,
                          real_t *b, real_t *A) {
    int_t i, j, k;
    real_t *can_vm;
    real_t *rhs;
    int_t *perm;

    d_zeros(&can_vm, num_stages, num_stages);
    d_zeros(&rhs, num_stages, num_stages);
    int_zeros(&perm, num_stages, 1);

    for (j = 0; j < num_stages; j++) {
        for (i = 0; i < num_stages; i++)
            can_vm[i + j * num_stages] = pow(nodes[i], j);
    }

    for (i = 0; i < num_stages * num_stages; i++) rhs[i] = 0.0;
    for (i = 0; i < num_stages; i++) rhs[i * (num_stages + 1)] = 1.0;
    //        print_matrix("stdout", rhs, num_stages, num_stages);

    LU_system_solve(can_vm, rhs, perm, num_stages, num_stages);
    //        print_matrix("stdout", rhs, num_stages, num_stages);

    for (k = 0; k < num_stages; k++) {
        for (i = 0; i < num_stages; i++) {
            A[i * num_stages + k] = 0.0;
            for (j = 0; j < num_stages; j++) {
                A[i * num_stages + k] =
                    A[i * num_stages + k] +
                    pow(nodes[k], j + 1) / (j + 1) * rhs[i * num_stages + j];
            }
        }
    }

    for (i = 0; i < num_stages; i++) {
        b[i] = 0.0;
        for (j = 0; j < num_stages; j++) {
            b[i] = b[i] + 1.0 / (j + 1) * rhs[i * num_stages + j];
        }
    }

    free(can_vm);
    free(rhs);
    free(perm);
}

// TODO(rien): replace these LU codes with blasfeo
real_t LU_system_solve(real_t *const A, real_t *const b, int *const perm,
                       int dim, int dim2) {
    real_t det;
    real_t swap;
    real_t valueMax;

    int i, j, k;
    int index1;
    int indexMax;
    int intSwap;
    int DIM = dim;
    int DIM_RHS = dim2;
    real_t *bPerm;
    bPerm = (real_t *) calloc(DIM*DIM_RHS, sizeof(real_t));
    real_t tmp_var;

    for (i = 0; i < DIM; ++i) {
        perm[i] = i;
    }
    det = 1.0000000000000000e+00;
    for (i = 0; i < (DIM - 1); i++) {
        indexMax = i;
        valueMax = fabs(A[i * DIM + i]);
        for (j = (i + 1); j < DIM; j++) {
            swap = fabs(A[i * DIM + j]);
            if (swap > valueMax) {
                indexMax = j;
                valueMax = swap;
            }
        }
        if (indexMax > i) {
            for (k = 0; k < DIM; ++k) {
                swap = A[k * DIM + i];
                A[k * DIM + i] = A[k * DIM + indexMax];
                A[k * DIM + indexMax] = swap;
            }
            intSwap = perm[i];
            perm[i] = perm[indexMax];
            perm[indexMax] = intSwap;
        }
        for (j = i + 1; j < DIM; j++) {
            A[i * DIM + j] = -A[i * DIM + j] / A[i * DIM + i];
            for (k = i + 1; k < DIM; k++) {
                A[k * DIM + j] += A[i * DIM + j] * A[k * DIM + i];
            }
        }
    }

    for (i = 0; i < DIM; ++i) {
        index1 = perm[i];
        for (j = 0; j < DIM_RHS; ++j) {
            bPerm[j * DIM + i] = b[j * DIM + index1];
        }
    }
    for (j = 1; j < DIM; ++j) {
        for (i = 0; i < j; ++i) {
            tmp_var = A[i * DIM + j];
            for (k = 0; k < DIM_RHS; ++k) {
                bPerm[k * DIM + j] += tmp_var * bPerm[k * DIM + i];
            }
        }
    }
    for (i = DIM - 1; - 1 < i; --i) {
        for (j = DIM - 1; i < j; --j) {
            tmp_var = A[j * DIM + i];
            for (k = 0; k < DIM_RHS; ++k) {
                bPerm[k * DIM + i] -= tmp_var * bPerm[k * DIM + j];
            }
        }
        tmp_var = 1.0 / A[i * (DIM + 1)];
        for (k = 0; k < DIM_RHS; ++k) {
            bPerm[k * DIM + i] = tmp_var * bPerm[k * DIM + i];
        }
    }
    for (k = 0; k < DIM * DIM_RHS; ++k) {
        b[k] = bPerm[k];
    }

    return det;
}
