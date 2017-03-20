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

#include "acados/ocp_qp/condensing.h"

#include "acados/ocp_qp/condensing_helper_functions.c"

void condensing_N2_fixed_initial_state(condensing_in *input, condensing_out *output,
    condensing_workspace *workspace) {

    const real_t *x0 = input->qp_input->lb[0];
    int_t offset = 0;

    calculate_transition_vector(input->qp_input, workspace, x0);
    calculate_transition_matrix(input->qp_input, workspace);

    calculate_gradient(input->qp_input, output, workspace, offset, x0);
    calculate_hessian(input->qp_input, output, workspace, offset);

    calculate_simple_bounds(input->qp_input, output);
    calculate_constraint_bounds(input->qp_input, output, workspace, x0);
    calculate_constraint_matrix(input->qp_input, output, workspace);
}

// static void propagateCX(real_t *C_, real_t *A_) {
//     int i, j, k;
//     for ( j = 0; j < NX; j++ ) {
//         for ( i = 0; i < NX; i++ ) {
//             for ( k = 0; k < NX; k++ ) {
//                 C_[j*NNN*NX+i] += A_[k*NX+i]*C_[j*NNN*NX-NX+k];
//             }
//         }
//     }
// }
//
// void propagate_x0_to_G(condensing_in **in, condensing_out *out,
//     condensing_workspace *ws) {
//     int_t i, j;
//     /* propagate x0: */
//     for ( j = 0; j < NX; j++ ) {
//         for ( i = 0; i < NX; i++ ) {
//             ws->G[j*NNN*NX+i] = data.A[j*NX+i]; /* A_0 */
//         }
//     }
//     for ( i = 1; i < NNN; i++ ) {
//         propagateCX(&ws->G[i*NX], &data.A[i*NX*NX]); /* G_{i,0} <- A_{i-1}*G_{i-1,0}  */
//     }
// }
//
// static void computeWx(real_t *Q_, real_t *C_, real_t *A_) {
//     int i, j , k;
//     for ( i = 0; i < NX*NX; i++ ) data.W1_x[i] = data.W2_x[i];
//     for ( j = 0; j < NX; j++ ) {
//         for ( i = 0; i < NX; i++ ) {
//             data.W2_x[j*NX+i] = 0.0;
//             for ( k = 0; k < NX; k++ ) {
//                 data.W2_x[j*NX+i] += Q_[k*NX+i]*C_[j*NNN*NX+k];
//                 data.W2_x[j*NX+i] += A_[i*NX+k]*data.W1_x[j*NX+k];
//             }
//         }
//     }
// }
//
// static void computeH_DX() {
//     int i, j , k;
//     for ( j = 0; j < NX; j++ ) {
//         for ( i = 0; i < NU; i++ ) {
//             data.Hc[j*NVC+NX+i] = data.S[j*NU+i];
//             for ( k = 0; k < NX; k++ ) {
//                 data.Hc[j*NVC+NX+i] += data.B[i*NX+k]*data.W2_x[j*NX+k];
//             }
//         }
//     }
//     for ( j = 0; j < NX; j++ ) {
//         for ( i = 0; i < NX; i++ ) {
//             data.Hc[j*NVC+i] = data.Q[j*NX+i];
//             for ( k = 0; k < NX; k++ ) {
//                 data.Hc[j*NVC+i] += data.A[i*NX+k]*data.W2_x[j*NX+k];
//             }
//         }
//     }
// }
//
// static void computeH_offDX(real_t *Hc_, real_t *S_, real_t *C_, real_t *B_) {
//     int i, j , k;
//     for ( j = 0; j < NX; j++ ) {
//         for ( i = 0; i < NU; i++ ) {
//             Hc_[j*NVC+i] = 0.0;
//             for ( k = 0; k < NX; k++ ) {
//                 Hc_[j*NVC+i] += S_[k*NU+i]*C_[j*NNN*NX+k];
//                 Hc_[j*NVC+i] += B_[i*NX+k]*data.W2_x[j*NX+k];
//             }
//         }
//     }
// }
//
// void propagate_x0_to_H(condensing_in *in, condensing_out *out,
//     condensing_workspace *ws) {
//     /* propagate x0: */
//     computeWx(&data.Q[NNN*NX*NX], &ws->G[(NNN-1)*NX], &data.A[0]);
//     for ( int_t i = NNN-1; i > 0; i-- ) {
//         computeH_offDX(&data.Hc[NX+i*NU], &data.S[i*NX*NU], &ws->G[(i-1)*NX], &data.B[i*NX*NU]);
//         computeWx(&data.Q[i*NX*NX], &ws->G[(i-1)*NX], &data.A[i*NX*NX]);
//     }
//     computeH_DX();
// }
//
// static void compute_G_free_initial_state() {
//     int i, j;
//     /* x0 */
//     for ( i = 0; i < NX; i++ ) data.gc[i] = data.f[i];
//     for ( j = 0; j < NX; j++ ) {
//         for ( i = 0; i < NX; i++ ) {
//             data.gc[i] += data.A[i*NX+j]*data.w2[j];
//         }
//     }
//     /* first control */
//     for ( i = 0; i < NU; i++ ) data.gc[NX+i] = data.f[NX+i];
//     for ( j = 0; j < NX; j++ ) {
//         for ( i = 0; i < NU; i++ ) {
//             data.gc[NX+i] += data.B[i*NX+j]*data.w2[j];
//         }
//     }
// }
//
// void condensingN2_free_initial_state(condensing_in *input, condensing_out *output,
//     condensing_workspace *workspace) {
//     int_t i, j, offset;
//
//     /* Copy bound values */
//     offset = NX;
//     for ( i = 0; i < NX; i++ ) data.lbU[i] = data.lb[i];
//     for ( i = 0; i < NX; i++ ) data.ubU[i] = data.ub[i];
//     for ( i = 0; i < NNN; i++ ) {
//         for ( j = 0; j < NU; j++ ) {
//             data.lbU[NX+i*NU+j] = data.lb[i*(NX+NU)+NX+j];
//             data.ubU[NX+i*NU+j] = data.ub[i*(NX+NU)+NX+j];
//         }
//     }
//
//     /* Create matrix G, NOTE: this is a sparse matrix which is currently
//     stored as a dense one! */
//     propagate_x0_to_G(input, output, workspace);
//     /* propagate controls: */
//     for ( j = 0; j < NNN; j++ ) {
//         // calculate_G_column(j, offset);
//         if ( j > 0 ) {
//             // calculate_g_row(&data.g[j*NX], &data.A[j*NX*NX], &data.b[j*NX]);
//         } else {
//             for ( i = 0; i < NX; i++ ) data.g[i] = data.b[i];
//         }
//     }
//
//     calculate_constraint_bounds(input, output, workspace, 0);
//
//     /* !! Hessian propagation !! */
//     propagate_x0_to_H(input, output, workspace);
//     calculate_hessian(input, output, workspace, offset);
//
//     /* !! gradient propagation !! */
//     calculate_gradient(input, output, workspace, offset, 0);
//     compute_G_free_initial_state();
// }
