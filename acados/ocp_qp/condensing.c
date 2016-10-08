#include "acados/ocp_qp/condensing.h"

static void calculate_transition_vector(ocp_qp_in *in,
    condensing_workspace *ws, const real_t *x0) {

    for (int_t k = 0; k < in->nx[0]; k++) {
        ws->g[0][k] = in->b[0][k];
        for (int_t i = 0; i < in->nx[0]; i++) {
            ws->g[0][k] += in->A[0][k+i*in->nx[0]]*x0[i];
        }
    }
    for (int_t k = 1; k < in->N; k++) {
        for (int_t j = 0; j < in->nx[0]; j++) {
            ws->g[k][j] = in->b[k][j];
            for (int_t i = 0; i < in->nx[0]; i++) {
                ws->g[k][j] += in->A[k][j+i*in->nx[0]]*ws->g[k-1][i];
            }
        }
    }
}

static void diag_trans_blk(const real_t *B, real_t *G, int_t NX, int_t NU) {
    for (int_t col = 0; col < NU; col++) {
        for (int_t row = 0; row < NX; row++) {
            G[col*NX+row] = B[col*NX+row];
        }
    }
}

static void offdiag_trans_blk(const real_t *A, const real_t *G_prev,
        real_t *G, int_t NX, int_t NU) {
    for (int_t j = 0; j < NU; j++) {
        for (int_t i = 0; i < NX; i++) {
            for (int_t k = 0; k < NX; k++) {
                G[j*NX+i] += A[k*NX+i]*G_prev[j*NX+k];
            }
        }
    }
}

static void calculate_transition_matrix(ocp_qp_in *in, condensing_workspace *ws) {
    for (int_t j = 0; j < in->N; j++) {
        diag_trans_blk(in->B[j], ws->G[j][j], in->nx[j], in->nu[j]);
        for (int_t i = j+1; i < in->N; i++) {
            offdiag_trans_blk(in->A[i], ws->G[i-1][j], ws->G[i][j], in->nx[i], in->nu[i]);
        }
    }
}

static void update_w(condensing_workspace *ws, const real_t *q, const real_t *Q,
    const real_t *g, const real_t *A, int_t NX) {

    for (int_t i = 0; i < NX; i++) {
        ws->w1[i] = ws->w2[i];
        ws->w2[i] = q[i];
    }
    for (int_t j = 0; j < NX; j++) {
        for (int_t i = 0; i < NX; i++) {
            ws->w2[i] += Q[j*NX+i]*g[j];
            ws->w2[i] += A[i*NX+j]*ws->w1[j];
        }
    }
}

static void calc_gradient_blk(condensing_workspace *ws, real_t *h, const real_t *r,
    const real_t *S, const real_t *g, const real_t *B, int_t NX, int_t NU) {

    for (int_t i = 0; i < NU; i++) h[i] = r[i];
    for (int_t j = 0; j < NX; j++) {
        for (int_t i = 0; i < NU; i++) {
            h[i] += S[j*NU+i]*g[j];
            h[i] += B[i*NX+j]*ws->w2[j];
        }
    }
}

static void corr_grad_fixd_init_state(ocp_qp_in *in, condensing_out *out,
    condensing_workspace *ws, const real_t *x0) {

    for (int_t j = 0; j < in->nu[0]; j++) ws->Sx0[j] = 0;
    for (int_t j = 0; j < in->nu[0]; j++) {
        for (int_t i = 0; i < in->nx[0]; i++) ws->Sx0[j] += in->S[0][i+j*in->nx[0]]*x0[i];
    }
    /* first control */
    for (int_t i = 0; i < in->nu[0]; i++) out->h[i] = in->r[0][i];
    for (int_t j = 0; j < in->nx[0]; j++) {
        for (int_t i = 0; i < in->nu[0]; i++) {
            out->h[i] += in->B[0][i*in->nx[0]+j]*ws->w2[j] + ws->Sx0[i];
        }
    }
}

static void calculate_gradient(ocp_qp_in *in, condensing_out *out, condensing_workspace *ws,
    int_t offset, const real_t *x0) {

    update_w(ws, in->q[in->N], in->Q[in->N], ws->g[in->N-1], in->A[0], in->nx[0]);
    for (int_t i = in->N-1; i > 0; i--) {
        calc_gradient_blk(ws, &out->h[offset+i*in->nu[0]], in->r[i], in->S[i],
                    ws->g[i-1], in->B[i], in->nx[i], in->nu[i]);
        update_w(ws, in->q[i], in->Q[i], ws->g[i-1], in->A[i], in->nx[i]);
    }
    corr_grad_fixd_init_state(in, out, ws, x0);
}

static void update_W(condensing_workspace *ws, const real_t *Q, const real_t *G,
    const real_t *A, int_t NX, int_t NU) {

    for (int_t i = 0; i < NX*NU; i++) ws->W1_u[i] = ws->W2_u[i];
    for (int_t j = 0; j < NU; j++) {
        for (int_t i = 0; i < NX; i++) {
            ws->W2_u[j*NX+i] = 0.0;
            for (int_t k = 0; k < NX; k++) {
                ws->W2_u[j*NX+i] += Q[k*NX+i]*G[j*NX+k];
                ws->W2_u[j*NX+i] += A[i*NX+k]*ws->W1_u[j*NX+k];
            }
        }
    }
}

static void offdiag_hess_blk(condensing_workspace *ws, real_t *H,
    const real_t *S, const real_t *G, const real_t *B, int_t NX, int_t NU) {

    int_t ncv = ws->nconvars;
    for (int_t j = 0; j < NU; j++) {
        for (int_t i = 0; i < NU; i++) {
            H[j*ncv+i] = 0.0;
            for (int_t k = 0; k < NX; k++) {
                H[j*ncv+i] += S[k*NU+i]*G[j*NX+k];
                H[j*ncv+i] += B[i*NX+k]*ws->W2_u[j*NX+k];
            }
        }
    }
}

static void diag_hess_blk(condensing_workspace *ws, real_t *H,
    const real_t *R, const real_t *B, int_t NX, int_t NU) {

    int_t ncv = ws->nconvars;
    for (int_t j = 0; j < NU; j++) {
        for (int_t i = 0; i < NU; i++) {
            H[j*ncv+i] = R[j*NU+i];
            for (int_t k = 0; k < NX; k++) {
                H[j*ncv+i] += B[i*NX+k]*ws->W2_u[j*NX+k];
            }
        }
    }
}

static void calculate_hessian(ocp_qp_in *in, condensing_out *out,
    condensing_workspace *ws, int_t offset) {

    int_t ncv = ws->nconvars;
    for (int_t j = 0; j < in->N; j++) {
        for (int_t i = 0; i < in->nx[0]*in->nu[0]; i++) ws->W2_u[i] = 0.0;
        update_W(ws, in->Q[in->N], ws->G[in->N-1][j], in->A[0], in->nx[j], in->nu[j]);
        for (int_t i = in->N-1; i > j; i--) {
            offdiag_hess_blk(ws, &out->H[(offset+j*in->nu[0])*ncv+offset+i*in->nu[0]],
                in->S[i], ws->G[i-1][j], in->B[i], in->nx[i], in->nu[i]);
            update_W(ws, in->Q[i], ws->G[i-1][j], in->A[i], in->nx[i], in->nu[i]);
        }
        diag_hess_blk(ws, &out->H[(offset+j*in->nu[0])*ncv+offset+j*in->nu[0]],
                in->R[j], in->B[j], in->nx[j], in->nu[j]);
        // Symmetrize H
        for (int_t i = 0; i < ncv; i++) {
            for (int_t j = 0; j < i; j++) {
                out->H[i*ncv+j] = out->H[j*ncv+i];
            }
        }
    }
}

static void calculate_simple_bounds(ocp_qp_in *in, condensing_out *out) {
    for (int_t i = 0; i < in->N; i++) {
        for (int_t j = 0; j < in->nb[i]; j++) {
            if (in->nx[0] <= in->idxb[i][j] && in->idxb[i][j] < in->nx[0]+in->nu[0]) {
                out->lb[i*in->nu[0]-in->nx[0]+in->idxb[i][j]] = in->lb[i][j];
                out->ub[i*in->nu[0]-in->nx[0]+in->idxb[i][j]] = in->ub[i][j];
            }
        }
    }
}

static void calculate_constraint_bounds(ocp_qp_in *in, condensing_out *out,
    condensing_workspace *ws, const real_t *x0) {

    // State simple bounds
    int_t idx;
    int_t ctr = in->nc[0];
    for (int_t i = 0; i < in->N; i++) {
        for (int_t j = 0; j < ws->nstate_bounds[i+1]; j++) {
            idx = in->idxb[i+1][j];
            if (idx < in->nx[i+1]) {
                out->lbA[ctr+idx] = in->lb[i+1][j] - ws->g[i][idx];
                out->ubA[ctr+idx] = in->ub[i+1][j] - ws->g[i][idx];
            }
        }
        ctr += ws->nstate_bounds[i+1] + in->nc[i+1];
    }
    // State polytopic constraints
    for (int_t i = 0; i < in->nc[0]; i++) {
        out->lbA[i] = in->lc[0][i];
        out->ubA[i] = in->uc[0][i];
        for (int_t j = 0; j < in->nx[0]; j++) {
            out->lbA[i] -= in->Cx[0][j*in->nc[0]+i]*x0[j];
            out->ubA[i] -= in->Cx[0][j*in->nc[0]+i]*x0[j];
        }
    }
    ctr = 0;
    for (int_t i = 1; i <= in->N; i++) {
        ctr += ws->nstate_bounds[i] + in->nc[i-1];
        for (int_t j = 0; j < in->nc[i]; j++) {
            out->lbA[ctr+j] = in->lc[i][j];
            out->ubA[ctr+j] = in->uc[i][j];
            for (int_t k = 0; k < in->nx[0]; k++) {
                out->lbA[ctr+j] -= in->Cx[i][k*in->nc[i]+j]*ws->g[i-1][k];
                out->ubA[ctr+j] -= in->Cx[i][k*in->nc[i]+j]*ws->g[i-1][k];
            }
        }
    }
}

static void offdiag_D_blk(int_t nc, const real_t *Cx, const real_t *G,
        real_t *D, int_t NX, int_t NU) {
    for (int_t j = 0; j < NU; j++) {
        for (int_t i = 0; i < nc; i++) {
            for (int_t k = 0; k < NX; k++) {
                D[j*nc+i] += Cx[k*nc+i]*G[j*NX+k];
            }
        }
    }
}

static void calculate_D(ocp_qp_in *in, condensing_workspace *ws) {
    for (int_t k = 0; k < in->N; k++) {
        for (int_t j = 0; j < in->nu[0]; j++) {
            for (int_t i = 0; i < in->nc[k]; i++) {
                ws->D[k][k][j*in->nc[k]+i] = in->Cu[k][j*in->nc[k]+i];
            }
        }
    }
    for (int_t i = 1; i < in->N+1; i++) {
        if (in->nc[i] > 0) {
            for (int_t j = 0; j < i; j++) {
                offdiag_D_blk(in->nc[i], in->Cx[i], ws->G[i-1][j], ws->D[i][j],
                        in->nx[i], in->nu[i]);
            }
        }
    }
}

static void calculate_constraint_matrix(ocp_qp_in *in, condensing_out *out,
    condensing_workspace *ws) {

    int_t ldA = ws->nconstraints;

    if (ldA) {
    int_t ctr = 0, ctr2 = 0;
    for (int_t i = 0; i < in->N; i++) {
        ctr2 = ctr;
        for (int_t j = 0; j <= i; j++) {
            for (int_t k = 0; k < in->nu[0]; k++) {
                for (int_t l = 0; l < in->nc[i]; l++) {
                    out->A[ctr2+k*ldA+l] = ws->D[i][j][k*in->nc[i]+l];
                }
            }
            ctr2 += ldA*in->nu[0];
        }
        ctr += in->nc[i] + ws->nstate_bounds[i+1];
    }
    ctr = 0;
    ctr2 = 0;
    for (int_t i = 0; i < in->N; i++) {
        ctr += in->nc[i];
        ctr2 = ctr;
        for (int_t j = 0; j <= i; j++) {
            for (int_t k = 0; k < in->nu[0]; k++) {
                for (int_t l = 0; l < in->nx[0]; l++) {
                    out->A[ctr2+k*ldA+l] = ws->G[i][j][k*in->nx[0]+l];
                }
            }
            ctr2 += ldA*in->nu[0];
        }
        ctr += ws->nstate_bounds[i+1];
    }
    for (int_t j = 0; j < in->N; j++) {
        for (int_t k = 0; k < in->nu[0]; k++) {
            for (int_t l = 0; l < in->nx[0]; l++) {
                out->A[ctr+j*ldA*in->nu[0]+k*ldA+l] = ws->D[in->N][j][k*in->nc[in->N]+l];
            }
        }
    }
    }
}

void condensing_N2_fixed_initial_state(condensing_in *input, condensing_out *output,
    condensing_workspace *workspace) {

    const real_t *x0 = input->qp_input->lb[0];
    int_t offset = 0;

    calculate_transition_vector(input->qp_input, workspace, x0);
    calculate_transition_matrix(input->qp_input, workspace);

    calculate_gradient(input->qp_input, output, workspace, offset, x0);
    calculate_hessian(input->qp_input, output, workspace, offset);

    calculate_simple_bounds(input->qp_input, output);
    calculate_D(input->qp_input, workspace);
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
