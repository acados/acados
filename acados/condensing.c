#include "condensing.h"

static void calculate_transition_vector(const condensing_in *in,
    condensing_workspace *ws, real_t *x0) {

    for (int_t k = 0; k < NX; k++) {
        ws->g[0][k] = in->b[0][k];
        for (int_t i = 0; i < NX; i++) {
            ws->g[0][k] += in->A[0][k+i*NX]*x0[i];
        }
    }
    for (int_t k = 1; k < NNN; k++) {
        for (int_t j = 0; j < NX; j++) {
            ws->g[k][j] = in->b[k][j];
            for (int_t i = 0; i < NX; i++) {
                ws->g[k][j] += in->A[k][j+i*NX]*ws->g[k-1][i];
            }
        }
    }
}

static void diag_trans_blk(const real_t *B, real_t *G) {
    for (int_t col = 0; col < NU; col++) {
        for (int_t row = 0; row < NX; row++) {
            G[col*NX+row] = B[col*NX+row];
        }
    }
}

static void offdiag_trans_blk(real_t *A, real_t *G_prev, real_t *G) {
    for (int_t j = 0; j < NU; j++) {
        for (int_t i = 0; i < NX; i++) {
            for (int_t k = 0; k < NX; k++) {
                G[j*NX+i] += A[k*NX+i]*G_prev[j*NX+k];
            }
        }
    }
}

static void calculate_transition_matrix(condensing_in *in, condensing_workspace *ws) {
    for (int_t j = 0; j < NNN; j++) {
        diag_trans_blk(in->B[j], ws->G[j][j]);
        for (int_t i = j+1; i < NNN; i++) {
            offdiag_trans_blk(in->A[i], ws->G[i-1][j], ws->G[i][j]);
        }
    }
}

static void update_w(condensing_workspace *ws, real_t *q, real_t *Q,
    real_t *g, real_t *A) {

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

static void calc_gradient_blk(condensing_workspace *ws, real_t *h, real_t *r,
    real_t *S, real_t *g, const real_t *B) {

    for (int_t i = 0; i < NU; i++) h[i] = r[i];
    for (int_t j = 0; j < NX; j++) {
        for (int_t i = 0; i < NU; i++) {
            h[i] += S[j*NU+i]*g[j];
            h[i] += B[i*NX+j]*ws->w2[j];
        }
    }
}

static void corr_grad_fixd_init_state(condensing_in *in, condensing_out *out,
    condensing_workspace *ws, real_t *x0) {

    real_t Sx0[NU] = {0};
    for (int_t j = 0; j < NU; j++) {
        for (int_t i = 0; i < NX; i++) Sx0[j] += in->S[0][i+j*NX]*x0[i];
    }
    /* first control */
    for (int_t i = 0; i < NU; i++) out->h[i] = in->r[0][i];
    for (int_t j = 0; j < NX; j++) {
        for (int_t i = 0; i < NU; i++) {
            out->h[i] += in->B[0][i*NX+j]*ws->w2[j] + Sx0[i];
        }
    }
}

static void calculate_gradient(condensing_in *in, condensing_out *out, condensing_workspace *ws,
    int_t offset, real_t *x0) {

    update_w(ws, in->q[NNN], in->Q[NNN], ws->g[NNN-1], in->A[0]);
    for (int_t i = NNN-1; i > 0; i--) {
        calc_gradient_blk(ws, &out->h[offset+i*NU], in->r[i], in->S[i],
                    ws->g[i-1], in->B[i]);
        update_w(ws, in->q[i], in->Q[i], ws->g[i-1], in->A[i]);
    }
    corr_grad_fixd_init_state(in, out, ws, x0);
}

static void update_W(condensing_workspace *ws, real_t *Q, real_t *G, real_t *A) {
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

static void offdiag_hess_blk(condensing_workspace *ws, real_t *H, real_t *S,
    real_t *G, const real_t *B) {

    for (int_t j = 0; j < NU; j++) {
        for (int_t i = 0; i < NU; i++) {
            H[j*NVC+i] = 0.0;
            for (int_t k = 0; k < NX; k++) {
                H[j*NVC+i] += S[k*NU+i]*G[j*NX+k];
                H[j*NVC+i] += B[i*NX+k]*ws->W2_u[j*NX+k];
            }
        }
    }
}

static void diag_hess_blk(condensing_workspace *ws, real_t *H, real_t *R, const real_t *B) {
    for (int_t j = 0; j < NU; j++) {
        for (int_t i = 0; i < NU; i++) {
            H[j*NVC+i] = R[j*NU+i];
            for (int_t k = 0; k < NX; k++) {
                H[j*NVC+i] += B[i*NX+k]*ws->W2_u[j*NX+k];
            }
        }
    }
}

static void calculate_hessian(condensing_in *in, condensing_out *out,
    condensing_workspace *ws, int_t offset) {

    for (int_t j = 0; j < NNN; j++) {
        for (int_t i = 0; i < NX*NU; i++) ws->W2_u[i] = 0.0;
        update_W(ws, in->Q[NNN], ws->G[NNN-1][j], in->A[0]);
        for (int_t i = NNN-1; i > j; i--) {
            offdiag_hess_blk(ws, &out->H[(offset+j*NU)*NVC+offset+i*NU], in->S[i],
                    ws->G[i-1][j], in->B[i]);
            update_W(ws, in->Q[i], ws->G[i-1][j], in->A[i]);
        }
        diag_hess_blk(ws, &out->H[(offset+j*NU)*NVC+offset+j*NU], in->R[j], in->B[j]);
        // Symmetrize H
        for (int_t i = 0; i < NVC; i++) {
            for (int_t j = 0; j < i; j++) {
                out->H[i*NVC+j] = out->H[j*NVC+i];
            }
        }
    }
}

static void calculate_simple_bounds(condensing_in *in, condensing_out *out) {
    for (int_t i = 0; i < NNN; i++) {
        for (int_t j = 0; j < in->nb[i]; j++) {
            if (NX <= in->idxb[i][j] && in->idxb[i][j] < NX+NU) {
                out->lb[i*NU-NX+in->idxb[i][j]] = in->lb[i][j];
                out->ub[i*NU-NX+in->idxb[i][j]] = in->ub[i][j];
            }
        }
    }
}

static void calculate_constraint_bounds(condensing_in *in, condensing_out *out,
    condensing_workspace *ws, real_t *x0) {

    // State simple bounds
    int ctr = in->nc[0];
    for (int_t i = 0; i < NNN; i++) {
        for (int_t j = 0; j < NX; j++) {
            out->lbA[ctr+j] = in->lb[i+1][j] - ws->g[i][j];
            out->ubA[ctr+j] = in->ub[i+1][j] - ws->g[i][j];
        }
        ctr += NX + in->nc[i+1];
    }
    // State polytopic constraints
    for (int_t i = 0; i < in->nc[0]; i++) {
        out->lbA[i] = in->lc[0][i];
        out->ubA[i] = in->uc[0][i];
        for (int_t j = 0; j < NX; j++) {
            out->lbA[i] -= in->Cx[0][j*in->nc[0]+i]*x0[j];
            out->ubA[i] -= in->Cx[0][j*in->nc[0]+i]*x0[j];
        }
    }
    ctr = NX + in->nc[0];
    for (int_t i = 1; i <= NNN; i++) {
        for (int_t j = 0; j < in->nc[i]; j++) {
            out->lbA[ctr+j] = in->lc[i][j];
            out->ubA[ctr+j] = in->uc[i][j];
            for (int_t k = 0; k < NX; k++) {
                out->lbA[ctr+j] -= in->Cx[i][k*in->nc[i]+j]*ws->g[i-1][k];
                out->ubA[ctr+j] -= in->Cx[i][k*in->nc[i]+j]*ws->g[i-1][k];
            }
        }
        ctr += NX + in->nc[i];
    }
}

static void offdiag_D_blk(int_t nc, real_t *Cx, real_t *G, real_t *D) {
    for (int_t j = 0; j < NU; j++) {
        for (int_t i = 0; i < nc; i++) {
            for (int_t k = 0; k < NX; k++) {
                D[j*nc+i] += Cx[k*nc+i]*G[j*NX+k];
            }
        }
    }
}

static void calculate_D(condensing_in *in, condensing_workspace *ws) {
    for (int_t k = 0; k < NNN; k++) {
        for (int_t j = 0; j < NU; j++) {
            for (int_t i = 0; i < in->nc[k]; i++) {
                ws->D[k][k][j*in->nc[k]+i] = in->Cu[k][j*in->nc[k]+i];
            }
        }
    }
    for (int_t i = 1; i < NNN+1; i++) {
        for (int_t j = 0; j < i; j++) {
            offdiag_D_blk(in->nc[i], in->Cx[i], ws->G[i-1][j], ws->D[i][j]);
        }
    }
}

static void calculate_constraint_matrix(condensing_in *in, condensing_out *out,
    condensing_workspace *ws) {

    int_t ldA = 0;
    for (int_t i = 0; i < NNN; i++) {
        ldA += in->nc[i] + NX;
    }
    ldA += in->nc[NNN];

    int_t ctr = 0, ctr2 = 0;
    for (int_t i = 0; i < NNN; i++) {
        ctr2 = ctr;
        for (int_t j = 0; j <= i; j++) {
            for (int_t k = 0; k < NU; k++) {
                for (int_t l = 0; l < NX; l++) {
                    out->A[ctr2+k*ldA+l] = ws->D[i][j][k*in->nc[i]+l];
                }
            }
            ctr2 += ldA*NU;
        }
        ctr += in->nc[i] + NX;
    }
    ctr = 0;
    ctr2 = 0;
    for (int_t i = 0; i < NNN; i++) {
        ctr += in->nc[i];
        ctr2 = ctr;
        for (int_t j = 0; j <= i; j++) {
            for (int_t k = 0; k < NU; k++) {
                for (int_t l = 0; l < NX; l++) {
                    out->A[ctr2+k*ldA+l] = ws->G[i][j][k*NX+l];
                }
            }
            ctr2 += ldA*NU;
        }
        ctr += NX;
    }
    for (int_t j = 0; j < NNN; j++) {
        for (int_t k = 0; k < NU; k++) {
            for (int_t l = 0; l < NX; l++) {
                out->A[ctr+j*ldA*NU+k*ldA+l] = ws->D[NNN][j][k*in->nc[NNN]+l];
            }
        }
    }
}

void condensingN2_fixed_initial_state(condensing_in *input, condensing_out *output,
    condensing_workspace *workspace) {

    real_t *x0 = input->lb[0];
    int_t offset = 0;

    calculate_transition_vector(input, workspace, x0);
    calculate_transition_matrix(input, workspace);

    calculate_gradient(input, output, workspace, offset, x0);
    calculate_hessian(input, output, workspace, offset);

    calculate_simple_bounds(input, output);
    calculate_D(input, workspace);
    calculate_constraint_bounds(input, output, workspace, x0);
    calculate_constraint_matrix(input, output, workspace);
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
