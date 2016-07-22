#include "condensing.h"

extern data_struct data;

static void construct_g_row(real_t* g_, real_t* A_, real_t* c_) {
    int i, j;
    for ( j = 0; j < NX; j++ ) {
        for ( i = 0; i < NX; i++ ) {
            g_[i] += A_[j*NX+i]*g_[j-NX];
        }
        g_[j] += c_[j];
    }
}

static void propagateCU(real_t* C_, real_t* A_) {
    int i, j, k;
    for ( j = 0; j < NU; j++ ) {
        for ( i = 0; i < NX; i++ ) {
            for ( k = 0; k < NX; k++ ) {
                C_[j*NNN*NX+i] += A_[k*NX+i]*C_[j*NNN*NX-NX+k];
            }
        }
    }
}

static void propagateCX(real_t* C_, real_t* A_) {
    int i, j, k;
    for ( j = 0; j < NX; j++ ) {
        for ( i = 0; i < NX; i++ ) {
            for ( k = 0; k < NX; k++ ) {
                C_[j*NNN*NX+i] += A_[k*NX+i]*C_[j*NNN*NX-NX+k];
            }
        }
    }
}

static void computeWx(real_t* Q_, real_t* C_, real_t* A_) {
    int i, j , k;
    for ( i = 0; i < NX*NX; i++ ) data.W1_x[i] = data.W2_x[i];
    for ( j = 0; j < NX; j++ ) {
        for ( i = 0; i < NX; i++ ) {
            data.W2_x[j*NX+i] = 0.0;
            for ( k = 0; k < NX; k++ ) {
                data.W2_x[j*NX+i] += Q_[k*NX+i]*C_[j*NNN*NX+k];
                data.W2_x[j*NX+i] += A_[i*NX+k]*data.W1_x[j*NX+k];
            }
        }
    }
}

static void computeWu(real_t* Q_, real_t* C_, real_t* A_) {
    int i, j , k;
    for ( i = 0; i < NX*NU; i++ ) data.W1_u[i] = data.W2_u[i];
    for ( j = 0; j < NU; j++ ) {
        for ( i = 0; i < NX; i++ ) {
            data.W2_u[j*NX+i] = 0.0;
            for ( k = 0; k < NX; k++ ) {
                data.W2_u[j*NX+i] += Q_[k*NX+i]*C_[j*NNN*NX+k];
                data.W2_u[j*NX+i] += A_[i*NX+k]*data.W1_u[j*NX+k];
            }
        }
    }
}

static void computeH_offDX(real_t* Hc_, real_t* S_, real_t* C_, real_t* B_) {
    int i, j , k;
    for ( j = 0; j < NX; j++ ) {
        for ( i = 0; i < NU; i++ ) {
            Hc_[j*NVC+i] = 0.0;
            for ( k = 0; k < NX; k++ ) {
                Hc_[j*NVC+i] += S_[k*NU+i]*C_[j*NNN*NX+k];
                Hc_[j*NVC+i] += B_[i*NX+k]*data.W2_x[j*NX+k];
            }
        }
    }
}

static void computeH_offDU(real_t* Hc_, real_t* S_, real_t* C_, real_t* B_) {
    int i, j , k;
    for ( j = 0; j < NU; j++ ) {
        for ( i = 0; i < NU; i++ ) {
            Hc_[j*NVC+i] = 0.0;
            for ( k = 0; k < NX; k++ ) {
                Hc_[j*NVC+i] += S_[k*NU+i]*C_[j*NNN*NX+k];
                Hc_[j*NVC+i] += B_[i*NX+k]*data.W2_u[j*NX+k];
            }
        }
    }
}

static void computeH_DX() {
    int i, j , k;
    for ( j = 0; j < NX; j++ ) {
        for ( i = 0; i < NU; i++ ) {
            data.Hc[j*NVC+NX+i] = data.S[j*NU+i];
            for ( k = 0; k < NX; k++ ) {
                data.Hc[j*NVC+NX+i] += data.B[i*NX+k]*data.W2_x[j*NX+k];
            }
        }
    }
    for ( j = 0; j < NX; j++ ) {
        for ( i = 0; i < NX; i++ ) {
            data.Hc[j*NVC+i] = data.Q[j*NX+i];
            for ( k = 0; k < NX; k++ ) {
                data.Hc[j*NVC+i] += data.A[i*NX+k]*data.W2_x[j*NX+k];
            }
        }
    }
}

static void computeH_DU(real_t* Hc_, real_t* R_, real_t* B_) {
    int i, j , k;
    for ( j = 0; j < NU; j++ ) {
        for ( i = 0; i < NU; i++ ) {
            Hc_[j*NVC+i] = R_[j*NU+i];
            for ( k = 0; k < NX; k++ ) {
                Hc_[j*NVC+i] += B_[i*NX+k]*data.W2_u[j*NX+k];
            }
        }
    }
}


static void computeWg(real_t* q_, real_t* Q_, real_t* d_, real_t* A_) {
    int i, j;
    for ( i = 0; i < NX; i++ ) {
        data.w1[i] = data.w2[i];
        data.w2[i] = q_[i];
    }
    for ( j = 0; j < NX; j++ ) {
        for ( i = 0; i < NX; i++ ) {
            data.w2[i] += Q_[j*NX+i]*d_[j];
            data.w2[i] += A_[i*NX+j]*data.w1[j];
        }
    }
}

static void computeG_off(real_t* gc_, real_t* r_, real_t* S_, real_t* d_, real_t* B_) {
    int i, j;
    for ( i = 0; i < NU; i++ ) gc_[i] = r_[i];
    for ( j = 0; j < NX; j++ ) {
        for ( i = 0; i < NU; i++ ) {
            gc_[i] += S_[j*NU+i]*d_[j];
            gc_[i] += B_[i*NX+j]*data.w2[j];
        }
    }
}

static void computeG_fixed_initial_state(real_t* x0) {
    int i, j;
    real_t Sx0[NU] = {0};
    for ( j = 0; j < NU; j++ ) {
        for ( i = 0; i < NX; i++ ) Sx0[j] += data.S[i+j*NX]*x0[i];
    }
    /* first control */
    for ( i = 0; i < NU; i++ ) data.gc[i] = data.g[NX+i];
    for ( j = 0; j < NX; j++ ) {
        for ( i = 0; i < NU; i++ ) {
            data.gc[i] += data.B[i*NX+j]*data.w2[j] + Sx0[i];
        }
    }
}

void construct_G_column(int_t j, int_t offset) {
    for ( int_t c = 0; c < NU; c++ ) {
        for ( int_t r = 0; r < NX; r++ ) {
            data.G[(offset+j*NU+c)*NNN*NX+j*NX+r] = data.B[j*NX*NU+c*NX+r];
        }
    }
    for ( int_t i = j+1; i < NNN; i++ ) {
        /* G_{i,0} <- A_{i-1}*G_{i-1,0}  */
        propagateCU(&data.G[(offset+j*NU)*NNN*NX+i*NX], &data.A[i*NX*NX]);
    }
}

void set_bounds() {
    /* correct lbA and ubA with d values */
    int_t NA = 0;
    for ( int_t i = 0; i < NNN; i++ ) {
        for ( int_t j = 0; j < NX; j++ ) {
            // State simple bounds
            data.lbA[i*(NX+NA)+j] = data.lb[(i+1)*(NX+NU)+j] - data.d[i*NX+j];
            data.ubA[i*(NX+NA)+j] = data.ub[(i+1)*(NX+NU)+j] - data.d[i*NX+j];
        }
        // Sample code for polytopic constraints
        // for ( int_t j = 0; j < NA; j++ ) {
        //     // State simple bounds
        //     data.lbA[NX+i*(NX+NA)+j] = data.lbA[(i+1)*(NX+NU)+j] - D*data.d[i*NX+j];
        //     data.ubA[NX+i*(NX+NA)+j] = data.ubA[(i+1)*(NX+NU)+j] - D*data.d[i*NX+j];
        // }
    }
}

void propagate_controls_to_H(int_t offset) {
    int_t i, j;
    /* propagate controls: */
    for ( j = 0; j < NNN; j++ ) {
        for ( i = 0; i < NX*NU; i++ ) data.W2_u[i] = 0.0;
        /* A is unused here because W is zero */
        computeWu(&data.Q[NNN*NX*NX], &data.G[(offset+j*NU)*NNN*NX+(NNN-1)*NX], &data.A[0]);
        for ( i = NNN-1; i > j; i-- ) {
            computeH_offDU(&data.Hc[(offset+j*NU)*NVC+offset+i*NU], &data.S[i*NX*NU], \
                    &data.G[(offset+j*NU)*NNN*NX+(i-1)*NX], &data.B[i*NX*NU]);
            computeWu(&data.Q[i*NX*NX], &data.G[(offset+j*NU)*NNN*NX+(i-1)*NX], &data.A[i*NX*NX]);
        }
        computeH_DU(&data.Hc[(offset+j*NU)*NVC+offset+j*NU], &data.R[j*NU*NU], &data.B[j*NX*NU]);
    }
}

void propagate_to_G(int_t offset) {
    computeWg(&data.g[NNN*(NX+NU)], &data.Q[NNN*NX*NX], &data.d[(NNN-1)*NX], &data.A[0]);
    for ( int_t i = NNN-1; i > 0; i-- ) {
        computeG_off(&data.gc[offset+i*NU], &data.g[i*(NX+NU)+NX], &data.S[i*NX*NU], \
                    &data.d[(i-1)*NX], &data.B[i*NX*NU]);
        computeWg(&data.g[i*(NX+NU)], &data.Q[i*NX*NX], &data.d[(i-1)*NX], &data.A[i*NX*NX]);
    }
}

void calculate_G(int_t offset) {
    for ( int_t j = 0; j < NNN; j++ ) {
        construct_G_column(j, offset);
    }
}

void calculate_g() {
    for ( int_t k = 0; k < NX; k++ ) data.d[k] = data.b[k];
    for ( int_t j = 1; j < NNN; j++ ) {
        construct_g_row(&data.d[j*NX], &data.A[j*NX*NX], &data.b[j*NX]);
        // TODO(robin): for ( i = 0; i < NX; i++ ) data.d[k] += data.A[i*NX+k]*x0[i];
    }
}

void condensingN2_fixed_initial_state() {
    int_t i, j;
    int_t offset = 0;
    real_t x0[NX] = {0};
    // Fix x0 to lower bound (== upper bound)
    for ( i = 0; i < NX; i++ ) x0[i] = data.lb[i];
    for ( i = 0; i < NNN; i++ ) {
        for ( j = 0; j < NU; j++ ) {
            data.lbU[i*NU+j] = data.lb[i*(NX+NU)+NX+j];
            data.ubU[i*NU+j] = data.ub[i*(NX+NU)+NX+j];
        }
    }

    calculate_G(offset);
    calculate_g();

    set_bounds();

    /* !! Hessian propagation !! */
    propagate_controls_to_H(offset);

    /* !! gradient propagation !! */
    /* A is unused for this operation because the W's are zero */
    propagate_to_G(offset);
    computeG_fixed_initial_state(&(x0[0]));
}

void propagate_x0_to_G() {
    int_t i, j;
    /* propagate x0: */
    for ( j = 0; j < NX; j++ ) {
        for ( i = 0; i < NX; i++ ) {
            data.G[j*NNN*NX+i] = data.A[j*NX+i]; /* A_0 */
        }
    }
    for ( i = 1; i < NNN; i++ ) {
        propagateCX(&data.G[i*NX], &data.A[i*NX*NX]); /* G_{i,0} <- A_{i-1}*G_{i-1,0}  */
    }
}

void propagate_x0_to_H() {
    /* propagate x0: */
    /* A is unused for this operation because the W's are zero */
    computeWx(&data.Q[NNN*NX*NX], &data.G[(NNN-1)*NX], &data.A[0]);
    for ( int_t i = NNN-1; i > 0; i-- ) {
        computeH_offDX(&data.Hc[NX+i*NU], &data.S[i*NX*NU], &data.G[(i-1)*NX], &data.B[i*NX*NU]);
        computeWx(&data.Q[i*NX*NX], &data.G[(i-1)*NX], &data.A[i*NX*NX]);
    }
    computeH_DX();
}

static void compute_G_free_initial_state() {
    int i, j;
    /* x0 */
    for ( i = 0; i < NX; i++ ) data.gc[i] = data.g[i];
    for ( j = 0; j < NX; j++ ) {
        for ( i = 0; i < NX; i++ ) {
            data.gc[i] += data.A[i*NX+j]*data.w2[j];
        }
    }
    /* first control */
    for ( i = 0; i < NU; i++ ) data.gc[NX+i] = data.g[NX+i];
    for ( j = 0; j < NX; j++ ) {
        for ( i = 0; i < NU; i++ ) {
            data.gc[NX+i] += data.B[i*NX+j]*data.w2[j];
        }
    }
}

void condensingN2_free_initial_state() {
    int_t i, j, offset;

    /* Copy bound values */
    offset = NX;
    for ( i = 0; i < NX; i++ ) data.lbU[i] = data.lb[i];
    for ( i = 0; i < NX; i++ ) data.ubU[i] = data.ub[i];
    for ( i = 0; i < NNN; i++ ) {
        for ( j = 0; j < NU; j++ ) {
            data.lbU[NX+i*NU+j] = data.lb[i*(NX+NU)+NX+j];
            data.ubU[NX+i*NU+j] = data.ub[i*(NX+NU)+NX+j];
        }
    }

    /* Create matrix G, NOTE: this is a sparse matrix which is currently stored as a dense one! */
    propagate_x0_to_G();
    /* propagate controls: */
    for ( j = 0; j < NNN; j++ ) {
        construct_G_column(j, offset);
        if ( j > 0 ) {
            construct_g_row(&data.d[j*NX], &data.A[j*NX*NX], &data.b[j*NX]);
        } else {
            for ( i = 0; i < NX; i++ ) data.d[i] = data.b[i];
        }
    }

    set_bounds();

    /* !! Hessian propagation !! */
    propagate_x0_to_H();
    propagate_controls_to_H(offset);

    /* !! gradient propagation !! */
    /* A is unused for this operation because the W's are zero */
    propagate_to_G(offset);
    compute_G_free_initial_state();
}
