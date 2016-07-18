#include "condensing.h"

extern data_struct data;

extern void printMatrix(const char* name, real_t* mat, unsigned nRows, unsigned nCols);

static void propagatec(real_t* d_, real_t* A_, real_t* b_) {
    int i, j;
    for ( j = 0; j < NX; j++ ) {
        for ( i = 0; i < NX; i++ ) {
            d_[i] += A_[j*NX+i]*d_[j-NX];
        }
        d_[j] += b_[j];
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

#if FIXED_INITIAL_STATE == 0
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
#endif

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

#if FIXED_INITIAL_STATE == 0
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
#endif

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

#if FIXED_INITIAL_STATE == 0
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
#endif

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

#if FIXED_INITIAL_STATE == 0
static void computeG() {
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
#endif

#if FIXED_INITIAL_STATE == 1
static void computeG_fixed_initial_state( real_t* Sx0 ) {
    int i, j;
    /* first control */
    for ( i = 0; i < NU; i++ ) data.gc[i] = data.g[NX+i];
    for ( j = 0; j < NX; j++ ) {
        for ( i = 0; i < NU; i++ ) {
            data.gc[i] += data.B[i*NX+j]*data.w2[j] + Sx0[i];
        }
    }
}
#endif

/* This contains an implementation of our block condensing algorithm. */
void block_condensing() {
    int_t i, j, r, c, offset;

    /* Copy bound values */
    #if FIXED_INITIAL_STATE == 1
    int_t k;
    offset = 0;
    real_t x0[NX] = {0};
    // Read x0 from lower bound (= upper bound)
    for ( i = 0; i < NX; i++ ) x0[i] = data.lb[i];
    for ( i = 0; i < NNN; i++ ) {
        for ( j = 0; j < NU; j++ ) {
            data.lbU[i*NU+j] = data.lb[i*(NX+NU)+NX+j];
            data.ubU[i*NU+j] = data.ub[i*(NX+NU)+NX+j];
        }
    }
    #else
    offset = NX;
    for ( i = 0; i < NX; i++ ) data.lbU[i] = data.lb[i];
    for ( i = 0; i < NX; i++ ) data.ubU[i] = data.ub[i];
    for ( i = 0; i < NNN; i++ ) {
        for ( j = 0; j < NU; j++ ) {
            data.lbU[NX+i*NU+j] = data.lb[i*(NX+NU)+NX+j];
            data.ubU[NX+i*NU+j] = data.ub[i*(NX+NU)+NX+j];
        }
    }
    #endif

    /* Create matrix G, NOTE: this is a sparse matrix which is currently stored as a dense one! */
    #if FIXED_INITIAL_STATE == 0
    /* propagate x0: */
    for ( j = 0; j < NX; j++ ) {
        for ( i = 0; i < NX; i++ ) {
            data.C[j*NNN*NX+i] = data.A[j*NX+i]; /* A_0 */
        }
    }
    for ( i = 1; i < NNN; i++ ) {
        propagateCX(&data.C[i*NX], &data.A[i*NX*NX]); /* G_{i,0} <- A_{i-1}*G_{i-1,0}  */
    }
    #endif

    /* propagate controls: */
    for ( j = 0; j < NNN; j++ ) {
        for ( c = 0; c < NU; c++ ) {
            for ( r = 0; r < NX; r++ ) {
                data.C[(offset+j*NU+c)*NNN*NX+j*NX+r] = data.B[j*NX*NU+c*NX+r]; /* B_j */
            }
        }
        for ( i = j+1; i < NNN; i++ ) {
            /* G_{i,0} <- A_{i-1}*G_{i-1,0}  */
            propagateCU(&data.C[(offset+j*NU)*NNN*NX+i*NX], &data.A[i*NX*NX]);
        }
        #if FIXED_INITIAL_STATE == 1
        if ( j > 0 ) {
            propagatec(&data.d[j*NX], &data.A[j*NX*NX], &data.b[j*NX]);
        } else {
            for ( k = 0; k < NX; k++ ) {
                data.d[k] = data.b[k];
                for ( i = 0; i < NX; i++ ) data.d[k] += data.A[i*NX+k]*x0[i];
            }
        }
        #else
        if ( j > 0 ) {
            propagatec(&data.d[j*NX], &data.A[j*NX*NX], &data.b[j*NX]);
        } else {
            for ( i = 0; i < NX; i++ ) data.d[i] = data.b[i];
        }
        #endif
    }

    /* correct lbA and ubA with d values */
    for ( i = 0; i < NNN; i++ ) {
        for ( j = 0; j < NX; j++ ) {
            data.lbA[i*NX+j] = data.lb[(i+1)*(NX+NU)+j] - data.d[i*NX+j];
            data.ubA[i*NX+j] = data.ub[(i+1)*(NX+NU)+j] - data.d[i*NX+j];
        }
    }

    /* !! Hessian propagation !! */
    #if FIXED_INITIAL_STATE == 0
    /* propagate x0: */
    /* A is unused for this operation because the W's are zero */
    computeWx(&data.Q[NNN*NX*NX], &data.C[(NNN-1)*NX], &data.A[0]);
    for ( i = NNN-1; i > 0; i-- ) {
        computeH_offDX(&data.Hc[NX+i*NU], &data.S[i*NX*NU], &data.C[(i-1)*NX], &data.B[i*NX*NU]);
        computeWx(&data.Q[i*NX*NX], &data.C[(i-1)*NX], &data.A[i*NX*NX]);
    }
    computeH_DX();
    #endif

    /* propagate controls: */
    for ( j = 0; j < NNN; j++ ) {
        for ( i = 0; i < NX*NU; i++ ) data.W2_u[i] = 0.0;
        /* A is unused here because W is zero */
        computeWu(&data.Q[NNN*NX*NX], &data.C[(offset+j*NU)*NNN*NX+(NNN-1)*NX], &data.A[0]);
        for ( i = NNN-1; i > j; i-- ) {
            computeH_offDU(&data.Hc[(offset+j*NU)*NVC+offset+i*NU], &data.S[i*NX*NU], \
                    &data.C[(offset+j*NU)*NNN*NX+(i-1)*NX], &data.B[i*NX*NU]);
            computeWu(&data.Q[i*NX*NX], &data.C[(offset+j*NU)*NNN*NX+(i-1)*NX], &data.A[i*NX*NX]);
        }
        computeH_DU(&data.Hc[(offset+j*NU)*NVC+offset+j*NU], &data.R[j*NU*NU], &data.B[j*NX*NU]);
    }

    /* !! gradient propagation !! */
    /* A is unused for this operation because the W's are zero */
    computeWg(&data.g[NNN*(NX+NU)], &data.Q[NNN*NX*NX], &data.d[(NNN-1)*NX], &data.A[0]);
    for ( i = NNN-1; i > 0; i-- ) {
        computeG_off(&data.gc[offset+i*NU], &data.g[i*(NX+NU)+NX], &data.S[i*NX*NU], \
                    &data.d[(i-1)*NX], &data.B[i*NX*NU]);
        computeWg(&data.g[i*(NX+NU)], &data.Q[i*NX*NX], &data.d[(i-1)*NX], &data.A[i*NX*NX]);
    }
    #if FIXED_INITIAL_STATE == 1
    real_t Sx0[NU] = {0};
    for ( j = 0; j < NU; j++ ) {
        for ( i = 0; i < NX; i++ ) Sx0[j] += data.S[i+j*NX]*x0[i];
    }
    computeG_fixed_initial_state(&Sx0[0]);
    #else
    computeG();
    #endif
}
