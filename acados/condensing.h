#ifndef CONDENSING_H
#define CONDENSING_H

#include "acados_types.h"

#define FIXED_INITIAL_STATE 1
#if FIXED_INITIAL_STATE == 1
#define NVC 60
#else
#define NVC 68
#endif

#define NX 8
#define NU 3
#define NNN 20

typedef struct data_struct_ {
    /* INPUT DATA */
    real_t Q[(NNN+1)*NX*NX];
    real_t R[NNN*NU*NU];
    real_t S[NNN*NX*NU];
    real_t A[NNN*NX*NX];
    real_t B[NNN*NX*NU];
    real_t g[NNN*(NX+NU)+NX];
    real_t b[NNN*NX];
    real_t lb[NNN*(NX+NU)+NX];
    real_t ub[NNN*(NX+NU)+NX];
    /* OUTPUT DATA */
    real_t Hc[NVC*NVC];
    real_t gc[NVC];
    /*real_t Ac[NNN*NX*NVC];*/
    real_t lbA[NNN*NX];
    real_t ubA[NNN*NX];
    real_t lbU[NVC];
    real_t ubU[NVC];
    real_t C[NNN*NX*NVC];
    real_t d[NNN*NX];
    real_t W1_x[NX*NX];
    real_t W2_x[NX*NX];
    real_t W1_u[NX*NU];
    real_t W2_u[NX*NU];
    real_t w1[NX];
    real_t w2[NX];
} data_struct;

void condensingN2_fixed_initial_state();

void condensingN2_free_initial_state();

#endif  // CONDENSING_H
