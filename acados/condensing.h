#ifndef ACADOS_CONDENSING_H_
#define ACADOS_CONDENSING_H_

#include "acados_types.h"
#include "ocp_qp_common.h"

#define FIXED_INITIAL_STATE 1
#if FIXED_INITIAL_STATE == 1
#define NVC 60
#else
#define NVC 68
#endif

#define NX 8
#define NU 3
#define NNN 20
#define NA 11
#define NCONSTRAINTS 169

typedef struct condensing_in_ {
    ocp_qp_input *qp_input;
} condensing_in;

typedef struct condensing_out_ {
    real_t *H;
    real_t *h;
    real_t *lb;
    real_t *ub;
    real_t *A;
    real_t *lbA;
    real_t *ubA;
} condensing_out;

typedef struct condensing_memory_ {
    real_t dummy;
} condensing_memory;

typedef struct condensing_workspace_ {
    real_t ***G;
    real_t **g;
    real_t ***D;
    real_t *W1_x;
    real_t *W2_x;
    real_t *W1_u;
    real_t *W2_u;
    real_t *w1;
    real_t *w2;
} condensing_workspace;

void condensingN2_fixed_initial_state(condensing_in *input, condensing_out *output,
    condensing_workspace *workspace);

void condensingN2_free_initial_state();

#endif  // ACADOS_CONDENSING_H_
