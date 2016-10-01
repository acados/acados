#ifndef ACADOS_OCP_QP_CONDENSING_H_
#define ACADOS_OCP_QP_CONDENSING_H_

#include "acados/utils/types.h"
#include "acados/ocp_qp/ocp_qp_common.h"

typedef struct condensing_in_ {
    ocp_qp_input *qp_input;
} condensing_input;

typedef struct condensing_out_ {
    real_t *H;
    real_t *h;
    real_t *lb;
    real_t *ub;
    real_t *A;
    real_t *lbA;
    real_t *ubA;
} condensing_output;

typedef struct condensing_memory_ {
    real_t dummy;
} condensing_memory;

typedef struct condensing_workspace_ {
    int_t nconvars;
    int_t nconstraints;
    int_t *nstate_bounds;
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

void condensingN2_fixed_initial_state(condensing_input *in, condensing_output *out,
    condensing_workspace *work);

void condensingN2_free_initial_state();

#endif  // ACADOS_OCP_QP_CONDENSING_H_
