#ifndef OCP_QP_CONDENSING_QPOASES_H_
#define OCP_QP_CONDENSING_QPOASES_H_

#include "acados_types.h"

// OCP QP interface
// struct of arguments to the solver
struct ocp_qp_condensing_qpoases_args{
    real_t tol;
    int_t max_iter;
    real_t min_step;
    real_t mu0;
    real_t sigma_min;
};

int_t ocp_qp_condensing_qpoases(int_t N, int_t *nx, int_t *nu, int_t *nb, int_t *ng, \
    real_t **A, real_t **B, real_t **b, \
    real_t **Q, real_t **S, real_t **R, real_t **q, real_t **r, \
    int_t **idxb, real_t **lb, real_t **ub, \
    real_t **C, real_t **D, real_t **lg, real_t **ug, \
    real_t **x, real_t **u, \
    struct ocp_qp_condensing_qpoases_args *args, real_t *work);

int_t ocp_qp_condensing_qpoases_workspace_size(int_t N, int_t *nxx, int_t *nuu,
    int_t *nbb, int_t *ngg, struct ocp_qp_condensing_qpoases_args *args);

void initialise_qpoases();

#endif  // OCP_QP_CONDENSING_QPOASES_H_
