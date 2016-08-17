#ifndef OCP_QP_CONDENSING_QPOASES_H_
#define OCP_QP_CONDENSING_QPOASES_H_

#include "acados_types.h"

// OCP QP interface
// struct of arguments to the solver
typedef struct ocp_qp_condensing_qpoases_args_ {
    real_t dummy;
} ocp_qp_condensing_qpoases_args;

int_t ocp_qp_condensing_qpoases(int_t N, int_t *nx, int_t *nu, int_t *nb, int_t *nc,
    double **A, double **B, double **b,
    double **Q, double **S, double **R, double **q, double **r,
    int_t **idxb, double **lb, double **ub,
    double **Cx, double **Cu, double **lc, double **uc,
    double **x, double **u,
    ocp_qp_condensing_qpoases_args *args, double *work);

int_t ocp_qp_condensing_qpoases_workspace_size(int_t N, int_t *nxx, int_t *nuu,
    int_t *nbb, int_t *ngg, ocp_qp_condensing_qpoases_args *args);

void initialise_qpoases();

#endif  // OCP_QP_CONDENSING_QPOASES_H_
