#ifndef ACADOS_OCP_QP_COMMON_H_
#define ACADOS_OCP_QP_COMMON_H_

typedef struct ocp_qp_input_ {
    int_t N;
    int_t *nx;
    int_t *nu;
    int_t *nb;
    int_t *nc;
    real_t **A;
    real_t **B;
    real_t **b;
    real_t **Q;
    real_t **S;
    real_t **R;
    real_t **q;
    real_t **r;
    int_t **idxb;
    real_t **lb;
    real_t **ub;
    real_t **Cx;
    real_t **Cu;
    real_t **lc;
    real_t **uc;
} ocp_qp_input;

typedef struct ocp_qp_output_ {
    real_t **x;
    real_t **u;
} ocp_qp_output;

#endif  // ACADOS_OCP_QP_COMMON_H_
