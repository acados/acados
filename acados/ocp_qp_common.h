#ifndef ACADOS_OCP_QP_COMMON_H_
#define ACADOS_OCP_QP_COMMON_H_

typedef struct ocp_qp_input_ {
    int_t N;
    const int_t *nx;
    const int_t *nu;
    const int_t *nb;
    const int_t *nc;
    const real_t **A;
    const real_t **B;
    const real_t **b;
    const real_t **Q;
    const real_t **S;
    const real_t **R;
    const real_t **q;
    const real_t **r;
    const int_t **idxb;
    const real_t **lb;
    const real_t **ub;
    const real_t **Cx;
    const real_t **Cu;
    const real_t **lc;
    const real_t **uc;
} ocp_qp_input;

typedef struct ocp_qp_output_ {
    real_t **x;
    real_t **u;
} ocp_qp_output;

#endif  // ACADOS_OCP_QP_COMMON_H_
