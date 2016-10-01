#ifndef ACADOS_OCP_QP_HPMPC_H_
#define ACADOS_OCP_QP_HPMPC_H_

#include "acados/utils/types.h"

// OCP QP interface
// struct of arguments to the solver
struct ocp_qp_hpmpc_args{
    double tol;
    int max_iter;
    double min_step;
    double mu0;
    double sigma_min;
};

// struct ocp_qpStruct {}

int ocp_qp_hpmpc(int N, int *nx, int *nu, int *nb, int *ng, \
    double **A, double **B, double **b, \
    double **Q, double **S, double **R, double **q, double **r, \
    int **idxb, double **lb, double **ub, \
    double **C, double **D, double **lg, double **ug, \
    double **x, double **u, \
    struct ocp_qp_hpmpc_args *args, double *work);

// int ocp_qp_hpmpc(qpStruct *qp, args *args,double *work)

int ocp_qp_hpmpc_workspace_size(int N, int *nxx, int *nuu, int *nbb, int *ngg,
                                struct ocp_qp_hpmpc_args *args);


#endif  // ACADOS_OCP_QP_HPMPC_H_
