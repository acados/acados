#ifndef ACADOS_OCP_QP_OCP_QP_HPMPC_H_
#define ACADOS_OCP_QP_OCP_QP_HPMPC_H_

#include "acados/utils/types.h"
#include "acados/ocp_qp/ocp_qp_common.h"

// OCP QP interface
// struct of arguments to the solver
typedef struct ocp_qp_hpmpc_args_
	{
	double tol;
	int max_iter;
//	double min_step;
	double mu0;
//	double sigma_min;
	int warm_start;
	int N2; // horizion length of the partially condensed problem
	} ocp_qp_hpmpc_args;

int ocp_qp_hpmpc(ocp_qp_in *qp_in, ocp_qp_out *qp_out, ocp_qp_hpmpc_args *qp_args, void *workspace);

int ocp_qp_hpmpc_workspace_size_bytes(int N, int *nx, int *nu, int *nb, int *ng, int **hidxb, ocp_qp_hpmpc_args *hpmpc_args);


#endif  // ACADOS_OCP_QP_OCP_QP_HPMPC_H_
