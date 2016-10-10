#include <stdlib.h>

#include "acados/utils/types.h"
#include "acados/ocp_qp/ocp_qp_common.h"
#include "acados/ocp_qp/ocp_qp_hpmpc.h"

#include "hpmpc/include/c_interface.h"



// work space size
int ocp_qp_hpmpc_workspace_size_bytes(int N, int *nx, int *nu, int *nb, int *ng, int **hidxb, \
    ocp_qp_hpmpc_args *hpmpc_args) {

    int ii;

    int N2 = hpmpc_args->N2;

    int k_max = hpmpc_args->max_iter;

    int workspace_size = 8 + 5*k_max*sizeof(double);  // alignment to single-word (4-byte)

    for (ii=0; ii <= N; ii++)
        workspace_size += nb[ii]*sizeof(int);

    workspace_size += hpmpc_d_ip_ocp_hard_tv_work_space_size_bytes(N, nx, nu, nb, hidxb, ng, N2);

    return workspace_size;
}



int ocp_qp_hpmpc(ocp_qp_in *qp_in, ocp_qp_out *qp_out, ocp_qp_hpmpc_args *hpmpc_args, \
    void *workspace) {

    // initialize return code
    int acados_status = ACADOS_SUCCESS;

    // loop index
    int ii, jj;

    // extract input struct members
    int N = qp_in->N;
    int *nx = (int *) qp_in->nx;
    int *nu = (int *) qp_in->nu;
    int *nb = (int *) qp_in->nb;
    int **hidxb_swp = (int **) qp_in->idxb;
    int *ng = (int *) qp_in->nc;
    double **hA = (double **) qp_in->A;
    double **hB = (double **) qp_in->B;
    double **hb = (double **) qp_in->b;
    double **hQ = (double **) qp_in->Q;
    double **hS = (double **) qp_in->S;
    double **hR = (double **) qp_in->R;
    double **hq = (double **) qp_in->q;
    double **hr = (double **) qp_in->r;
    double **hlb = (double **) qp_in->lb;
    double **hub = (double **) qp_in->ub;
    double **hC = (double **) qp_in->Cx;
    double **hD = (double **) qp_in->Cu;
    double **hlg = (double **) qp_in->lc;
    double **hug = (double **) qp_in->uc;

    // extract output struct members
    double **hx = qp_out->x;
    double **hu = qp_out->u;
    double **hpi = qp_out->pi;
    double **hlam = qp_out->lam;

    // extract args struct members
    double mu_tol = hpmpc_args->tol;
    int k_max = hpmpc_args->max_iter;
//  double alpha = hpmpc_args->min_step; // fixed in the solver
    double mu0 = hpmpc_args->mu0;
    int warm_start = hpmpc_args->warm_start;
    int N2 = hpmpc_args->N2;  // horizon length of the partially condensed problem

    //  other solver arguments
    int kk = -1;  // actual number of iterations
    double inf_norm_res[4];  // inf norm of residuals
    for (ii = 0; ii < 4; ii++) inf_norm_res[ii] = 0.0;  // zero

    // memory for stat
    size_t addr = ( ( (size_t) workspace) + 7) / 8 * 8;  // align to 8-byte boundaries
    double *ptr_double = (double *) addr;
    double *stat = ptr_double;
    ptr_double += 5*k_max;
    for (ii = 0; ii < 5*k_max; ii++) stat[ii] = 0.0;  // zero

    // memory for idxb
    int *hidxb[N+1];
//  size_t addr = (( (size_t) workspace ) + 3 ) / 4 * 4;  // align to 4-byte boundaries
    int *ptr_int = (int *) ptr_double;
    for (ii = 0; ii <= N; ii++) {
        hidxb[ii] = ptr_int;
        ptr_int += nb[ii];
    }
    workspace = (void *) ptr_int;

    //  swap x and u in bounds (by updating their indeces)
    for (ii = 0; ii <= N; ii++) {
        jj = 0;
        for (; jj < nb[ii]; jj++) {
            if (hidxb_swp[ii][jj] < nx[ii]) {  // state
                hidxb[ii][jj] = hidxb_swp[ii][jj]+nu[ii];
            } else {  // input
                hidxb[ii][jj] = hidxb_swp[ii][jj]-nx[ii];
            }
        }
    }

    int hpmpc_status = fortran_order_d_ip_ocp_hard_tv(&kk, k_max, mu0, mu_tol, N, nx, nu, nb, \
        hidxb, ng, N2, warm_start, hA, hB, hb, hQ, hS, hR, hq, hr, hlb, hub, hC, hD, hlg, hug, \
        hx, hu, hpi, hlam, inf_norm_res, workspace, stat);

    if (hpmpc_status == 1) acados_status = ACADOS_MAXITER;

    if (hpmpc_status == 2) acados_status = ACADOS_MINSTEP;

    // return
    return acados_status;
}
