
#include <stddef.h>

#include "acados/ocp_qp/ocp_qp_common.h"
#include "acados/utils/print.h"
#include "acados_c/ocp_qp.h"

// QP data printed e.g. from matlab
#include "./ocp_qp_bugs/ocp_qp_data.c"

int main() {

    /************************************************
    * choose solver
    ************************************************/

    ocp_qp_solver_plan plan;
    plan.qp_solver = PARTIAL_CONDENSING_QPDUNES;

    /************************************************
    * set up dims struct
    ************************************************/

    ocp_qp_dims dims;

    dims.N = N;
    dims.nx = nx;
    dims.nu = nu;
    dims.nb = nb;
    dims.ng = ng;
    dims.ns = ns;
    dims.nbx = nbx;
    dims.nbu = nbu;

    /************************************************
    * set up dynamics
    ************************************************/

    double **hA = malloc(N*sizeof(double *));
    double **hB = malloc(N*sizeof(double *));
    double **hb = malloc(N*sizeof(double *));

    int sum_A = 0;
    int sum_B = 0;
    int sum_b = 0;

    for (int ii = 0; ii < N; ii++)
    {
        // TODO(dimitris): test for varying dimensions
        hA[ii] = &Av[sum_A];
        sum_A += dims.nx[ii+1]*dims.nx[ii];

        hB[ii] = &Bv[sum_B];
        sum_B += dims.nx[ii+1]*dims.nu[ii];

        hb[ii] = &b[sum_b];
        sum_b += dims.nx[ii+1];
    }

    /************************************************
    * set up objective
    ************************************************/

    int sum_Q = 0;
    int sum_R = 0;
    int sum_S = 0;
    int sum_q = 0;
    int sum_r = 0;

    double **hQ = malloc((N+1)*sizeof(double *));
    double **hR = malloc((N+1)*sizeof(double *));
    double **hS = malloc((N+1)*sizeof(double *));
    double **hq = malloc((N+1)*sizeof(double *));
    double **hr = malloc((N+1)*sizeof(double *));

    for (int ii = 0; ii < N+1; ii++)
    {
        hQ[ii] = &Qv[sum_Q];
        sum_Q += dims.nx[ii]*dims.nx[ii];

        hR[ii] = &Rv[sum_R];
        sum_R += dims.nu[ii]*dims.nu[ii];

        hS[ii] = &Sv[sum_S];
        sum_S += dims.nx[ii]*dims.nu[ii];

        hq[ii] = &q[sum_q];
        sum_q += dims.nx[ii];

        hr[ii] = &r[sum_r];
        sum_r += dims.nu[ii];
    }

    /************************************************
    * set up constraints
    ************************************************/

    int **hidxb = malloc((N+1)*sizeof(int *));
    double **hlb = malloc((N+1)*sizeof(double *));
    double **hub = malloc((N+1)*sizeof(double *));
    double **hC = malloc((N+1)*sizeof(double *));
    double **hD = malloc((N+1)*sizeof(double *));
    double **hlg = malloc((N+1)*sizeof(double *));
    double **hug = malloc((N+1)*sizeof(double *));

    int sum_nb = 0;

    for (int ii = 0; ii < N+1; ii++)
    {
        hidxb[ii] = &idxb[sum_nb];
        hlb[ii] = &lb[sum_nb];
        hub[ii] = &ub[sum_nb];
        sum_nb += dims.nb[ii];
    }

    ocp_qp_in *qp_in = create_ocp_qp_in(&dims);

    d_cvt_colmaj_to_ocp_qp(hA, hB, hb, hQ, hS, hR, hq, hr, hidxb, hlb, hub, hC, hD, hlg, hug, NULL, NULL, NULL, NULL, NULL, qp_in);

    // print_ocp_qp_in(qp_in);

    void *args = ocp_qp_create_args(&plan, &dims);

    ocp_qp_out *qp_out = create_ocp_qp_out(&dims);

    ocp_qp_solver *qp_solver = ocp_qp_create(&plan, &dims, args);

    int acados_return = ocp_qp_solve(qp_solver, qp_in, qp_out);

    print_ocp_qp_out(qp_out);

    if (acados_return != ACADOS_SUCCESS)
    {
        printf("QP SOLVER FAILED WITH FLAG %d\n", acados_return);
        return -1;
    }

    free(hA);
    free(hB);
    free(hb);

    free(hQ);
    free(hR);
    free(hS);
    free(hq);
    free(hr);


    free(hidxb);
    free(hlb);
    free(hub);
    free(hC);
    free(hD);
    free(hlg);
    free(hug);

    return 0;
}
