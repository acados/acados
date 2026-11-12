#ifdef __cplusplus
extern "C" {
#endif

#include "acados/utils/external_function_generic.h"
#include "blasfeo_d_blas.h"
#include "blasfeo_d_aux.h"

enum cost_stage
{
    COST_INITIAL,
    COST_PATH,
    COST_TERMINAL
};

static void cost_value(void **in, double *fun, void *params, enum cost_stage stage)
{
    const int nx = 8;
    const int nu = stage == COST_TERMINAL ? 0 : 3;
    const int include_x = stage != COST_INITIAL;
    const double Qdiag[] = {1, 1, 1, 1, 1, 1, 1, 1};
    const double Rdiag[] = {2, 2, 2};
    struct blasfeo_dvec_args *x_args = in[0];
    struct blasfeo_dvec *x = x_args->x;
    int xi = x_args->xi;

    *fun = 0.0;
    if (include_x)
    {
        for (int ii = 0; ii < nx; ii++)
            *fun += BLASFEO_DVECEL(x, xi + ii) * BLASFEO_DVECEL(x, xi + ii) * Qdiag[ii];
    }
    if (nu > 0)
    {
        struct blasfeo_dvec_args *u_args = in[1];
        struct blasfeo_dvec *u = u_args->x;
        int ui = u_args->xi;
        for (int ii = 0; ii < nu; ii++)
            *fun += Rdiag[ii] * BLASFEO_DVECEL(u, ui + ii) * BLASFEO_DVECEL(u, ui + ii);
    }
    *fun *= 0.5;
}

static void cost_jac(void **in, void **out, void *params, enum cost_stage stage)
{
    const int nx = 8;
    const int nu = stage == COST_TERMINAL ? 0 : 3;
    const int include_x = stage != COST_INITIAL;
    const double Qdiag[] = {1, 1, 1, 1, 1, 1, 1, 1};
    const double Rdiag[] = {2, 2, 2};
    struct blasfeo_dvec_args *x_args = in[0];
    struct blasfeo_dvec *x = x_args->x;
    int xi = x_args->xi;
    struct blasfeo_dvec *grad = out[1];

    blasfeo_dvecse(nu + nx, 0.0, grad, 0);
    if (nu > 0)
    {
        struct blasfeo_dvec_args *u_args = in[1];
        struct blasfeo_dvec *u = u_args->x;
        int ui = u_args->xi;
        for (int ii = 0; ii < nu; ii++)
            BLASFEO_DVECEL(grad, ii) = Rdiag[ii] * BLASFEO_DVECEL(u, ui + ii);
    }
    if (include_x)
    {
        for (int ii = 0; ii < nx; ii++)
            BLASFEO_DVECEL(grad, nu + ii) = Qdiag[ii] * BLASFEO_DVECEL(x, xi + ii);
    }
    cost_value(in, out[0], params, stage);
}

static void cost_jac_hess(void **in, void **out, void *params, enum cost_stage stage)
{
    const int nx = 8;
    const int nu = stage == COST_TERMINAL ? 0 : 3;
    const int include_x = stage != COST_INITIAL;
    const double Qdiag[] = {1, 1, 1, 1, 1, 1, 1, 1};
    const double Rdiag[] = {2, 2, 2};
    struct blasfeo_dmat *hess = out[2];

    blasfeo_dgese(nu + nx, nu + nx, 0.0, hess, 0, 0);
    for (int ii = 0; ii < nu; ii++)
        BLASFEO_DMATEL(hess, ii, ii) = Rdiag[ii];
    if (include_x)
    {
        for (int ii = 0; ii < nx; ii++)
            BLASFEO_DMATEL(hess, nu + ii, nu + ii) = Qdiag[ii];
    }
    cost_jac(in, out, params, stage);
}

int ext_cost_0(void **in, void **out, void *params) { cost_jac_hess(in, out, params, COST_INITIAL); return 0; }
int ext_cost_0_jac(void **in, void **out, void *params) { cost_jac(in, out, params, COST_INITIAL); return 0; }
int ext_cost_0_jac_hess(void **in, void **out, void *params) { cost_jac_hess(in, out, params, COST_INITIAL); return 0; }
int ext_cost(void **in, void **out, void *params) { cost_jac_hess(in, out, params, COST_PATH); return 0; }
int ext_cost_jac(void **in, void **out, void *params) { cost_jac(in, out, params, COST_PATH); return 0; }
int ext_cost_jac_hess(void **in, void **out, void *params) { cost_jac_hess(in, out, params, COST_PATH); return 0; }
int ext_costN(void **in, void **out, void *params) { cost_jac_hess(in, out, params, COST_TERMINAL); return 0; }
int ext_costN_jac(void **in, void **out, void *params) { cost_jac(in, out, params, COST_TERMINAL); return 0; }
int ext_costN_jac_hess(void **in, void **out, void *params) { cost_jac_hess(in, out, params, COST_TERMINAL); return 0; }

#ifdef __cplusplus
} /* extern "C" */
#endif
