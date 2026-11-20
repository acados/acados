/*
 * Copyright (c) The acados authors.
 *
 * This file is part of acados.
 *
 * Licensed under the 2-Clause BSD License.
 */

#include "acados/utils/external_function_generic.h"
#include "blasfeo_d_aux.h"

static void set_constraint_jacobian(struct blasfeo_dmat *jac, int ai, int aj)
{
    blasfeo_dgese(3, 3, 0.0, jac, ai, aj);
    BLASFEO_DMATEL(jac, ai, aj) = 1.0;
    BLASFEO_DMATEL(jac, ai + 1, aj + 1) = 1.0;
    BLASFEO_DMATEL(jac, ai + 2, aj + 2) = 1.0;
}

int generic_constraints_fun(void **in, void **out, void *params)
{
    struct blasfeo_dvec_args *x_args = in[0];
    struct blasfeo_dvec_args *u_args = in[1];
    struct blasfeo_dvec_args *fun_args = out[0];
    struct blasfeo_dvec *value = fun_args->x;

    blasfeo_dveccp(1, u_args->x, u_args->xi, value, fun_args->xi);
    blasfeo_dveccp(2, x_args->x, x_args->xi, value, fun_args->xi + 1);
    return 0;
}

int generic_constraints_fun_jac(void **in, void **out, void *params)
{
    generic_constraints_fun(in, out, params);

    struct blasfeo_dmat_args *jac_args = out[1];
    set_constraint_jacobian(jac_args->A, jac_args->ai, jac_args->aj);
    return 0;
}

int generic_constraints_fun_jac_hess(void **in, void **out, void *params)
{
    generic_constraints_fun_jac(in, out, params);

    struct blasfeo_dmat_args *hess_args = out[2];
    blasfeo_dgese(3, 3, 0.0, hess_args->A, hess_args->ai, hess_args->aj);
    return 0;
}

int generic_constraints_fun_e(void **in, void **out, void *params)
{
    struct blasfeo_dvec_args *x_args = in[0];
    struct blasfeo_dvec *x = x_args->x;
    (void) in[1];
    struct blasfeo_dvec_args *fun_args = out[0];
    struct blasfeo_dvec *value = fun_args->x;

    blasfeo_dveccp(2, x, x_args->xi, value, fun_args->xi);
    return 0;
}

int generic_constraints_fun_jac_e(void **in, void **out, void *params)
{
    generic_constraints_fun_e(in, out, params);

    struct blasfeo_dmat_args *jac_args = out[1];
    blasfeo_dgese(2, 2, 0.0, jac_args->A, jac_args->ai, jac_args->aj);
    BLASFEO_DMATEL(jac_args->A, jac_args->ai, jac_args->aj) = 1.0;
    BLASFEO_DMATEL(jac_args->A, jac_args->ai + 1, jac_args->aj + 1) = 1.0;
    return 0;
}

int generic_constraints_fun_jac_hess_e(void **in, void **out, void *params)
{
    generic_constraints_fun_jac_e(in, out, params);

    struct blasfeo_dmat_args *hess_args = out[2];
    blasfeo_dgese(2, 2, 0.0, hess_args->A, hess_args->ai, hess_args->aj);
    return 0;
}
