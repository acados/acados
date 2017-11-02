#include <stddef.h>

#include "acados/sim/sim_casadi_wrapper.h"

void vde_fun(const int_t nx, const int_t nu, const real_t *in, real_t *out, casadi_function_t vde) {

    const double *x = in;
    const double *Sx = in + nx;
    const double *Su = in + nx + nx * nx;
    const double *u = in + nx + nx * (nx + nu);

    double *x_out = out;
    double *Sx_out = out + nx;
    double *Su_out = out + nx + nx * nx;

    int casadi_mem = 0;
    int *casadi_iw = NULL;
    double *casadi_w = NULL;

    const double *casadi_arg[4];
    double *casadi_res[3];

    casadi_arg[0] = x;
    casadi_arg[1] = Sx;
    casadi_arg[2] = Su;
    casadi_arg[3] = u;

    casadi_res[0] = x_out;
    casadi_res[1] = Sx_out;
    casadi_res[2] = Su_out;

    vde(casadi_arg, casadi_res, casadi_iw, casadi_w, casadi_mem);
}

void jac_fun(const int_t nx, const real_t *in, real_t *out, casadi_function_t jac) {

    const double *x = in;
    const double *u = in + nx;

    double *x_out = out;
    double *jac_out = out + nx;

    int casadi_mem = 0;
    int *casadi_iw = NULL;
    double *casadi_w = NULL;

    const double *casadi_arg[2];
    double *casadi_res[2];

    casadi_arg[0] = x;
    casadi_arg[1] = u;

    casadi_res[0] = x_out;
    casadi_res[1] = jac_out;

    jac(casadi_arg, casadi_res, casadi_iw, casadi_w, casadi_mem);
}

void discrete_model_fun(const int_t nx, const int_t nu, const real_t *in, real_t *out,
                        casadi_function_t discrete_model) {

    const double *x = in;
    const double *u = in + nx;

    double *x_out = out;
    double *jac_out = out + nx;

    int casadi_mem = 0;
    int *casadi_iw = NULL;
    double *casadi_w = NULL;

    const double *casadi_arg[2];
    double *casadi_res[2];

    casadi_arg[0] = x;
    casadi_arg[1] = u;

    casadi_res[0] = x_out;
    casadi_res[1] = jac_out;

    discrete_model(casadi_arg, casadi_res, casadi_iw, casadi_w, casadi_mem);
}
