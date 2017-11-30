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

void vde_hess_fun(const int_t nx, const int_t nu, const real_t* in, real_t* out,
    casadi_function_t vde_hess) {

    const double* x = in;
    const double* Sx = in + nx;
    const double* Su = in + nx + nx*nx;
    const double* lambdaX = in + nx*(1+nx+nu);
    const double* u  = in + nx*(2+nx+nu);

    double* adj_out = out;
    double* hess_out = out+nx+nu;

    int casadi_mem = 0;
    int *casadi_iw = 0;
    double *casadi_w = 0;

    const double *casadi_arg[5];
    double *casadi_res[2];

    casadi_arg[0] = x;
    casadi_arg[1] = Sx;
    casadi_arg[2] = Su;
    casadi_arg[3] = lambdaX;
    casadi_arg[4] = u;

    casadi_res[0] = adj_out;
    casadi_res[1] = hess_out;

    vde_hess(casadi_arg, casadi_res, casadi_iw, casadi_w, casadi_mem);
}

void vde_adj_fun(const int_t nx, const int_t nu, const real_t* in, real_t* out,
                 casadi_function_t vde_adj) {

    const real_t* x = in;
    const real_t* lambdaX = in + nx;
    const real_t* u  = in + 2*nx;

    real_t* adj_out = out;

    int casadi_mem = 0;
    int *casadi_iw = 0;
    real_t *casadi_w = 0;

    const real_t *casadi_arg[3];
    real_t *casadi_res[1];

    casadi_arg[0] = x;
    casadi_arg[1] = lambdaX;
    casadi_arg[2] = u;

    casadi_res[0] = adj_out;

    vde_adj(casadi_arg, casadi_res, casadi_iw, casadi_w, casadi_mem);
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

void impl_ode_fun(const int_t nx, const int_t nu, const real_t *in, real_t *out, casadi_function_t impl_ode){
    const double *x = in;
    const double *xdot = in + nx;
    const double *u = in + 2 * nx;

    double *res = out;

    int casadi_mem = 0;
    int *casadi_iw = NULL;
    double *casadi_w = NULL;

    const double *casadi_arg[3];
    double *casadi_res[1];


    casadi_arg[0] = x;
    casadi_arg[1] = xdot;
    casadi_arg[2] = u;

    casadi_res[0] = res;  

    impl_ode(casadi_arg, casadi_res, casadi_iw, casadi_w, casadi_mem);
}

void impl_jac_x_fun(const int_t nx, const int_t nu, const real_t *in, real_t *out, casadi_function_t impl_jac_x){
    const double *x = in;
    const double *xdot = in + nx;
    const double *u = in + 2 * nx;

    double *jac_x = out;

    int casadi_mem = 0;
    int *casadi_iw = NULL;
    double *casadi_w = NULL;

    const double *casadi_arg[3];
    double *casadi_res[1];

    casadi_arg[0] = x;
    casadi_arg[1] = xdot;
    casadi_arg[2] = u;

    casadi_res[0] = jac_x;

    impl_jac_x(casadi_arg, casadi_res, casadi_iw, casadi_w, casadi_mem);
}

void impl_jac_xdot_fun(const int_t nx, const int_t nu, const real_t *in, real_t *out, casadi_function_t impl_jac_xdot){
    const double *x = in;
    const double *xdot = in + nx;
    const double *u = in + 2 * nx;

    double *jac_xdot = out;

    int casadi_mem = 0;
    int *casadi_iw = NULL;
    double *casadi_w = NULL;

    const double *casadi_arg[3];
    double *casadi_res[1];

    casadi_arg[0] = x;
    casadi_arg[1] = xdot;
    casadi_arg[2] = u;

    casadi_res[0] = jac_xdot;

    impl_jac_xdot(casadi_arg, casadi_res, casadi_iw, casadi_w, casadi_mem);
}

void impl_jac_u_fun(const int_t nx, const int_t nu, const real_t *in, real_t *out, casadi_function_t impl_jac_u){
    const double *x = in;
    const double *xdot = in + nx;
    const double *u = in + 2 * nx;

    double *jac_u = out;

    int casadi_mem = 0;
    int *casadi_iw = NULL;
    double *casadi_w = NULL;

    const double *casadi_arg[3];
    double *casadi_res[1];

    casadi_arg[0] = x;
    casadi_arg[1] = xdot;
    casadi_arg[2] = u;

    casadi_res[0] = jac_u;

    impl_jac_u(casadi_arg, casadi_res, casadi_iw, casadi_w, casadi_mem);
}