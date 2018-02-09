#include <stddef.h>

#include "acados/sim/sim_gnsf_casadi_wrapper.h"

void res_inc_Jff_wrapped(const int nx1, const int nu, const int n_out, const int num_stages, const real_t *in, real_t *out, casadi_function_t res_inc_Jff_fun) {

    const double *ff = in;
    const double *x1 = in + n_out * num_stages;
    const double *u  = in + n_out * num_stages + nx1;

    // double *res_inc_Jff_out = out;

    int casadi_mem = 0;
    int *casadi_iw = NULL;
    double *casadi_w = NULL;

    const double *casadi_arg[3];
    double *casadi_res[1];

    casadi_arg[0] = ff;
    casadi_arg[1] = x1;
    casadi_arg[2] = u;

    // casadi_res[0] = res_inc_Jff_out;
    casadi_res[0] = out;


    res_inc_Jff_fun(casadi_arg, casadi_res, casadi_iw, casadi_w, casadi_mem);
}



void jac_res_ffx1u_wrapped(const int nx1, const int nu, const int n_out, const int num_stages, const real_t *in, real_t *out, casadi_function_t jac_res_ffx1u_fun) {

    const double *ff = in;
    const double *x1 = in + n_out * num_stages;
    const double *u  = in + n_out * num_stages + nx1;

    double *res_inc_Jff_out = out;

    int casadi_mem = 0;
    int *casadi_iw = NULL;
    double *casadi_w = NULL;

    const double *casadi_arg[3];
    double *casadi_res[1];

    casadi_arg[0] = ff;
    casadi_arg[1] = x1;
    casadi_arg[2] = u;

    casadi_res[0] = res_inc_Jff_out;

    jac_res_ffx1u_fun(casadi_arg, casadi_res, casadi_iw, casadi_w, casadi_mem);
}


void export_from_ML_wrapped(const real_t *in, real_t *out, casadi_function_t get_ints_fun){
    
    int casadi_mem = 0;
    int *casadi_iw = NULL;
    double *casadi_w = NULL;

    double *ints_out = out;

    const double *casadi_arg[1];
    double *casadi_res[1];

    casadi_arg[0] = in;

    casadi_res[0] = ints_out;

    get_ints_fun(casadi_arg, casadi_res, casadi_iw, casadi_w, casadi_mem);
}



// template:
// void vde_fun(const int_t nx, const int_t nu, const real_t *in, real_t *out, casadi_function_t vde) {

//     const double *x = in;
//     const double *Sx = in + nx;
//     const double *Su = in + nx + nx * nx;
//     const double *u = in + nx + nx * (nx + nu);

//     double *x_out = out;
//     double *Sx_out = out + nx;
//     double *Su_out = out + nx + nx * nx;

//     int casadi_mem = 0;
//     int *casadi_iw = NULL;
//     double *casadi_w = NULL;

//     const double *casadi_arg[4];
//     double *casadi_res[3];

//     casadi_arg[0] = x;
//     casadi_arg[1] = Sx;
//     casadi_arg[2] = Su;
//     casadi_arg[3] = u;

//     casadi_res[0] = x_out;
//     casadi_res[1] = Sx_out;
//     casadi_res[2] = Su_out;

//     vde(casadi_arg, casadi_res, casadi_iw, casadi_w, casadi_mem);
// }