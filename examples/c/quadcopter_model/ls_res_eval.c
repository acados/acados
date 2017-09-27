#pragma GCC diagnostic ignored "-Wunused-parameter"
#pragma GCC diagnostic ignored "-Wunused-function"
#include "ls_res_eval.h"

// wrapper to casadi autogen code for residual evaluation

void ls_res_eval_quadcopter(const real_t* in, real_t* out,
                      int (*ls_res)(const real_t**, real_t**, int*, real_t*,
                                 int)) {

    const int nx = 11;
    const int nu = 4;
    const int nr = 15;
    const int np = 0;

    const real_t* x = in;
    const real_t* u = in + nx;
    const real_t* p = in + nx + nu;

    real_t* r  = out;
    real_t* drdw = out + nr;
    real_t* rref = out + nr + nr*(nx+nu);

    // TODO(Andrea): need to remove this and allocate memory somewhere else
    int casadi_mem = 0;
    int *casadi_iw = 0;
    real_t *casadi_w = 0;

    // this stuff is hard coded because dimensions won't change (I assume?)
    const real_t *casadi_arg[4];
    real_t *casadi_res[3];

    casadi_arg[0] = x;
    casadi_arg[1] = u;
    casadi_arg[2] = p;

    casadi_res[0] = r;
    casadi_res[1] = drdw;

    ls_res(casadi_arg, casadi_res, casadi_iw, casadi_w, casadi_mem);

}

void ls_res_eval_quadcopter_end(const real_t* in, real_t* out,
                      int (*ls_res)(const real_t**, real_t**, int*, real_t*,
                                 int)) {

    const int nx = 11;
    const int nu = 0;
    const int nr = 11;
    const int np = 0;

    const real_t* x = in;
    const real_t* u = in + nx;
    const real_t* p = in + nx + nu;

    real_t* r  = out;
    real_t* drdw = out + nr;
    real_t* rref = out + nr + nr*(nx+nu);

    // TODO(Andrea): need to remove this and allocate memory somewhere else
    int casadi_mem = 0;
    int *casadi_iw = 0;
    real_t *casadi_w = 0;

    // this stuff is hard coded because dimensions won't change (I assume?)
    const real_t *casadi_arg[4];
    real_t *casadi_res[3];

    casadi_arg[0] = x;
    casadi_arg[1] = u;
    casadi_arg[2] = p;

    casadi_res[0] = r;
    casadi_res[1] = drdw;

    ls_res_end(casadi_arg, casadi_res, casadi_iw, casadi_w, casadi_mem);

}
