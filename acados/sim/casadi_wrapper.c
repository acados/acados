#include "acados/sim/casadi_wrapper.h"

void vde_fun(const int_t nx, const int_t nu, const real_t* in, real_t* out,
             int (*vde)(const real_t**, real_t**, int*, real_t*, int)) {

    const real_t* x = in;
    const real_t* Sx = in + nx;
    const real_t* Su = in + nx + nx * nx;
    const real_t* u = in + nx + nx * (nx + nu);

    real_t* x_out = out;
    real_t* Sx_out = out + nx;
    real_t* Su_out = out + nx + nx * nx;

    int_t casadi_mem = 0;
    int_t* casadi_iw = 0;
    real_t* casadi_w = 0;

    const real_t* casadi_arg[4];
    real_t* casadi_res[3];

    casadi_arg[0] = x;
    casadi_arg[1] = Sx;
    casadi_arg[2] = Su;
    casadi_arg[3] = u;

    casadi_res[0] = x_out;
    casadi_res[1] = Sx_out;
    casadi_res[2] = Su_out;

    vde(casadi_arg, casadi_res, casadi_iw, casadi_w, casadi_mem);
}
