/* The model comes from \cite{Wirsching2006} */
#include "examples/c/quadcopter_model/quadcopter_model.h"

// extern int jacFun(const real_t** arg, real_t** res);

// The auto-generated VDE functions from CasADi:

void VDE_fun_quadcopter(const real_t* in, real_t* out,
                      int (*vde)(const real_t**, real_t**, int*, real_t*,
                                 int)) {
    int_t NX = 11;
    int_t NU = 4;
    const real_t* x = in;
    const real_t* Sx = in + NX;
    const real_t* Su = in + NX + NX * NX;
    const real_t* u = in + NX + NX * (NX + NU);

    real_t* x_out = out;
    real_t* Sx_out = out + NX;
    real_t* Su_out = out + NX + NX * NX;

    const real_t* casadi_arg[4];
    real_t* casadi_res[3];

    casadi_arg[0] = x;
    casadi_arg[1] = Sx;
    casadi_arg[2] = Su;
    casadi_arg[3] = u;

    casadi_res[0] = x_out;
    casadi_res[1] = Sx_out;
    casadi_res[2] = Su_out;

    int_t* iw = 0;
    real_t* w = 0;
    int_t mem = 0;

    vde(casadi_arg, casadi_res, iw, w, mem);
}
