#include "test/sim/pendulum/casadi/casadi_pendulum.h"

extern int vde_forw_pendulum(const real_t** arg, real_t** res, int* iw, real_t* w, int mem);
extern int vde_adj_pendulum(const real_t** arg, real_t** res, int* iw, real_t* w, int mem);
extern int vde_hess_pendulum(const real_t** arg, real_t** res, int* iw, real_t* w, int mem);
extern int jac_pendulum(const real_t** arg, real_t** res, int* iw, real_t* w, int mem);

// The auto-generated VDE functions from CasADi:
void VDE_forw_pendulum(const real_t* in, real_t* out, int (*vde)(const real_t**, real_t**, int*, real_t*, int)) {
    int_t NX = 4;
    int_t NU = 1;
    const real_t* x = in;
    const real_t* Sx = in + NX;
    const real_t* Su = in + NX + NX*NX;
    const real_t* u  = in + NX + NX*(NX+NU);

    real_t* x_out = out;
    real_t* Sx_out = out + NX;
    real_t* Su_out = out + NX + NX*NX;

    int casadi_mem = 0;
    int *casadi_iw = 0;
    real_t *casadi_w = 0;

    const real_t *casadi_arg[4];
    real_t *casadi_res[3];

    casadi_arg[0] = x;
    casadi_arg[1] = Sx;
    casadi_arg[2] = Su;
    casadi_arg[3] = u;

    casadi_res[0] = x_out;
    casadi_res[1] = Sx_out;
    casadi_res[2] = Su_out;
    vde(casadi_arg, casadi_res, casadi_iw, casadi_w, casadi_mem);
    vde_forw_pendulum(casadi_arg, casadi_res, casadi_iw, casadi_w, casadi_mem);
}


// The auto-generated VDE functions from CasADi:
void VDE_hess_pendulum(const real_t* in, real_t* out) {
    int_t NX = 4;
    int_t NU = 1;
    const real_t* x = in;
    const real_t* Sx = in + NX;
    const real_t* Su = in + NX + NX*NX;
    const real_t* lambdaX = in + NX*(1+NX+NU);
    const real_t* u  = in + NX*(2+NX+NU);

    real_t* adj_out = out;
    real_t* hess_out = out+NX+NU;

    int casadi_mem = 0;
    int *casadi_iw = 0;
    real_t *casadi_w = 0;

    const real_t *casadi_arg[5];
    real_t *casadi_res[2];

    casadi_arg[0] = x;
    casadi_arg[1] = Sx;
    casadi_arg[2] = Su;
    casadi_arg[3] = lambdaX;
    casadi_arg[4] = u;

    casadi_res[0] = adj_out;
    casadi_res[1] = hess_out;

    vde_hess_pendulum(casadi_arg, casadi_res, casadi_iw, casadi_w, casadi_mem);
}


// The auto-generated VDE functions from CasADi:
void VDE_adj_pendulum(const real_t* in, real_t* out) {
    int_t NX = 4;
    const real_t* x = in;
    const real_t* lambdaX = in + NX;
    const real_t* u  = in + 2*NX;

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

    vde_adj_pendulum(casadi_arg, casadi_res, casadi_iw, casadi_w, casadi_mem);
}


// The auto-generated Jacobian functions from CasADi:
void jac_fun_pendulum(const real_t* in, real_t* out) {
    int_t NX = 4;
    const real_t* x = in;
    const real_t* u  = in + NX;

    real_t* x_out = out;
    real_t* jac_out = out + NX;

    int casadi_mem = 0;
    int *casadi_iw = 0;
    real_t *casadi_w = 0;

    const real_t *casadi_arg[4];
    real_t *casadi_res[3];

    casadi_arg[0] = x;
    casadi_arg[1] = u;

    casadi_res[0] = x_out;
    casadi_res[1] = jac_out;

    jac_pendulum(casadi_arg, casadi_res, casadi_iw, casadi_w, casadi_mem);
}
