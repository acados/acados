/* The model comes from \cite{Wirsching2006} */
#pragma GCC diagnostic ignored "-Wunused-parameter"
#pragma GCC diagnostic ignored "-Wunused-function"
#include "residual_u_eval_wrapper.h"
#include "acados_kitty_interface.h"

// wrapper to casadi autogen code for residual evaluation

void residual_u_eval_wrapper(nmpc_data* nmpc_data, const real_t* in, real_t* out, void *internal_mem) {

    const int nu = nmpc_data->NU;
    const int nr = nmpc_data->NU;  // !!!! NU residuals
    const int np = nmpc_data->NP;

    const real_t* u = in;
    const real_t* p = in + nu;
    const real_t* t0 = in + nu + np;
    const real_t* h = in + nu + np + 1;


    real_t* f  = out;
    real_t* dfdu = out + nr;
    real_t* fref = out + nr + nr*nu;

    // TODO(Andrea): need to remove this and allocate memory somewhere else
    int casadi_mem = 0;
    int *casadi_iw = 0;
    real_t *casadi_w = internal_mem;

    // this stuff is hard coded because dimensions won't change (I assume?)
    const real_t *casadi_arg[4];
    real_t *casadi_res[3];

    casadi_arg[0] = u;
    casadi_arg[1] = p;
    casadi_arg[2] = t0;
    casadi_arg[3] = h;

    casadi_res[0] = f;
    casadi_res[1] = dfdu;
    casadi_res[2] = fref;

    nmpc_data->res_u(casadi_arg, casadi_res, casadi_iw, casadi_w, casadi_mem);

}
