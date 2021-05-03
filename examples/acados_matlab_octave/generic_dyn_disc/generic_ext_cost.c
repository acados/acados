#ifdef __cplusplus
extern "C" {
#endif

#include <math.h>
#include "acados/utils/external_function_generic.h"
#include "blasfeo_d_blas.h"
#include "blasfeo_d_aux.h"

//void ext_cost(void *ext_fun, ext_fun_arg_t *type_in, void **in, ext_fun_arg_t *type_out, void **out)
int ext_cost(void **in, void **out, void *params)
{
    int ii, jj;

    const int nx = 8;
    const int nu = 3;

    // Inspired from https://github.com/acados/acados/blob/master/examples/c/no_interface_examples/mass_spring_nmpc_example.c#L158-L262
    // Manual implementation of weighted least squares
    // fun = 0.5 * (x'*Q*x + u'*R*u)
    // grad_x = Q*x
    // grad_u = R*u
    // hess_xx = Q
    // hess_xu = 0
    // hess_uu = R
    // hess_ux = 0

    // Where:
    // Q = diag(Qdiag)
    // R = R
//     const double *Qdiag = (double *)params;
    const double Qdiag[] = {1, 1, 1, 1, 1, 1, 1, 1};
    const double Rdiag[] = {2, 2, 2};

    // extract inputs
    // 0: [x], size: nx, type: BLASFEO_DVEC_ARGS
    struct blasfeo_dvec_args *x_args = in[0];
    struct blasfeo_dvec *x = x_args->x;
    int xi = x_args->xi; // offset in vector
    // 1: [u], size: nu, type: BLASFEO_DVEC_ARGS
    struct blasfeo_dvec_args *u_args = in[1];
    struct blasfeo_dvec *u = u_args->x;
    int ui = u_args->xi; // offset in vector
    
    //printf("q0 = %f, q1 = %f, q2 = %f, q3 = %f\n", Qdiag[0], Qdiag[1], Qdiag[2], Qdiag[3]);
    
    // extract outputs
    // 0: fun: COLMAJ
    double *fun = out[0];
    // 1: [grad_u; grad_x], size: nu+nx, type: BLASFEO_DVEC
    struct blasfeo_dvec *grad = out[1];
    // 2: [hess_uu, hess_ux; hess_xu, hess_xx], size: (nu+nx)*(nu+nx), type: BLASFEO_DMAT
    struct blasfeo_dmat *hess = out[2];

    // Hessian
    blasfeo_dgese(nu+nx, nu+nx, 0.0, hess, 0, 0);
    for(ii=0; ii<nu; ii++)
        BLASFEO_DMATEL(hess, ii, ii) = Rdiag[ii]; // R
    for(jj=0; ii<nu+nx; ii++, jj++)
        BLASFEO_DMATEL(hess, ii, ii) = Qdiag[jj]; // Q

    // gradient
    for(ii=0; ii<nu; ii++)
        BLASFEO_DVECEL(grad, ii) = BLASFEO_DMATEL(hess, ii, ii) * BLASFEO_DVECEL(u, ui+ii); // r
    for(ii=0; ii<nx; ii++)
        BLASFEO_DVECEL(grad, nu+ii) = BLASFEO_DMATEL(hess, nu+ii, nu+ii) * BLASFEO_DVECEL(x, xi+ii); // q

    // function
    *fun = 0.0;
    for(ii=0; ii<nu; ii++)
        *fun += BLASFEO_DVECEL(grad, ii) * BLASFEO_DVECEL(u, ui+ii); // r
    for(ii=0; ii<nx; ii++)
        *fun += BLASFEO_DVECEL(grad, nu+ii) * BLASFEO_DVECEL(x, xi+ii); // q

    *fun *= 0.5;

    return 0;

}

//void ext_costN(void *ext_fun, ext_fun_arg_t *type_in, void **in, ext_fun_arg_t *type_out, void **out)
int ext_costN(void **in, void **out, void *params)
{
    int ii, jj;

    const int nx = 8;
    const int nu = 0; // because of terminal cost

    // Inspired from https://github.com/acados/acados/blob/master/examples/c/no_interface_examples/mass_spring_nmpc_example.c#L158-L262
    // Manual implementation of weighted least squares
    // fun = 0.5 * (x'*Q*x + u'*R*u)
    // grad_x = Q*x
    // grad_u = R*u
    // hess_xx = Q
    // hess_xu = 0
    // hess_uu = R
    // hess_ux = 0

    // Where:
    // Q = diag(Qdiag)
//     const double *Qdiag = (double *)params;
    const double Qdiag[] = {1, 1, 1, 1, 1, 1, 1, 1};

    // extract inputs
    // 0: [x], size: nx, type: BLASFEO_DVEC_ARGS
    struct blasfeo_dvec_args *x_args = in[0];
    struct blasfeo_dvec *x = x_args->x;
    int xi = x_args->xi;
    // 1: [u], size: nu, type: BLASFEO_DVEC_ARGS
    struct blasfeo_dvec_args *u_args = in[1];
    struct blasfeo_dvec *u = u_args->x;
    int ui = u_args->xi;

    // extract outputs
    // 0: fun: COLMAJ
    double *fun = out[0];
    // 1: [grad_u; grad_x], size: nu+nx, type: BLASFEO_DVEC
    struct blasfeo_dvec *grad = out[1];
    // 2: [hess_uu, hess_ux; hess_xu, hess_xx], size: (nu+nx)*(nu+nx), type: BLASFEO_DMAT
    struct blasfeo_dmat *hess = out[2];

    // Hessian
    blasfeo_dgese(nx, nx, 0.0, hess, 0, 0);
    for(ii=0; ii<nx; ii++)
        BLASFEO_DMATEL(hess, ii, ii) = Qdiag[ii]; // Q


    // gradient
    for(ii=0; ii<nx; ii++)
        BLASFEO_DVECEL(grad, ii) = BLASFEO_DMATEL(hess, ii, ii) * BLASFEO_DVECEL(x, xi+ii); // q

    // function
    *fun = 0.0;
    for(ii=0; ii<nx; ii++)
        *fun += BLASFEO_DVECEL(grad, ii) * BLASFEO_DVECEL(x, xi+ii); // q

    *fun *= 0.5;

    return 0;

}


#ifdef __cplusplus
} /* extern "C" */
#endif
