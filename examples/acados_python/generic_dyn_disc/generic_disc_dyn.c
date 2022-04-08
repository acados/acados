#ifdef __cplusplus
extern "C" {
#endif

#include "acados/utils/external_function_generic.h"
#include "acados/utils/math.h"
#include "blasfeo_d_blas.h"
#include "blasfeo_d_aux.h"    
#include "blasfeo_d_aux_ext_dep.h"

// auxilary function to compute discrete dynamics, see below
static void mass_spring_system(double Ts, int nx, int nu, double *A, double *B, double *b)
{
    int nx2 = nx * nx;

    int info = 0;

    int pp = nx / 2;  // number of masses

    /***********************************************
    * build the continuous time system
    ***********************************************/

    double *T;
    d_zeros(&T, pp, pp);

    for (int ii = 0; ii < pp; ii++) T[ii * (pp + 1)] = -3;
    for (int ii = 0; ii < pp - 1; ii++) T[ii * (pp + 1) + 1] = 1;
    for (int ii = 1; ii < pp; ii++) T[ii * (pp + 1) - 1] = 1;

    double *Z;
    d_zeros(&Z, pp, pp);
    double *I;
    d_zeros(&I, pp, pp);
    for (int ii = 0; ii < pp; ii++) I[ii * (pp + 1)] = 1.0;  // I = eye(pp);
    double *Ac;
    d_zeros(&Ac, nx, nx);
    dmcopy(pp, pp, Z, pp, Ac, nx);
    dmcopy(pp, pp, T, pp, Ac + pp, nx);
    dmcopy(pp, pp, I, pp, Ac + pp * nx, nx);
    dmcopy(pp, pp, Z, pp, Ac + pp * (nx + 1), nx);
    free(T);
    free(Z);
    free(I);

    d_zeros(&I, nu, nu);
    for (int ii = 0; ii < nu; ii++) I[ii * (nu + 1)] = 1.0;  // I = eye(nu);
    double *Bc;
    d_zeros(&Bc, nx, nu);
    dmcopy(nu, nu, I, nu, Bc + pp, nx);
    free(I);

    /************************************************
    * compute the discrete time system
    ************************************************/

    double *bb;
    d_zeros(&bb, nx, 1);
    dmcopy(nx, 1, bb, nx, b, nx);

    dmcopy(nx, nx, Ac, nx, A, nx);
    dscal_3l(nx2, Ts, A);
    expm(nx, A);

    d_zeros(&T, nx, nx);
    d_zeros(&I, nx, nx);
    for (int ii = 0; ii < nx; ii++) I[ii * (nx + 1)] = 1.0;  // I = eye(nx);
    dmcopy(nx, nx, A, nx, T, nx);
    daxpy_3l(nx2, -1.0, I, T);
    dgemm_nn_3l(nx, nu, nx, T, nx, Bc, nx, B, nx);
    free(T);
    free(I);

    int *ipiv = (int *)malloc(nx * sizeof(int));
    dgesv_3l(nx, nu, Ac, nx, ipiv, B, nx, &info);
    free(ipiv);

    free(Ac);
    free(Bc);
    free(bb);
}


int disc_dyn_fun_jac_hess(void **in, void **out, void *params)
{
    int ii;

    int nu = 3;
    int nx = 8;

    // compute mass sping dynamics
    double *A;
    d_zeros(&A, nx, nx);  // states update matrix
    double *B;
    d_zeros(&B, nx, nu);  // inputs matrix
    double *b;
    d_zeros(&b, nx, 1);  // states offset
    double Ts = 0.5;

    mass_spring_system(Ts, nx, nu, A, B, b);

    for (ii=0; ii<nx; ii++)
        b[ii] = 0.0;

    // extract inputs
    // 0: [x], size: nx, type: BLASFEO_DVEC_ARGS
    struct blasfeo_dvec_args *x_args = in[0];
    struct blasfeo_dvec *x = x_args->x;
    int xi = x_args->xi; // offset in vector
    // 1: [u], size: nu, type: BLASFEO_DVEC_ARGS
    struct blasfeo_dvec_args *u_args = in[1];
    struct blasfeo_dvec *u = u_args->x;
    int ui = u_args->xi; // offset in vector

    // extract outputs
    // 0: [fun], size: nx1, type: BLASFEO_DVEC_ARGS
    struct blasfeo_dvec_args *f_args = out[0];
    struct blasfeo_dvec *fun = f_args->x;
    // 1: [jac_u'; jac_x'], size: (nu+nx)*nx1, type: BLASFEO_DMAT_ARGS
    struct blasfeo_dmat_args *j_args = out[1];
    struct blasfeo_dmat *jac = j_args->A;
    // 2: [hess_ux], size: (nu+nx)*(nu+nx)
    struct blasfeo_dmat_args *h_args = out[2];
    struct blasfeo_dmat *hess = h_args->A;
    // hess
    blasfeo_dgese(nu+nx, nu+nx, 0.0, hess, 0, 0);

    // jac
    blasfeo_pack_tran_dmat(nx, nu, B, nx, jac, 0, 0);
    blasfeo_pack_tran_dmat(nx, nx, A, nx, jac, nu, 0);

    // fun
    blasfeo_dgemv_t(nu, nx, 1.0, jac, 0, 0, u, ui, 0.0, fun, 0, fun, 0);
    blasfeo_dgemv_t(nx, nx, 1.0, jac, nu, 0, x, xi, 1.0, fun, 0, fun, 0);

    // free memory
    free(A);
    free(B);
    free(b);

    return 0;
}



int disc_dyn_fun_jac(void **in, void **out, void *params)
{
    int ii;

    int nu = 3;
    int nx = 8;

    // compute mass sping dynamics
    double *A;
    d_zeros(&A, nx, nx);  // states update matrix
    double *B;
    d_zeros(&B, nx, nu);  // inputs matrix
    double *b;
    d_zeros(&b, nx, 1);  // states offset
    double Ts = 0.5;

    mass_spring_system(Ts, nx, nu, A, B, b);

    for (ii=0; ii<nx; ii++)
        b[ii] = 0.0;

    // extract inputs
    // 0: [x], size: nx, type: BLASFEO_DVEC_ARGS
    struct blasfeo_dvec_args *x_args = in[0];
    struct blasfeo_dvec *x = x_args->x;
    int xi = x_args->xi; // offset in vector
    // 1: [u], size: nu, type: BLASFEO_DVEC_ARGS
    struct blasfeo_dvec_args *u_args = in[1];
    struct blasfeo_dvec *u = u_args->x;
    int ui = u_args->xi; // offset in vector

    // extract outputs
    // 0: [fun], size: nx1, type: BLASFEO_DVEC_ARGS
    struct blasfeo_dvec_args *f_args = out[0];
    struct blasfeo_dvec *fun = f_args->x;
    // 1: [jac_u'; jac_x'], size: (nu+nx)*nx1, type: BLASFEO_DMAT_ARGS
    struct blasfeo_dmat_args *j_args = out[1];
    struct blasfeo_dmat *jac = j_args->A;

    // jac
    blasfeo_pack_tran_dmat(nx, nu, B, nx, jac, 0, 0);
    blasfeo_pack_tran_dmat(nx, nx, A, nx, jac, nu, 0);

    // fun
    blasfeo_dgemv_t(nu, nx, 1.0, jac, 0, 0, u, ui, 0.0, fun, 0, fun, 0);
    blasfeo_dgemv_t(nx, nx, 1.0, jac, nu, 0, x, xi, 1.0, fun, 0, fun, 0);

    // free memory
    free(A);
    free(B);
    free(b);

    return 0;
}


int disc_dyn_fun(void **in, void **out, void *params)
{
    int ii;

    int nu = 3;
    int nx = 8;

    // compute mass sping dynamics
    double *A;
    d_zeros(&A, nx, nx);  // states update matrix
    double *B;
    d_zeros(&B, nx, nu);  // inputs matrix
    double *b;
    d_zeros(&b, nx, 1);  // states offset
    double Ts = 0.5;

    mass_spring_system(Ts, nx, nu, A, B, b);

    for (ii=0; ii<nx; ii++)
        b[ii] = 0.0;

    // extract inputs
    // 0: [x], size: nx, type: BLASFEO_DVEC_ARGS
    struct blasfeo_dvec_args *x_args = in[0];
    struct blasfeo_dvec *x = x_args->x;
    int xi = x_args->xi; // offset in vector
    // 1: [u], size: nu, type: BLASFEO_DVEC_ARGS
    struct blasfeo_dvec_args *u_args = in[1];
    struct blasfeo_dvec *u = u_args->x;
    int ui = u_args->xi; // offset in vector

    // extract outputs
    // 0: [fun], size: nx1, type: BLASFEO_DVEC_ARGS
    struct blasfeo_dvec_args *f_args = out[0];
    struct blasfeo_dvec *fun = f_args->x;

    // tmp_jac
    struct blasfeo_dmat tmp_jac;            // matrix structure
    blasfeo_allocate_dmat(nx, nx+nu, &tmp_jac);  // allocate and assign memory
    blasfeo_pack_tran_dmat(nx, nu, B, nx, &tmp_jac, 0, 0);
    blasfeo_pack_tran_dmat(nx, nx, A, nx, &tmp_jac, nu, 0);

    // fun
    blasfeo_dgemv_t(nu, nx, 1.0, &tmp_jac, 0, 0, u, ui, 0.0, fun, 0, fun, 0);
    blasfeo_dgemv_t(nx, nx, 1.0, &tmp_jac, nu, 0, x, xi, 1.0, fun, 0, fun, 0);

    // free memory
    blasfeo_free_dmat(&tmp_jac);
    free(A);
    free(B);
    free(b);

    return 0;
}

#ifdef __cplusplus
} /* extern "C" */
#endif
