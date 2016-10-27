#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "blasfeo/include/blasfeo_target.h"
#include "blasfeo/include/blasfeo_common.h"
#include "blasfeo/include/blasfeo_i_aux.h"
#include "blasfeo/include/blasfeo_d_aux.h"
#include "blasfeo/include/blasfeo_d_kernel.h"
#include "blasfeo/include/blasfeo_d_blas.h"

#include "acados/sim/sim_lifted_irk_integrator.h"
#include "acados/utils/print.h"

#if TRIPLE_LOOP

#if CODE_GENERATION
#define DIM 6       // num_stages*NX
#define DIM_RHS 9   // NX+NU
#endif

real_t LU_system_ACADO(real_t* const A, int* const perm, int dim, int* nswaps) {
    real_t det;
    real_t swap;
    real_t valueMax;

#if !CODE_GENERATION
    int DIM = dim;
#else
    dim += 0;
#endif

    int i, j, k;
    int indexMax;
    int intSwap;

    for (i = 0; i < DIM; ++i) {
        perm[i] = i;
    }
    det = 1.0000000000000000e+00;
    for (i=0; i < (DIM-1); i++) {
        indexMax = i;
        valueMax = fabs(A[i*DIM+i]);
        for (j=(i+1); j < DIM; j++) {
            swap = fabs(A[i*DIM+j]);
            if (swap > valueMax) {
                indexMax = j;
                valueMax = swap;
            }
        }
        if (indexMax > i) {
            nswaps[0] += 1;
            for (k = 0; k < DIM; ++k) {
                swap = A[k*DIM+i];
                A[k*DIM+i] = A[k*DIM+indexMax];
                A[k*DIM+indexMax] = swap;
            }
            intSwap = perm[i];
            perm[i] = perm[indexMax];
            perm[indexMax] = intSwap;
        }
//        det *= A[i*DIM+i];
        for (j=i+1; j < DIM; j++) {
            A[i*DIM+j] = -A[i*DIM+j]/A[i*DIM+i];
            for (k=i+1; k < DIM; k++) {
                A[k*DIM+j] += A[i*DIM+j] * A[k*DIM+i];
            }
        }
    }
//    det *= A[DIM*DIM-1];
//    det = fabs(det);
    return det;
}

real_t solve_system_ACADO(real_t* const A, real_t* const b, int* const perm, int dim, int dim2) {
    int i, j, k;
    int index1;

#if !CODE_GENERATION
    int DIM = dim;
    int DIM_RHS = dim2;
#else
    dim += 0;
    dim2 += 0;
#endif
    real_t bPerm[DIM*DIM_RHS];
    real_t tmp_var;

    for (i = 0; i < DIM; ++i) {
        index1 = perm[i];
        for (j = 0; j < DIM_RHS; ++j) {
            bPerm[j*DIM+i] = b[j*DIM+index1];
        }
    }
    for (j = 1; j < DIM; ++j) {
        for (i = 0; i < j; ++i) {
            tmp_var = A[i*DIM+j];
            for (k = 0; k < DIM_RHS; ++k) {
                bPerm[k*DIM+j] += tmp_var*bPerm[k*DIM+i];
            }
        }
    }
    for (i = DIM-1; -1 < i; --i) {
        for (j = DIM-1; i < j; --j) {
            tmp_var = A[j*DIM+i];
            for (k = 0; k < DIM_RHS; ++k) {
                bPerm[k*DIM+i] -= tmp_var*bPerm[k*DIM+j];
            }
        }
        tmp_var = 1.0/A[i*(DIM+1)];
        for (k = 0; k < DIM_RHS; ++k) {
            bPerm[k*DIM+i] = tmp_var*bPerm[k*DIM+i];
        }
    }
    for (k = 0; k < DIM*DIM_RHS; ++k) {
        b[k] = bPerm[k];
    }
    return 0;
}


#endif

#if FIXED_STEP_SIZE == 1
// Integration step size.
real_t H_INT = 1.0/100;
// Number of integration steps.
#define NSTEPS  10
#endif

void sim_lifted_irk(const sim_in *in, sim_out *out, const sim_RK_opts *opts,
        sim_lifted_irk_memory *mem, sim_lifted_irk_workspace *work ) {
    int_t nx = in->nx;
    int_t nu = in->nu;
    int_t num_stages = opts->num_stages;
    int_t i, s1, s2, j, istep;
#if WARM_SWAP
    int_t iswap;
    int *ipiv_old = mem->ipiv;  // pivoting vector
    int_t ipiv_tmp[num_stages*nx];
#endif
#if FIXED_STEP_SIZE == 0
    real_t H_INT = in->step;
    int_t NSTEPS = in->nSteps;
#endif

    real_t *A_mat = opts->A_mat;
    real_t *b_vec = opts->b_vec;
//    real_t *c_vec = opts->c_vec;

//    print_matrix("stdout", A_mat, num_stages, num_stages);

    real_t *VDE_tmp = work->VDE_tmp;
    real_t *out_tmp = work->out_tmp;
    real_t *rhs_in = work->rhs_in;

    real_t *K_traj = mem->K_traj;
    real_t *DK_traj = mem->DK_traj;

    int *ipiv = work->ipiv;  // pivoting vector
    real_t *sys_mat = work->sys_mat;
    real_t *sys_sol = work->sys_sol;
#if !TRIPLE_LOOP
    struct d_strmat *str_mat = work->str_mat;
#if !TRANSPOSED
    struct d_strmat *str_sol = work->str_sol;
#else
    struct d_strmat *str_sol_t = work->str_sol_t;
#endif
#endif

    acado_timer timer, timer_la, timer_ad;
    real_t timing_la = 0.0;
    real_t timing_ad = 0.0;

    acado_tic(&timer);

    for (i = 0; i < nx; i++) out_tmp[i] = in->x[i];
    for (i = 0; i < nx*(nx+nu); i++) out_tmp[nx+i] = 0.0;  // sensitivities
    for (i = 0; i < nx; i++) out_tmp[nx+i*nx+i] = 1.0;     // sensitivities wrt x

    for (i = 0; i < nu; i++) rhs_in[nx*(1+nx+nu)+i] = in->u[i];

    // Newton step of the collocation variables with respect to the inputs:
    for (istep = 0; istep < NSTEPS; istep++) {
        for (s1 = 0; s1 < num_stages; s1++) {
            for (j = 0; j < nx; j++) {  // step in X
                for (i = 0; i < nx; i++) {
                    K_traj[(istep*num_stages+s1)*nx+i] +=
                            DK_traj[(istep*num_stages+s1)*nx*(nx+nu)+j*nx+i]*
                            (in->x[j]-mem->x[j]);  // RK step
                }
                mem->x[j] = in->x[j];
            }
            for (j = 0; j < nu; j++) {  // step in U
                for (i = 0; i < nx; i++) {
                    K_traj[(istep*num_stages+s1)*nx+i] +=
                            DK_traj[(istep*num_stages+s1)*nx*(nx+nu)+(nx+j)*nx+i]*
                            (in->u[j]-mem->u[j]);  // RK step
                }
                mem->u[j] = in->u[j];
            }
        }
    }

    mem->nswaps += 0;
    for (istep = 0; istep < NSTEPS; istep++) {
        // form exact linear system matrix (explicit ODE case):
        for (i = 0; i < num_stages*nx*num_stages*nx; i++ ) sys_mat[i] = 0.0;
#if WARM_SWAP
        for (i = 0; i < num_stages*nx; i++ ) {
            iswap = ipiv_old[i];
            sys_mat[i*num_stages*nx+iswap] = 1.0;  // identity matrix
        }
#else
        for (i = 0; i < num_stages*nx; i++ ) sys_mat[i*(num_stages*nx+1)] = 1.0;  // identity matrix
#endif

        for (s1 = 0; s1 < num_stages; s1++) {
            for (i = 0; i < nx; i++) {
                rhs_in[i] = out_tmp[i];
            }
            for (i = 0; i < nu; i++) rhs_in[nx+i] = in->u[i];
            for (s2 = 0; s2 < num_stages; s2++) {
                for (i = 0; i < nx; i++) {
                    rhs_in[i] += H_INT*A_mat[s2*num_stages+s1]*K_traj[istep*num_stages*nx+s2*nx+i];
                }
            }
            acado_tic(&timer_ad);
            in->jac_fun(rhs_in, VDE_tmp);  // k evaluation
            timing_ad += acado_toc(&timer_ad);

            // put VDE_tmp in sys_mat:
            for (s2 = 0; s2 < num_stages; s2++) {
#if WARM_SWAP
                for (i = 0; i < nx; i++) {
                    iswap = ipiv_old[s1*nx+i];
                    for (j = 0; j < nx; j++) {
                        sys_mat[(s2*nx+j)*num_stages*nx+iswap] -=
                                H_INT*A_mat[s2*num_stages+s1]*VDE_tmp[nx+j*nx+i];
                    }
                }
#else
                for (j = 0; j < nx; j++) {
                    for (i = 0; i < nx; i++) {
                        sys_mat[(s2*nx+j)*num_stages*nx+s1*nx+i] -=
                                H_INT*A_mat[s2*num_stages+s1]*VDE_tmp[nx+j*nx+i];
                    }
                }
#endif
            }
        }

        acado_tic(&timer_la);
#if TRIPLE_LOOP
        LU_system_ACADO(sys_mat, ipiv, num_stages*nx, &mem->nswaps);
#else
        // ---- BLASFEO: LU factorization ----
        d_cvt_mat2strmat(num_stages*nx, num_stages*nx, sys_mat, num_stages*nx,
                str_mat, 0, 0);  // mat2strmat
        dgetrf_libstr(num_stages*nx, num_stages*nx, str_mat, 0, 0, str_mat, 0,
                0, ipiv);  // Gauss elimination
        // ---- BLASFEO: LU factorization ----
#endif
        timing_la += acado_toc(&timer_la);

        for (s1 = 0; s1 < num_stages; s1++) {
            for (i = 0; i < nx*(1+nx+nu); i++) {
                rhs_in[i] = out_tmp[i];
            }
            for (s2 = 0; s2 < num_stages; s2++) {
                for (i = 0; i < nx; i++) {
                    rhs_in[i] += H_INT*A_mat[s2*num_stages+s1]*K_traj[istep*num_stages*nx+s2*nx+i];
                }
            }
            acado_tic(&timer_ad);
            in->VDE_fun(rhs_in, VDE_tmp);  // k evaluation
            timing_ad += acado_toc(&timer_ad);

            // put VDE_tmp in sys_sol:
#if WARM_SWAP
            for (i = 0; i < nx; i++) {
                iswap = ipiv_old[s1*nx+i];
                for (j = 0; j < nx+nu+1; j++) {
                    sys_sol[j*num_stages*nx+iswap] = VDE_tmp[j*nx+i];
                }
                sys_sol[iswap] -= K_traj[istep*num_stages*nx+s1*nx+i];
            }
#else
            for (j = 0; j < nx+nu+1; j++) {
                for (i = 0; i < nx; i++) {
                    sys_sol[j*num_stages*nx+s1*nx+i] = VDE_tmp[j*nx+i];
                }
            }
            for (i = 0; i < nx; i++) {
                sys_sol[s1*nx+i] -= K_traj[istep*num_stages*nx+s1*nx+i];
            }
#endif
        }

        acado_tic(&timer_la);
#if TRIPLE_LOOP
        solve_system_ACADO(sys_mat, sys_sol, ipiv, num_stages*nx, nx+nu+1);
#else
#if TRANSPOSED
        // ---- BLASFEO: row transformations + backsolve ----
        d_cvt_tran_mat2strmat(num_stages*nx, nx+nu+1, sys_sol, num_stages*nx,
                str_sol_t, 0, 0);  // mat2strmat
        dcolpe_libstr(num_stages*nx, ipiv, str_sol_t);  // col permutations
        dtrsm_rltu_libstr(nx+nu+1, num_stages*nx, 1.0, str_mat, 0, 0, str_sol_t,
                0, 0, str_sol_t, 0, 0);  // L backsolve
        dtrsm_rutn_libstr(nx+nu+1, num_stages*nx, 1.0, str_mat, 0, 0, str_sol_t,
                0, 0, str_sol_t, 0, 0);  // U backsolve
        d_cvt_tran_strmat2mat(nx+nu+1, num_stages*nx, str_sol_t, 0, 0, sys_sol,
                num_stages*nx);  // strmat2mat
        // ---- BLASFEO: row transformations + backsolve ----
#else
        // ---- BLASFEO: row transformations + backsolve ----
        d_cvt_mat2strmat(num_stages*nx, nx+nu+1, sys_sol, num_stages*nx,
                str_sol, 0, 0);  // mat2strmat
        drowpe_libstr(num_stages*nx, ipiv, str_sol);  // row permutations
        dtrsm_llnu_libstr(num_stages*nx, nx+nu+1, 1.0, str_mat, 0, 0, str_sol,
                0, 0, str_sol, 0, 0);  // L backsolve
        dtrsm_lunn_libstr(num_stages*nx, nx+nu+1, 1.0, str_mat, 0, 0, str_sol,
                0, 0, str_sol, 0, 0);  // U backsolve
        d_cvt_strmat2mat(num_stages*nx, nx+nu+1, str_sol, 0, 0, sys_sol,
                num_stages*nx);  // strmat2mat
        // ---- BLASFEO: row transformations + backsolve ----
#endif
#endif
#if WARM_SWAP
#if TRIPLE_LOOP
        for (i = 0; i < num_stages*nx; i++) ipiv_tmp[i] = ipiv[i];
#else
        for (i = 0; i < num_stages*nx; i++) ipiv_tmp[i] = i;
        for (i = 0; i < num_stages*nx; i++) {
            j = ipiv[i];
            if (j != i) {
                iswap = ipiv_tmp[i];
                ipiv_tmp[i] = ipiv_tmp[j];
                ipiv_tmp[j] = iswap;
            }
        }
#endif

//        fprintf(stdout, "ipiv_old: ");
//        for (i = 0; i < num_stages*nx; i++) fprintf(stdout, "%d ", ipiv_old[i]);
//        fprintf(stdout, "\n");
//        fprintf(stdout, "ipiv: ");
//        for (i = 0; i < num_stages*nx; i++) fprintf(stdout, "%d ", ipiv_tmp[i]);
//        fprintf(stdout, "\n");
        for (i = 0; i < num_stages*nx; i++) {
            ipiv_tmp[i] = ipiv_old[ipiv_tmp[i]];
        }
//        fprintf(stdout, "ipiv_tmp: ");
//        for (i = 0; i < num_stages*nx; i++) fprintf(stdout, "%d ", ipiv_tmp[i]);
//        fprintf(stdout, "\n");
        for (i = 0; i < num_stages*nx; i++) {
            ipiv_old[i] = ipiv_tmp[i];
        }
#endif
        timing_la += acado_toc(&timer_la);

        // Newton step of the collocation variables
        for (i = 0; i < num_stages*nx; i++) {
            K_traj[istep*num_stages*nx+i] += sys_sol[i];
        }
        for (s1 = 0; s1 < num_stages; s1++) {
            for (i = 0; i < nx; i++) {
                out_tmp[i] += H_INT*b_vec[s1]*K_traj[istep*num_stages*nx+s1*nx+i];  // RK step
            }
        }

        // Sensitivities collocation variables
        for (s1 = 0; s1 < num_stages; s1++) {
            for (j = 0; j < (nx+nu); j++) {
                for (i = 0; i < nx; i++) {
                    DK_traj[(istep*num_stages+s1)*nx*(nx+nu)+j*nx+i] =
                            sys_sol[(j+1)*num_stages*nx+s1*nx+i];
                }
            }
        }
        for (s1 = 0; s1 < num_stages; s1++) {
            for (i = 0; i < nx*(nx+nu); i++) {
                out_tmp[nx+i] += H_INT*b_vec[s1]*
                        DK_traj[(istep*num_stages+s1)*nx*(nx+nu)+i];  // RK step
            }
        }
    }
    for (i = 0; i < nx; i++)    out->xn[i] = out_tmp[i];
    for (i = 0; i < nx*nx; i++) out->Sx[i] = out_tmp[nx+i];
    for (i = 0; i < nx*nu; i++) out->Su[i] = out_tmp[nx+nx*nx+i];

    out->info->CPUtime = acado_toc(&timer);
    out->info->LAtime = timing_la;
    out->info->ADtime = timing_ad;
}


void sim_lifted_irk_create_workspace(const sim_in *in, sim_RK_opts *opts,
        sim_lifted_irk_workspace *work) {
    int_t nx = in->nx;
    int_t nu = in->nu;
    int_t num_stages = opts->num_stages;

    work->rhs_in = malloc(sizeof(*work->rhs_in) * (nx*(1+nx+nu)+nu));
    work->VDE_tmp = malloc(sizeof(*work->VDE_tmp) * (nx*(1+nx+nu)));
    work->out_tmp = malloc(sizeof(*work->out_tmp) * (nx*(1+nx+nu)));
    work->ipiv = malloc(sizeof(*work->ipiv) * (num_stages*nx));
    work->sys_mat = malloc(sizeof(*work->sys_mat) * (num_stages*nx)*(num_stages*nx));
    work->sys_sol = malloc(sizeof(*work->sys_sol) * (num_stages*nx)*(1+nx+nu));

#if !TRIPLE_LOOP
    // matrices in matrix struct format:
    int size_strmat = d_size_strmat(num_stages*nx, num_stages*nx)
            + d_size_strmat(num_stages*nx, 1+nx+nu)
            + d_size_strmat(1+nx+nu, num_stages*nx);
    void *memory_strmat;
    v_zeros_align(&memory_strmat, size_strmat);
    char *ptr_memory_strmat = (char *) memory_strmat;

    d_create_strmat(num_stages*nx, num_stages*nx, work->str_mat, ptr_memory_strmat);
    ptr_memory_strmat += work->str_mat->memory_size;

    d_create_strmat(num_stages*nx, 1+nx+nu, work->str_sol, ptr_memory_strmat);
    ptr_memory_strmat += work->str_sol->memory_size;

    d_create_strmat(1+nx+nu, num_stages*nx, work->str_sol_t, ptr_memory_strmat);
    ptr_memory_strmat += work->str_sol_t->memory_size;
#endif
}


void sim_lifted_irk_create_memory(const sim_in *in, sim_RK_opts *opts,
        sim_lifted_irk_memory *mem) {
    int_t i;
    int_t nx = in->nx;
    int_t nu = in->nu;
    int_t nSteps = in->nSteps;
    int_t num_stages = opts->num_stages;

    mem->K_traj = malloc(sizeof(*mem->K_traj) * (nSteps*num_stages*nx));
    mem->DK_traj = malloc(sizeof(*mem->DK_traj) * (nSteps*num_stages*nx*(nx+nu)));
    mem->x = malloc(sizeof(*mem->x) * nx);
    mem->u = malloc(sizeof(*mem->u) * nu);

#if WARM_SWAP
    mem->ipiv = malloc(sizeof(*mem->ipiv) * (num_stages*nx));
    for ( i = 0; i < num_stages*nx; i++ ) mem->ipiv[i] = i;
#endif

    for ( i = 0; i < nSteps*num_stages*nx; i++ ) mem->K_traj[i] = 0.0;
    for ( i = 0; i < nSteps*num_stages*nx*(nx+nu); i++ ) mem->DK_traj[i] = 0.0;
}


void sim_irk_create_opts(const int_t num_stages, const char* name, sim_RK_opts *opts) {
    if ( strcmp(name, "Gauss") == 0 ) {  // GAUSS METHODS
        if ( num_stages == 1 ) {
            opts->num_stages = 1;       // GL2
            opts->A_mat = malloc(sizeof(*opts->A_mat) * (num_stages*num_stages));
            opts->b_vec = malloc(sizeof(*opts->b_vec) * (num_stages));
            opts->c_vec = malloc(sizeof(*opts->c_vec) * (num_stages));
            opts->A_mat[0] = 1.0/2.0;
            opts->b_vec[0] = 1.0;
            opts->c_vec[0] = 1.0/2.0;
        } else if ( num_stages == 2 ) {
            opts->num_stages = 2;       // GL4
            opts->A_mat = malloc(sizeof(*opts->A_mat) * (num_stages*num_stages));
            opts->b_vec = malloc(sizeof(*opts->b_vec) * (num_stages));
            opts->c_vec = malloc(sizeof(*opts->c_vec) * (num_stages));

            memcpy(opts->A_mat,
                    ((real_t[]) {1.0/4.0, (1.0/4.0-sqrt(3.0)/6.0),
                    (1.0/4.0+sqrt(3.0)/6.0), 1.0/4.0}),
                    sizeof(*opts->A_mat) * (num_stages*num_stages));
            memcpy(opts->b_vec,
                    ((real_t[]) {1.0/2.0, 1.0/2.0}), sizeof(*opts->b_vec) * (num_stages));
            memcpy(opts->c_vec,
                    ((real_t[]) {1.0/2.0+sqrt(3.0)/6.0, 1.0/2.0-sqrt(3.0)/6.0}),
                    sizeof(*opts->c_vec) * (num_stages));
        } else if ( num_stages == 3 ) {
            opts->num_stages = 3;       // GL6
            opts->A_mat = malloc(sizeof(*opts->A_mat) * (num_stages*num_stages));
            opts->b_vec = malloc(sizeof(*opts->b_vec) * (num_stages));
            opts->c_vec = malloc(sizeof(*opts->c_vec) * (num_stages));

            memcpy(opts->A_mat,
                    ((real_t[]) {5.0/36.0, 5.0/36.0+1.0/24.0*sqrt(15.0),
                    5.0/36.0+1.0/30.0*sqrt(15.0), 2.0/9.0-1.0/15.0*sqrt(15.0),
                    2.0/9.0, 2.0/9.0+1.0/15.0*sqrt(15.0), 5.0/36.0-1.0/30.0*sqrt(15.0),
                    5.0/36.0-1.0/24.0*sqrt(15.0), 5.0/36.0}),
                    sizeof(*opts->A_mat) * (num_stages*num_stages));
            memcpy(opts->b_vec,
                    ((real_t[]) {5.0/18.0, 4.0/9.0, 5.0/18.0}),
                    sizeof(*opts->b_vec) * (num_stages));
            memcpy(opts->c_vec,
                    ((real_t[]) {1.0/2.0-sqrt(15.0)/10.0, 1.0/2.0, 1.0/2.0+sqrt(15.0)/10.0}),
                    sizeof(*opts->c_vec) * (num_stages));
        } else {
            // throw error somehow?
        }
    } else if ( strcmp(name, "Radau") == 0 ) {
        // TODO(rien): add Radau IIA collocation schemes
    } else {
        // throw error somehow?
    }
}
