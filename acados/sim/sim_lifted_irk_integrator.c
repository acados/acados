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

real_t LU_system_ACADO(real_t* const A, int* const perm, int dim) {
    real_t det;
    real_t swap;
    real_t valueMax;
//    printf("LU_system_ACADO, dim: %d \n", dim);

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
//    printf("solve_system_ACADO, dim: %d, dim2: %d \n", dim, dim2);

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

real_t solve_system_trans_ACADO(real_t* const A, real_t* const b,
        int* const perm, int dim, int dim2) {
    int i, j, k;
    int index1;
//    printf("solve_system_trans_ACADO, dim: %d, dim2: %d \n", dim, dim2);

#if !CODE_GENERATION
    int DIM = dim;
    int DIM_RHS = dim2;
#else
    dim += 0;
    dim2 += 0;
#endif
    real_t bPerm[DIM*DIM_RHS];
    real_t tmp_var;


    for (k = 0; k < DIM*DIM_RHS; ++k) {
        bPerm[k] = b[k];
    }
    for (i = 0; i < DIM; i++) {
        for (j = 0; j < i; j++) {
            tmp_var = A[i*DIM+j];
            for (k = 0; k < DIM_RHS; ++k) {
                bPerm[k*DIM+i] -= tmp_var*bPerm[k*DIM+j];
            }
        }
        tmp_var = 1.0/A[i*(DIM+1)];
        for (k = 0; k < DIM_RHS; ++k) {
            bPerm[k*DIM+i] = tmp_var*bPerm[k*DIM+i];
        }
    }
    for (j = DIM-1; j > -1; --j) {
        for (i = DIM-1; i > j; --i) {
            tmp_var = A[j*DIM+i];
            for (k = 0; k < DIM_RHS; ++k) {
                bPerm[k*DIM+j] += tmp_var*bPerm[k*DIM+i];
            }
        }
    }
    for (i = 0; i < DIM; ++i) {
        index1 = perm[i];
        for (j = 0; j < DIM_RHS; ++j) {
            b[j*DIM+index1] = bPerm[j*DIM+i];
        }
    }
    return 0;
}


#endif

void transform_mat(real_t *mat, real_t *trans, real_t *mat_trans,
        const int_t stages, const int_t  n, const int_t m) {
    for (int_t s1 = 0; s1 < stages; s1++) {
        for (int_t i = 0; i < n; i++) {
            for (int_t j = 0; j < m; j++) {
                mat_trans[j*(stages*n)+s1*n+i] = 0.0;
                for (int_t s2 = 0; s2 < stages; s2++) {
                    mat_trans[j*(stages*n)+s1*n+i] +=
                            trans[s2*stages+s1]*mat[j*(stages*n)+s2*n+i];
                }
            }
        }
    }
}

void construct_subsystems(real_t *mat, real_t **mat2,
        const int_t stages, const int_t n, const int_t m) {
    int_t idx = 0;
    for (int_t s1 = 0; s1 < stages; s1++) {
        if ((s1+1) < stages) {  // complex conjugate pair of eigenvalues
            for (int_t j = 0; j < m; j++) {
                for (int_t i = 0; i < n; i++) {
                    mat2[idx][j*2*n+i] = mat[j*(stages*n)+s1*n+i];
                }
            }
            s1++;
            for (int_t j = 0; j < m; j++) {
                for (int_t i = 0; i < n; i++) {
                    mat2[idx][j*2*n+n+i] = mat[j*(stages*n)+s1*n+i];
                }
            }
        } else {  // real eigenvalue
            for (int_t j = 0; j < m; j++) {
                for (int_t i = 0; i < n; i++) {
                    mat2[idx][j*n+i] = mat[j*(stages*n)+s1*n+i];
                }
            }
        }
        idx++;
    }
}

void construct_subsystems_combined(real_t *mat, real_t **mat2,
        const int_t stages, const int_t n, const int_t m, real_t **transf,
        const sim_in *in, bool trans) {
    int_t idx = 0;
    real_t H_INT = in->step;
    sim_RK_opts *opts = in->opts;
    real_t tmp_eig3;
    for (int_t s1 = 0; s1 < stages; s1++) {
        if ((s1+1) < stages) {  // complex conjugate pair of eigenvalues
            tmp_eig3 = H_INT/opts->scheme.eig[s1+1];
            for (int_t j = 0; j < m; j++) {
                for (int_t i = 0; i < n; i++) {
                    if (trans)  mat2[idx+1][j*n+i] = tmp_eig3*mat[j*(stages*n)+s1*n+i];
                    else        mat2[idx+1][j*n+i] = -tmp_eig3*mat[j*(stages*n)+s1*n+i];
                }
            }
//            print_matrix_name("stdout", "transf_mat", transf[idx+1], n, n);
            for (int_t j = 0; j < m; j++) {
                for (int_t i = 0; i < n; i++) {
                    mat2[idx][j*n+i] = 0.0;
                    for (int_t k = 0; k < n; k++) {
                        if (trans) {
                            mat2[idx][j*n+i] += transf[idx+1][i*n+k]*mat[j*(stages*n)+s1*n+k];
                        } else {
                            mat2[idx][j*n+i] += transf[idx+1][k*n+i]*mat[j*(stages*n)+s1*n+k];
                        }
                    }
                    if (trans)  mat2[idx][j*n+i] -= mat[j*(stages*n)+(s1+1)*n+i];
                    else        mat2[idx][j*n+i] += mat[j*(stages*n)+(s1+1)*n+i];
                }
            }
            s1++;
            idx++;
        } else {  // real eigenvalue
            for (int_t j = 0; j < m; j++) {
                for (int_t i = 0; i < n; i++) {
                    mat2[idx][j*n+i] = mat[j*(stages*n)+s1*n+i];
                }
            }
        }
        idx++;
    }
}

void destruct_subsystems(real_t *mat, real_t **mat2,
        const int_t stages, const int_t n, const int_t m) {
    int_t idx = 0;
    for (int_t s1 = 0; s1 < stages; s1++) {
        if ((s1+1) < stages) {  // complex conjugate pair of eigenvalues
            for (int_t j = 0; j < m; j++) {
                for (int_t i = 0; i < n; i++) {
                    mat[j*(stages*n)+s1*n+i] = mat2[idx][j*2*n+i];
                }
            }
            s1++;
            for (int_t j = 0; j < m; j++) {
                for (int_t i = 0; i < n; i++) {
                    mat[j*(stages*n)+s1*n+i] = mat2[idx][j*2*n+n+i];
                }
            }
        } else {  // real eigenvalue
            for (int_t j = 0; j < m; j++) {
                for (int_t i = 0; i < n; i++) {
                    mat[j*(stages*n)+s1*n+i] = mat2[idx][j*n+i];
                }
            }
        }
        idx++;
    }
}

void destruct_subsystems_combined(real_t *mat, real_t **mat2,
        const int_t stages, const int_t n, const int_t m,
        real_t **transf, bool trans) {
    int_t idx = 0;
    for (int_t s1 = 0; s1 < stages; s1++) {
        if ((s1+1) < stages) {  // complex conjugate pair of eigenvalues
            for (int_t j = 0; j < m; j++) {
                for (int_t i = 0; i < n; i++) {
                    mat[j*(stages*n)+s1*n+i] = mat2[idx][j*n+i];
                }
            }
//            print_matrix_name("stdout", "transf_mat", transf[idx+1], n, n);
            for (int_t j = 0; j < m; j++) {
                for (int_t i = 0; i < n; i++) {
                    mat[j*(stages*n)+(s1+1)*n+i] = mat2[idx+1][j*n+i];
                    for (int_t k = 0; k < n; k++) {
                        if (trans)  {
                            mat[j*(stages*n)+(s1+1)*n+i] -=
                                    transf[idx+1][i*n+k]*mat2[idx][j*n+k];
                        } else {
                            mat[j*(stages*n)+(s1+1)*n+i] +=
                                    transf[idx+1][k*n+i]*mat2[idx][j*n+k];
                        }
                    }
                }
            }
            s1++;
            idx++;
        } else {  // real eigenvalue
            for (int_t j = 0; j < m; j++) {
                for (int_t i = 0; i < n; i++) {
                    mat[j*(stages*n)+s1*n+i] = mat2[idx][j*n+i];
                }
            }
        }
        idx++;
    }
}

void form_linear_system_matrix(int_t istep, const sim_in *in, sim_lifted_irk_memory *mem,
        sim_lifted_irk_workspace *work, real_t *sys_mat, real_t **sys_mat2, real_t timing_ad) {
    int_t nx = in->nx;
    int_t nu = in->nu;
    real_t H_INT = in->step;
    sim_RK_opts *opts = in->opts;
    int_t num_stages = opts->num_stages;
    real_t *A_mat = opts->A_mat;

    real_t tmp_eig, tmp_eig2, tmp_eig3;
    acado_timer timer_ad;

    real_t *tmp_mat;
    real_t *out_tmp = work->out_tmp;
    real_t *rhs_in = work->rhs_in;
    real_t *jac_tmp = work->jac_tmp;

    real_t **jac_traj = mem->jac_traj;
    real_t *K_traj = mem->K_traj;

    int_t i, j, k, s1, s2;
    if ((opts->scheme.type == simplified_in || opts->scheme.type == simplified_inis
            || opts->scheme.type == simplified_inis2) && istep == 0) {
        int idx = 0;
        for (s1 = 0; s1 < num_stages; s1++) {
            if ((s1+1) == num_stages) {  // real eigenvalue
                for (i = 0; i < nx*nx; i++) sys_mat2[idx][i] = 0.0;
                tmp_eig = 1.0/H_INT*opts->scheme.eig[s1];
                for (i = 0; i < nx; i++) {
                    sys_mat2[idx][i*(nx+1)] = tmp_eig;
                }
            } else if (opts->scheme.type != simplified_inis2) {  // complex conjugate pair
                for (i = 0; i < 4*nx*nx; i++) sys_mat2[idx][i] = 0.0;
                tmp_eig = 1.0/H_INT*opts->scheme.eig[s1];
                tmp_eig2 = 1.0/H_INT*opts->scheme.eig[s1+1];
                for (i = 0; i < nx; i++) {
                    sys_mat2[idx][i*(2*nx+1)] = tmp_eig;
                    sys_mat2[idx][(nx+i)*(2*nx+1)] = tmp_eig;

                    sys_mat2[idx][i*2*nx+(nx+i)] = tmp_eig2;
                    sys_mat2[idx][(nx+i)*2*nx+i] = -tmp_eig2;
                }
                s1++;  // skip the complex conjugate eigenvalue
            }
            idx++;
        }
    } else if (opts->scheme.type == exact || (opts->scheme.type == approx && istep == 0)) {
        for (i = 0; i < num_stages*nx*num_stages*nx; i++) sys_mat[i] = 0.0;
        for (i = 0; i < num_stages*nx; i++ ) sys_mat[i*(num_stages*nx+1)] = 1.0;  // identity
    }

    for (s1 = 0; s1 < num_stages; s1++) {
        //                if (opts->scheme.type == exact || s1 == 0) {
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
        in->jac_fun(rhs_in, jac_tmp);  // k evaluation
        timing_ad += acado_toc(&timer_ad);
        //                }
        if (opts->scheme.type == simplified_in || opts->scheme.type == simplified_inis
                || opts->scheme.type == approx || opts->scheme.type == simplified_inis2) {
            for (i = 0; i < nx*nx; i++) jac_traj[istep*num_stages+s1][i] = jac_tmp[nx+i];
        }

        // put jac_tmp in sys_mat:
        if (opts->scheme.type == exact) {
            for (s2 = 0; s2 < num_stages; s2++) {
                for (j = 0; j < nx; j++) {
                    for (i = 0; i < nx; i++) {
                        sys_mat[(s2*nx+j)*num_stages*nx+s1*nx+i] -=
                                H_INT*A_mat[s2*num_stages+s1]*jac_tmp[nx+j*nx+i];
                    }
                }
            }
        }
    }

    int idx = 0;
    for (s1 = 0; s1 < num_stages; s1++) {
        // put jac_traj[0] in sys_mat:
        if ((opts->scheme.type == simplified_in || opts->scheme.type == simplified_inis
                || opts->scheme.type == simplified_inis2) && istep == 0) {
            if ((s1+1) == num_stages) {  // real eigenvalue
                for (j = 0; j < nx; j++) {
                    for (i = 0; i < nx; i++) {
                        sys_mat2[idx][j*nx+i] -= jac_traj[0][j*nx+i];
                    }
                }
            } else if (opts->scheme.type != simplified_inis2) {  // complex conjugate pair
                for (j = 0; j < nx; j++) {
                    for (i = 0; i < nx; i++) {
                        sys_mat2[idx][j*2*nx+i] -= jac_traj[0][j*nx+i];
                        sys_mat2[idx][(nx+j)*2*nx+nx+i] -= jac_traj[0][j*nx+i];
                    }
                }
                s1++;  // skip the complex conjugate eigenvalue
            } else if (opts->scheme.type == simplified_inis2) {  // complex conjugate pair
                tmp_mat = work->sys_mat;
                for (i = 0; i < nx*nx; i++) {
                    tmp_mat[i] = 0.0;
                    sys_mat2[idx+1][i] = 0.0;
                }
                tmp_eig = 1.0/H_INT*opts->scheme.eig[s1];
                tmp_eig2 = 1.0/H_INT*opts->scheme.eig[s1+1];
                tmp_eig3 = 1.0/tmp_eig2;
                for (i = 0; i < nx; i++) {
                    tmp_mat[i*(nx+1)] = tmp_eig;
                    sys_mat2[idx+1][i*(nx+1)] = tmp_eig*tmp_eig3;  // alpha/beta
                }

                for (j = 0; j < nx; j++) {
                    for (i = 0; i < nx; i++) {
                        tmp_mat[j*nx+i] -= jac_traj[0][j*nx+i];
                        sys_mat2[idx+1][j*nx+i] -= jac_traj[0][j*nx+i]*tmp_eig3;
                    }
                }

                // Combine second subsystem into the first one:
                for (j = 0; j < nx; j++) {
                    for (i = 0; i < nx; i++) {
                        sys_mat2[idx][j*nx+i] = 0.0;
                        for (k = 0; k < nx; k++) {
                            sys_mat2[idx][j*nx+i] += tmp_mat[k*nx+i]*sys_mat2[idx+1][j*nx+k];
                        }
                    }
                }
                for (i = 0; i < nx; i++) sys_mat2[idx][i*(nx+1)] += tmp_eig2;

                idx++;
                s1++;  // skip the complex conjugate eigenvalue
            }
            idx++;
        } else if (opts->scheme.type == approx && istep == 0) {
            for (s2 = 0; s2 < num_stages; s2++) {
                for (j = 0; j < nx; j++) {
                    for (i = 0; i < nx; i++) {
                        sys_mat[(s2*nx+j)*num_stages*nx+s1*nx+i] -=
                                H_INT*A_mat[s2*num_stages+s1]*jac_traj[0][j*nx+i];
                    }
                }
            }
        }
    }
}


void sim_lifted_irk(const sim_in *in, sim_out *out,
        void *mem_, void *work_ ) {
    int_t nx = in->nx;
    int_t nu = in->nu;
    sim_RK_opts *opts = in->opts;
    int_t num_stages = opts->num_stages;
    int_t dim_sys = num_stages*nx;
    int_t i, s1, s2, j, istep;
    sim_lifted_irk_memory *mem = (sim_lifted_irk_memory*) mem_;
    sim_lifted_irk_workspace *work = (sim_lifted_irk_workspace*) work_;
    real_t H_INT = in->step;
    int_t NSTEPS = in->nSteps;
    int_t NF = in->nsens_forw;

    real_t *A_mat = opts->A_mat;
    real_t *b_vec = opts->b_vec;
//    real_t *c_vec = opts->c_vec;

//    print_matrix("stdout", A_mat, num_stages, num_stages);

    real_t **VDE_tmp = work->VDE_tmp;
    real_t *out_tmp = work->out_tmp;
    real_t *rhs_in = work->rhs_in;

    real_t **jac_traj = mem->jac_traj;

    real_t *K_traj = mem->K_traj;
    real_t *DK_traj = mem->DK_traj;
    real_t *mu_traj = mem->mu_traj;
    real_t *adj_traj = mem->adj_traj;

    real_t *adj_tmp = work->out_adj_tmp;

    int *ipiv = work->ipiv;  // pivoting vector
    int **ipiv2 = mem->ipiv2;  // pivoting vector
    real_t *sys_mat = work->sys_mat;
    real_t **sys_mat2 = mem->sys_mat2;
    real_t *sys_sol = work->sys_sol;
    real_t **sys_sol2 = mem->sys_sol2;
    real_t *sys_sol_trans = work->sys_sol_trans;
#if !TRIPLE_LOOP
    struct d_strmat *str_mat = work->str_mat;
    struct d_strmat **str_mat2 = mem->str_mat2;

    struct d_strmat *str_sol = work->str_sol;
    struct d_strmat **str_sol2 = mem->str_sol2;
#endif  // !TRIPLE_LOOP

    acado_timer timer, timer_la, timer_ad;
    real_t timing_la = 0.0;
    real_t timing_ad = 0.0;

    acado_tic(&timer);

    for (i = 0; i < nx; i++) out_tmp[i] = in->x[i];
    for (i = 0; i < nx*NF; i++) out_tmp[nx+i] = in->S_forw[i];  // sensitivities

    for (i = 0; i < nu; i++) rhs_in[nx*(1+NF)+i] = in->u[i];

    if (opts->scheme.type == simplified_in || opts->scheme.type == simplified_inis
            || opts->scheme.type == approx || opts->scheme.type == simplified_inis2) {
        for (i = 0; i < nx; i++) adj_tmp[i] = in->S_adj[i];
    }

    // Newton step of the collocation variables with respect to the inputs:
    if (NF == (nx+nu) && in->sens_adj) {
        for (istep = NSTEPS-1; istep > -1; istep--) {  // ADJOINT update
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

            // Newton step of the Lagrange multipliers mu, based on adj_traj:
            if (opts->scheme.type == simplified_in || opts->scheme.type == simplified_inis
                    || opts->scheme.type == approx || opts->scheme.type == simplified_inis2) {
                for (s1 = 0; s1 < num_stages; s1++) {
                    for (i = 0; i < nx; i++) {
                        sys_sol[s1*nx+i] = -mu_traj[istep*num_stages*nx+s1*nx+i];
                        for (s2 = 0; s2 < num_stages; s2++) {
                            sys_sol[s1*nx+i] += H_INT*A_mat[s1*num_stages+s2]*
                                    adj_traj[istep*num_stages*nx+s2*nx+i];
                        }
                        sys_sol[s1*nx+i] -= H_INT*b_vec[s1]*adj_tmp[i];
                    }
                }
//                print_matrix("stdout", sys_sol, 1, num_stages*nx);
//                print_matrix("stdout", sys_sol, 1, 1);

                // TRANSFORM using transf1_T:
                if (opts->scheme.type == simplified_in || opts->scheme.type == simplified_inis
                        || opts->scheme.type == simplified_inis2) {
                    // apply the transf1 operation:
                    for (s1 = 0; s1 < num_stages; s1++) {
                        for (s2 = 0; s2 < num_stages; s2++) {
                            work->trans[s2*num_stages+s1] = 1.0/H_INT*
                                    opts->scheme.transf1_T[s2*num_stages+s1];
                        }
                    }
                    transform_mat(sys_sol, work->trans, sys_sol_trans, num_stages, nx, 1);

                    // construct sys_sol2 from sys_sol_trans:
                    if (opts->scheme.type != simplified_inis2) {
                        construct_subsystems(sys_sol_trans, sys_sol2, num_stages, nx, 1);
                    } else {
                        construct_subsystems_combined(sys_sol_trans, sys_sol2,
                                num_stages, nx, 1, sys_mat2, in, true);
                    }
                }

                acado_tic(&timer_la);
                int_t idx = 0;
                for (s1 = 0; s1 < num_stages; s1++) {
                    // THIS LOOP IS PARALLELIZABLE BECAUSE OF DECOMPOSABLE LINEAR SUBSYSTEMS
                    if (opts->scheme.type == simplified_in
                            || opts->scheme.type == simplified_inis
                            || opts->scheme.type == simplified_inis2) {
                        sys_mat = sys_mat2[idx];
                        ipiv = ipiv2[idx];
                        sys_sol = sys_sol2[idx];
#if !TRIPLE_LOOP
                        str_mat = str_mat2[idx];
                        str_sol = str_sol2[idx];
#endif
                        idx++;
                        if ((s1+1) < num_stages && opts->scheme.type != simplified_inis2) {
                            dim_sys = 2*nx;
                            s1++;  // complex conjugate pair of eigenvalues
                        } else if ((s1+1) < num_stages && opts->scheme.type == simplified_inis2) {
                            dim_sys = nx;
                            idx++;
                            s1++;  // complex conjugate pair of eigenvalues
                        } else {
                            dim_sys = nx;
                        }
                    } else {
                        dim_sys = num_stages*nx;
                        s1 = num_stages;  // break out of for-loop
                    }
#if TRIPLE_LOOP
                    solve_system_trans_ACADO(sys_mat, sys_sol, ipiv, dim_sys, 1);
#else  // TRIPLE_LOOP
#error: NOT YET IMPLEMENTED
#endif  // TRIPLE_LOOP
                }
                timing_la += acado_toc(&timer_la);

                // TRANSFORM using transf2_T:
                if (opts->scheme.type == simplified_in || opts->scheme.type == simplified_inis
                        || opts->scheme.type == simplified_inis2) {
                    // construct sys_sol_trans from sys_sol2:
                    if (opts->scheme.type != simplified_inis2) {
                        destruct_subsystems(sys_sol_trans, sys_sol2, num_stages, nx, 1);
                    } else {
                        destruct_subsystems_combined(sys_sol_trans, sys_sol2,
                                num_stages, nx, 1, sys_mat2, true);
                    }

                    // apply the transf2 operation:
                    sys_sol = work->sys_sol;
                    transform_mat(sys_sol_trans, opts->scheme.transf2_T, sys_sol,
                            num_stages, nx, 1);
                }

                // update mu_traj
                for (s1 = 0; s1 < num_stages; s1++) {
                    for (i = 0; i < nx; i++) {
                        mu_traj[istep*num_stages*nx+s1*nx+i] += sys_sol[s1*nx+i];
                    }
                }

                // update adj_tmp:
                // TODO(rien): USE ADJOINT DIFFERENTIATION HERE INSTEAD !!:
                for (j = 0; j < nx; j++) {
                    for (s1 = 0; s1 < num_stages; s1++) {
                        for (i = 0; i < nx; i++) {
                            adj_tmp[j] += mu_traj[istep*num_stages*nx+s1*nx+i]*
                                    jac_traj[istep*num_stages+s1][j*nx+i];
                        }
                    }
                }
            }
        }
    }

    if (opts->scheme.type == simplified_in || opts->scheme.type == simplified_inis
            || opts->scheme.type == approx || opts->scheme.type == simplified_inis2) {
        for (i = 0; i < NF; i++) out->grad[i] = 0.0;
    }

    for (istep = 0; istep < NSTEPS; istep++) {
        // form linear system matrix (explicit ODE case):
        form_linear_system_matrix(istep, in, mem, work, sys_mat, sys_mat2, timing_ad);

        int_t idx;
        if (opts->scheme.type == exact || istep == 0) {
            acado_tic(&timer_la);
            idx = 0;
            for (s1 = 0; s1 < num_stages; s1++) {
                // THIS LOOP IS PARALLELIZABLE BECAUSE OF DECOMPOSABLE LINEAR SUBSYSTEMS
                if (opts->scheme.type == simplified_in || opts->scheme.type == simplified_inis
                        || opts->scheme.type == simplified_inis2) {
                    sys_mat = sys_mat2[idx];
                    ipiv = ipiv2[idx];
#if !TRIPLE_LOOP
                    str_mat = str_mat2[idx];
#endif
                    idx++;
                    if ((s1+1) < num_stages && opts->scheme.type != simplified_inis2) {
                        dim_sys = 2*nx;
                        s1++;  // complex conjugate pair of eigenvalues
                    } else if ((s1+1) < num_stages && opts->scheme.type == simplified_inis2) {
                        dim_sys = nx;
                        idx++;
                        s1++;  // complex conjugate pair of eigenvalues
                    } else {
                        dim_sys = nx;
                    }
                } else {
                    dim_sys = num_stages*nx;
                    s1 = num_stages;  // break out of for-loop
                }
#if TRIPLE_LOOP
                LU_system_ACADO(sys_mat, ipiv, dim_sys);
#else   // TRIPLE_LOOP
                // ---- BLASFEO: LU factorization ----
#if defined(LA_HIGH_PERFORMANCE)
                d_cvt_mat2strmat(dim_sys, dim_sys, sys_mat, dim_sys,
                        str_mat, 0, 0);  // mat2strmat
#endif  // LA_BLAS | LA_REFERENCE
                dgetrf_libstr(dim_sys, dim_sys, str_mat, 0, 0, str_mat, 0,
                        0, ipiv);  // Gauss elimination
                // ---- BLASFEO: LU factorization ----
#endif   // TRIPLE_LOOP
            }
            timing_la += acado_toc(&timer_la);
        }

        if (opts->scheme.type == simplified_in || opts->scheme.type == simplified_inis
                || opts->scheme.type == simplified_inis2) {
            sys_sol = work->sys_sol;
            sys_sol_trans = work->sys_sol_trans;
        }
        for (s1 = 0; s1 < num_stages; s1++) {
            for (i = 0; i < nx*(1+NF); i++) {
                rhs_in[i] = out_tmp[i];
            }
            for (s2 = 0; s2 < num_stages; s2++) {
                for (i = 0; i < nx; i++) {
                    rhs_in[i] += H_INT*A_mat[s2*num_stages+s1]*K_traj[istep*num_stages*nx+s2*nx+i];
                }
                if (opts->scheme.type == simplified_inis || opts->scheme.type == simplified_inis2) {
                    for (j = 0; j < NF; j++) {
                        for (i = 0; i < nx; i++) {
                            rhs_in[(j+1)*nx+i] += H_INT*A_mat[s2*num_stages+s1]*
                                    DK_traj[(istep*num_stages+s2)*nx*NF+j*nx+i];
                        }
                    }
                }
            }
            acado_tic(&timer_ad);
            in->VDE_forw(rhs_in, VDE_tmp[s1]);  // k evaluation
            timing_ad += acado_toc(&timer_ad);

            // put VDE_tmp in sys_sol:
            for (j = 0; j < 1+NF; j++) {
                for (i = 0; i < nx; i++) {
                    sys_sol[j*num_stages*nx+s1*nx+i] = VDE_tmp[s1][j*nx+i];
                }
            }
            for (i = 0; i < nx; i++) {
                sys_sol[s1*nx+i] -= K_traj[istep*num_stages*nx+s1*nx+i];
            }
            if (opts->scheme.type == simplified_inis || opts->scheme.type == simplified_inis2) {
                for (j = 0; j < NF; j++) {
                    for (i = 0; i < nx; i++) {
                        sys_sol[(j+1)*num_stages*nx+s1*nx+i] -=
                                DK_traj[(istep*num_stages+s1)*nx*NF+j*nx+i];
                    }
                }
            }
        }

        // Inexact Newton with Iterated Sensitivities (INIS) based gradient correction:
        if (opts->scheme.type == simplified_inis || opts->scheme.type == simplified_inis2) {
            for (j = 0; j < NF; j++) {
                for (i = 0; i < num_stages*nx; i++) {
                    out->grad[j] += mu_traj[istep*num_stages*nx+i]*sys_sol[(j+1)*num_stages*nx+i];
                }
            }
        }

        if (opts->scheme.type == simplified_in || opts->scheme.type == simplified_inis
                || opts->scheme.type == simplified_inis2) {
            // apply the transf1 operation:
            for (s1 = 0; s1 < num_stages; s1++) {
                for (s2 = 0; s2 < num_stages; s2++) {
                    work->trans[s2*num_stages+s1] = 1.0/H_INT*
                            opts->scheme.transf1[s2*num_stages+s1];
                }
            }
            transform_mat(sys_sol, work->trans, sys_sol_trans, num_stages, nx, 1+NF);

            // construct sys_sol2 from sys_sol_trans:
            if (opts->scheme.type != simplified_inis2) {
                construct_subsystems(sys_sol_trans, sys_sol2, num_stages, nx, 1+NF);
            } else {
                construct_subsystems_combined(sys_sol_trans, sys_sol2,
                        num_stages, nx, 1+NF, sys_mat2, in, false);
            }
        }

        acado_tic(&timer_la);
        idx = 0;
        for (s1 = 0; s1 < num_stages; s1++) {
            // THIS LOOP IS PARALLELIZABLE BECAUSE OF DECOMPOSABLE LINEAR SUBSYSTEMS
            if (opts->scheme.type == simplified_in || opts->scheme.type == simplified_inis
                    || opts->scheme.type == simplified_inis2) {
                sys_mat = sys_mat2[idx];
                ipiv = ipiv2[idx];
                sys_sol = sys_sol2[idx];
#if !TRIPLE_LOOP
                str_mat = str_mat2[idx];
                str_sol = str_sol2[idx];
#endif
                idx++;
                if ((s1+1) < num_stages && opts->scheme.type != simplified_inis2) {
                    dim_sys = 2*nx;
                    s1++;  // complex conjugate pair of eigenvalues
                } else if ((s1+1) < num_stages && opts->scheme.type == simplified_inis2) {
                    dim_sys = nx;
                    idx++;
                    s1++;  // complex conjugate pair of eigenvalues
                } else {
                    dim_sys = nx;
                }
            } else {
                dim_sys = num_stages*nx;
                s1 = num_stages;  // break out of for-loop
            }
#if TRIPLE_LOOP
            solve_system_ACADO(sys_mat, sys_sol, ipiv, dim_sys, 1+NF);
#else  // TRIPLE_LOOP
#if defined(LA_HIGH_PERFORMANCE)
            // ---- BLASFEO: row transformations + backsolve ----
            d_cvt_mat2strmat(dim_sys, 1+NF, sys_sol, dim_sys,
                    str_sol, 0, 0);  // mat2strmat
            drowpe_libstr(dim_sys, ipiv, str_sol);  // row permutations
            dtrsm_llnu_libstr(dim_sys, 1+NF, 1.0, str_mat, 0, 0, str_sol,
                    0, 0, str_sol, 0, 0);  // L backsolve
            dtrsm_lunn_libstr(dim_sys, 1+NF, 1.0, str_mat, 0, 0, str_sol,
                    0, 0, str_sol, 0, 0);  // U backsolve
            d_cvt_strmat2mat(dim_sys, 1+NF, str_sol, 0, 0, sys_sol,
                    dim_sys);  // strmat2mat
            // ---- BLASFEO: row transformations + backsolve ----
#else   // LA_BLAS | LA_REFERENCE
            // ---- BLASFEO: row transformations + backsolve ----
            drowpe_libstr(dim_sys, ipiv, str_sol);  // row permutations
            dtrsm_llnu_libstr(dim_sys, 1+NF, 1.0, str_mat, 0, 0, str_sol,
                    0, 0, str_sol, 0, 0);  // L backsolve
            dtrsm_lunn_libstr(dim_sys, 1+NF, 1.0, str_mat, 0, 0, str_sol,
                    0, 0, str_sol, 0, 0);  // U backsolve
            // ---- BLASFEO: row transformations + backsolve ----
#endif  // LA_BLAFEO
#endif  // TRIPLE_LOOP
        }
        timing_la += acado_toc(&timer_la);

        if (opts->scheme.type == simplified_in || opts->scheme.type == simplified_inis
                || opts->scheme.type == simplified_inis2) {
            // construct sys_sol_trans from sys_sol2:
            if (opts->scheme.type != simplified_inis2) {
                destruct_subsystems(sys_sol_trans, sys_sol2, num_stages, nx, 1+NF);
            } else {
                destruct_subsystems_combined(sys_sol_trans, sys_sol2,
                        num_stages, nx, 1+NF, sys_mat2, false);
            }

            // apply the transf2 operation:
            sys_sol = work->sys_sol;
            transform_mat(sys_sol_trans, opts->scheme.transf2, sys_sol, num_stages, nx, 1+NF);
        }

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
            for (j = 0; j < NF; j++) {
                for (i = 0; i < nx; i++) {
                    if (opts->scheme.type == simplified_inis
                            || opts->scheme.type == simplified_inis2) {
                        DK_traj[(istep*num_stages+s1)*nx*NF+j*nx+i] +=
                            sys_sol[(j+1)*num_stages*nx+s1*nx+i];
                    } else {
                        DK_traj[(istep*num_stages+s1)*nx*NF+j*nx+i] =
                            sys_sol[(j+1)*num_stages*nx+s1*nx+i];
                    }
                }
            }
        }
        for (s1 = 0; s1 < num_stages; s1++) {
            for (i = 0; i < nx*NF; i++) {
                out_tmp[nx+i] += H_INT*b_vec[s1]*
                        DK_traj[(istep*num_stages+s1)*nx*NF+i];  // RK step
            }
        }
        if (opts->scheme.type == simplified_inis || opts->scheme.type == simplified_inis2
                || opts->scheme.type == simplified_in
                || opts->scheme.type == approx) {  // Adjoint derivatives:
            for (s1 = 0; s1 < num_stages; s1++) {
                for (j = 0; j < nx; j++) {
                    adj_traj[istep*num_stages*nx+s1*nx+j] = 0.0;
                    for (i = 0; i < nx; i++) {
                        adj_traj[istep*num_stages*nx+s1*nx+j] +=
                                mu_traj[istep*num_stages*nx+s1*nx+i]*
                                jac_traj[istep*num_stages+s1][j*nx+i];
                    }
                }
            }
        }
//        print_matrix_name("stdout", "adj_traj", &adj_traj[0], 1, num_stages*nx);
//        print_matrix_name("stdout", "mu_traj", &mu_traj[0], 1, num_stages*nx);
//        print_matrix_name("stdout", "DK_traj", &DK_traj[0], 1, nx*NF);
//        print_matrix_name("stdout", "VDE_tmp[0]", &VDE_tmp[0][0], 1, nx*(NF+1));
//        print_matrix_name("stdout", "VDE_tmp[1]", &VDE_tmp[1][0], 1, nx*(NF+1));
        if (opts->scheme.type == simplified_in || opts->scheme.type == approx) {
            // Standard Inexact Newton based gradient correction:
            for (j = 0; j < NF; j++) {
                for (s1 = 0; s1 < num_stages; s1++) {
                    for (i = 0; i < nx; i++) {
                        out->grad[j] += mu_traj[istep*num_stages*nx+s1*nx+i]*
                                VDE_tmp[s1][(j+1)*nx+i];
                    }
                }
            }
            for (j = 0; j < NF; j++) {
                for (s1 = 0; s1 < num_stages; s1++) {
                    for (i = 0; i < nx; i++) {
                        out->grad[j] -= mu_traj[istep*num_stages*nx+s1*nx+i]*
                                DK_traj[(istep*num_stages+s1)*nx*NF+j*nx+i];
                    }
                }
                for (s2 = 0; s2 < num_stages; s2++) {
                    for (s1 = 0; s1 < num_stages; s1++) {
                        for (i = 0; i < nx; i++) {
                            out->grad[j] += H_INT*A_mat[s2*num_stages+s1]*
                                    adj_traj[istep*num_stages*nx+s1*nx+i]*
                                    DK_traj[(istep*num_stages+s2)*nx*NF+j*nx+i];
                        }
                    }
                }
            }
        }
//        print_matrix_name("stdout", "grad", &out->grad[0], 1, NF);
    }
    for (i = 0; i < nx; i++)    out->xn[i] = out_tmp[i];
    for (i = 0; i < nx*NF; i++) out->S_forw[i] = out_tmp[nx+i];

    out->info->CPUtime = acado_toc(&timer);
    out->info->LAtime = timing_la;
    out->info->ADtime = timing_ad;
}


void sim_lifted_irk_create_workspace(const sim_in *in,
        sim_lifted_irk_workspace *work) {
    int_t nx = in->nx;
    int_t nu = in->nu;
    sim_RK_opts *opts = in->opts;
    int_t num_stages = opts->num_stages;
    int_t NF = in->nsens_forw;
//    int_t num_sys = ceil(num_stages/2.0);
    int_t dim_sys = num_stages*nx;
    if (opts->scheme.type == simplified_in || opts->scheme.type == simplified_inis) {
        dim_sys = nx;
        if (num_stages > 1) dim_sys = 2*nx;
    } else if (opts->scheme.type == simplified_inis2) {
        dim_sys = nx;
    }

    work->rhs_in = malloc(sizeof(*work->rhs_in) * (nx*(1+NF)+nu));
    work->out_tmp = malloc(sizeof(*work->out_tmp) * (nx*(1+NF)));
    if (opts->scheme.type == approx || opts->scheme.type == exact) {
        work->ipiv = malloc(sizeof(*work->ipiv) * (dim_sys));
        for (int_t i = 0; i < dim_sys; i++) work->ipiv[i] = i;
        work->sys_mat = malloc(sizeof(*work->sys_mat) * (dim_sys*dim_sys));
    }
    work->sys_sol = malloc(sizeof(*work->sys_sol) * (num_stages*nx)*(1+NF));
    work->VDE_tmp = malloc(sizeof(*work->VDE_tmp) * num_stages);
    for (int_t i = 0; i < num_stages; i++) {
        work->VDE_tmp[i] = malloc(sizeof(*work->VDE_tmp[i]) * (nx*(1+NF)));
    }
    work->jac_tmp = malloc(sizeof(*work->jac_tmp) * nx*(nx+1));

    if (opts->scheme.type == simplified_in || opts->scheme.type == simplified_inis
            || opts->scheme.type == simplified_inis2) {
        work->sys_sol_trans = malloc(sizeof(*work->sys_sol_trans) * (num_stages*nx)*(1+NF));
        work->trans = malloc(sizeof(*work->trans) * (num_stages*num_stages));
    }

    if (opts->scheme.type == simplified_in || opts->scheme.type == simplified_inis
            || opts->scheme.type == simplified_inis2
            || opts->scheme.type == approx) {
        work->out_adj_tmp = malloc(sizeof(*work->out_adj_tmp) * (nx));
    }

#if !TRIPLE_LOOP
    if (opts->scheme.type == simplified_in || opts->scheme.type == simplified_inis
            || opts->scheme.type == simplified_inis2) {
        work->str_mat2 = malloc(sizeof(*work->str_mat2) * num_sys);
        work->str_sol2 = malloc(sizeof(*work->str_sol2) * num_sys);
    }
#if !defined(LA_HIGH_PERFORMANCE)
    real_t *sys_mat, *sys_sol;
#endif
    for (int_t i = 0; i < num_sys; i++) {
        if (opts->scheme.type == simplified_in || opts->scheme.type == simplified_inis
                || opts->scheme.type == simplified_inis2) {
#if !defined(LA_HIGH_PERFORMANCE)
            sys_mat = work->sys_mat2[i];
            sys_sol = work->sys_sol2[i];
#endif
            if ((i+1) == num_sys && num_sys != floor(num_stages/2.0)) {  // odd number of stages
                dim_sys = nx;
            } else {
                dim_sys = 2*nx;
            }
        } else {
#if !defined(LA_HIGH_PERFORMANCE)
            sys_mat = work->sys_mat;
            sys_sol = work->sys_sol;
#endif
            dim_sys = num_stages*nx;
            i = num_sys;  // break out of for-loop
        }

#if defined(LA_HIGH_PERFORMANCE)
        // matrices in matrix struct format:
        int size_strmat = 0;
        size_strmat += d_size_strmat(dim_sys, dim_sys);
        size_strmat += d_size_strmat(dim_sys, 1+NF);

        // accocate memory
        void *memory_strmat;
        v_zeros_align(&memory_strmat, size_strmat);
        char *ptr_memory_strmat = (char *) memory_strmat;

        d_create_strmat(dim_sys, dim_sys, work->str_mat, ptr_memory_strmat);
        ptr_memory_strmat += work->str_mat->memory_size;
        d_create_strmat(dim_sys, 1+NF, work->str_sol, ptr_memory_strmat);
        ptr_memory_strmat += work->str_sol->memory_size;

#elif defined(LA_REFERENCE)

        //  pointer to column-major matrix
        d_create_strmat(dim_sys, dim_sys, work->str_mat, sys_mat);
        d_create_strmat(dim_sys, 1+NF, work->str_sol, sys_sol);

        // allocate new memory only for the diagonal
        int size_strmat = 0;
        size_strmat += d_size_diag_strmat(dim_sys, dim_sys);

        void *memory_strmat = malloc(size_strmat);
        //    void *memory_strmat;
        //    v_zeros_align(&memory_strmat, size_strmat);
        char *ptr_memory_strmat = (char *) memory_strmat;

        d_cast_diag_mat2strmat((double *) ptr_memory_strmat, work->str_mat);
        //    ptr_memory_strmat += d_size_diag_strmat(dim_sys, dim_sys);

#else  // LA_BLAS

        // not allocate new memory: point to column-major matrix
        d_create_strmat(dim_sys, dim_sys, work->str_mat, sys_mat);
        d_create_strmat(dim_sys, 1+NF, work->str_sol, sys_sol);

#endif  // LA_HIGH_PERFORMANCE

        if (opts->scheme.type == simplified_in || opts->scheme.type == simplified_inis
                || opts->scheme.type == simplified_inis2) {
            work->str_mat2[i] = work->str_mat;
            work->str_sol2[i] = work->str_sol;
        }
    }
#endif  // !TRIPLE_LOOP
}


void sim_lifted_irk_create_memory(const sim_in *in,
        sim_lifted_irk_memory *mem) {
    int_t i;
    int_t nx = in->nx;
    int_t nu = in->nu;
    int_t nSteps = in->nSteps;
    sim_RK_opts *opts = in->opts;
    int_t num_stages = opts->num_stages;
    int_t NF = in->nsens_forw;
    int_t num_sys = ceil(num_stages/2.0);
    if (opts->scheme.type == simplified_inis2) num_sys = num_stages;
//    printf("num_stages: %d \n", num_stages);
//    printf("ceil(num_stages/2.0): %d \n", (int)ceil(num_stages/2.0));
//    printf("floor(num_stages/2.0): %d \n", (int)floor(num_stages/2.0));

    mem->K_traj = malloc(sizeof(*mem->K_traj) * (nSteps*num_stages*nx));
    mem->DK_traj = malloc(sizeof(*mem->DK_traj) * (nSteps*num_stages*nx*NF));
    mem->mu_traj = malloc(sizeof(*mem->mu_traj) * (nSteps*num_stages*nx));
    mem->x = malloc(sizeof(*mem->x) * nx);
    mem->u = malloc(sizeof(*mem->u) * nu);

    if (opts->scheme.type == simplified_in || opts->scheme.type == simplified_inis
            || opts->scheme.type == simplified_inis2
            || opts->scheme.type == approx) {
        mem->adj_traj = malloc(sizeof(*mem->adj_traj) * (nSteps*num_stages*nx));
        for ( i = 0; i < nSteps*num_stages*nx; i++ ) mem->adj_traj[i] = 0.0;

        mem->jac_traj = malloc(sizeof(*mem->jac_traj) * nSteps*num_stages);
        for (int_t i = 0; i < nSteps*num_stages; i++) {
            mem->jac_traj[i] = malloc(sizeof(*mem->jac_traj[i]) * (nx*nx));
        }
    }

    for ( i = 0; i < nSteps*num_stages*nx; i++ ) mem->K_traj[i] = 0.0;
    for ( i = 0; i < nSteps*num_stages*nx*NF; i++ ) mem->DK_traj[i] = 0.0;
    for ( i = 0; i < nSteps*num_stages*nx; i++ ) mem->mu_traj[i] = 0.0;

    if (opts->scheme.type == simplified_in || opts->scheme.type == simplified_inis
            || opts->scheme.type == simplified_inis2) {
        mem->sys_mat2 = malloc(sizeof(*mem->sys_mat2) * num_sys);
        mem->ipiv2 = malloc(sizeof(*mem->ipiv2) * num_sys);
        mem->sys_sol2 = malloc(sizeof(*mem->sys_sol2) * num_sys);
        for (int_t i = 0; i < num_sys; i++) {
            if ((i+1) == num_sys && num_sys != floor(num_stages/2.0)
                    && opts->scheme.type != simplified_inis2) {  // odd number of stages
                mem->sys_mat2[i] = malloc(sizeof(*(mem->sys_mat2[i])) * (nx*nx));
                mem->ipiv2[i] = malloc(sizeof(*(mem->ipiv2[i])) * (nx));
                mem->sys_sol2[i] = malloc(sizeof(*(mem->sys_sol2[i])) * (nx*(1+NF)));

                for (int_t j = 0; j < nx*nx; j++) mem->sys_mat2[i][j] = 0.0;
                for (int_t j = 0; j < nx; j++) mem->sys_mat2[i][j*(nx+1)] = 1.0;
            } else if (opts->scheme.type != simplified_inis2) {
                mem->sys_mat2[i] = malloc(sizeof(*(mem->sys_mat2[i])) * (4*nx*nx));
                mem->ipiv2[i] = malloc(sizeof(*(mem->ipiv2[i])) * (2*nx));
                mem->sys_sol2[i] = malloc(sizeof(*(mem->sys_sol2[i])) * (2*nx*(1+NF)));

                for (int_t j = 0; j < 4*nx*nx; j++) mem->sys_mat2[i][j] = 0.0;
                for (int_t j = 0; j < 2*nx; j++) mem->sys_mat2[i][j*(2*nx+1)] = 1.0;
            } else if (opts->scheme.type == simplified_inis2) {
                mem->sys_mat2[i] = malloc(sizeof(*(mem->sys_mat2[i])) * (nx*nx));
                mem->ipiv2[i] = malloc(sizeof(*(mem->ipiv2[i])) * (nx));
                mem->sys_sol2[i] = malloc(sizeof(*(mem->sys_sol2[i])) * (nx*(1+NF)));

                for (int_t j = 0; j < nx*nx; j++) mem->sys_mat2[i][j] = 0.0;
                for (int_t j = 0; j < nx; j++) mem->sys_mat2[i][j*(nx+1)] = 1.0;
            }
        }
    }
}


void sim_irk_create_opts(const int_t num_stages, const char* name, sim_RK_opts *opts) {
    opts->num_stages = num_stages;
    opts->A_mat = malloc(sizeof(*opts->A_mat) * (num_stages*num_stages));
    opts->b_vec = malloc(sizeof(*opts->b_vec) * (num_stages));
    opts->c_vec = malloc(sizeof(*opts->c_vec) * (num_stages));
    opts->scheme.type = exact;

    if ( strcmp(name, "Gauss") == 0 ) {  // GAUSS METHODS
        get_Gauss_nodes(opts->num_stages, opts->c_vec);
        create_Butcher_table(opts->num_stages, opts->c_vec, opts->b_vec, opts->A_mat);
    } else if ( strcmp(name, "Radau") == 0 ) {
        // TODO(rien): add Radau IIA collocation schemes
//        get_Radau_nodes(opts->num_stages, opts->c_vec);
        create_Butcher_table(opts->num_stages, opts->c_vec, opts->b_vec, opts->A_mat);
    } else {
        // throw error somehow?
    }
//    print_matrix("stdout", opts->c_vec, num_stages, 1);
//    print_matrix("stdout", opts->A_mat, num_stages, num_stages);
//    print_matrix("stdout", opts->b_vec, num_stages, 1);
}


void sim_irk_create_Newton_scheme(const int_t num_stages, const char* name,
        sim_RK_opts *opts, enum Newton_type_collocation type) {
    opts->scheme.type = type;
    if ( strcmp(name, "Gauss") == 0 ) {  // GAUSS METHODS
        if (num_stages <= 15 && (type == simplified_in || type == simplified_inis
                || type == simplified_inis2)) {
            opts->scheme.eig = malloc(sizeof(*opts->scheme.eig) * (num_stages));
            opts->scheme.transf1 = malloc(sizeof(*opts->scheme.transf1) * (num_stages*num_stages));
            opts->scheme.transf2 = malloc(sizeof(*opts->scheme.transf2) * (num_stages*num_stages));
            opts->scheme.transf1_T =
                    malloc(sizeof(*opts->scheme.transf1_T) * (num_stages*num_stages));
            opts->scheme.transf2_T =
                    malloc(sizeof(*opts->scheme.transf2_T) * (num_stages*num_stages));
            read_Gauss_simplified(opts->num_stages, &opts->scheme);
        }
    } else if ( strcmp(name, "Radau") == 0 ) {
        // TODO(rien): add Radau IIA collocation schemes
    } else {
        // throw error somehow?
    }
}
