/*
 *    This file is part of acados.
 *
 *    acados is free software; you can redistribute it and/or
 *    modify it under the terms of the GNU Lesser General Public
 *    License as published by the Free Software Foundation; either
 *    version 3 of the License, or (at your option) any later version.
 *
 *    acados is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *    Lesser General Public License for more details.
 *
 *    You should have received a copy of the GNU Lesser General Public
 *    License along with acados; if not, write to the Free Software Foundation,
 *    Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 *
 */

#include "acados/sim/sim_lifted_irk_integrator.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "blasfeo/include/blasfeo_target.h"
#include "blasfeo/include/blasfeo_common.h"
#include "blasfeo/include/blasfeo_d_aux_ext_dep.h"
#include "blasfeo/include/blasfeo_d_blas.h"
#include "blasfeo/include/blasfeo_d_kernel.h"
#include "blasfeo/include/blasfeo_i_aux_ext_dep.h"

#include "acados/utils/mem.h"
#include "acados/utils/print.h"


int sim_lifted_irk_opts_calculate_size(sim_dims *dims)
{

    int size = sizeof(sim_rk_opts);

    int ns = dims->num_stages;
    size += ns * ns * sizeof(double);  // A_mat
    size += ns * sizeof(double);  // b_vec
    size += ns * sizeof(double);  // c_vec

    size += sizeof(Newton_scheme);

    make_int_multiple_of(8, &size);

    size += ns * sizeof(double);  // eig

    size += ns*ns * sizeof(double);  // transf1
    size += ns*ns * sizeof(double);  // transf2
    size += ns*ns * sizeof(double);  // transf1_T
    size += ns*ns * sizeof(double);  // transf2_T

    make_int_multiple_of(8, &size);
    size += 2 * 8;

    return size;
}



void *sim_lifted_irk_assign_opts(sim_dims *dims, void *raw_memory)
{
    char *c_ptr = (char *) raw_memory;

    sim_rk_opts *opts = (sim_rk_opts *) c_ptr;
    c_ptr += sizeof(sim_rk_opts);

    int ns = dims->num_stages;
    opts->num_stages = ns;

    align_char_to(8, &c_ptr);

    assign_double(ns*ns, &opts->A_mat, &c_ptr);
    assign_double(ns, &opts->b_vec, &c_ptr);
    assign_double(ns, &opts->c_vec, &c_ptr);

    opts->scheme = (Newton_scheme *) c_ptr;
    c_ptr += sizeof(Newton_scheme);

    align_char_to(8, &c_ptr);

    assign_double(ns, &opts->scheme->eig, &c_ptr);

    assign_double(ns*ns, &opts->scheme->transf1, &c_ptr);
    assign_double(ns*ns, &opts->scheme->transf2, &c_ptr);
    assign_double(ns*ns, &opts->scheme->transf1_T, &c_ptr);
    assign_double(ns*ns, &opts->scheme->transf2_T, &c_ptr);

    assert((char*)raw_memory + sim_lifted_irk_opts_calculate_size(dims) >= c_ptr);

    return (void *)opts;
}


void sim_lifted_irk_initialize_default_args(sim_dims *dims, void *opts_) {
    sim_rk_opts *opts = (sim_rk_opts*) opts_;
    enum Newton_type_collocation type = exact;
    opts->scheme->type = type;
    opts->scheme->freeze = false;
    get_Gauss_nodes(opts->num_stages, opts->c_vec);
    create_Butcher_table(opts->num_stages, opts->c_vec, opts->b_vec, opts->A_mat);
    if (dims->num_stages <= 15 && (type == simplified_in || type == simplified_inis)) {
        read_Gauss_simplified(opts->num_stages, opts->scheme);
    } else {
        opts->scheme->type = exact;
    }

    opts->num_steps = 2;
    opts->num_forw_sens = dims->nx + dims->nu;
    opts->sens_forw = true;
    opts->sens_adj = false;
    opts->sens_hess = false;
}


int sim_lifted_irk_calculate_memory_size(sim_dims *dims, void *opts_) {
    sim_rk_opts *opts = (sim_rk_opts *) opts_;

    int nx = dims->nx;
    int nu = dims->nu;
    int num_steps = opts->num_steps;
    int num_stages = opts->num_stages;
    int nf = opts->num_forw_sens;
    int num_sys = (int) ceil(num_stages/2.0);

    int size = sizeof(sim_lifted_irk_memory);

    size += num_steps * num_stages * sizeof(double *);  // jac_traj
    size += num_sys * sizeof(double *);  // sys_mat2
    size += num_sys * sizeof(int *);  // ipiv2
    size += num_sys * sizeof(double *);  // sys_sol2
    size += num_sys * sizeof(struct d_strmat *);  // str_mat2
    size += num_sys * sizeof(struct d_strmat *);  // str_sol2

    make_int_multiple_of(8, &size);

    size += nf * sizeof(double);  // grad_correction
    size += nx * num_stages * sizeof(double);  // grad_K

    size += num_steps * num_stages * nx * sizeof(double);  // K_traj
    size += num_steps * num_stages * nx * nf * sizeof(double);  // DK_traj
    size += num_steps * num_stages * nx * sizeof(double);  // mu_traj
    size += nx * sizeof(double);  // x
    size += nu * sizeof(double);  // u
    if (opts->scheme->type == simplified_inis)
        size += num_steps * num_stages * nx * nf * sizeof(double);  // delta_DK_traj

    if (opts->scheme->type == simplified_in || opts->scheme->type == simplified_inis) {
        size += num_steps * num_stages * nx * sizeof(double);  // adj_traj
        size += num_steps * num_stages * nx * nx * sizeof(double);  // jac_traj
    }

    if (opts->scheme->type == simplified_in || opts->scheme->type == simplified_inis) {
        int dim_sys = 2 * nx;
        size += (num_sys - 1) * dim_sys * dim_sys * sizeof(double);  // sys_mat2
        size += (num_sys - 1) * dim_sys * sizeof(double);  // ipiv2
        size += (num_sys - 1) * dim_sys * (1 + nf) * sizeof(double);  // sys_sol2

        if (num_sys != floor(num_stages / 2.0)) // odd number of stages
            dim_sys = nx;
        size += dim_sys * dim_sys * sizeof(double);  // sys_mat2
        size += dim_sys * sizeof(int);  // ipiv2
        size += dim_sys * (1 + nf) * sizeof(double);  // sys_sol2
    } else {
        num_sys = 1;
    }

#if !TRIPLE_LOOP
    if (opts->scheme->type == simplified_in || opts->scheme->type == simplified_inis) {

        int dim_sys = 2 * nx;
        for (int i = 0; i < num_sys; i++) {
            if ((i + 1) == num_sys && num_sys != floor(num_stages / 2.0))  // odd number of stages
                dim_sys = nx;

#if defined(LA_HIGH_PERFORMANCE)
            size += d_size_strmat(dim_sys, dim_sys);  // str_mat2
            size += d_size_strmat(dim_sys, 1 + NF);  // str_sol2
#elif defined(LA_REFERENCE)
            size += d_size_diag_strmat(dim_sys, dim_sys);
#else  // LA_BLAS
            size += 0;
#endif  // LA_HIGH_PERFORMANCE
        }
    }
#endif  // !TRIPLE_LOOP
    size += 2 * 8;
    return size;
}

void *sim_lifted_irk_assign_memory(sim_dims *dims, void *opts_, void *raw_memory) {
    sim_rk_opts *opts = (sim_rk_opts *) opts_;

    int nx = dims->nx;
    int nu = dims->nu;
    int num_steps = opts->num_steps;
    int num_stages = opts->num_stages;
    int nf = opts->num_forw_sens;
    int num_sys = (int) ceil(num_stages/2.0);

    char *c_ptr = raw_memory;

    sim_lifted_irk_memory *memory = raw_memory;
    c_ptr += sizeof(sim_lifted_irk_memory);

    align_char_to(8, &c_ptr);

    assign_double_ptrs(num_steps * num_stages, &memory->jac_traj, &c_ptr);
    assign_double_ptrs(num_sys, &memory->sys_mat2, &c_ptr);
    assign_int_ptrs(num_sys, &memory->ipiv2, &c_ptr);
    assign_double_ptrs(num_sys, &memory->sys_sol2, &c_ptr);
    assign_strmat_ptrs_to_ptrs(num_sys, &memory->str_mat2, &c_ptr);
    assign_strmat_ptrs_to_ptrs(num_sys, &memory->str_sol2, &c_ptr);

    align_char_to(8, &c_ptr);

    assign_double(nf, &memory->grad_correction, &c_ptr);
    assign_double(nx * num_stages, &memory->grad_K, &c_ptr);

    assign_double(num_steps * num_stages * nx, &memory->K_traj, &c_ptr);
    assign_double(num_steps * num_stages * nx * nf, &memory->DK_traj, &c_ptr);
    assign_double(num_steps * num_stages * nx, &memory->mu_traj, &c_ptr);
    assign_double(nx, &memory->x, &c_ptr);
    assign_double(nu, &memory->u, &c_ptr);

    if (opts->scheme->type == simplified_inis)
        assign_double(num_steps * num_stages * nx * nf, &memory->delta_DK_traj, &c_ptr);

    if (opts->scheme->type == simplified_in || opts->scheme->type == simplified_inis) {
        assign_double(num_steps * num_stages * nx, &memory->adj_traj, &c_ptr);
        for (int i = 0; i < num_steps * num_stages; ++i)
            assign_double(nx * nx, &memory->jac_traj[i], &c_ptr);
    }

    if (opts->scheme->type == simplified_in || opts->scheme->type == simplified_inis) {
        int dim_sys = 2 * nx;
        for (int i = 0; i < num_sys; ++i) {
            assign_double(dim_sys * dim_sys, &memory->sys_mat2[i], &c_ptr);
            assign_double(dim_sys, &memory->sys_mat2[i], &c_ptr);
            assign_double(dim_sys * (1 + nf), &memory->sys_mat2[i], &c_ptr);

            for (int j = 0; j < dim_sys * dim_sys; ++j)
                memory->sys_mat2[i][j] = 0.0;
            for (int j = 0; j < dim_sys; j++)
                memory->sys_mat2[i][j * (dim_sys + 1)] = 1.0;

        }
        if (num_sys != floor(num_stages / 2.0)) // odd number of stages
            dim_sys = nx;
        assign_double(dim_sys * dim_sys, &memory->sys_mat2[num_sys], &c_ptr);
        assign_double(dim_sys, &memory->sys_mat2[num_sys], &c_ptr);
        assign_double(dim_sys * (1 + nf), &memory->sys_mat2[num_sys], &c_ptr);

        for (int j = 0; j < dim_sys * dim_sys; ++j)
            memory->sys_mat2[num_sys][j] = 0.0;
        for (int j = 0; j < dim_sys; j++)
            memory->sys_mat2[num_sys][j * (dim_sys + 1)] = 1.0;
    } else {
        num_sys = 1;
    }

#if !TRIPLE_LOOP
    if (opts->scheme->type == simplified_in || opts->scheme->type == simplified_inis) {

    int dim_sys = 2 * nx;
    for (int i = 0; i < num_sys; i++) {
        if ((i + 1) == num_sys && num_sys != floor(num_stages / 2.0))  // odd number of stages
            dim_sys = nx;

#if defined(LA_HIGH_PERFORMANCE)
        assign_strmat(dim_sys, dim_sys, memory->str_mat2[i], &c_ptr);
        assign_strmat(dim_sys, 1 + nf, memory->str_sol2[i], &c_ptr);
#elif defined(LA_REFERENCE)
        assign_strmat(dim_sys, dim_sys, memory->str_mat2[i], memory->sys_mat2[i]);
        assign_strmat(dim_sys, 1 + nf, memory->str_sol2[i], memory->sys_sol2[i]);
        d_cast_diag_mat2strmat((double *) c_ptr, memory->str_mat2[i]);
        c_ptr += dim_sys * sizeof(double);
#else  // LA_BLAS
        assign_strmat(dim_sys, dim_sys, memory->str_mat2[i], memory->sys_mat2[i]);
        assign_strmat(dim_sys, 1 + nf, memory->str_sol2[i], memory->sys_sol2[i]);
#endif  // LA_HIGH_PERFORMANCE
        dgesc_libstr(dim_sys, dim_sys, 0.0, memory->str_mat2[i], 0, 0);
        dgesc_libstr(dim_sys, 1 + nf, 0.0, memory->str_sol2[i], 0, 0);
        }
    }
#endif  // !TRIPLE_LOOP

    assert((char*)raw_memory + sim_lifted_irk_calculate_memory_size(dims, opts) >= c_ptr);

    // initialize
    for (int i = 0; i < num_steps * num_stages * nx; ++i)
        memory->K_traj[i] = 0.0;
    for (int i = 0; i < num_steps * num_stages * nx * nf; ++i)
        memory->DK_traj[i] = 0.0;
    for (int i = 0; i < num_steps * num_stages * nx; ++i)
        memory->mu_traj[i] = 0.0;

    if (opts->scheme->type == simplified_inis)
        for (int i = 0; i < num_steps * num_stages * nx * nf; ++i)
            memory->delta_DK_traj[i] = 0.0;
    if (opts->scheme->type == simplified_in || opts->scheme->type == simplified_inis) {
        for (int i = 0; i < num_steps * num_stages * nx; ++i)
            memory->adj_traj[i] = 0.0;
        for (int i = 0; i < num_steps * num_stages; ++i)
            for (int j = 0; j < nx * nx; ++j)
                memory->jac_traj[i][j] = 0.0;
    }

    return memory;
}



int sim_lifted_irk_calculate_workspace_size(sim_dims *dims, void *args) {
    int nx = dims->nx;
    int nu = dims->nu;
    sim_rk_opts *opts = (sim_rk_opts *)args;
    int num_stages = opts->num_stages;
    int NF = opts->num_forw_sens;
    //    int num_sys = ceil(num_stages/2.0);
    int dim_sys = num_stages * nx;
    if (opts->scheme->type == simplified_in ||
        opts->scheme->type == simplified_inis) {
        dim_sys = nx;
        if (num_stages > 1) dim_sys = 2 * nx;
    }

    int size = sizeof(sim_lifted_irk_workspace);
    size += (nx * (1 + NF) + nu + 1) * sizeof(real_t);  // rhs_in
    size += (nx * (1 + NF)) * sizeof(real_t);           // out_tmp
    if (opts->scheme->type == exact) {
        size += (dim_sys) * sizeof(int);             // ipiv
        size += (dim_sys * dim_sys) * sizeof(real_t);  // sys_mat
    }
    size += ((num_stages * nx) * (1 + NF)) * sizeof(real_t);  // sys_sol
    size += (num_stages) * sizeof(real_t *);                  // VDE_tmp
    size += (num_stages * nx * (1 + NF)) * sizeof(real_t);    // VDE_tmp[...]
    size += (nx * (nx + 1)) * sizeof(real_t);                 // jac_tmp

    if (opts->scheme->type == simplified_in ||
        opts->scheme->type == simplified_inis) {
        size +=
            ((num_stages * nx) * (1 + NF)) * sizeof(real_t);  // sys_sol_trans
        size += (num_stages * num_stages) * sizeof(real_t);   // trans
    }

    if (opts->scheme->type == simplified_in ||
        opts->scheme->type == simplified_inis) {
        size += (nx) * sizeof(real_t);  // out_adj_tmp
    }

#if !TRIPLE_LOOP
    size += sizeof(struct d_strmat);  // str_mat
    size += sizeof(struct d_strmat);  // str_sol

    int size_strmat = 0;
#if defined(LA_HIGH_PERFORMANCE)
    // matrices in matrix struct format:
    size_strmat += d_size_strmat(dim_sys, dim_sys);
    size_strmat += d_size_strmat(dim_sys, 1 + NF);

#elif defined(LA_REFERENCE)
    // allocate new memory only for the diagonal
    size_strmat += d_size_diag_strmat(dim_sys, dim_sys);

#endif  // LA_HIGH_PERFORMANCE
    size += size_strmat;
#endif  // !TRIPLE_LOOP

    return size;
}


static void sim_lifted_irk_cast_workspace(sim_lifted_irk_workspace *work,
                                          const sim_in *in, void *args) {
    int nx = in->nx;
    int nu = in->nu;
    sim_rk_opts *opts = (sim_rk_opts *)args;
    int num_stages = opts->num_stages;
    int NF = opts->num_forw_sens;
    //    int num_sys = ceil(num_stages/2.0);
    int dim_sys = num_stages * nx;
    if (opts->scheme->type == simplified_in ||
        opts->scheme->type == simplified_inis) {
        dim_sys = nx;
        if (num_stages > 1) dim_sys = 2 * nx;
    }

    char *ptr = (char *)work;
    ptr += sizeof(sim_lifted_irk_workspace);
    work->rhs_in = (real_t *)ptr;
    ptr += (nx * (1 + NF) + nu + 1) * sizeof(real_t);  // rhs_in
    work->out_tmp = (real_t *)ptr;
    ptr += (nx * (1 + NF)) * sizeof(real_t);  // out_tmp
    if (opts->scheme->type == exact) {
        work->ipiv = (int *)ptr;
        ptr += (dim_sys) * sizeof(int);  // ipiv
        work->sys_mat = (real_t *)ptr;
        ptr += (dim_sys * dim_sys) * sizeof(real_t);  // sys_mat
    }
    work->sys_sol = (real_t *)ptr;
    ptr += ((num_stages * nx) * (1 + NF)) * sizeof(real_t);  // sys_sol
    work->VDE_tmp = (real_t **)ptr;
    ptr += (num_stages) * sizeof(real_t *);  // VDE_tmp
    for (int i = 0; i < num_stages; i++) {
        work->VDE_tmp[i] = (real_t *)ptr;
        ptr += (nx * (1 + NF)) * sizeof(real_t);  // VDE_tmp[i]
    }
    work->jac_tmp = (real_t *)ptr;
    ptr += (nx * (nx + 1)) * sizeof(real_t);  // jac_tmp

    if (opts->scheme->type == simplified_in ||
        opts->scheme->type == simplified_inis) {
        work->sys_sol_trans = (real_t *)ptr;
        ptr +=
            ((num_stages * nx) * (1 + NF)) * sizeof(real_t);  // sys_sol_trans
        work->trans = (real_t *)ptr;
        ptr += (num_stages * num_stages) * sizeof(real_t);  // trans
    }

    if (opts->scheme->type == simplified_in ||
        opts->scheme->type == simplified_inis) {
        work->out_adj_tmp = (real_t *)ptr;
        ptr += (nx) * sizeof(real_t);  // out_adj_tmp
    }

#if !TRIPLE_LOOP
    work->str_mat = (struct d_strmat *)ptr;
    ptr += sizeof(struct d_strmat);
    work->str_sol = (struct d_strmat *)ptr;
    ptr += sizeof(struct d_strmat);
#if defined(LA_HIGH_PERFORMANCE)
    // matrices in matrix struct format:
    int size_strmat = 0;
    size_strmat += d_size_strmat(dim_sys, dim_sys);
    size_strmat += d_size_strmat(dim_sys, 1 + NF);

    d_create_strmat(dim_sys, dim_sys, work->str_mat, ptr);
    ptr += work->str_mat->memory_size;
    d_create_strmat(dim_sys, 1 + NF, work->str_sol, ptr);
    ptr += work->str_sol->memory_size;

#elif defined(LA_REFERENCE)

    //  pointer to column-major matrix
    d_create_strmat(dim_sys, dim_sys, work->str_mat, work->sys_mat);
    d_create_strmat(dim_sys, 1 + NF, work->str_sol, work->sys_sol);

    d_cast_diag_mat2strmat((double *)ptr, work->str_mat);
    ptr += d_size_diag_strmat(dim_sys, dim_sys);

#else  // LA_BLAS

    // not allocate new memory: point to column-major matrix
    d_create_strmat(dim_sys, dim_sys, work->str_mat, work->sys_mat);
    d_create_strmat(dim_sys, 1 + NF, work->str_sol, work->sys_sol);

#endif  // LA_HIGH_PERFORMANCE
#endif  // !TRIPLE_LOOP
}

#if TRIPLE_LOOP

#if CODE_GENERATION
#define DIM 6      // num_stages*NX
#define DIM_RHS 9  // NX+NU
#endif

real_t LU_system_ACADO(real_t *const A, int *const perm, int dim) {
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
    for (i = 0; i < (DIM - 1); i++) {
        indexMax = i;
        valueMax = fabs(A[i * DIM + i]);
        for (j = (i + 1); j < DIM; j++) {
            swap = fabs(A[i * DIM + j]);
            if (swap > valueMax) {
                indexMax = j;
                valueMax = swap;
            }
        }
        if (indexMax > i) {
            for (k = 0; k < DIM; ++k) {
                swap = A[k * DIM + i];
                A[k * DIM + i] = A[k * DIM + indexMax];
                A[k * DIM + indexMax] = swap;
            }
            intSwap = perm[i];
            perm[i] = perm[indexMax];
            perm[indexMax] = intSwap;
        }
        //        det *= A[i*DIM+i];
        for (j = i + 1; j < DIM; j++) {
            A[i * DIM + j] = -A[i * DIM + j] / A[i * DIM + i];
            for (k = i + 1; k < DIM; k++) {
                A[k * DIM + j] += A[i * DIM + j] * A[k * DIM + i];
            }
        }
    }
    //    det *= A[DIM*DIM-1];
    //    det = fabs(det);
    return det;
}

real_t solve_system_ACADO(real_t *const A, real_t *const b, int *const perm,
                          int dim, int dim2) {
    int i, j, k;
    int index1;

#if !CODE_GENERATION
    int DIM = dim;
    int DIM_RHS = dim2;
#else
    dim += 0;
    dim2 += 0;
#endif
    real_t *bPerm;
    bPerm = (real_t *) calloc(DIM*DIM_RHS, sizeof(real_t));
    real_t tmp_var;

    for (i = 0; i < DIM; ++i) {
        index1 = perm[i];
        for (j = 0; j < DIM_RHS; ++j) {
            bPerm[j * DIM + i] = b[j * DIM + index1];
        }
    }
    for (j = 1; j < DIM; ++j) {
        for (i = 0; i < j; ++i) {
            tmp_var = A[i * DIM + j];
            for (k = 0; k < DIM_RHS; ++k) {
                bPerm[k * DIM + j] += tmp_var * bPerm[k * DIM + i];
            }
        }
    }
    for (i = DIM - 1; - 1 < i; --i) {
        for (j = DIM - 1; i < j; --j) {
            tmp_var = A[j * DIM + i];
            for (k = 0; k < DIM_RHS; ++k) {
                bPerm[k * DIM + i] -= tmp_var * bPerm[k * DIM + j];
            }
        }
        tmp_var = 1.0 / A[i * (DIM + 1)];
        for (k = 0; k < DIM_RHS; ++k) {
            bPerm[k * DIM + i] = tmp_var * bPerm[k * DIM + i];
        }
    }
    for (k = 0; k < DIM * DIM_RHS; ++k) {
        b[k] = bPerm[k];
    }
    return 0;
}

real_t solve_system_trans_ACADO(real_t *const A, real_t *const b,
                                int *const perm, int dim, int dim2) {
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
    real_t *bPerm;
    bPerm = (real_t *) calloc(DIM*DIM_RHS, sizeof(real_t));
    real_t tmp_var;

    for (k = 0; k < DIM * DIM_RHS; ++k) {
        bPerm[k] = b[k];
    }
    for (i = 0; i < DIM; i++) {
        for (j = 0; j < i; j++) {
            tmp_var = A[i * DIM + j];
            for (k = 0; k < DIM_RHS; ++k) {
                bPerm[k * DIM + i] -= tmp_var * bPerm[k * DIM + j];
            }
        }
        tmp_var = 1.0 / A[i * (DIM + 1)];
        for (k = 0; k < DIM_RHS; ++k) {
            bPerm[k * DIM + i] = tmp_var * bPerm[k * DIM + i];
        }
    }
    for (j = DIM - 1; j > -1; --j) {
        for (i = DIM - 1; i > j; --i) {
            tmp_var = A[j * DIM + i];
            for (k = 0; k < DIM_RHS; ++k) {
                bPerm[k * DIM + j] += tmp_var * bPerm[k * DIM + i];
            }
        }
    }
    for (i = 0; i < DIM; ++i) {
        index1 = perm[i];
        for (j = 0; j < DIM_RHS; ++j) {
            b[j * DIM + index1] = bPerm[j * DIM + i];
        }
    }
    return 0;
}

#endif

void transform_mat(real_t *mat, real_t *trans, real_t *mat_trans,
                   const int stages, const int n, const int m) {
    //    print_matrix_name("stdout", "trans", trans, stages, stages);
    for (int i = 0; i < stages * n * m; i++) mat_trans[i] = 0.0;
    for (int j = 0; j < m; j++) {
        for (int s2 = 0; s2 < stages; s2++) {
            for (int s1 = 0; s1 < stages; s1++) {
                if (trans[s2 * stages + s1] != 0) {
                    for (int i = 0; i < n; i++) {
                        mat_trans[j * (stages * n) + s1 * n + i] +=
                            trans[s2 * stages + s1] *
                            mat[j * (stages * n) + s2 * n + i];
                    }
                }
            }
        }
    }
}

void transform_vec(real_t *mat, real_t *trans, real_t *mat_trans,
                   const int stages, const int n) {
    for (int i = 0; i < stages * n; i++) mat_trans[i] = 0.0;
    for (int s2 = 0; s2 < stages; s2++) {
        for (int s1 = 0; s1 < stages; s1++) {
            for (int i = 0; i < n; i++) {
                mat_trans[s1 * n + i] +=
                    trans[s2 * stages + s1] * mat[s2 * n + i];
            }
        }
    }
}

void construct_subsystems(real_t *mat, real_t **mat2, const int stages,
                          const int n, const int m) {
    int idx = 0;
    for (int s1 = 0; s1 < stages; s1++) {
        if ((s1 + 1) < stages) {  // complex conjugate pair of eigenvalues
            for (int j = 0; j < m; j++) {
                for (int i = 0; i < n; i++) {
                    mat2[idx][j * 2 * n + i] =
                        mat[j * (stages * n) + s1 * n + i];
                }
            }
            s1++;
            for (int j = 0; j < m; j++) {
                for (int i = 0; i < n; i++) {
                    mat2[idx][j * 2 * n + n + i] =
                        mat[j * (stages * n) + s1 * n + i];
                }
            }
        } else {  // real eigenvalue
            for (int j = 0; j < m; j++) {
                for (int i = 0; i < n; i++) {
                    mat2[idx][j * n + i] = mat[j * (stages * n) + s1 * n + i];
                }
            }
        }
        idx++;
    }
}

void destruct_subsystems(real_t *mat, real_t **mat2, const int stages,
                         const int n, const int m) {
    int idx = 0;
    for (int s1 = 0; s1 < stages; s1++) {
        if ((s1 + 1) < stages) {  // complex conjugate pair of eigenvalues
            for (int j = 0; j < m; j++) {
                for (int i = 0; i < n; i++) {
                    mat[j * (stages * n) + s1 * n + i] =
                        mat2[idx][j * 2 * n + i];
                }
            }
            s1++;
            for (int j = 0; j < m; j++) {
                for (int i = 0; i < n; i++) {
                    mat[j * (stages * n) + s1 * n + i] =
                        mat2[idx][j * 2 * n + n + i];
                }
            }
        } else {  // real eigenvalue
            for (int j = 0; j < m; j++) {
                for (int i = 0; i < n; i++) {
                    mat[j * (stages * n) + s1 * n + i] = mat2[idx][j * n + i];
                }
            }
        }
        idx++;
    }
}

void form_linear_system_matrix(int istep, const sim_in *in, void *args,
                               sim_lifted_irk_memory *mem,
                               sim_lifted_irk_workspace *work, real_t *sys_mat,
                               real_t **sys_mat2, real_t timing_ad) {
    int nx = in->nx;
    int nu = in->nu;
    real_t H_INT = in->step;
    sim_rk_opts *opts = (sim_rk_opts *)args;
    int num_stages = opts->num_stages;
    real_t *A_mat = opts->A_mat;
    real_t *c_vec = opts->c_vec;

    real_t tmp_eig, tmp_eig2;
    acados_timer timer_ad;
    real_t *out_tmp = work->out_tmp;
    real_t *rhs_in = work->rhs_in;
    real_t *jac_tmp = work->jac_tmp;

    real_t **jac_traj = mem->jac_traj;
    real_t *K_traj = mem->K_traj;

    int i, j, s1, s2;
    if ((opts->scheme->type == simplified_in ||
         opts->scheme->type == simplified_inis) &&
        istep == 0 && !opts->scheme->freeze) {
        int idx = 0;
        for (s1 = 0; s1 < num_stages; s1++) {
            if ((s1 + 1) == num_stages) {  // real eigenvalue
                for (i = 0; i < nx * nx; i++) sys_mat2[idx][i] = 0.0;
                tmp_eig = 1.0 / H_INT * opts->scheme->eig[s1];
                for (i = 0; i < nx; i++) {
                    sys_mat2[idx][i * (nx + 1)] = tmp_eig;
                }
            } else {  // complex conjugate pair
                for (i = 0; i < 4 * nx * nx; i++) sys_mat2[idx][i] = 0.0;
                tmp_eig = 1.0 / H_INT * opts->scheme->eig[s1];
                tmp_eig2 = 1.0 / H_INT * opts->scheme->eig[s1 + 1];
                for (i = 0; i < nx; i++) {
                    sys_mat2[idx][i * (2 * nx + 1)] = tmp_eig;
                    sys_mat2[idx][(nx + i) * (2 * nx + 1)] = tmp_eig;

                    sys_mat2[idx][i * 2 * nx + (nx + i)] = tmp_eig2;
                    sys_mat2[idx][(nx + i) * 2 * nx + i] = -tmp_eig2;
                }
                s1++;  // skip the complex conjugate eigenvalue
            }
            idx++;
        }
    } else if (opts->scheme->type == exact) {
        for (i = 0; i < num_stages * nx * num_stages * nx; i++)
            sys_mat[i] = 0.0;
        for (i = 0; i < num_stages * nx; i++)
            sys_mat[i * (num_stages * nx + 1)] = 1.0;  // identity
    }

    for (s1 = 0; s1 < num_stages; s1++) {
        //                if (opts->scheme->type == exact || s1 == 0) {
        for (i = 0; i < nx; i++) {
            rhs_in[i] = out_tmp[i];
        }
        for (i = 0; i < nu; i++) rhs_in[nx + i] = in->u[i];
        rhs_in[nx + nu] = (istep + c_vec[s1]) / opts->num_steps;  // time
        for (s2 = 0; s2 < num_stages; s2++) {
            for (i = 0; i < nx; i++) {
                rhs_in[i] += H_INT * A_mat[s2 * num_stages + s1] *
                             K_traj[istep * num_stages * nx + s2 * nx + i];
            }
        }
        acados_tic(&timer_ad);
        in->jacobian_wrapper(nx, rhs_in, jac_tmp, in->jac);  // k evaluation
        timing_ad += acados_toc(&timer_ad);
        //                }
        if (opts->scheme->type == simplified_in ||
            opts->scheme->type == simplified_inis) {
            for (i = 0; i < nx * nx; i++)
                jac_traj[istep * num_stages + s1][i] = jac_tmp[nx + i];
        }

        // put jac_tmp in sys_mat:
        if (opts->scheme->type == exact) {
            for (s2 = 0; s2 < num_stages; s2++) {
                for (j = 0; j < nx; j++) {
                    for (i = 0; i < nx; i++) {
                        sys_mat[(s2 * nx + j) * num_stages * nx + s1 * nx +
                                i] -= H_INT * A_mat[s2 * num_stages + s1] *
                                      jac_tmp[nx + j * nx + i];
                    }
                }
            }
        }
    }

    int idx = 0;
    for (s1 = 0; s1 < num_stages; s1++) {
        // put jac_traj[0] in sys_mat:
        if ((opts->scheme->type == simplified_in ||
             opts->scheme->type == simplified_inis) &&
            istep == 0 && !opts->scheme->freeze) {
            if ((s1 + 1) == num_stages) {  // real eigenvalue
                for (j = 0; j < nx; j++) {
                    for (i = 0; i < nx; i++) {
                        sys_mat2[idx][j * nx + i] -= jac_traj[0][j * nx + i];
                    }
                }
            } else {  // complex conjugate pair
                for (j = 0; j < nx; j++) {
                    for (i = 0; i < nx; i++) {
                        sys_mat2[idx][j * 2 * nx + i] -=
                            jac_traj[0][j * nx + i];
                        sys_mat2[idx][(nx + j) * 2 * nx + nx + i] -=
                            jac_traj[0][j * nx + i];
                    }
                }
                s1++;  // skip the complex conjugate eigenvalue
            }
            idx++;
        }
    }
}

int sim_lifted_irk(sim_in *in, sim_out *out, void *args, void *mem_, void *work_) {

    int nx = in->nx;
    int nu = in->nu;
    sim_rk_opts *opts = (sim_rk_opts *)args;
    int num_stages = opts->num_stages;
    int dim_sys = num_stages * nx;
    int i, s1, s2, j, istep;
    sim_lifted_irk_memory *mem = (sim_lifted_irk_memory *)mem_;
    sim_lifted_irk_workspace *work = (sim_lifted_irk_workspace *)work_;
    sim_lifted_irk_cast_workspace(work, in, args);
    real_t H_INT = in->step;
    int NF = opts->num_forw_sens;

    real_t *A_mat = opts->A_mat;
    real_t *b_vec = opts->b_vec;
    real_t *c_vec = opts->c_vec;

    real_t **VDE_tmp = work->VDE_tmp;
    real_t *out_tmp = work->out_tmp;
    real_t *rhs_in = work->rhs_in;

    real_t **jac_traj = mem->jac_traj;

    real_t *K_traj = mem->K_traj;
    real_t *DK_traj = mem->DK_traj;
    real_t *delta_DK_traj = mem->delta_DK_traj;
    real_t *mu_traj = mem->mu_traj;
    real_t *adj_traj = mem->adj_traj;

    real_t *adj_tmp = work->out_adj_tmp;

    int *ipiv = work->ipiv;    // pivoting vector
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

    acados_timer timer, timer_la, timer_ad;
    real_t timing_la = 0.0;
    real_t timing_ad = 0.0;

    assert(NF == nx + nu && "Not implemented yet for other num_forw_sens");

    acados_tic(&timer);
    for (i = 0; i < nx; i++) out_tmp[i] = in->x[i];
    for (i = 0; i < nx * NF; i++)
        out_tmp[nx + i] = in->S_forw[i];  // sensitivities

    for (i = 0; i < nu; i++) rhs_in[nx * (1 + NF) + i] = in->u[i];

    if (opts->scheme->type == simplified_in ||
        opts->scheme->type == simplified_inis) {
        for (i = 0; i < nx; i++) adj_tmp[i] = in->S_adj[i];
    }

    // Newton step of the collocation variables with respect to the inputs:
    if (NF == (nx + nu) && opts->sens_adj) {
        for (istep = opts->num_steps - 1; istep > -1; istep--) {  // ADJOINT update
            for (s1 = 0; s1 < num_stages; s1++) {
                for (j = 0; j < nx; j++) {  // step in X
                    for (i = 0; i < nx; i++) {
                        K_traj[(istep * num_stages + s1) * nx + i] +=
                            DK_traj[(istep * num_stages + s1) * nx * (nx + nu) +
                                    j * nx + i] *
                            (in->x[j] - mem->x[j]);  // RK step
                    }
                    mem->x[j] = in->x[j];
                }
                for (j = 0; j < nu; j++) {  // step in U
                    for (i = 0; i < nx; i++) {
                        K_traj[(istep * num_stages + s1) * nx + i] +=
                            DK_traj[(istep * num_stages + s1) * nx * (nx + nu) +
                                    (nx + j) * nx + i] *
                            (in->u[j] - mem->u[j]);  // RK step
                    }
                    mem->u[j] = in->u[j];
                }
            }
            if (opts->scheme->type == simplified_inis && !opts->scheme->freeze) {
                for (s1 = 0; s1 < num_stages; s1++) {
                    for (j = 0; j < NF; j++) {
                        for (i = 0; i < nx; i++) {
                            DK_traj[(istep * num_stages + s1) * nx * NF +
                                    j * nx + i] +=
                                delta_DK_traj[(istep * num_stages + s1) * nx *
                                                  NF +
                                              j * nx + i];
                        }
                    }
                }
            }

            // Newton step of the Lagrange multipliers mu, based on adj_traj:
            if (opts->scheme->type == simplified_in ||
                opts->scheme->type == simplified_inis) {
                for (s1 = 0; s1 < num_stages; s1++) {
                    for (i = 0; i < nx; i++) {
                        sys_sol[s1 * nx + i] =
                            -mu_traj[istep * num_stages * nx + s1 * nx + i];
                        for (s2 = 0; s2 < num_stages; s2++) {
                            sys_sol[s1 * nx + i] +=
                                H_INT * A_mat[s1 * num_stages + s2] *
                                adj_traj[istep * num_stages * nx + s2 * nx + i];
                        }
                        sys_sol[s1 * nx + i] -= H_INT * b_vec[s1] * adj_tmp[i];
                        sys_sol[s1 * nx + i] += mem->grad_K[s1 * nx + i];
                    }
                }
                //                print_matrix("stdout", sys_sol, 1,
                //                num_stages*nx); print_matrix("stdout",
                //                sys_sol, 1, 1);

                // TRANSFORM using transf1_T:
                if (opts->scheme->type == simplified_in ||
                    opts->scheme->type == simplified_inis) {
                    // apply the transf1 operation:
                    for (s1 = 0; s1 < num_stages; s1++) {
                        for (s2 = 0; s2 < num_stages; s2++) {
                            work->trans[s2 * num_stages + s1] =
                                1.0 / H_INT *
                                opts->scheme->transf1_T[s2 * num_stages + s1];
                        }
                    }
                    transform_vec(sys_sol, work->trans, sys_sol_trans,
                                  num_stages, nx);

                    // construct sys_sol2 from sys_sol_trans:
                    construct_subsystems(sys_sol_trans, sys_sol2, num_stages,
                                         nx, 1);
                }
                acados_tic(&timer_la);
                int idx = 0;
                for (s1 = 0; s1 < num_stages; s1++) {
                    // THIS LOOP IS PARALLELIZABLE BECAUSE OF DECOMPOSABLE
                    // LINEAR SUBSYSTEMS
                    if (opts->scheme->type == simplified_in ||
                        opts->scheme->type == simplified_inis) {
                        sys_mat = sys_mat2[idx];
                        ipiv = ipiv2[idx];
                        sys_sol = sys_sol2[idx];
#if !TRIPLE_LOOP
                        str_mat = str_mat2[idx];
                        str_sol = str_sol2[idx];
#endif
                        idx++;
                        if ((s1 + 1) < num_stages) {
                            dim_sys = 2 * nx;
                            s1++;  // complex conjugate pair of eigenvalues
                        } else {
                            dim_sys = nx;
                        }
                    } else {
                        dim_sys = num_stages * nx;
                        s1 = num_stages;  // break out of for-loop
                    }
#if TRIPLE_LOOP
                    solve_system_trans_ACADO(sys_mat, sys_sol, ipiv, dim_sys,
                                             1);
#else  // TRIPLE_LOOP
#error : NOT YET IMPLEMENTED
#endif  // TRIPLE_LOOP
                }
                timing_la += acados_toc(&timer_la);
                // TRANSFORM using transf2_T:
                if (opts->scheme->type == simplified_in ||
                    opts->scheme->type == simplified_inis) {
                    // construct sys_sol_trans from sys_sol2:
                    destruct_subsystems(sys_sol_trans, sys_sol2, num_stages, nx,
                                        1);

                    // apply the transf2 operation:
                    sys_sol = work->sys_sol;
                    transform_vec(sys_sol_trans, opts->scheme->transf2_T,
                                  sys_sol, num_stages, nx);
                }

                // update mu_traj
                for (s1 = 0; s1 < num_stages; s1++) {
                    for (i = 0; i < nx; i++) {
                        mu_traj[istep * num_stages * nx + s1 * nx + i] +=
                            sys_sol[s1 * nx + i];
                    }
                }

                // update adj_tmp:
                // TODO(rien): USE ADJOINT DIFFERENTIATION HERE INSTEAD !!:
                for (j = 0; j < nx; j++) {
                    for (s1 = 0; s1 < num_stages; s1++) {
                        for (i = 0; i < nx; i++) {
                            adj_tmp[j] -=
                                mu_traj[istep * num_stages * nx + s1 * nx + i] *
                                jac_traj[istep * num_stages + s1][j * nx + i];
                        }
                    }
                }
            }
        }
    }

    if (opts->scheme->type == simplified_in ||
        opts->scheme->type == simplified_inis) {
        for (i = 0; i < NF; i++) out->grad[i] = 0.0;
    }

    for (istep = 0; istep < opts->num_steps; istep++) {
        // form linear system matrix (explicit ODE case):

        form_linear_system_matrix(istep, in, args, mem, work, sys_mat, sys_mat2, timing_ad);

        int idx;
        if (opts->scheme->type == exact || (istep == 0 && !opts->scheme->freeze)) {
            acados_tic(&timer_la);
            idx = 0;
            for (s1 = 0; s1 < num_stages; s1++) {
                // THIS LOOP IS PARALLELIZABLE BECAUSE OF DECOMPOSABLE LINEAR
                // SUBSYSTEMS
                if (opts->scheme->type == simplified_in ||
                    opts->scheme->type == simplified_inis) {
                    sys_mat = sys_mat2[idx];
                    ipiv = ipiv2[idx];
#if !TRIPLE_LOOP
                    str_mat = str_mat2[idx];
#endif
                    idx++;
                    if ((s1 + 1) < num_stages) {
                        dim_sys = 2 * nx;
                        s1++;  // complex conjugate pair of eigenvalues
                    } else {
                        dim_sys = nx;
                    }
                } else {
                    dim_sys = num_stages * nx;
                    s1 = num_stages;  // break out of for-loop
                }
#if TRIPLE_LOOP
                LU_system_ACADO(sys_mat, ipiv, dim_sys);
#else  // TRIPLE_LOOP
// ---- BLASFEO: LU factorization ----
#if defined(LA_HIGH_PERFORMANCE)
                d_cvt_mat2strmat(dim_sys, dim_sys, sys_mat, dim_sys, str_mat, 0,
                                 0);  // mat2strmat
#endif  // LA_BLAS | LA_REFERENCE
                dgetrf_libstr(dim_sys, dim_sys, str_mat, 0, 0, str_mat, 0, 0,
                              ipiv);  // Gauss elimination
// ---- BLASFEO: LU factorization ----
#endif  // TRIPLE_LOOP
            }
            timing_la += acados_toc(&timer_la);
        }

        if (opts->scheme->type == simplified_in ||
            opts->scheme->type == simplified_inis) {
            sys_sol = work->sys_sol;
            sys_sol_trans = work->sys_sol_trans;
        }
        for (s1 = 0; s1 < num_stages; s1++) {
            for (i = 0; i < nx * (1 + NF); i++) {
                rhs_in[i] = out_tmp[i];
            }
            for (s2 = 0; s2 < num_stages; s2++) {
                for (i = 0; i < nx; i++) {
                    rhs_in[i] += H_INT * A_mat[s2 * num_stages + s1] *
                                 K_traj[istep * num_stages * nx + s2 * nx + i];
                }
                if (opts->scheme->type == simplified_inis) {
                    for (j = 0; j < NF; j++) {
                        for (i = 0; i < nx; i++) {
                            rhs_in[(j + 1) * nx + i] +=
                                H_INT * A_mat[s2 * num_stages + s1] *
                                DK_traj[(istep * num_stages + s2) * nx * NF +
                                        j * nx + i];
                        }
                    }
                }
            }
            rhs_in[nx*(1+NF)+nu] = ((real_t) istep+c_vec[s1])/((real_t) opts->num_steps);  // time

            acados_tic(&timer_ad);
            in->forward_vde_wrapper(nx, nu, rhs_in, VDE_tmp[s1], in->vde);  // k evaluation
            timing_ad += acados_toc(&timer_ad);

            // put VDE_tmp in sys_sol:
            for (j = 0; j < 1 + NF; j++) {
                for (i = 0; i < nx; i++) {
                    sys_sol[j * num_stages * nx + s1 * nx + i] =
                        VDE_tmp[s1][j * nx + i];
                }
            }
            for (i = 0; i < nx; i++) {
                sys_sol[s1 * nx + i] -=
                    K_traj[istep * num_stages * nx + s1 * nx + i];
            }
            if (opts->scheme->type == simplified_inis) {
                for (j = 0; j < NF; j++) {
                    for (i = 0; i < nx; i++) {
                        sys_sol[(j + 1) * num_stages * nx + s1 * nx + i] -=
                            DK_traj[(istep * num_stages + s1) * nx * NF +
                                    j * nx + i];
                    }
                }
            }
        }

        // Inexact Newton with Iterated Sensitivities (INIS) based gradient
        // correction:
        if (opts->scheme->type == simplified_inis) {
            for (j = 0; j < NF; j++) {
                for (i = 0; i < num_stages * nx; i++) {
                    out->grad[j] += mu_traj[istep * num_stages * nx + i] *
                                    sys_sol[(j + 1) * num_stages * nx + i];
                }
            }
        }

        if (opts->scheme->type == simplified_in ||
            opts->scheme->type == simplified_inis) {
            // apply the transf1 operation:
            for (s1 = 0; s1 < num_stages; s1++) {
                for (s2 = 0; s2 < num_stages; s2++) {
                    work->trans[s2 * num_stages + s1] =
                        1.0 / H_INT *
                        opts->scheme->transf1[s2 * num_stages + s1];
                }
            }
            if (!opts->scheme->freeze) {
                transform_mat(sys_sol, work->trans, sys_sol_trans, num_stages,
                              nx, 1 + NF);
            } else {
                transform_mat(sys_sol, work->trans, sys_sol_trans, num_stages,
                              nx, 1);
            }

            // construct sys_sol2 from sys_sol_trans:
            if (!opts->scheme->freeze) {
                construct_subsystems(sys_sol_trans, sys_sol2, num_stages, nx,
                                     1 + NF);
            } else {
                construct_subsystems(sys_sol_trans, sys_sol2, num_stages, nx,
                                     1);
            }
        }

        acados_tic(&timer_la);
        idx = 0;
        for (s1 = 0; s1 < num_stages; s1++) {
            // THIS LOOP IS PARALLELIZABLE BECAUSE OF DECOMPOSABLE LINEAR
            // SUBSYSTEMS
            if (opts->scheme->type == simplified_in ||
                opts->scheme->type == simplified_inis) {
                sys_mat = sys_mat2[idx];
                ipiv = ipiv2[idx];
                sys_sol = sys_sol2[idx];
#if !TRIPLE_LOOP
                str_mat = str_mat2[idx];
                str_sol = str_sol2[idx];
#endif
                idx++;
                if ((s1 + 1) < num_stages) {
                    dim_sys = 2 * nx;
                    s1++;  // complex conjugate pair of eigenvalues
                } else {
                    dim_sys = nx;
                }
            } else {
                dim_sys = num_stages * nx;
                s1 = num_stages;  // break out of for-loop
            }
#if TRIPLE_LOOP
            if (!opts->scheme->freeze) {
                solve_system_ACADO(sys_mat, sys_sol, ipiv, dim_sys, 1 + NF);
            } else {
                solve_system_ACADO(sys_mat, sys_sol, ipiv, dim_sys, 1);
            }
#else  // TRIPLE_LOOP
#if defined(LA_HIGH_PERFORMANCE)
            // ---- BLASFEO: row transformations + backsolve ----
            d_cvt_mat2strmat(dim_sys, 1 + NF, sys_sol, dim_sys, str_sol, 0,
                             0);                    // mat2strmat
            drowpe_libstr(dim_sys, ipiv, str_sol);  // row permutations
            dtrsm_llnu_libstr(dim_sys, 1 + NF, 1.0, str_mat, 0, 0, str_sol, 0,
                              0, str_sol, 0, 0);  // L backsolve
            dtrsm_lunn_libstr(dim_sys, 1 + NF, 1.0, str_mat, 0, 0, str_sol, 0,
                              0, str_sol, 0, 0);  // U backsolve
            d_cvt_strmat2mat(dim_sys, 1 + NF, str_sol, 0, 0, sys_sol,
                             dim_sys);  // strmat2mat
                                        // BLASFEO: row transformations + backsolve
#else   // LA_BLAS | LA_REFERENCE
            // ---- BLASFEO: row transformations + backsolve ----
            drowpe_libstr(dim_sys, ipiv, str_sol);  // row permutations
            dtrsm_llnu_libstr(dim_sys, 1 + NF, 1.0, str_mat, 0, 0, str_sol, 0,
                              0, str_sol, 0, 0);  // L backsolve
            dtrsm_lunn_libstr(dim_sys, 1 + NF, 1.0, str_mat, 0, 0, str_sol, 0,
                              0, str_sol, 0, 0);  // U backsolve
                                                  // BLASFEO: row transformations + backsolve
#endif  // LA_BLAFEO
#endif  // TRIPLE_LOOP
        }
        timing_la += acados_toc(&timer_la);
        if (opts->scheme->type == simplified_in || opts->scheme->type == simplified_inis) {
            // construct sys_sol_trans from sys_sol2:
            if (!opts->scheme->freeze) {
                destruct_subsystems(sys_sol_trans, sys_sol2, num_stages, nx,
                                    1 + NF);
            } else {
                destruct_subsystems(sys_sol_trans, sys_sol2, num_stages, nx, 1);
            }

            // apply the transf2 operation:
            sys_sol = work->sys_sol;
            if (!opts->scheme->freeze) {
                transform_mat(sys_sol_trans, opts->scheme->transf2, sys_sol,
                              num_stages, nx, 1 + NF);
            } else {
                transform_mat(sys_sol_trans, opts->scheme->transf2, sys_sol,
                              num_stages, nx, 1);
            }
        }

        // Newton step of the collocation variables
        for (i = 0; i < num_stages * nx; i++) {
            K_traj[istep * num_stages * nx + i] += sys_sol[i];
        }
        for (s1 = 0; s1 < num_stages; s1++) {
            for (i = 0; i < nx; i++) {
                out_tmp[i] +=
                    H_INT * b_vec[s1] *
                    K_traj[istep * num_stages * nx + s1 * nx + i];  // RK step
            }
        }

        // Sensitivities collocation variables
        for (s1 = 0; s1 < num_stages; s1++) {
            for (j = 0; j < NF; j++) {
                for (i = 0; i < nx; i++) {
                    if (opts->scheme->type == simplified_inis && opts->sens_adj &&
                        !opts->scheme->freeze) {
                        delta_DK_traj[(istep * num_stages + s1) * nx * NF +
                                      j * nx + i] =
                            sys_sol[(j + 1) * num_stages * nx + s1 * nx + i];
                    } else if (opts->scheme->type == simplified_inis &&
                               !opts->scheme->freeze) {
                        DK_traj[(istep * num_stages + s1) * nx * NF + j * nx +
                                i] +=
                            sys_sol[(j + 1) * num_stages * nx + s1 * nx + i];
                    } else if (!opts->scheme->freeze) {
                        DK_traj[(istep * num_stages + s1) * nx * NF + j * nx +
                                i] =
                            sys_sol[(j + 1) * num_stages * nx + s1 * nx + i];
                    }
                }
            }
        }
        for (s1 = 0; s1 < num_stages; s1++) {
            for (i = 0; i < nx * NF; i++) {
                out_tmp[nx + i] += H_INT * b_vec[s1] *
                                   DK_traj[(istep * num_stages + s1) * nx * NF +
                                           i];  // RK step
            }
        }
        if (opts->scheme->type == simplified_inis ||
            opts->scheme->type == simplified_in) {
            // Adjoint derivatives:
            for (s1 = 0; s1 < num_stages; s1++) {
                for (j = 0; j < nx; j++) {
                    adj_traj[istep * num_stages * nx + s1 * nx + j] = 0.0;
                    for (i = 0; i < nx; i++) {
                        adj_traj[istep * num_stages * nx + s1 * nx + j] +=
                            mu_traj[istep * num_stages * nx + s1 * nx + i] *
                            jac_traj[istep * num_stages + s1][j * nx + i];
                    }
                }
            }
        }
        //        print_matrix_name("stdout", "adj_traj", &adj_traj[0], 1,
        //        num_stages*nx); print_matrix_name("stdout", "mu_traj",
        //        &mu_traj[0], 1, num_stages*nx); print_matrix_name("stdout",
        //        "DK_traj", &DK_traj[0], 1, nx*NF); print_matrix_name("stdout",
        //        "VDE_tmp[0]", &VDE_tmp[0][0], 1, nx*(NF+1));
        //        print_matrix_name("stdout", "VDE_tmp[1]", &VDE_tmp[1][0], 1,
        //        nx*(NF+1));
        if (opts->scheme->type == simplified_in) {
            // Standard Inexact Newton based gradient correction:
            for (j = 0; j < NF; j++) {
                for (s1 = 0; s1 < num_stages; s1++) {
                    for (i = 0; i < nx; i++) {
                        out->grad[j] +=
                            mu_traj[istep * num_stages * nx + s1 * nx + i] *
                            VDE_tmp[s1][(j + 1) * nx + i];
                    }
                }
            }
            for (j = 0; j < NF; j++) {
                for (s1 = 0; s1 < num_stages; s1++) {
                    for (i = 0; i < nx; i++) {
                        out->grad[j] -=
                            mu_traj[istep * num_stages * nx + s1 * nx + i] *
                            DK_traj[(istep * num_stages + s1) * nx * NF +
                                    j * nx + i];
                    }
                }
                for (s2 = 0; s2 < num_stages; s2++) {
                    for (s1 = 0; s1 < num_stages; s1++) {
                        for (i = 0; i < nx; i++) {
                            out->grad[j] +=
                                H_INT * A_mat[s2 * num_stages + s1] *
                                adj_traj[istep * num_stages * nx + s1 * nx +
                                         i] *
                                DK_traj[(istep * num_stages + s2) * nx * NF +
                                        j * nx + i];
                        }
                    }
                }
            }
        }
        //        print_matrix_name("stdout", "grad", &out->grad[0], 1, NF);
    }
    for (i = 0; i < nx; i++) out->xn[i] = out_tmp[i];
    for (i = 0; i < nx * NF; i++) out->S_forw[i] = out_tmp[nx + i];

    out->info->CPUtime = acados_toc(&timer);
    out->info->LAtime = timing_la;
    out->info->ADtime = timing_ad;

    return 0;  // success
}
