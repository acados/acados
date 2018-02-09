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



#define FW_VDE_NUMIN 5
#define FW_VDE_NUMOUT 3
#define JAC_ODE_NUMIN 3
#define JAC_ODE_NUMOUT 2



int sim_lifted_irk_forward_vde_input_dims[FW_VDE_NUMIN] = {0};
int sim_lifted_irk_forward_vde_output_dims[FW_VDE_NUMOUT] = {0};
external_function_dims sim_lifted_irk_forward_vde_dims = {
    FW_VDE_NUMIN,
    FW_VDE_NUMOUT,
    sim_lifted_irk_forward_vde_input_dims,
    sim_lifted_irk_forward_vde_output_dims
};



int sim_lifted_irk_jacobian_ode_input_dims[JAC_ODE_NUMIN] = {0};
int sim_lifted_irk_jacobian_ode_output_dims[JAC_ODE_NUMOUT] = {0};
external_function_dims sim_lifted_irk_jacobian_ode_dims = {
    JAC_ODE_NUMIN,
    JAC_ODE_NUMOUT,
    sim_lifted_irk_jacobian_ode_input_dims,
    sim_lifted_irk_jacobian_ode_output_dims
};



int sim_lifted_irk_integrator_calculate_args_size(sim_dims *dims, void *submodules_)
{
    sim_lifted_irk_integrator_submodules *submodules = (sim_lifted_irk_integrator_submodules *) submodules_;

    int size = sizeof(sim_lifted_irk_integrator_args);

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

    size += 2*sizeof(external_function_fcn_ptrs);

    if (submodules->forward_vde != NULL) {
        size += submodules->forward_vde->calculate_args_size(&sim_lifted_irk_forward_vde_dims, submodules->forward_vde->submodules);
    }

    if (submodules->jacobian_ode != NULL) {
        size += submodules->jacobian_ode->calculate_args_size(&sim_lifted_irk_jacobian_ode_dims, submodules->jacobian_ode->submodules);
    }

    make_int_multiple_of(8, &size);
    size += 2 * 8;

    return size;
}



void *sim_lifted_irk_integrator_assign_args(sim_dims *dims, void **submodules_, void *raw_memory)
{
    sim_lifted_irk_integrator_submodules *submodules = (sim_lifted_irk_integrator_submodules *) *submodules_;
    
    char *c_ptr = (char *)raw_memory;

    sim_lifted_irk_integrator_args *args = (sim_lifted_irk_integrator_args *)c_ptr;
    c_ptr += sizeof(sim_lifted_irk_integrator_args);

    int ns = dims->num_stages;
    args->num_stages = ns;

    align_char_to(8, &c_ptr);

    assign_double(ns * ns, &args->A_mat, &c_ptr);
    assign_double(ns, &args->b_vec, &c_ptr);
    assign_double(ns, &args->c_vec, &c_ptr);

    args->scheme = (Newton_scheme *)c_ptr;
    c_ptr += sizeof(Newton_scheme);

    align_char_to(8, &c_ptr);

    assign_double(ns, &args->scheme->eig, &c_ptr);

    assign_double(ns * ns, &args->scheme->transf1, &c_ptr);
    assign_double(ns * ns, &args->scheme->transf2, &c_ptr);
    assign_double(ns * ns, &args->scheme->transf1_T, &c_ptr);
    assign_double(ns * ns, &args->scheme->transf2_T, &c_ptr);

    if (submodules->forward_vde != NULL) {
        args->submodules.forward_vde = (external_function_fcn_ptrs *)c_ptr;
        c_ptr += sizeof(external_function_fcn_ptrs);

        void *forward_vde_submodules = submodules->forward_vde->submodules;
        args->forward_vde_args = submodules->forward_vde->assign_args(&sim_lifted_irk_forward_vde_dims, &(forward_vde_submodules), c_ptr);
        c_ptr += submodules->forward_vde->calculate_args_size(&sim_lifted_irk_forward_vde_dims, submodules->forward_vde->submodules);

        *(args->submodules.forward_vde) = *(submodules->forward_vde);
        args->submodules.forward_vde->submodules = forward_vde_submodules;
    } else {
        args->submodules.forward_vde = NULL;
    }

    if (submodules->jacobian_ode != NULL) {
        args->submodules.jacobian_ode = (external_function_fcn_ptrs *)c_ptr;
        c_ptr += sizeof(external_function_fcn_ptrs);

        void *jacobian_ode_submodules = submodules->jacobian_ode->submodules;
        args->jacobian_ode_args = submodules->jacobian_ode->assign_args(&sim_lifted_irk_jacobian_ode_dims, &(jacobian_ode_submodules), c_ptr);
        c_ptr += submodules->jacobian_ode->calculate_args_size(&sim_lifted_irk_jacobian_ode_dims, submodules->jacobian_ode->submodules);

        *(args->submodules.jacobian_ode) = *(submodules->jacobian_ode);
        args->submodules.jacobian_ode->submodules = jacobian_ode_submodules;
    } else {
        args->submodules.jacobian_ode = NULL;
    }

    assert((char *)raw_memory + sim_lifted_irk_integrator_calculate_args_size(dims, *submodules_) >= c_ptr);

    // Update submodules pointer
    *submodules_ = (void *) &(args->submodules);

    return (void *)args;
}



void *sim_lifted_irk_integrator_copy_args(sim_dims *dims, void *raw_memory, void *source_)
{
    sim_lifted_irk_integrator_args *source = (sim_lifted_irk_integrator_args *) source_;
    sim_lifted_irk_integrator_args *dest;

    sim_lifted_irk_integrator_submodules *submodules = &source->submodules;

    dest = sim_lifted_irk_integrator_assign_args(dims, (void **) &submodules, raw_memory);

    dest->interval = source->interval;
    dest->num_stages = source->num_stages;
    dest->num_steps = source->num_steps;
    dest->num_forw_sens = source->num_forw_sens;
    dest->sens_forw = source->sens_forw;
    dest->sens_adj = source->sens_adj;
    dest->sens_hess = source->sens_hess;
    dest->newton_iter = source->newton_iter;

    int ns = dims->num_stages;

    memcpy(dest->A_mat, source->A_mat, ns*ns*sizeof(double));
    memcpy(dest->c_vec, source->c_vec, ns*sizeof(double));
    memcpy(dest->b_vec, source->b_vec, ns*sizeof(double));

    dest->scheme->type = source->scheme->type;
    dest->scheme->single = source->scheme->single;
    dest->scheme->freeze = source->scheme->freeze;
    memcpy(dest->scheme->eig, source->scheme->eig, ns*sizeof(double));
    memcpy(dest->scheme->transf1, source->scheme->transf1, ns*ns*sizeof(double));
    memcpy(dest->scheme->transf2, source->scheme->transf2, ns*ns*sizeof(double));
    memcpy(dest->scheme->transf1_T, source->scheme->transf1_T, ns*ns*sizeof(double));
    memcpy(dest->scheme->transf2_T, source->scheme->transf2_T, ns*ns*sizeof(double));
    dest->scheme->low_tria = source->scheme->low_tria;

    extern external_function_dims sim_lifted_irk_forward_vde_dims;
    extern external_function_dims sim_lifted_irk_jacobian_ode_dims;
    source->submodules.forward_vde->copy_args(&sim_lifted_irk_forward_vde_dims, dest->forward_vde_args, source->forward_vde_args);
    source->submodules.jacobian_ode->copy_args(&sim_lifted_irk_jacobian_ode_dims, dest->jacobian_ode_args, source->jacobian_ode_args);

    return (void *)dest;
}



void sim_lifted_irk_integrator_initialize_default_args(sim_dims *dims, void *args_)
{
    sim_lifted_irk_integrator_args *args = (sim_lifted_irk_integrator_args *) args_;
    enum Newton_type_collocation type = exact;
    args->scheme->type = type;
    args->scheme->freeze = false;
    get_Gauss_nodes(args->num_stages, args->c_vec);
    create_Butcher_table(args->num_stages, args->c_vec, args->b_vec, args->A_mat);
    if (dims->num_stages <= 15 && (type == simplified_in || type == simplified_inis)) {
        read_Gauss_simplified(args->num_stages, args->scheme);
    } else {
        args->scheme->type = exact;
    }

    args->num_steps = 2;
    args->num_forw_sens = dims->nx + dims->nu;
    args->sens_forw = true;
    args->sens_adj = false;
    args->sens_hess = false;

    args->submodules.forward_vde->initialize_default_args(args->forward_vde_args);
    args->submodules.jacobian_ode->initialize_default_args(args->jacobian_ode_args);
}



int sim_lifted_irk_integrator_calculate_memory_size(sim_dims *dims, void *args_)
{
    sim_lifted_irk_integrator_args *args = (sim_lifted_irk_integrator_args *) args_;

    int nx = dims->nx;
    int nu = dims->nu;
    int num_steps = args->num_steps;
    int num_stages = args->num_stages;
    int nf = args->num_forw_sens;
    int num_sys = (int) ceil(num_stages/2.0);

    int size = sizeof(sim_lifted_irk_integrator_memory);

    size += num_steps * num_stages * sizeof(double *);  // jac_traj
    size += num_sys * sizeof(double *);  // sys_mat2
    size += num_sys * sizeof(int *);  // ipiv2
    size += num_sys * sizeof(double *);  // sys_sol2
    size += num_sys * sizeof(struct blasfeo_dmat *);  // str_mat2
    size += num_sys * sizeof(struct blasfeo_dmat *);  // str_sol2

    make_int_multiple_of(8, &size);

    size += nf * sizeof(double);  // grad_correction
    size += nx * num_stages * sizeof(double);  // grad_K

    size += num_steps * num_stages * nx * sizeof(double);  // K_traj
    size += num_steps * num_stages * nx * nf * sizeof(double);  // DK_traj
    size += num_steps * num_stages * nx * sizeof(double);  // mu_traj
    size += nx * sizeof(double);  // x
    size += nu * sizeof(double);  // u
    if (args->scheme->type == simplified_inis)
        size += num_steps * num_stages * nx * nf * sizeof(double);  // delta_DK_traj

    if (args->scheme->type == simplified_in || args->scheme->type == simplified_inis) {
        size += num_steps * num_stages * nx * sizeof(double);  // adj_traj
        size += num_steps * num_stages * nx * nx * sizeof(double);  // jac_traj
    }

    if (args->scheme->type == simplified_in || args->scheme->type == simplified_inis) {
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
    if (args->scheme->type == simplified_in || args->scheme->type == simplified_inis) {

        int dim_sys = 2 * nx;
        for (int i = 0; i < num_sys; i++) {
            if ((i + 1) == num_sys && num_sys != floor(num_stages / 2.0))  // odd number of stages
                dim_sys = nx;

#if defined(LA_HIGH_PERFORMANCE)
            size += blasfeo_memsize_dmat(dim_sys, dim_sys);  // str_mat2
            size += blasfeo_memsize_dmat(dim_sys, 1 + NF);  // str_sol2
#elif defined(LA_REFERENCE)
            size += blasfeo_memsize_diag_dmat(dim_sys, dim_sys);
#else  // LA_BLAS
            size += 0;
#endif  // LA_HIGH_PERFORMANCE
        }
    }
#endif  // !TRIPLE_LOOP

    size += args->submodules.forward_vde->calculate_memory_size(&sim_lifted_irk_forward_vde_dims, args->forward_vde_args);
    
    size += args->submodules.jacobian_ode->calculate_memory_size(&sim_lifted_irk_jacobian_ode_dims, args->jacobian_ode_args);

    size += 2 * 8;
    return size;
}



void *sim_lifted_irk_integrator_assign_memory(sim_dims *dims, void *args_, void *raw_memory)
{
    sim_lifted_irk_integrator_args *args = (sim_lifted_irk_integrator_args *) args_;

    int nx = dims->nx;
    int nu = dims->nu;
    int num_steps = args->num_steps;
    int num_stages = args->num_stages;
    int nf = args->num_forw_sens;
    int num_sys = (int) ceil(num_stages/2.0);

    char *c_ptr = raw_memory;

    sim_lifted_irk_integrator_memory *memory = raw_memory;
    c_ptr += sizeof(sim_lifted_irk_integrator_memory);

    align_char_to(8, &c_ptr);

    assign_double_ptrs(num_steps * num_stages, &memory->jac_traj, &c_ptr);
    assign_double_ptrs(num_sys, &memory->sys_mat2, &c_ptr);
    assign_int_ptrs(num_sys, &memory->ipiv2, &c_ptr);
    assign_double_ptrs(num_sys, &memory->sys_sol2, &c_ptr);
    assign_blasfeo_dmat_ptrs(num_sys, &memory->str_mat2, &c_ptr);
    assign_blasfeo_dmat_ptrs(num_sys, &memory->str_sol2, &c_ptr);

    align_char_to(8, &c_ptr);

    assign_double(nf, &memory->grad_correction, &c_ptr);
    assign_double(nx * num_stages, &memory->grad_K, &c_ptr);

    assign_double(num_steps * num_stages * nx, &memory->K_traj, &c_ptr);
    assign_double(num_steps * num_stages * nx * nf, &memory->DK_traj, &c_ptr);
    assign_double(num_steps * num_stages * nx, &memory->mu_traj, &c_ptr);
    assign_double(nx, &memory->x, &c_ptr);
    assign_double(nu, &memory->u, &c_ptr);

    if (args->scheme->type == simplified_inis)
        assign_double(num_steps * num_stages * nx * nf, &memory->delta_DK_traj, &c_ptr);

    if (args->scheme->type == simplified_in || args->scheme->type == simplified_inis) {
        assign_double(num_steps * num_stages * nx, &memory->adj_traj, &c_ptr);
        for (int i = 0; i < num_steps * num_stages; ++i)
            assign_double(nx * nx, &memory->jac_traj[i], &c_ptr);
    }

    if (args->scheme->type == simplified_in || args->scheme->type == simplified_inis) {
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
    if (args->scheme->type == simplified_in || args->scheme->type == simplified_inis) {

    int dim_sys = 2 * nx;
    for (int i = 0; i < num_sys; i++) {
        if ((i + 1) == num_sys && num_sys != floor(num_stages / 2.0))  // odd number of stages
            dim_sys = nx;

#if defined(LA_HIGH_PERFORMANCE)
        assign_blasfeo_dmat_mem(dim_sys, dim_sys, memory->str_mat2[i], &c_ptr);
        assign_blasfeo_dmat_mem(dim_sys, 1 + nf, memory->str_sol2[i], &c_ptr);
#elif defined(LA_REFERENCE)
        assign_blasfeo_dmat_mem(dim_sys, dim_sys, memory->str_mat2[i], memory->sys_mat2[i]);
        assign_blasfeo_dmat_mem(dim_sys, 1 + nf, memory->str_sol2[i], memory->sys_sol2[i]);
        d_cast_diag_mat2strmat((double *) c_ptr, memory->str_mat2[i]);
        c_ptr += dim_sys * sizeof(double);
#else  // LA_BLAS
        assign_blasfeo_dmat_mem(dim_sys, dim_sys, memory->str_mat2[i], memory->sys_mat2[i]);
        assign_blasfeo_dmat_mem(dim_sys, 1 + nf, memory->str_sol2[i], memory->sys_sol2[i]);
#endif  // LA_HIGH_PERFORMANCE
        blasfeo_dgesc(dim_sys, dim_sys, 0.0, memory->str_mat2[i], 0, 0);
        blasfeo_dgesc(dim_sys, 1 + nf, 0.0, memory->str_sol2[i], 0, 0);
        }
    }
#endif  // !TRIPLE_LOOP

    memory->forward_vde_mem = args->submodules.forward_vde->assign_memory(&sim_lifted_irk_forward_vde_dims, args->forward_vde_args, c_ptr);
    c_ptr += args->submodules.forward_vde->calculate_memory_size(&sim_lifted_irk_forward_vde_dims, args->forward_vde_args);

    memory->jacobian_ode_mem = args->submodules.jacobian_ode->assign_memory(&sim_lifted_irk_jacobian_ode_dims, args->jacobian_ode_args, c_ptr);
    c_ptr += args->submodules.jacobian_ode->calculate_memory_size(&sim_lifted_irk_jacobian_ode_dims, args->jacobian_ode_args);

    assert((char*)raw_memory + sim_lifted_irk_integrator_calculate_memory_size(dims, args) >= c_ptr);

    // initialize
    for (int i = 0; i < num_steps * num_stages * nx; ++i)
        memory->K_traj[i] = 0.0;
    for (int i = 0; i < num_steps * num_stages * nx * nf; ++i)
        memory->DK_traj[i] = 0.0;
    for (int i = 0; i < num_steps * num_stages * nx; ++i)
        memory->mu_traj[i] = 0.0;

    if (args->scheme->type == simplified_inis)
        for (int i = 0; i < num_steps * num_stages * nx * nf; ++i)
            memory->delta_DK_traj[i] = 0.0;
    if (args->scheme->type == simplified_in || args->scheme->type == simplified_inis) {
        for (int i = 0; i < num_steps * num_stages * nx; ++i)
            memory->adj_traj[i] = 0.0;
        for (int i = 0; i < num_steps * num_stages; ++i)
            for (int j = 0; j < nx * nx; ++j)
                memory->jac_traj[i][j] = 0.0;
    }

    return memory;
}



int sim_lifted_irk_integrator_calculate_workspace_size(sim_dims *dims, void *args_)
{
    int nx = dims->nx;
    int nu = dims->nu;
    int np = dims->np;
    sim_lifted_irk_integrator_args *args = (sim_lifted_irk_integrator_args *) args_;
    int num_stages = args->num_stages;
    int NF = args->num_forw_sens;
    //    int num_sys = ceil(num_stages/2.0);
    int dim_sys = num_stages * nx;
    if (args->scheme->type == simplified_in ||
        args->scheme->type == simplified_inis) {
        dim_sys = nx;
        if (num_stages > 1) dim_sys = 2 * nx;
    }

    int size = sizeof(sim_lifted_irk_integrator_workspace);
    size += (nx * (1 + NF) + nu + np + 1) * sizeof(double);  // rhs_in
    size += (nx * (1 + NF)) * sizeof(double);           // out_tmp
    if (args->scheme->type == exact) {
        size += (dim_sys) * sizeof(int);             // ipiv
        size += (dim_sys * dim_sys) * sizeof(double);  // sys_mat
    }
    size += ((num_stages * nx) * (1 + NF)) * sizeof(double);  // sys_sol
    size += (num_stages) * sizeof(double *);                  // VDE_tmp
    size += (num_stages * nx * (1 + NF)) * sizeof(double);    // VDE_tmp[...]
    size += (nx * (nx + 1)) * sizeof(double);                 // jac_tmp

    if (args->scheme->type == simplified_in ||
        args->scheme->type == simplified_inis) {
        size +=
            ((num_stages * nx) * (1 + NF)) * sizeof(double);  // sys_sol_trans
        size += (num_stages * num_stages) * sizeof(double);   // trans
    }

    if (args->scheme->type == simplified_in ||
        args->scheme->type == simplified_inis) {
        size += (nx) * sizeof(double);  // out_adj_tmp
    }

#if !TRIPLE_LOOP
    size += sizeof(struct blasfeo_dmat);  // str_mat
    size += sizeof(struct blasfeo_dmat);  // str_sol

    int size_strmat = 0;
#if defined(LA_HIGH_PERFORMANCE)
    // matrices in matrix struct format:
    size_strmat += blasfeo_memsize_dmat(dim_sys, dim_sys);
    size_strmat += blasfeo_memsize_dmat(dim_sys, 1 + NF);

#elif defined(LA_REFERENCE)
    // allocate new memory only for the diagonal
    size_strmat += blasfeo_memsize_diag_dmat(dim_sys, dim_sys);

#endif  // LA_HIGH_PERFORMANCE
    size += size_strmat;
#endif  // !TRIPLE_LOOP

    size += args->submodules.forward_vde->calculate_workspace_size(&sim_lifted_irk_forward_vde_dims, args->forward_vde_args);
    
    size += args->submodules.jacobian_ode->calculate_workspace_size(&sim_lifted_irk_jacobian_ode_dims, args->jacobian_ode_args);

    return size;
}



static void *cast_workspace(sim_dims *dims, void *args_, void *raw_memory)
{
    sim_lifted_irk_integrator_args *args = (sim_lifted_irk_integrator_args *) args_;
    
    int nx = dims->nx;
    int nu = dims->nu;
    int np = dims->np;

    int num_stages = args->num_stages;
    int NF = args->num_forw_sens;
    //    int num_sys = ceil(num_stages/2.0);
    int dim_sys = num_stages * nx;
    if (args->scheme->type == simplified_in ||
        args->scheme->type == simplified_inis) {
        dim_sys = nx;
        if (num_stages > 1) dim_sys = 2 * nx;
    }

    char *c_ptr = (char *)raw_memory;

    sim_lifted_irk_integrator_workspace *work = (sim_lifted_irk_integrator_workspace *) c_ptr;
    c_ptr += sizeof(sim_lifted_irk_integrator_workspace);

    work->rhs_in = (double *)c_ptr;
    c_ptr += (nx * (1 + NF) + nu + np + 1) * sizeof(double);  // rhs_in
    work->out_tmp = (double *)c_ptr;
    c_ptr += (nx * (1 + NF)) * sizeof(double);  // out_tmp
    if (args->scheme->type == exact) {
        work->ipiv = (int *)c_ptr;
        c_ptr += (dim_sys) * sizeof(int);  // ipiv
        work->sys_mat = (double *)c_ptr;
        c_ptr += (dim_sys * dim_sys) * sizeof(double);  // sys_mat
    }
    work->sys_sol = (double *)c_ptr;
    c_ptr += ((num_stages * nx) * (1 + NF)) * sizeof(double);  // sys_sol
    work->VDE_tmp = (double **)c_ptr;
    c_ptr += (num_stages) * sizeof(double *);  // VDE_tmp
    for (int i = 0; i < num_stages; i++) {
        work->VDE_tmp[i] = (double *)c_ptr;
        c_ptr += (nx * (1 + NF)) * sizeof(double);  // VDE_tmp[i]
    }
    work->jac_tmp = (double *)c_ptr;
    c_ptr += (nx * (nx + 1)) * sizeof(double);  // jac_tmp

    if (args->scheme->type == simplified_in ||
        args->scheme->type == simplified_inis) {
        work->sys_sol_trans = (double *)c_ptr;
        c_ptr +=
            ((num_stages * nx) * (1 + NF)) * sizeof(double);  // sys_sol_trans
        work->trans = (double *)c_ptr;
        c_ptr += (num_stages * num_stages) * sizeof(double);  // trans
    }

    if (args->scheme->type == simplified_in ||
        args->scheme->type == simplified_inis) {
        work->out_adj_tmp = (double *)c_ptr;
        c_ptr += (nx) * sizeof(double);  // out_adj_tmp
    }

#if !TRIPLE_LOOP
    work->str_mat = (struct blasfeo_dmat *)c_ptr;
    c_ptr += sizeof(struct blasfeo_dmat);
    work->str_sol = (struct blasfeo_dmat *)c_ptr;
    c_ptr += sizeof(struct blasfeo_dmat);
#if defined(LA_HIGH_PERFORMANCE)
    // matrices in matrix struct format:
    int size_strmat = 0;
    size_strmat += blasfeo_memsize_dmat(dim_sys, dim_sys);
    size_strmat += blasfeo_memsize_dmat(dim_sys, 1 + NF);

    blasfeo_create_dmat(dim_sys, dim_sys, work->str_mat, c_ptr);
    c_ptr += work->str_mat->memory_size;
    blasfeo_create_dmat(dim_sys, 1 + NF, work->str_sol, c_ptr);
    c_ptr += work->str_sol->memory_size;

#elif defined(LA_REFERENCE)

    //  pointer to column-major matrix
    blasfeo_create_dmat(dim_sys, dim_sys, work->str_mat, work->sys_mat);
    blasfeo_create_dmat(dim_sys, 1 + NF, work->str_sol, work->sys_sol);

    d_cast_diag_mat2strmat((double *)c_ptr, work->str_mat);
    c_ptr += blasfeo_memsize_diag_dmat(dim_sys, dim_sys);

#else  // LA_BLAS

    // not allocate new memory: point to column-major matrix
    blasfeo_create_dmat(dim_sys, dim_sys, work->str_mat, work->sys_mat);
    blasfeo_create_dmat(dim_sys, 1 + NF, work->str_sol, work->sys_sol);

#endif  // LA_HIGH_PERFORMANCE
#endif  // !TRIPLE_LOOP

    work->forward_vde_work = (void *) c_ptr;
    c_ptr += args->submodules.forward_vde->calculate_workspace_size(&sim_lifted_irk_forward_vde_dims, args->forward_vde_args);

    work->jacobian_ode_work = (void *) c_ptr;
    c_ptr += args->submodules.jacobian_ode->calculate_workspace_size(&sim_lifted_irk_jacobian_ode_dims, args->jacobian_ode_args);

    assert((char*)work + sim_lifted_irk_integrator_calculate_workspace_size(dims, args_) >= c_ptr);

    return (void *)work;
}



static void compute_forward_vde(const int nx, const int nu, const int np, double *in, double *out, 
                         sim_lifted_irk_integrator_args *args,
                         sim_lifted_irk_integrator_memory *mem,
                         sim_lifted_irk_integrator_workspace *work) 
{
    double *x = in;
    double *Sx = in + nx;
    double *Su = in + nx + nx * nx;
    double *u = in + nx + nx * (nx + nu);
    double *p = in + nx + nx * (nx + nu) + nu;

    double *x_out = out;
    double *Sx_out = out + nx;
    double *Su_out = out + nx + nx * nx;

    double *fw_in_inputs[FW_VDE_NUMIN];
    fw_in_inputs[0] = x;
    fw_in_inputs[1] = Sx;
    fw_in_inputs[2] = Su;
    fw_in_inputs[3] = u;
    fw_in_inputs[4] = p;

    bool fw_in_compute_output[FW_VDE_NUMOUT] = {true, true, true};

    double *fw_out_outputs[FW_VDE_NUMOUT];
    fw_out_outputs[0] = x_out;
    fw_out_outputs[1] = Sx_out;
    fw_out_outputs[2] = Su_out;

    external_function_in fw_vde_in;
    fw_vde_in.inputs = fw_in_inputs;
    fw_vde_in.compute_output = fw_in_compute_output;

    external_function_out fw_vde_out;
    fw_vde_out.outputs = fw_out_outputs;

    args->submodules.forward_vde->fun(&fw_vde_in, &fw_vde_out, args->forward_vde_args, mem->forward_vde_mem, work->forward_vde_work);
}



static void compute_jacobian_ode(const int nx, const int nu, const int np, double *in, double *out, 
                         sim_lifted_irk_integrator_args *args,
                         sim_lifted_irk_integrator_memory *mem,
                         sim_lifted_irk_integrator_workspace *work) 
{
    double *x = in;
    double *u = in + nx;
    double *p = in + nx + nu;

    double *x_out = out;
    double *jac_out = out + nx;

    double *jac_in_inputs[JAC_ODE_NUMIN];
    jac_in_inputs[0] = x;
    jac_in_inputs[1] = u;
    jac_in_inputs[2] = p;

    bool jac_in_compute_output[JAC_ODE_NUMOUT] = {true, true};

    double *jac_out_outputs[JAC_ODE_NUMOUT];
    jac_out_outputs[0] = x_out;
    jac_out_outputs[1] = jac_out;

    external_function_in jac_ode_in;
    jac_ode_in.inputs = jac_in_inputs;
    jac_ode_in.compute_output = jac_in_compute_output;

    external_function_out jac_ode_out;
    jac_ode_out.outputs = jac_out_outputs;

    args->submodules.forward_vde->fun(&jac_ode_in, &jac_ode_out, args->jacobian_ode_args, mem->jacobian_ode_mem, work->jacobian_ode_work);
}



#if TRIPLE_LOOP

#if CODE_GENERATION
#define DIM 6      // num_stages*NX
#define DIM_RHS 9  // NX+NU
#endif



double LU_system_ACADO(double *const A, int *const perm, int dim)
{
    double det;
    double swap;
    double valueMax;
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



double solve_system_ACADO(double *const A, double *const b, int *const perm,
                          int dim, int dim2)
{
    int i, j, k;
    int index1;

#if !CODE_GENERATION
    int DIM = dim;
    int DIM_RHS = dim2;
#else
    dim += 0;
    dim2 += 0;
#endif
    double *bPerm;
    bPerm = (double *) calloc(DIM*DIM_RHS, sizeof(double));
    double tmp_var;

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



double solve_system_trans_ACADO(double *const A, double *const b,
                                int *const perm, int dim, int dim2)
{
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
    double *bPerm;
    bPerm = (double *) calloc(DIM*DIM_RHS, sizeof(double));
    double tmp_var;

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

void transform_mat(double *mat, double *trans, double *mat_trans,
                   const int stages, const int n, const int m)
{
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

void transform_vec(double *mat, double *trans, double *mat_trans,
                   const int stages, const int n)
{
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




void construct_subsystems(double *mat, double **mat2, const int stages,
                          const int n, const int m)
{
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



void destruct_subsystems(double *mat, double **mat2, const int stages,
                         const int n, const int m)
{
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



void form_linear_system_matrix(int istep, const sim_in *in, void *args_,
                               sim_lifted_irk_integrator_memory *mem,
                               sim_lifted_irk_integrator_workspace *work, double *sys_mat,
                               double **sys_mat2, double timing_ad)
{
    int nx = in->nx;
    int nu = in->nu;
    int np = in->np;
    double H_INT = in->step;
    sim_lifted_irk_integrator_args *args = (sim_lifted_irk_integrator_args *) args_;
    int num_stages = args->num_stages;
    double *A_mat = args->A_mat;
    double *c_vec = args->c_vec;

    double tmp_eig, tmp_eig2;
    acados_timer timer_ad;
    double *out_tmp = work->out_tmp;
    double *rhs_in = work->rhs_in;
    double *jac_tmp = work->jac_tmp;

    double **jac_traj = mem->jac_traj;
    double *K_traj = mem->K_traj;

    int i, j, s1, s2;
    if ((args->scheme->type == simplified_in ||
         args->scheme->type == simplified_inis) &&
        istep == 0 && !args->scheme->freeze) {
        int idx = 0;
        for (s1 = 0; s1 < num_stages; s1++) {
            if ((s1 + 1) == num_stages) {  // real eigenvalue
                for (i = 0; i < nx * nx; i++) sys_mat2[idx][i] = 0.0;
                tmp_eig = 1.0 / H_INT * args->scheme->eig[s1];
                for (i = 0; i < nx; i++) {
                    sys_mat2[idx][i * (nx + 1)] = tmp_eig;
                }
            } else {  // complex conjugate pair
                for (i = 0; i < 4 * nx * nx; i++) sys_mat2[idx][i] = 0.0;
                tmp_eig = 1.0 / H_INT * args->scheme->eig[s1];
                tmp_eig2 = 1.0 / H_INT * args->scheme->eig[s1 + 1];
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
    } else if (args->scheme->type == exact) {
        for (i = 0; i < num_stages * nx * num_stages * nx; i++)
            sys_mat[i] = 0.0;
        for (i = 0; i < num_stages * nx; i++)
            sys_mat[i * (num_stages * nx + 1)] = 1.0;  // identity
    }

    for (s1 = 0; s1 < num_stages; s1++) {
        //                if (args->scheme->type == exact || s1 == 0) {
        for (i = 0; i < nx; i++) {
            rhs_in[i] = out_tmp[i];
        }
        for (i = 0; i < nu; i++) rhs_in[nx + i] = in->u[i];
        for (i = 0; i < np; i++) rhs_in[nx + nu + i] = in->p[i];
        rhs_in[nx + nu] = (istep + c_vec[s1]) / args->num_steps;  // time
        for (s2 = 0; s2 < num_stages; s2++) {
            for (i = 0; i < nx; i++) {
                rhs_in[i] += H_INT * A_mat[s2 * num_stages + s1] *
                             K_traj[istep * num_stages * nx + s2 * nx + i];
            }
        }
        acados_tic(&timer_ad);

        compute_jacobian_ode(nx, nu, np, rhs_in, jac_tmp, args, mem, work);  // k evaluation 

        timing_ad += acados_toc(&timer_ad);
        //                }
        if (args->scheme->type == simplified_in ||
            args->scheme->type == simplified_inis) {
            for (i = 0; i < nx * nx; i++)
                jac_traj[istep * num_stages + s1][i] = jac_tmp[nx + i];
        }

        // put jac_tmp in sys_mat:
        if (args->scheme->type == exact) {
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
        if ((args->scheme->type == simplified_in ||
             args->scheme->type == simplified_inis) &&
            istep == 0 && !args->scheme->freeze) {
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



int sim_lifted_irk_integrator(sim_in *in, sim_out *out, void *args_, void *mem_, void *work_)
{
    int nx = in->nx;
    int nu = in->nu;
    int np = in->np;
    sim_lifted_irk_integrator_args *args = (sim_lifted_irk_integrator_args *) args_;
    int num_stages = args->num_stages;
    int dim_sys = num_stages * nx;
    int i, s1, s2, j, istep;

    sim_dims dims = {
        args->num_stages,
        in->nx,
        in->nu,
        in->np
    };

    sim_lifted_irk_integrator_memory *mem = (sim_lifted_irk_integrator_memory *) mem_;
    sim_lifted_irk_integrator_workspace *work = (sim_lifted_irk_integrator_workspace *) cast_workspace(&dims, args, work_);
    
    double H_INT = in->step;
    int NF = args->num_forw_sens;

    double *A_mat = args->A_mat;
    double *b_vec = args->b_vec;
    double *c_vec = args->c_vec;

    double **VDE_tmp = work->VDE_tmp;
    double *out_tmp = work->out_tmp;
    double *rhs_in = work->rhs_in;

    double **jac_traj = mem->jac_traj;

    double *K_traj = mem->K_traj;
    double *DK_traj = mem->DK_traj;
    double *delta_DK_traj = mem->delta_DK_traj;
    double *mu_traj = mem->mu_traj;
    double *adj_traj = mem->adj_traj;

    double *adj_tmp = work->out_adj_tmp;

    int *ipiv = work->ipiv;    // pivoting vector
    int **ipiv2 = mem->ipiv2;  // pivoting vector
    double *sys_mat = work->sys_mat;
    double **sys_mat2 = mem->sys_mat2;
    double *sys_sol = work->sys_sol;
    double **sys_sol2 = mem->sys_sol2;
    double *sys_sol_trans = work->sys_sol_trans;
#if !TRIPLE_LOOP
    struct blasfeo_dmat *str_mat = work->str_mat;
    struct blasfeo_dmat **str_mat2 = mem->str_mat2;

    struct blasfeo_dmat *str_sol = work->str_sol;
    struct blasfeo_dmat **str_sol2 = mem->str_sol2;
#endif  // !TRIPLE_LOOP

    acados_timer timer, timer_la, timer_ad;
    double timing_la = 0.0;
    double timing_ad = 0.0;

    assert(NF == nx + nu && "Not implemented yet for other num_forw_sens");

    acados_tic(&timer);
    for (i = 0; i < nx; i++) out_tmp[i] = in->x[i];
    for (i = 0; i < nx * NF; i++)
        out_tmp[nx + i] = in->S_forw[i];  // sensitivities

    for (i = 0; i < nu; i++) rhs_in[nx * (1 + NF) + i] = in->u[i];
    for (i = 0; i < np; i++) rhs_in[nx * (1 + NF) + nu + i] = in->p[i];

    if (args->scheme->type == simplified_in ||
        args->scheme->type == simplified_inis) {
        for (i = 0; i < nx; i++) adj_tmp[i] = in->S_adj[i];
    }

    // Newton step of the collocation variables with respect to the inputs:
    if (NF == (nx + nu) && args->sens_adj) {
        for (istep = args->num_steps - 1; istep > -1;
             istep--) {  // ADJOINT update
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
            if (args->scheme->type == simplified_inis &&
                !args->scheme->freeze) {
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
            if (args->scheme->type == simplified_in ||
                args->scheme->type == simplified_inis) {
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
                if (args->scheme->type == simplified_in ||
                    args->scheme->type == simplified_inis) {
                    // apply the transf1 operation:
                    for (s1 = 0; s1 < num_stages; s1++) {
                        for (s2 = 0; s2 < num_stages; s2++) {
                            work->trans[s2 * num_stages + s1] =
                                1.0 / H_INT *
                                args->scheme->transf1_T[s2 * num_stages + s1];
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
                    if (args->scheme->type == simplified_in ||
                        args->scheme->type == simplified_inis) {
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
                if (args->scheme->type == simplified_in ||
                    args->scheme->type == simplified_inis) {
                    // construct sys_sol_trans from sys_sol2:
                    destruct_subsystems(sys_sol_trans, sys_sol2, num_stages, nx,
                                        1);

                    // apply the transf2 operation:
                    sys_sol = work->sys_sol;
                    transform_vec(sys_sol_trans, args->scheme->transf2_T,
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

    if (args->scheme->type == simplified_in ||
        args->scheme->type == simplified_inis) {
        for (i = 0; i < NF; i++) out->grad[i] = 0.0;
    }

    for (istep = 0; istep < args->num_steps; istep++) {
        // form linear system matrix (explicit ODE case):

        form_linear_system_matrix(istep, in, args, mem, work, sys_mat, sys_mat2,
                                  timing_ad);

        int idx;
        if (args->scheme->type == exact ||
            (istep == 0 && !args->scheme->freeze)) {
            acados_tic(&timer_la);
            idx = 0;
            for (s1 = 0; s1 < num_stages; s1++) {
                // THIS LOOP IS PARALLELIZABLE BECAUSE OF DECOMPOSABLE LINEAR
                // SUBSYSTEMS
                if (args->scheme->type == simplified_in ||
                    args->scheme->type == simplified_inis) {
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
                blasfeo_pack_dmat(dim_sys, dim_sys, sys_mat, dim_sys, str_mat,
                                  0,
                                  0);  // mat2strmat
#endif  // LA_BLAS | LA_REFERENCE
                blasfeo_dgetrf_rowpivot(dim_sys, dim_sys, str_mat, 0, 0,
                                        str_mat, 0, 0,
                                        ipiv);  // Gauss elimination
// ---- BLASFEO: LU factorization ----
#endif  // TRIPLE_LOOP
            }
            timing_la += acados_toc(&timer_la);
        }

        if (args->scheme->type == simplified_in ||
            args->scheme->type == simplified_inis) {
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
                if (args->scheme->type == simplified_inis) {
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
            rhs_in[nx * (1 + NF) + nu + np] = ((double)istep + c_vec[s1]) /
                                              ((double)args->num_steps);  // time

            acados_tic(&timer_ad);

            compute_forward_vde(nx, nu, np, rhs_in, VDE_tmp[s1], args, mem, work); // k evaluation

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
            if (args->scheme->type == simplified_inis) {
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
        if (args->scheme->type == simplified_inis) {
            for (j = 0; j < NF; j++) {
                for (i = 0; i < num_stages * nx; i++) {
                    out->grad[j] += mu_traj[istep * num_stages * nx + i] *
                                    sys_sol[(j + 1) * num_stages * nx + i];
                }
            }
        }

        if (args->scheme->type == simplified_in ||
            args->scheme->type == simplified_inis) {
            // apply the transf1 operation:
            for (s1 = 0; s1 < num_stages; s1++) {
                for (s2 = 0; s2 < num_stages; s2++) {
                    work->trans[s2 * num_stages + s1] =
                        1.0 / H_INT *
                        args->scheme->transf1[s2 * num_stages + s1];
                }
            }
            if (!args->scheme->freeze) {
                transform_mat(sys_sol, work->trans, sys_sol_trans, num_stages,
                              nx, 1 + NF);
            } else {
                transform_mat(sys_sol, work->trans, sys_sol_trans, num_stages,
                              nx, 1);
            }

            // construct sys_sol2 from sys_sol_trans:
            if (!args->scheme->freeze) {
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
            if (args->scheme->type == simplified_in ||
                args->scheme->type == simplified_inis) {
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
            if (!args->scheme->freeze) {
                solve_system_ACADO(sys_mat, sys_sol, ipiv, dim_sys, 1 + NF);
            } else {
                solve_system_ACADO(sys_mat, sys_sol, ipiv, dim_sys, 1);
            }
#else  // TRIPLE_LOOP
#if defined(LA_HIGH_PERFORMANCE)
            // ---- BLASFEO: row transformations + backsolve ----
            blasfeo_pack_dmat(dim_sys, 1 + NF, sys_sol, dim_sys, str_sol, 0,
                              0);                    // mat2strmat
            blasfeo_drowpe(dim_sys, ipiv, str_sol);  // row permutations
            blasfeo_dtrsm_llnu(dim_sys, 1 + NF, 1.0, str_mat, 0, 0, str_sol, 0,
                               0, str_sol, 0, 0);  // L backsolve
            blasfeo_dtrsm_lunn(dim_sys, 1 + NF, 1.0, str_mat, 0, 0, str_sol, 0,
                               0, str_sol, 0, 0);  // U backsolve
            blasfeo_unpack_dmat(
                dim_sys, 1 + NF, str_sol, 0, 0, sys_sol,
                dim_sys);  // strmat2mat
                           // BLASFEO: row transformations + backsolve
#else   // LA_BLAS | LA_REFERENCE
        // ---- BLASFEO: row transformations + backsolve ----
            blasfeo_drowpe(dim_sys, ipiv, str_sol);  // row permutations
            blasfeo_dtrsm_llnu(dim_sys, 1 + NF, 1.0, str_mat, 0, 0, str_sol, 0,
                               0, str_sol, 0, 0);  // L backsolve
            blasfeo_dtrsm_lunn(dim_sys, 1 + NF, 1.0, str_mat, 0, 0, str_sol, 0,
                               0, str_sol, 0,
                               0);  // U backsolve
                                    // BLASFEO: row transformations + backsolve
#endif  // LA_BLAFEO
#endif  // TRIPLE_LOOP
        }
        timing_la += acados_toc(&timer_la);
        if (args->scheme->type == simplified_in ||
            args->scheme->type == simplified_inis) {
            // construct sys_sol_trans from sys_sol2:
            if (!args->scheme->freeze) {
                destruct_subsystems(sys_sol_trans, sys_sol2, num_stages, nx,
                                    1 + NF);
            } else {
                destruct_subsystems(sys_sol_trans, sys_sol2, num_stages, nx, 1);
            }

            // apply the transf2 operation:
            sys_sol = work->sys_sol;
            if (!args->scheme->freeze) {
                transform_mat(sys_sol_trans, args->scheme->transf2, sys_sol,
                              num_stages, nx, 1 + NF);
            } else {
                transform_mat(sys_sol_trans, args->scheme->transf2, sys_sol,
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
                    if (args->scheme->type == simplified_inis &&
                        args->sens_adj && !args->scheme->freeze) {
                        delta_DK_traj[(istep * num_stages + s1) * nx * NF +
                                      j * nx + i] =
                            sys_sol[(j + 1) * num_stages * nx + s1 * nx + i];
                    } else if (args->scheme->type == simplified_inis &&
                               !args->scheme->freeze) {
                        DK_traj[(istep * num_stages + s1) * nx * NF + j * nx +
                                i] +=
                            sys_sol[(j + 1) * num_stages * nx + s1 * nx + i];
                    } else if (!args->scheme->freeze) {
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
        if (args->scheme->type == simplified_inis ||
            args->scheme->type == simplified_in) {
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
        if (args->scheme->type == simplified_in) {
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
