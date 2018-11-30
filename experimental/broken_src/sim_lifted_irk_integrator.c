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
/* OLD implementation removed from acados kernel, since there is a new implementation,
    which is easier to maintain - FreyJo */
#include "acados/sim/sim_lifted_irk_integrator.h"

#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "blasfeo/include/blasfeo_d_aux_ext_dep.h"
#include "blasfeo/include/blasfeo_d_blas.h"
#include "blasfeo/include/blasfeo_d_kernel.h"
#include "blasfeo/include/blasfeo_i_aux_ext_dep.h"

#include "acados/utils/mem.h"
#include "acados/utils/print.h"

/************************************************
 * dims
 ************************************************/

int sim_lifted_irk_dims_calculate_size()
{
    int size = sizeof(sim_lifted_irk_dims);

    return size;
}

void *sim_lifted_irk_dims_assign(void *config_, void *raw_memory)
{
    char *c_ptr = raw_memory;

    sim_lifted_irk_dims *dims = (sim_lifted_irk_dims *) c_ptr;
    c_ptr += sizeof(sim_lifted_irk_dims);

    dims->nx = 0;
    dims->nu = 0;
    dims->nz = 0;

    assert((char *) raw_memory + sim_lifted_irk_dims_calculate_size() >= c_ptr);

    return dims;
}

void sim_lifted_irk_set_nx(void *dims_, int nx)
{
    sim_lifted_irk_dims *dims = (sim_lifted_irk_dims *) dims_;
    dims->nx = nx;
}

void sim_lifted_irk_set_nu(void *dims_, int nu)
{
    sim_lifted_irk_dims *dims = (sim_lifted_irk_dims *) dims_;
    dims->nu = nu;
}

void sim_lifted_irk_set_nz(void *dims_, int nz)
{
    sim_lifted_irk_dims *dims = (sim_lifted_irk_dims *) dims_;
    dims->nz = nz;
}

void sim_lifted_irk_get_nx(void *dims_, int *nx)
{
    sim_lifted_irk_dims *dims = (sim_lifted_irk_dims *) dims_;
    *nx = dims->nx;
}

void sim_lifted_irk_get_nu(void *dims_, int *nu)
{
    sim_lifted_irk_dims *dims = (sim_lifted_irk_dims *) dims_;
    *nu = dims->nu;
}

void sim_lifted_irk_get_nz(void *dims_, int *nz)
{
    sim_lifted_irk_dims *dims = (sim_lifted_irk_dims *) dims_;
    *nz = dims->nz;
}

/************************************************
 * model
 ************************************************/

int sim_lifted_irk_model_calculate_size(void *config, void *dims)
{
    int size = 0;

    size += sizeof(lifted_irk_model);

    return size;
}

void *sim_lifted_irk_model_assign(void *config, void *dims, void *raw_memory)
{
    char *c_ptr = (char *) raw_memory;

    lifted_irk_model *data = (lifted_irk_model *) c_ptr;
    c_ptr += sizeof(lifted_irk_model);

    return data;
}

int sim_lifted_irk_model_set_function(void *model_, sim_function_t fun_type, void *fun)
{
    lifted_irk_model *model = model_;

    switch (fun_type)
    {
        case EXPL_ODE_JAC:
            model->expl_ode_jac = (external_function_generic *) fun;
            break;
        case EXPL_VDE_FOR:
            model->expl_vde_for = (external_function_generic *) fun;
            break;
        default:
            return ACADOS_FAILURE;
    }
    return ACADOS_SUCCESS;
}

/************************************************
 * opts
 ************************************************/

int sim_lifted_irk_opts_calculate_size(void *config_, void *dims)
{
    int ns_max = 15;

    int size = 0;

    size += sizeof(sim_rk_opts);

    size += ns_max * ns_max * sizeof(double);  // A_mat
    size += ns_max * sizeof(double);           // b_vec
    size += ns_max * sizeof(double);           // c_vec

    size += sizeof(Newton_scheme);

    make_int_multiple_of(8, &size);

    size += ns_max * sizeof(double);  // eig

    size += ns_max * ns_max * sizeof(double);  // transf1
    size += ns_max * ns_max * sizeof(double);  // transf2
    size += ns_max * ns_max * sizeof(double);  // transf1_T
    size += ns_max * ns_max * sizeof(double);  // transf2_T

    int tmp0 = gauss_nodes_work_calculate_size(ns_max);
    int tmp1 = butcher_table_work_calculate_size(ns_max);
    int tmp2 = gauss_simplified_work_calculate_size(ns_max);
    int work_size = tmp0 > tmp1 ? tmp0 : tmp1;
    work_size = tmp2 > work_size ? tmp2 : work_size;
    size += work_size;  // work

    make_int_multiple_of(8, &size);
    size += 2 * 8;

    return size;
}

// TODO(all): return pointer to sim_rk_opts instead
void *sim_lifted_irk_opts_assign(void *config_, void *dims, void *raw_memory)
{
    int ns_max = 15;

    char *c_ptr = (char *) raw_memory;

    sim_rk_opts *opts = (sim_rk_opts *) c_ptr;
    c_ptr += sizeof(sim_rk_opts);

    align_char_to(8, &c_ptr);

    assign_and_advance_double(ns_max * ns_max, &opts->A_mat, &c_ptr);
    assign_and_advance_double(ns_max, &opts->b_vec, &c_ptr);
    assign_and_advance_double(ns_max, &opts->c_vec, &c_ptr);

    opts->scheme = (Newton_scheme *) c_ptr;
    c_ptr += sizeof(Newton_scheme);

    align_char_to(8, &c_ptr);

    assign_and_advance_double(ns_max, &opts->scheme->eig, &c_ptr);

    assign_and_advance_double(ns_max * ns_max, &opts->scheme->transf1, &c_ptr);
    assign_and_advance_double(ns_max * ns_max, &opts->scheme->transf2, &c_ptr);
    assign_and_advance_double(ns_max * ns_max, &opts->scheme->transf1_T, &c_ptr);
    assign_and_advance_double(ns_max * ns_max, &opts->scheme->transf2_T, &c_ptr);

    // work
    int tmp0 = gauss_nodes_work_calculate_size(ns_max);
    int tmp1 = butcher_table_work_calculate_size(ns_max);
    int tmp2 = gauss_simplified_work_calculate_size(ns_max);
    int work_size = tmp0 > tmp1 ? tmp0 : tmp1;
    work_size = tmp2 > work_size ? tmp2 : work_size;
    opts->work = c_ptr;
    c_ptr += work_size;

    assert((char *) raw_memory + sim_lifted_irk_opts_calculate_size(config_, dims) >= c_ptr);

    return (void *) opts;
}

void sim_lifted_irk_opts_initialize_default(void *config_, void *dims_, void *opts_)
{
    sim_rk_opts *opts = opts_;
    sim_lifted_irk_dims *dims = (sim_lifted_irk_dims *) dims_;

    opts->ns = 3;  // GL 3
    int ns = opts->ns;

    assert(ns <= NS_MAX && "ns > NS_MAX!");

    // set tableau size
    opts->tableau_size = opts->ns;

    enum Newton_type_collocation type = exact;
    opts->scheme->type = type;
    opts->scheme->freeze = false;

    // gauss collocation nodes
    gauss_nodes(ns, opts->c_vec, opts->work);

    // butcher tableau
    butcher_table(ns, opts->c_vec, opts->b_vec, opts->A_mat, opts->work);

    // ???
    if (ns <= 15 && (type == simplified_in || type == simplified_inis))
    {
        gauss_simplified(ns, opts->scheme, opts->work);
    }
    else
    {
        opts->scheme->type = exact;
    }

    opts->num_steps = 2;
    opts->num_forw_sens = dims->nx + dims->nu;
    opts->sens_forw = true;
    opts->sens_adj = false;  // TODO(andrea): this effectively disables the expansion step in lifted
                             // integrators !!!!!
    opts->sens_hess = false;

    opts->output_z = false;
    opts->sens_algebraic = false;

    return;
}

void sim_lifted_irk_opts_update(void *config_, void *dims, void *opts_)
{
    sim_rk_opts *opts = opts_;

    int ns = opts->ns;

    assert(ns <= NS_MAX && "ns > NS_MAX!");

    // set tableau size
    opts->tableau_size = opts->ns;

    // gauss collocation nodes
    gauss_nodes(ns, opts->c_vec, opts->work);

    // butcher tableau
    butcher_table(ns, opts->c_vec, opts->b_vec, opts->A_mat, opts->work);

    // ???
    enum Newton_type_collocation type = exact;
    if (ns <= 15 && (type == simplified_in || type == simplified_inis))
    {
        gauss_simplified(ns, opts->scheme, opts->work);
    }
    else
    {
        opts->scheme->type = exact;
    }

    return;
}

/************************************************
 * memory
 ************************************************/

int sim_lifted_irk_memory_calculate_size(void *config_, void *dims_, void *opts_)
{
    sim_rk_opts *opts = opts_;
    sim_lifted_irk_dims *dims = (sim_lifted_irk_dims *) dims_;

    int ns = opts->ns;

    int nx = dims->nx;
    int nu = dims->nu;
    int num_steps = opts->num_steps;
    int nf = opts->num_forw_sens;
    int num_sys = (int) ceil(ns / 2.0);

    int size = sizeof(sim_lifted_irk_memory);

    size += num_steps * ns * sizeof(double *);        // jac_traj
    size += num_sys * sizeof(double *);               // sys_mat2
    size += num_sys * sizeof(int *);                  // ipiv2
    size += num_sys * sizeof(double *);               // sys_sol2
    size += num_sys * sizeof(struct blasfeo_dmat *);  // str_mat2
    size += num_sys * sizeof(struct blasfeo_dmat *);  // str_sol2

    make_int_multiple_of(8, &size);

    size += nf * sizeof(double);       // grad_correction
    size += nx * ns * sizeof(double);  // grad_K

    size += num_steps * ns * nx * sizeof(double);       // K_traj
    size += num_steps * ns * nx * nf * sizeof(double);  // DK_traj
    size += num_steps * ns * nx * sizeof(double);       // mu_traj
    size += nx * sizeof(double);                        // x
    size += nu * sizeof(double);                        // u
    if (opts->scheme->type == simplified_inis)
        size += num_steps * ns * nx * nf * sizeof(double);  // delta_DK_traj

    if (opts->scheme->type == simplified_in || opts->scheme->type == simplified_inis)
    {
        size += num_steps * ns * nx * sizeof(double);       // adj_traj
        size += num_steps * ns * nx * nx * sizeof(double);  // jac_traj
    }

    if (opts->scheme->type == simplified_in || opts->scheme->type == simplified_inis)
    {
        int dim_sys = 2 * nx;
        size += (num_sys - 1) * dim_sys * dim_sys * sizeof(double);   // sys_mat2
        size += (num_sys - 1) * dim_sys * sizeof(double);             // ipiv2
        size += (num_sys - 1) * dim_sys * (1 + nf) * sizeof(double);  // sys_sol2

        if (num_sys != floor(ns / 2.0))  // odd number of stages
            dim_sys = nx;
        size += dim_sys * dim_sys * sizeof(double);   // sys_mat2
        size += dim_sys * sizeof(int);                // ipiv2
        size += dim_sys * (1 + nf) * sizeof(double);  // sys_sol2
    }
    else
    {
        num_sys = 1;
    }

#if !TRIPLE_LOOP
    if (opts->scheme->type == simplified_in || opts->scheme->type == simplified_inis)
    {
        int dim_sys = 2 * nx;
        for (int i = 0; i < num_sys; i++)
        {
            if ((i + 1) == num_sys && num_sys != floor(ns / 2.0))  // odd number of stages
                dim_sys = nx;

#if defined(LA_HIGH_PERFORMANCE)
            size += blasfeo_memsize_dmat(dim_sys, dim_sys);  // str_mat2
            size += blasfeo_memsize_dmat(dim_sys, 1 + NF);   // str_sol2
#elif defined(LA_REFERENCE)
            size += blasfeo_memsize_diag_dmat(dim_sys, dim_sys);
#else   // LA_BLAS
            size += 0;
#endif  // LA_HIGH_PERFORMANCE
        }
    }
#endif  // !TRIPLE_LOOP
    size += 2 * 8;
    return size;
}

void *sim_lifted_irk_memory_assign(void *config_, void *dims_, void *opts_, void *raw_memory)
{
    sim_rk_opts *opts = opts_;
    sim_lifted_irk_dims *dims = (sim_lifted_irk_dims *) dims_;

    int ns = opts->ns;

    int nx = dims->nx;
    int nu = dims->nu;
    int num_steps = opts->num_steps;
    int nf = opts->num_forw_sens;
    int num_sys = (int) ceil(ns / 2.0);

    char *c_ptr = raw_memory;

    sim_lifted_irk_memory *memory = raw_memory;
    c_ptr += sizeof(sim_lifted_irk_memory);

    align_char_to(8, &c_ptr);

    assign_and_advance_double_ptrs(num_steps * ns, &memory->jac_traj, &c_ptr);
    assign_and_advance_double_ptrs(num_sys, &memory->sys_mat2, &c_ptr);
    assign_and_advance_int_ptrs(num_sys, &memory->ipiv2, &c_ptr);
    assign_and_advance_double_ptrs(num_sys, &memory->sys_sol2, &c_ptr);
    assign_and_advance_blasfeo_dmat_ptrs(num_sys, &memory->str_mat2, &c_ptr);
    assign_and_advance_blasfeo_dmat_ptrs(num_sys, &memory->str_sol2, &c_ptr);

    align_char_to(8, &c_ptr);

    assign_and_advance_double(nf, &memory->grad_correction, &c_ptr);
    assign_and_advance_double(nx * ns, &memory->grad_K, &c_ptr);

    assign_and_advance_double(num_steps * ns * nx, &memory->K_traj, &c_ptr);
    assign_and_advance_double(num_steps * ns * nx * nf, &memory->DK_traj, &c_ptr);
    assign_and_advance_double(num_steps * ns * nx, &memory->mu_traj, &c_ptr);
    assign_and_advance_double(nx, &memory->x, &c_ptr);
    assign_and_advance_double(nu, &memory->u, &c_ptr);

    if (opts->scheme->type == simplified_inis)
        assign_and_advance_double(num_steps * ns * nx * nf, &memory->delta_DK_traj, &c_ptr);

    if (opts->scheme->type == simplified_in || opts->scheme->type == simplified_inis)
    {
        assign_and_advance_double(num_steps * ns * nx, &memory->adj_traj, &c_ptr);
        for (int i = 0; i < num_steps * ns; ++i)
            assign_and_advance_double(nx * nx, &memory->jac_traj[i], &c_ptr);
    }

    if (opts->scheme->type == simplified_in || opts->scheme->type == simplified_inis)
    {
        int dim_sys = 2 * nx;
        for (int i = 0; i < num_sys; ++i)
        {
            assign_and_advance_double(dim_sys * dim_sys, &memory->sys_mat2[i], &c_ptr);
            assign_and_advance_double(dim_sys, &memory->sys_mat2[i], &c_ptr);
            assign_and_advance_double(dim_sys * (1 + nf), &memory->sys_mat2[i], &c_ptr);

            for (int j = 0; j < dim_sys * dim_sys; ++j) memory->sys_mat2[i][j] = 0.0;
            for (int j = 0; j < dim_sys; j++) memory->sys_mat2[i][j * (dim_sys + 1)] = 1.0;
        }
        if (num_sys != floor(ns / 2.0))  // odd number of stages
            dim_sys = nx;
        assign_and_advance_double(dim_sys * dim_sys, &memory->sys_mat2[num_sys], &c_ptr);
        assign_and_advance_double(dim_sys, &memory->sys_mat2[num_sys], &c_ptr);
        assign_and_advance_double(dim_sys * (1 + nf), &memory->sys_mat2[num_sys], &c_ptr);

        for (int j = 0; j < dim_sys * dim_sys; ++j) memory->sys_mat2[num_sys][j] = 0.0;
        for (int j = 0; j < dim_sys; j++) memory->sys_mat2[num_sys][j * (dim_sys + 1)] = 1.0;
    }
    else
    {
        num_sys = 1;
    }

#if !TRIPLE_LOOP
    if (opts->scheme->type == simplified_in || opts->scheme->type == simplified_inis)
    {
        int dim_sys = 2 * nx;
        for (int i = 0; i < num_sys; i++)
        {
            if ((i + 1) == num_sys && num_sys != floor(ns / 2.0))  // odd number of stages
                dim_sys = nx;

#if defined(LA_HIGH_PERFORMANCE)
            assign_and_advance_blasfeo_dmat_mem(dim_sys, dim_sys, memory->str_mat2[i], &c_ptr);
            assign_and_advance_blasfeo_dmat_mem(dim_sys, 1 + nf, memory->str_sol2[i], &c_ptr);
#elif defined(LA_REFERENCE)
            assign_and_advance_blasfeo_dmat_mem(dim_sys, dim_sys, memory->str_mat2[i],
                                                memory->sys_mat2[i]);
            assign_and_advance_blasfeo_dmat_mem(dim_sys, 1 + nf, memory->str_sol2[i],
                                                memory->sys_sol2[i]);
            d_cast_diag_mat2strmat((double *) c_ptr, memory->str_mat2[i]);
            c_ptr += dim_sys * sizeof(double);
#else   // LA_BLAS
            assign_and_advance_blasfeo_dmat_mem(dim_sys, dim_sys, memory->str_mat2[i],
                                                memory->sys_mat2[i]);
            assign_and_advance_blasfeo_dmat_mem(dim_sys, 1 + nf, memory->str_sol2[i],
                                                memory->sys_sol2[i]);
#endif  // LA_HIGH_PERFORMANCE
            blasfeo_dgesc(dim_sys, dim_sys, 0.0, memory->str_mat2[i], 0, 0);
            blasfeo_dgesc(dim_sys, 1 + nf, 0.0, memory->str_sol2[i], 0, 0);
        }
    }
#endif  // !TRIPLE_LOOP

    assert((char *) raw_memory + sim_lifted_irk_memory_calculate_size(config_, dims, opts) >=
           c_ptr);

    // initialize
    for (int i = 0; i < num_steps * ns * nx; ++i) memory->K_traj[i] = 0.0;
    for (int i = 0; i < num_steps * ns * nx * nf; ++i) memory->DK_traj[i] = 0.0;
    for (int i = 0; i < num_steps * ns * nx; ++i) memory->mu_traj[i] = 0.0;

    if (opts->scheme->type == simplified_inis)
    {
        for (int i = 0; i < num_steps * ns * nx * nf; ++i) memory->delta_DK_traj[i] = 0.0;
    }
    if (opts->scheme->type == simplified_in || opts->scheme->type == simplified_inis)
    {
        for (int i = 0; i < num_steps * ns * nx; ++i) memory->adj_traj[i] = 0.0;
        for (int i = 0; i < num_steps * ns; ++i)
            for (int j = 0; j < nx * nx; ++j) memory->jac_traj[i][j] = 0.0;
    }

    return memory;
}

/************************************************
 * workspace
 ************************************************/

int sim_lifted_irk_workspace_calculate_size(void *config_, void *dims_, void *opts_)
{
    sim_rk_opts *opts = opts_;
    sim_lifted_irk_dims *dims = (sim_lifted_irk_dims *) dims_;

    int ns = opts->ns;

    int nx = dims->nx;
    int nu = dims->nu;
    int NF = opts->num_forw_sens;
    //    int num_sys = ceil(ns/2.0);
    int dim_sys = ns * nx;
    if (opts->scheme->type == simplified_in || opts->scheme->type == simplified_inis)
    {
        dim_sys = nx;
        if (ns > 1) dim_sys = 2 * nx;
    }

    int size = sizeof(sim_lifted_irk_workspace);
    size += (nx * (1 + NF) + nu + 1) * sizeof(double);  // rhs_in
    size += (nx * (1 + NF)) * sizeof(double);           // out_tmp
    if (opts->scheme->type == exact)
    {
        size += (dim_sys) * sizeof(int);               // ipiv
        size += (dim_sys * dim_sys) * sizeof(double);  // sys_mat
    }
    size += ((ns * nx) * (1 + NF)) * sizeof(double);  // sys_sol
    size += (ns) * sizeof(double *);                  // VDE_tmp
    size += (ns * nx * (1 + NF)) * sizeof(double);    // VDE_tmp[...]
    size += (nx * (nx + 1)) * sizeof(double);         // jac_tmp

    if (opts->scheme->type == simplified_in || opts->scheme->type == simplified_inis)
    {
        size += ((ns * nx) * (1 + NF)) * sizeof(double);  // sys_sol_trans
        size += (ns * ns) * sizeof(double);               // trans
    }

    if (opts->scheme->type == simplified_in || opts->scheme->type == simplified_inis)
    {
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

    return size;
}

static void sim_lifted_irk_workspace_cast(void *config_, sim_lifted_irk_workspace *work,
                                          const sim_in *in, void *opts_)
{
    sim_rk_opts *opts = opts_;
    sim_lifted_irk_dims *dims = (sim_lifted_irk_dims *) in->dims;

    int ns = opts->ns;

    int nx = dims->nx;
    int nu = dims->nu;
    int NF = opts->num_forw_sens;
    //    int num_sys = ceil(ns/2.0);
    int dim_sys = ns * nx;
    if (opts->scheme->type == simplified_in || opts->scheme->type == simplified_inis)
    {
        dim_sys = nx;
        if (ns > 1) dim_sys = 2 * nx;
    }

    char *ptr = (char *) work;
    ptr += sizeof(sim_lifted_irk_workspace);
    work->rhs_in = (double *) ptr;
    ptr += (nx * (1 + NF) + nu + 1) * sizeof(double);  // rhs_in
    work->out_tmp = (double *) ptr;
    ptr += (nx * (1 + NF)) * sizeof(double);  // out_tmp
    if (opts->scheme->type == exact)
    {
        work->ipiv = (int *) ptr;
        ptr += (dim_sys) * sizeof(int);  // ipiv
        work->sys_mat = (double *) ptr;
        ptr += (dim_sys * dim_sys) * sizeof(double);  // sys_mat
    }
    work->sys_sol = (double *) ptr;
    ptr += ((ns * nx) * (1 + NF)) * sizeof(double);  // sys_sol
    work->VDE_tmp = (double **) ptr;
    ptr += (ns) * sizeof(double *);  // VDE_tmp
    for (int i = 0; i < ns; i++)
    {
        work->VDE_tmp[i] = (double *) ptr;
        ptr += (nx * (1 + NF)) * sizeof(double);  // VDE_tmp[i]
    }
    work->jac_tmp = (double *) ptr;
    ptr += (nx * (nx + 1)) * sizeof(double);  // jac_tmp

    if (opts->scheme->type == simplified_in || opts->scheme->type == simplified_inis)
    {
        work->sys_sol_trans = (double *) ptr;
        ptr += ((ns * nx) * (1 + NF)) * sizeof(double);  // sys_sol_trans
        work->trans = (double *) ptr;
        ptr += (ns * ns) * sizeof(double);  // trans
    }

    if (opts->scheme->type == simplified_in || opts->scheme->type == simplified_inis)
    {
        work->out_adj_tmp = (double *) ptr;
        ptr += (nx) * sizeof(double);  // out_adj_tmp
    }

#if !TRIPLE_LOOP
    work->str_mat = (struct blasfeo_dmat *) ptr;
    ptr += sizeof(struct blasfeo_dmat);
    work->str_sol = (struct blasfeo_dmat *) ptr;
    ptr += sizeof(struct blasfeo_dmat);
#if defined(LA_HIGH_PERFORMANCE)
    // matrices in matrix struct format:
    int size_strmat = 0;
    size_strmat += blasfeo_memsize_dmat(dim_sys, dim_sys);
    size_strmat += blasfeo_memsize_dmat(dim_sys, 1 + NF);

    blasfeo_create_dmat(dim_sys, dim_sys, work->str_mat, ptr);
    ptr += work->str_mat->memory_size;
    blasfeo_create_dmat(dim_sys, 1 + NF, work->str_sol, ptr);
    ptr += work->str_sol->memory_size;

#elif defined(LA_REFERENCE)

    //  pointer to column-major matrix
    blasfeo_create_dmat(dim_sys, dim_sys, work->str_mat, work->sys_mat);
    blasfeo_create_dmat(dim_sys, 1 + NF, work->str_sol, work->sys_sol);

    d_cast_diag_mat2strmat((double *) ptr, work->str_mat);
    ptr += blasfeo_memsize_diag_dmat(dim_sys, dim_sys);

#else  // LA_BLAS

    // not allocate new memory: point to column-major matrix
    blasfeo_create_dmat(dim_sys, dim_sys, work->str_mat, work->sys_mat);
    blasfeo_create_dmat(dim_sys, 1 + NF, work->str_sol, work->sys_sol);

#endif  // LA_HIGH_PERFORMANCE
#endif  // !TRIPLE_LOOP
}

#if TRIPLE_LOOP

#if CODE_GENERATION
#define DIM 6      // ns*NX
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

    for (i = 0; i < DIM; ++i)
    {
        perm[i] = i;
    }
    det = 1.0000000000000000e+00;
    for (i = 0; i < (DIM - 1); i++)
    {
        indexMax = i;
        valueMax = fabs(A[i * DIM + i]);
        for (j = (i + 1); j < DIM; j++)
        {
            swap = fabs(A[i * DIM + j]);
            if (swap > valueMax)
            {
                indexMax = j;
                valueMax = swap;
            }
        }
        if (indexMax > i)
        {
            for (k = 0; k < DIM; ++k)
            {
                swap = A[k * DIM + i];
                A[k * DIM + i] = A[k * DIM + indexMax];
                A[k * DIM + indexMax] = swap;
            }
            intSwap = perm[i];
            perm[i] = perm[indexMax];
            perm[indexMax] = intSwap;
        }
        //        det *= A[i*DIM+i];
        for (j = i + 1; j < DIM; j++)
        {
            A[i * DIM + j] = -A[i * DIM + j] / A[i * DIM + i];
            for (k = i + 1; k < DIM; k++)
            {
                A[k * DIM + j] += A[i * DIM + j] * A[k * DIM + i];
            }
        }
    }
    //    det *= A[DIM*DIM-1];
    //    det = fabs(det);
    return det;
}

/************************************************
 * functions
 ************************************************/

double solve_system_ACADO(double *const A, double *const b, int *const perm, int dim, int dim2)
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
    bPerm = (double *) calloc(DIM * DIM_RHS, sizeof(double));
    double tmp_var;

    for (i = 0; i < DIM; ++i)
    {
        index1 = perm[i];
        for (j = 0; j < DIM_RHS; ++j)
        {
            bPerm[j * DIM + i] = b[j * DIM + index1];
        }
    }
    for (j = 1; j < DIM; ++j)
    {
        for (i = 0; i < j; ++i)
        {
            tmp_var = A[i * DIM + j];
            for (k = 0; k < DIM_RHS; ++k)
            {
                bPerm[k * DIM + j] += tmp_var * bPerm[k * DIM + i];
            }
        }
    }
    for (i = DIM - 1; - 1 < i; --i)
    {
        for (j = DIM - 1; i < j; --j)
        {
            tmp_var = A[j * DIM + i];
            for (k = 0; k < DIM_RHS; ++k)
            {
                bPerm[k * DIM + i] -= tmp_var * bPerm[k * DIM + j];
            }
        }
        tmp_var = 1.0 / A[i * (DIM + 1)];
        for (k = 0; k < DIM_RHS; ++k)
        {
            bPerm[k * DIM + i] = tmp_var * bPerm[k * DIM + i];
        }
    }
    for (k = 0; k < DIM * DIM_RHS; ++k)
    {
        b[k] = bPerm[k];
    }
    free(bPerm);
    return 0;
}

double solve_system_trans_ACADO(double *const A, double *const b, int *const perm, int dim,
                                int dim2)
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
    bPerm = (double *) calloc(DIM * DIM_RHS, sizeof(double));
    double tmp_var;

    for (k = 0; k < DIM * DIM_RHS; ++k)
    {
        bPerm[k] = b[k];
    }
    for (i = 0; i < DIM; i++)
    {
        for (j = 0; j < i; j++)
        {
            tmp_var = A[i * DIM + j];
            for (k = 0; k < DIM_RHS; ++k)
            {
                bPerm[k * DIM + i] -= tmp_var * bPerm[k * DIM + j];
            }
        }
        tmp_var = 1.0 / A[i * (DIM + 1)];
        for (k = 0; k < DIM_RHS; ++k)
        {
            bPerm[k * DIM + i] = tmp_var * bPerm[k * DIM + i];
        }
    }
    for (j = DIM - 1; j > -1; --j)
    {
        for (i = DIM - 1; i > j; --i)
        {
            tmp_var = A[j * DIM + i];
            for (k = 0; k < DIM_RHS; ++k)
            {
                bPerm[k * DIM + j] += tmp_var * bPerm[k * DIM + i];
            }
        }
    }
    for (i = 0; i < DIM; ++i)
    {
        index1 = perm[i];
        for (j = 0; j < DIM_RHS; ++j)
        {
            b[j * DIM + index1] = bPerm[j * DIM + i];
        }
    }
    free(bPerm);
    return 0;
}

#endif

void transform_mat(double *mat, double *trans, double *mat_trans, const int stages, const int n,
                   const int m)
{
    //    print_matrix_name("stdout", "trans", trans, stages, stages);
    for (int i = 0; i < stages * n * m; i++) mat_trans[i] = 0.0;
    for (int j = 0; j < m; j++)
    {
        for (int s2 = 0; s2 < stages; s2++)
        {
            for (int s1 = 0; s1 < stages; s1++)
            {
                if (trans[s2 * stages + s1] != 0)
                {
                    for (int i = 0; i < n; i++)
                    {
                        mat_trans[j * (stages * n) + s1 * n + i] +=
                            trans[s2 * stages + s1] * mat[j * (stages * n) + s2 * n + i];
                    }
                }
            }
        }
    }
}

void transform_vec(double *mat, double *trans, double *mat_trans, const int stages, const int n)
{
    for (int i = 0; i < stages * n; i++) mat_trans[i] = 0.0;
    for (int s2 = 0; s2 < stages; s2++)
    {
        for (int s1 = 0; s1 < stages; s1++)
        {
            for (int i = 0; i < n; i++)
            {
                mat_trans[s1 * n + i] += trans[s2 * stages + s1] * mat[s2 * n + i];
            }
        }
    }
}

void construct_subsystems(double *mat, double **mat2, const int stages, const int n, const int m)
{
    int idx = 0;
    for (int s1 = 0; s1 < stages; s1++)
    {
        if ((s1 + 1) < stages)
        {  // complex conjugate pair of eigenvalues
            for (int j = 0; j < m; j++)
            {
                for (int i = 0; i < n; i++)
                {
                    mat2[idx][j * 2 * n + i] = mat[j * (stages * n) + s1 * n + i];
                }
            }
            s1++;
            for (int j = 0; j < m; j++)
            {
                for (int i = 0; i < n; i++)
                {
                    mat2[idx][j * 2 * n + n + i] = mat[j * (stages * n) + s1 * n + i];
                }
            }
        }
        else
        {  // real eigenvalue
            for (int j = 0; j < m; j++)
            {
                for (int i = 0; i < n; i++)
                {
                    mat2[idx][j * n + i] = mat[j * (stages * n) + s1 * n + i];
                }
            }
        }
        idx++;
    }
}

void destruct_subsystems(double *mat, double **mat2, const int stages, const int n, const int m)
{
    int idx = 0;
    for (int s1 = 0; s1 < stages; s1++)
    {
        if ((s1 + 1) < stages)
        {  // complex conjugate pair of eigenvalues
            for (int j = 0; j < m; j++)
            {
                for (int i = 0; i < n; i++)
                {
                    mat[j * (stages * n) + s1 * n + i] = mat2[idx][j * 2 * n + i];
                }
            }
            s1++;
            for (int j = 0; j < m; j++)
            {
                for (int i = 0; i < n; i++)
                {
                    mat[j * (stages * n) + s1 * n + i] = mat2[idx][j * 2 * n + n + i];
                }
            }
        }
        else
        {  // real eigenvalue
            for (int j = 0; j < m; j++)
            {
                for (int i = 0; i < n; i++)
                {
                    mat[j * (stages * n) + s1 * n + i] = mat2[idx][j * n + i];
                }
            }
        }
        idx++;
    }
}

static void form_linear_system_matrix(void *config_, int istep, const sim_in *in, void *opts_,
                                      sim_lifted_irk_memory *mem, sim_lifted_irk_workspace *work,
                                      double *sys_mat, double **sys_mat2, double timing_ad)
{
    sim_rk_opts *opts = opts_;
    int ns = opts->ns;

    sim_lifted_irk_dims *dims = (sim_lifted_irk_dims *) in->dims;

    int nx = dims->nx;
    int nu = dims->nu;

    double H_INT = in->T / opts->num_steps;
    double *A_mat = opts->A_mat;
    double *c_vec = opts->c_vec;

    double tmp_eig, tmp_eig2;
    acados_timer timer_ad;
    double *out_tmp = work->out_tmp;
    double *rhs_in = work->rhs_in;
    double *jac_tmp = work->jac_tmp;

    double **jac_traj = mem->jac_traj;
    double *K_traj = mem->K_traj;

    ext_fun_arg_t ext_fun_type_in[5];
    void *ext_fun_in[5];  // XXX large enough ?
    ext_fun_arg_t ext_fun_type_out[5];
    void *ext_fun_out[5];  // XXX large enough ?

    lifted_irk_model *model = in->model;

    int i, j, s1, s2;
    if ((opts->scheme->type == simplified_in || opts->scheme->type == simplified_inis) &&
        istep == 0 && !opts->scheme->freeze)
    {
        int idx = 0;
        for (s1 = 0; s1 < ns; s1++)
        {
            if ((s1 + 1) == ns)
            {  // real eigenvalue
                for (i = 0; i < nx * nx; i++) sys_mat2[idx][i] = 0.0;
                tmp_eig = 1.0 / H_INT * opts->scheme->eig[s1];
                for (i = 0; i < nx; i++)
                {
                    sys_mat2[idx][i * (nx + 1)] = tmp_eig;
                }
            }
            else
            {  // complex conjugate pair
                for (i = 0; i < 4 * nx * nx; i++) sys_mat2[idx][i] = 0.0;
                tmp_eig = 1.0 / H_INT * opts->scheme->eig[s1];
                tmp_eig2 = 1.0 / H_INT * opts->scheme->eig[s1 + 1];
                for (i = 0; i < nx; i++)
                {
                    sys_mat2[idx][i * (2 * nx + 1)] = tmp_eig;
                    sys_mat2[idx][(nx + i) * (2 * nx + 1)] = tmp_eig;

                    sys_mat2[idx][i * 2 * nx + (nx + i)] = tmp_eig2;
                    sys_mat2[idx][(nx + i) * 2 * nx + i] = -tmp_eig2;
                }
                s1++;  // skip the complex conjugate eigenvalue
            }
            idx++;
        }
    }
    else if (opts->scheme->type == exact)
    {
        for (i = 0; i < ns * nx * ns * nx; i++) sys_mat[i] = 0.0;
        for (i = 0; i < ns * nx; i++) sys_mat[i * (ns * nx + 1)] = 1.0;  // identity
    }

    for (s1 = 0; s1 < ns; s1++)
    {
        //                if (opts->scheme->type == exact || s1 == 0) {
        for (i = 0; i < nx; i++)
        {
            rhs_in[i] = out_tmp[i];
        }
        for (i = 0; i < nu; i++) rhs_in[nx + i] = in->u[i];
        rhs_in[nx + nu] = (istep + c_vec[s1]) / opts->num_steps;  // time
        for (s2 = 0; s2 < ns; s2++)
        {
            for (i = 0; i < nx; i++)
            {
                rhs_in[i] += H_INT * A_mat[s2 * ns + s1] * K_traj[istep * ns * nx + s2 * nx + i];
            }
        }
        acados_tic(&timer_ad);

        ext_fun_type_in[0] = COLMAJ;
        ext_fun_in[0] = rhs_in + 0;  // x: nx
        ext_fun_type_in[1] = COLMAJ;
        ext_fun_in[1] = rhs_in + nx;  // u: nu

        ext_fun_type_out[0] = COLMAJ;
        ext_fun_out[0] = jac_tmp + 0;  // fun: nx
        ext_fun_type_out[1] = COLMAJ;
        ext_fun_out[1] = jac_tmp + nx;  // jac_x: nx*nx

        model->expl_ode_jac->evaluate(model->expl_ode_jac, ext_fun_type_in, ext_fun_in,
                                      ext_fun_type_out, ext_fun_out);  // k evaluation

        timing_ad += acados_toc(&timer_ad);
        //                }
        if (opts->scheme->type == simplified_in || opts->scheme->type == simplified_inis)
        {
            for (i = 0; i < nx * nx; i++) jac_traj[istep * ns + s1][i] = jac_tmp[nx + i];
        }

        // put jac_tmp in sys_mat:
        if (opts->scheme->type == exact)
        {
            for (s2 = 0; s2 < ns; s2++)
            {
                for (j = 0; j < nx; j++)
                {
                    for (i = 0; i < nx; i++)
                    {
                        sys_mat[(s2 * nx + j) * ns * nx + s1 * nx + i] -=
                            H_INT * A_mat[s2 * ns + s1] * jac_tmp[nx + j * nx + i];
                    }
                }
            }
        }
    }

    int idx = 0;
    for (s1 = 0; s1 < ns; s1++)
    {
        // put jac_traj[0] in sys_mat:
        if ((opts->scheme->type == simplified_in || opts->scheme->type == simplified_inis) &&
            istep == 0 && !opts->scheme->freeze)
        {
            if ((s1 + 1) == ns)
            {  // real eigenvalue
                for (j = 0; j < nx; j++)
                {
                    for (i = 0; i < nx; i++)
                    {
                        sys_mat2[idx][j * nx + i] -= jac_traj[0][j * nx + i];
                    }
                }
            }
            else
            {  // complex conjugate pair
                for (j = 0; j < nx; j++)
                {
                    for (i = 0; i < nx; i++)
                    {
                        sys_mat2[idx][j * 2 * nx + i] -= jac_traj[0][j * nx + i];
                        sys_mat2[idx][(nx + j) * 2 * nx + nx + i] -= jac_traj[0][j * nx + i];
                    }
                }
                s1++;  // skip the complex conjugate eigenvalue
            }
            idx++;
        }
    }
}

int sim_lifted_irk(void *config_, sim_in *in, sim_out *out, void *opts_, void *mem_, void *work_)
{
    sim_solver_config *config = config_;
    sim_rk_opts *opts = opts_;

    if ( opts->ns != opts->tableau_size )
    {
        printf("Error in sim_lifted_irk: the Butcher tableau size does not match ns");
        return ACADOS_FAILURE;
    }

    int ns = opts->ns;

    sim_lifted_irk_dims *dims = (sim_lifted_irk_dims *) in->dims;

    int nx = dims->nx;
    int nu = dims->nu;
    int nz = dims->nz;

    // assert - only use supported features
    assert(nz == 0 && "nz should be zero - DAEs are not (yet) supported for this integrator");
    assert(opts->output_z == false &&
            "opts->output_z should be false - DAEs are not (yet) supported for this integrator");
    assert(opts->sens_algebraic == false &&
       "opts->sens_algebraic should be false - DAEs are not (yet) supported for this integrator");

    int dim_sys = ns * nx;
    int i, s1, s2, j, istep;
    sim_lifted_irk_memory *mem = (sim_lifted_irk_memory *) mem_;
    sim_lifted_irk_workspace *work = (sim_lifted_irk_workspace *) work_;
    sim_lifted_irk_workspace_cast(config, work, in, opts);
    double H_INT = in->T / opts->num_steps;
    int NF = opts->num_forw_sens;

    double *A_mat = opts->A_mat;
    double *b_vec = opts->b_vec;
    double *c_vec = opts->c_vec;

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

    ext_fun_arg_t ext_fun_type_in[5];
    void *ext_fun_in[5];  // XXX large enough ?
    ext_fun_arg_t ext_fun_type_out[5];
    void *ext_fun_out[5];  // XXX large enough ?

    lifted_irk_model *model = in->model;

    acados_timer timer, timer_la, timer_ad;
    double timing_la = 0.0;
    double timing_ad = 0.0;

    assert(NF == nx + nu && "Not implemented yet for other num_forw_sens");

    acados_tic(&timer);
    for (i = 0; i < nx; i++) out_tmp[i] = in->x[i];
    for (i = 0; i < nx * NF; i++) out_tmp[nx + i] = in->S_forw[i];  // sensitivities

    for (i = 0; i < nu; i++) rhs_in[nx * (1 + NF) + i] = in->u[i];

    if (opts->scheme->type == simplified_in || opts->scheme->type == simplified_inis)
    {
        for (i = 0; i < nx; i++) adj_tmp[i] = in->S_adj[i];
    }

    // Newton step of the collocation variables with respect to the inputs:
    if (NF == (nx + nu) && opts->sens_adj)
    {
        for (istep = opts->num_steps - 1; istep > -1; istep--)
        {  // ADJOINT update
            for (s1 = 0; s1 < ns; s1++)
            {
                for (j = 0; j < nx; j++)
                {  // step in X
                    for (i = 0; i < nx; i++)
                    {
                        K_traj[(istep * ns + s1) * nx + i] +=
                            DK_traj[(istep * ns + s1) * nx * (nx + nu) + j * nx + i] *
                            (in->x[j] - mem->x[j]);  // RK step
                    }
                    mem->x[j] = in->x[j];
                }
                for (j = 0; j < nu; j++)
                {  // step in U
                    for (i = 0; i < nx; i++)
                    {
                        K_traj[(istep * ns + s1) * nx + i] +=
                            DK_traj[(istep * ns + s1) * nx * (nx + nu) + (nx + j) * nx + i] *
                            (in->u[j] - mem->u[j]);  // RK step
                    }
                    mem->u[j] = in->u[j];
                }
            }
            if (opts->scheme->type == simplified_inis && !opts->scheme->freeze)
            {
                for (s1 = 0; s1 < ns; s1++)
                {
                    for (j = 0; j < NF; j++)
                    {
                        for (i = 0; i < nx; i++)
                        {
                            DK_traj[(istep * ns + s1) * nx * NF + j * nx + i] +=
                                delta_DK_traj[(istep * ns + s1) * nx * NF + j * nx + i];
                        }
                    }
                }
            }

            // Newton step of the Lagrange multipliers mu, based on adj_traj:
            if (opts->scheme->type == simplified_in || opts->scheme->type == simplified_inis)
            {
                for (s1 = 0; s1 < ns; s1++)
                {
                    for (i = 0; i < nx; i++)
                    {
                        sys_sol[s1 * nx + i] = -mu_traj[istep * ns * nx + s1 * nx + i];
                        for (s2 = 0; s2 < ns; s2++)
                        {
                            sys_sol[s1 * nx + i] += H_INT * A_mat[s1 * ns + s2] *
                                                    adj_traj[istep * ns * nx + s2 * nx + i];
                        }
                        sys_sol[s1 * nx + i] -= H_INT * b_vec[s1] * adj_tmp[i];
                        sys_sol[s1 * nx + i] += mem->grad_K[s1 * nx + i];
                    }
                }
                //                print_matrix("stdout", sys_sol, 1,
                //                ns*nx); print_matrix("stdout",
                //                sys_sol, 1, 1);

                // TRANSFORM using transf1_T:
                if (opts->scheme->type == simplified_in || opts->scheme->type == simplified_inis)
                {
                    // apply the transf1 operation:
                    for (s1 = 0; s1 < ns; s1++)
                    {
                        for (s2 = 0; s2 < ns; s2++)
                        {
                            work->trans[s2 * ns + s1] =
                                1.0 / H_INT * opts->scheme->transf1_T[s2 * ns + s1];
                        }
                    }
                    transform_vec(sys_sol, work->trans, sys_sol_trans, ns, nx);

                    // construct sys_sol2 from sys_sol_trans:
                    construct_subsystems(sys_sol_trans, sys_sol2, ns, nx, 1);
                }
                acados_tic(&timer_la);
                int idx = 0;
                for (s1 = 0; s1 < ns; s1++)
                {
                    // THIS LOOP IS PARALLELIZABLE BECAUSE OF DECOMPOSABLE
                    // LINEAR SUBSYSTEMS
                    if (opts->scheme->type == simplified_in ||
                        opts->scheme->type == simplified_inis)
                    {
                        sys_mat = sys_mat2[idx];
                        ipiv = ipiv2[idx];
                        sys_sol = sys_sol2[idx];
#if !TRIPLE_LOOP
                        str_mat = str_mat2[idx];
                        str_sol = str_sol2[idx];
#endif
                        idx++;
                        if ((s1 + 1) < ns)
                        {
                            dim_sys = 2 * nx;
                            s1++;  // complex conjugate pair of eigenvalues
                        }
                        else
                        {
                            dim_sys = nx;
                        }
                    }
                    else
                    {
                        dim_sys = ns * nx;
                        s1 = ns;  // break out of for-loop
                    }
#if TRIPLE_LOOP
                    solve_system_trans_ACADO(sys_mat, sys_sol, ipiv, dim_sys, 1);
#else  // TRIPLE_LOOP
#error : NOT YET IMPLEMENTED
#endif  // TRIPLE_LOOP
                }
                timing_la += acados_toc(&timer_la);
                // TRANSFORM using transf2_T:
                if (opts->scheme->type == simplified_in || opts->scheme->type == simplified_inis)
                {
                    // construct sys_sol_trans from sys_sol2:
                    destruct_subsystems(sys_sol_trans, sys_sol2, ns, nx, 1);

                    // apply the transf2 operation:
                    sys_sol = work->sys_sol;
                    transform_vec(sys_sol_trans, opts->scheme->transf2_T, sys_sol, ns, nx);
                }

                // update mu_traj
                for (s1 = 0; s1 < ns; s1++)
                {
                    for (i = 0; i < nx; i++)
                    {
                        mu_traj[istep * ns * nx + s1 * nx + i] += sys_sol[s1 * nx + i];
                    }
                }

                // update adj_tmp:
                // TODO(rien): USE ADJOINT DIFFERENTIATION HERE INSTEAD !!:
                for (j = 0; j < nx; j++)
                {
                    for (s1 = 0; s1 < ns; s1++)
                    {
                        for (i = 0; i < nx; i++)
                        {
                            adj_tmp[j] -= mu_traj[istep * ns * nx + s1 * nx + i] *
                                          jac_traj[istep * ns + s1][j * nx + i];
                        }
                    }
                }
            }
        }
    }

    if (opts->scheme->type == simplified_in || opts->scheme->type == simplified_inis)
    {
        for (i = 0; i < NF; i++) out->grad[i] = 0.0;
    }

    for (istep = 0; istep < opts->num_steps; istep++)
    {
        // form linear system matrix (explicit ODE case):

        form_linear_system_matrix(config_, istep, in, opts, mem, work, sys_mat, sys_mat2,
                                  timing_ad);

        int idx;
        if (opts->scheme->type == exact || (istep == 0 && !opts->scheme->freeze))
        {
            acados_tic(&timer_la);
            idx = 0;
            for (s1 = 0; s1 < ns; s1++)
            {
                // THIS LOOP IS PARALLELIZABLE BECAUSE OF DECOMPOSABLE LINEAR
                // SUBSYSTEMS
                if (opts->scheme->type == simplified_in || opts->scheme->type == simplified_inis)
                {
                    sys_mat = sys_mat2[idx];
                    ipiv = ipiv2[idx];
#if !TRIPLE_LOOP
                    str_mat = str_mat2[idx];
#endif
                    idx++;
                    if ((s1 + 1) < ns)
                    {
                        dim_sys = 2 * nx;
                        s1++;  // complex conjugate pair of eigenvalues
                    }
                    else
                    {
                        dim_sys = nx;
                    }
                }
                else
                {
                    dim_sys = ns * nx;
                    s1 = ns;  // break out of for-loop
                }
#if TRIPLE_LOOP
                LU_system_ACADO(sys_mat, ipiv, dim_sys);
#else  // TRIPLE_LOOP
// ---- BLASFEO: LU factorization ----
#if defined(LA_HIGH_PERFORMANCE)
                blasfeo_pack_dmat(dim_sys, dim_sys, sys_mat, dim_sys, str_mat, 0, 0);  // mat2strmat
#endif  // LA_BLAS | LA_REFERENCE
                blasfeo_dgetrf_rowpivot(dim_sys, dim_sys, str_mat, 0, 0, str_mat, 0, 0,
                                        ipiv);  // Gauss elimination
                                                // ---- BLASFEO: LU factorization ----
#endif  // TRIPLE_LOOP
            }
            timing_la += acados_toc(&timer_la);
        }

        if (opts->scheme->type == simplified_in || opts->scheme->type == simplified_inis)
        {
            sys_sol = work->sys_sol;
            sys_sol_trans = work->sys_sol_trans;
        }
        for (s1 = 0; s1 < ns; s1++)
        {
            for (i = 0; i < nx * (1 + NF); i++)
            {
                rhs_in[i] = out_tmp[i];
            }
            for (s2 = 0; s2 < ns; s2++)
            {
                for (i = 0; i < nx; i++)
                {
                    rhs_in[i] +=
                        H_INT * A_mat[s2 * ns + s1] * K_traj[istep * ns * nx + s2 * nx + i];
                }
                if (opts->scheme->type == simplified_inis)
                {
                    for (j = 0; j < NF; j++)
                    {
                        for (i = 0; i < nx; i++)
                        {
                            rhs_in[(j + 1) * nx + i] +=
                                H_INT * A_mat[s2 * ns + s1] *
                                DK_traj[(istep * ns + s2) * nx * NF + j * nx + i];
                        }
                    }
                }
            }
            rhs_in[nx * (1 + NF) + nu] =
                ((double) istep + c_vec[s1]) / ((double) opts->num_steps);  // time

            acados_tic(&timer_ad);

            ext_fun_type_in[0] = COLMAJ;
            ext_fun_in[0] = rhs_in + 0;  // x: nx
            ext_fun_type_in[1] = COLMAJ;
            ext_fun_in[1] = rhs_in + nx;  // Sx: nx*nx
            ext_fun_type_in[2] = COLMAJ;
            ext_fun_in[2] = rhs_in + nx + nx * nx;  // Su: nx*nu
            ext_fun_type_in[3] = COLMAJ;
            ext_fun_in[3] = rhs_in + nx + nx * nx + nx * nu;  // u: nu

            ext_fun_type_out[0] = COLMAJ;
            ext_fun_out[0] = VDE_tmp[s1] + 0;  // fun: nx
            ext_fun_type_out[1] = COLMAJ;
            ext_fun_out[1] = VDE_tmp[s1] + nx;  // Sx: nx*nx
            ext_fun_type_out[2] = COLMAJ;
            ext_fun_out[2] = VDE_tmp[s1] + nx + nx * nx;  // Su: nx*nu

            model->expl_vde_for->evaluate(model->expl_vde_for, ext_fun_type_in, ext_fun_in,
                                          ext_fun_type_out, ext_fun_out);  // k evaluation
            timing_ad += acados_toc(&timer_ad);

            // put VDE_tmp in sys_sol:
            for (j = 0; j < 1 + NF; j++)
            {
                for (i = 0; i < nx; i++)
                {
                    sys_sol[j * ns * nx + s1 * nx + i] = VDE_tmp[s1][j * nx + i];
                }
            }
            for (i = 0; i < nx; i++)
            {
                sys_sol[s1 * nx + i] -= K_traj[istep * ns * nx + s1 * nx + i];
            }
            if (opts->scheme->type == simplified_inis)
            {
                for (j = 0; j < NF; j++)
                {
                    for (i = 0; i < nx; i++)
                    {
                        sys_sol[(j + 1) * ns * nx + s1 * nx + i] -=
                            DK_traj[(istep * ns + s1) * nx * NF + j * nx + i];
                    }
                }
            }
        }

        // Inexact Newton with Iterated Sensitivities (INIS) based gradient
        // correction:
        if (opts->scheme->type == simplified_inis)
        {
            for (j = 0; j < NF; j++)
            {
                for (i = 0; i < ns * nx; i++)
                {
                    out->grad[j] += mu_traj[istep * ns * nx + i] * sys_sol[(j + 1) * ns * nx + i];
                }
            }
        }

        if (opts->scheme->type == simplified_in || opts->scheme->type == simplified_inis)
        {
            // apply the transf1 operation:
            for (s1 = 0; s1 < ns; s1++)
            {
                for (s2 = 0; s2 < ns; s2++)
                {
                    work->trans[s2 * ns + s1] = 1.0 / H_INT * opts->scheme->transf1[s2 * ns + s1];
                }
            }
            if (!opts->scheme->freeze)
            {
                transform_mat(sys_sol, work->trans, sys_sol_trans, ns, nx, 1 + NF);
            }
            else
            {
                transform_mat(sys_sol, work->trans, sys_sol_trans, ns, nx, 1);
            }

            // construct sys_sol2 from sys_sol_trans:
            if (!opts->scheme->freeze)
            {
                construct_subsystems(sys_sol_trans, sys_sol2, ns, nx, 1 + NF);
            }
            else
            {
                construct_subsystems(sys_sol_trans, sys_sol2, ns, nx, 1);
            }
        }

        acados_tic(&timer_la);
        idx = 0;
        for (s1 = 0; s1 < ns; s1++)
        {
            // THIS LOOP IS PARALLELIZABLE BECAUSE OF DECOMPOSABLE LINEAR
            // SUBSYSTEMS
            if (opts->scheme->type == simplified_in || opts->scheme->type == simplified_inis)
            {
                sys_mat = sys_mat2[idx];
                ipiv = ipiv2[idx];
                sys_sol = sys_sol2[idx];
#if !TRIPLE_LOOP
                str_mat = str_mat2[idx];
                str_sol = str_sol2[idx];
#endif
                idx++;
                if ((s1 + 1) < ns)
                {
                    dim_sys = 2 * nx;
                    s1++;  // complex conjugate pair of eigenvalues
                }
                else
                {
                    dim_sys = nx;
                }
            }
            else
            {
                dim_sys = ns * nx;
                s1 = ns;  // break out of for-loop
            }
#if TRIPLE_LOOP
            if (!opts->scheme->freeze)
            {
                solve_system_ACADO(sys_mat, sys_sol, ipiv, dim_sys, 1 + NF);
            }
            else
            {
                solve_system_ACADO(sys_mat, sys_sol, ipiv, dim_sys, 1);
            }
#else  // TRIPLE_LOOP
#if defined(LA_HIGH_PERFORMANCE)
            // ---- BLASFEO: row transformations + backsolve ----
            blasfeo_pack_dmat(dim_sys, 1 + NF, sys_sol, dim_sys, str_sol, 0, 0);  // mat2strmat
            blasfeo_drowpe(dim_sys, ipiv, str_sol);  // row permutations
            blasfeo_dtrsm_llnu(dim_sys, 1 + NF, 1.0, str_mat, 0, 0, str_sol, 0, 0, str_sol, 0,
                               0);  // L backsolve
            blasfeo_dtrsm_lunn(dim_sys, 1 + NF, 1.0, str_mat, 0, 0, str_sol, 0, 0, str_sol, 0,
                               0);                                                  // U backsolve
            blasfeo_unpack_dmat(dim_sys, 1 + NF, str_sol, 0, 0, sys_sol, dim_sys);  // strmat2mat

// BLASFEO: row transformations + backsolve
#else   // LA_BLAS | LA_REFERENCE
            // ---- BLASFEO: row transformations + backsolve ----
            blasfeo_drowpe(dim_sys, ipiv, str_sol);  // row permutations
            blasfeo_dtrsm_llnu(dim_sys, 1 + NF, 1.0, str_mat, 0, 0, str_sol, 0, 0, str_sol, 0,
                               0);  // L backsolve
            blasfeo_dtrsm_lunn(dim_sys, 1 + NF, 1.0, str_mat, 0, 0, str_sol, 0, 0, str_sol, 0,
                               0);  // U backsolve
                                    // BLASFEO: row transformations + backsolve
#endif  // LA_BLAFEO
#endif  // TRIPLE_LOOP
        }
        timing_la += acados_toc(&timer_la);
        if (opts->scheme->type == simplified_in || opts->scheme->type == simplified_inis)
        {
            // construct sys_sol_trans from sys_sol2:
            if (!opts->scheme->freeze)
            {
                destruct_subsystems(sys_sol_trans, sys_sol2, ns, nx, 1 + NF);
            }
            else
            {
                destruct_subsystems(sys_sol_trans, sys_sol2, ns, nx, 1);
            }

            // apply the transf2 operation:
            sys_sol = work->sys_sol;
            if (!opts->scheme->freeze)
            {
                transform_mat(sys_sol_trans, opts->scheme->transf2, sys_sol, ns, nx, 1 + NF);
            }
            else
            {
                transform_mat(sys_sol_trans, opts->scheme->transf2, sys_sol, ns, nx, 1);
            }
        }

        // Newton step of the collocation variables
        for (i = 0; i < ns * nx; i++)
        {
            K_traj[istep * ns * nx + i] += sys_sol[i];
        }
        for (s1 = 0; s1 < ns; s1++)
        {
            for (i = 0; i < nx; i++)
            {
                out_tmp[i] += H_INT * b_vec[s1] * K_traj[istep * ns * nx + s1 * nx + i];  // RK step
            }
        }

        // Sensitivities collocation variables
        for (s1 = 0; s1 < ns; s1++)
        {
            for (j = 0; j < NF; j++)
            {
                for (i = 0; i < nx; i++)
                {
                    if (opts->scheme->type == simplified_inis && opts->sens_adj &&
                        !opts->scheme->freeze)
                    {
                        delta_DK_traj[(istep * ns + s1) * nx * NF + j * nx + i] =
                            sys_sol[(j + 1) * ns * nx + s1 * nx + i];
                    }
                    else if (opts->scheme->type == simplified_inis && !opts->scheme->freeze)
                    {
                        DK_traj[(istep * ns + s1) * nx * NF + j * nx + i] +=
                            sys_sol[(j + 1) * ns * nx + s1 * nx + i];
                    }
                    else if (!opts->scheme->freeze)
                    {
                        DK_traj[(istep * ns + s1) * nx * NF + j * nx + i] =
                            sys_sol[(j + 1) * ns * nx + s1 * nx + i];
                    }
                }
            }
        }
        for (s1 = 0; s1 < ns; s1++)
        {
            for (i = 0; i < nx * NF; i++)
            {
                out_tmp[nx + i] +=
                    H_INT * b_vec[s1] * DK_traj[(istep * ns + s1) * nx * NF + i];  // RK step
            }
        }
        if (opts->scheme->type == simplified_inis || opts->scheme->type == simplified_in)
        {
            // Adjoint derivatives:
            for (s1 = 0; s1 < ns; s1++)
            {
                for (j = 0; j < nx; j++)
                {
                    adj_traj[istep * ns * nx + s1 * nx + j] = 0.0;
                    for (i = 0; i < nx; i++)
                    {
                        adj_traj[istep * ns * nx + s1 * nx + j] +=
                            mu_traj[istep * ns * nx + s1 * nx + i] *
                            jac_traj[istep * ns + s1][j * nx + i];
                    }
                }
            }
        }
        //        print_matrix_name("stdout", "adj_traj", &adj_traj[0], 1,
        //        ns*nx); print_matrix_name("stdout", "mu_traj",
        //        &mu_traj[0], 1, ns*nx); print_matrix_name("stdout",
        //        "DK_traj", &DK_traj[0], 1, nx*NF); print_matrix_name("stdout",
        //        "VDE_tmp[0]", &VDE_tmp[0][0], 1, nx*(NF+1));
        //        print_matrix_name("stdout", "VDE_tmp[1]", &VDE_tmp[1][0], 1,
        //        nx*(NF+1));
        if (opts->scheme->type == simplified_in)
        {
            // Standard Inexact Newton based gradient correction:
            for (j = 0; j < NF; j++)
            {
                for (s1 = 0; s1 < ns; s1++)
                {
                    for (i = 0; i < nx; i++)
                    {
                        out->grad[j] +=
                            mu_traj[istep * ns * nx + s1 * nx + i] * VDE_tmp[s1][(j + 1) * nx + i];
                    }
                }
            }
            for (j = 0; j < NF; j++)
            {
                for (s1 = 0; s1 < ns; s1++)
                {
                    for (i = 0; i < nx; i++)
                    {
                        out->grad[j] -= mu_traj[istep * ns * nx + s1 * nx + i] *
                                        DK_traj[(istep * ns + s1) * nx * NF + j * nx + i];
                    }
                }
                for (s2 = 0; s2 < ns; s2++)
                {
                    for (s1 = 0; s1 < ns; s1++)
                    {
                        for (i = 0; i < nx; i++)
                        {
                            out->grad[j] += H_INT * A_mat[s2 * ns + s1] *
                                            adj_traj[istep * ns * nx + s1 * nx + i] *
                                            DK_traj[(istep * ns + s2) * nx * NF + j * nx + i];
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

void sim_lifted_irk_config_initialize_default(void *config_)
{
    sim_solver_config *config = config_;

    config->evaluate = &sim_lifted_irk;
    config->opts_calculate_size = &sim_lifted_irk_opts_calculate_size;
    config->opts_assign = &sim_lifted_irk_opts_assign;
    config->opts_initialize_default = &sim_lifted_irk_opts_initialize_default;
    config->opts_update = &sim_lifted_irk_opts_update;
    config->memory_calculate_size = &sim_lifted_irk_memory_calculate_size;
    config->memory_assign = &sim_lifted_irk_memory_assign;
    config->workspace_calculate_size = &sim_lifted_irk_workspace_calculate_size;
    config->model_calculate_size = &sim_lifted_irk_model_calculate_size;
    config->model_assign = &sim_lifted_irk_model_assign;
    config->model_set_function = &sim_lifted_irk_model_set_function;
    config->dims_calculate_size = &sim_lifted_irk_dims_calculate_size;
    config->dims_assign = &sim_lifted_irk_dims_assign;

    config->get_nx = &sim_lifted_irk_get_nx;
    config->get_nu = &sim_lifted_irk_get_nu;
    config->get_nz = &sim_lifted_irk_get_nz;
    return;
}
