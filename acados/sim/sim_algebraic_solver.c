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

#include "acados/sim/sim_algebraic_solver.h"

// standard
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
// acados
#include "acados/utils/mem.h"
#include "acados/utils/print.h"
#include "acados/utils/math.h"

#include "acados/sim/sim_common.h"

#include "blasfeo/include/blasfeo_d_aux.h"
#include "blasfeo/include/blasfeo_d_blas.h"

/************************************************
 * dims
 ************************************************/

int sim_algebraic_solver_dims_calculate_size()
{
    int size = sizeof(sim_algebraic_solver_dims);

    return size;
}

void *sim_algebraic_solver_dims_assign(void *config_, void *raw_memory)
{
    char *c_ptr = raw_memory;

    sim_algebraic_solver_dims *dims = (sim_algebraic_solver_dims *) c_ptr;
    c_ptr += sizeof(sim_algebraic_solver_dims);

    dims->nx = 0;
    dims->nu = 0;
    dims->nz = 0;

    assert((char *) raw_memory + sim_algebraic_solver_dims_calculate_size() >= c_ptr);

    return dims;
}

static void sim_algebraic_solver_set_nu(void *config_, void *dims_, const int *nu)
{
    sim_algebraic_solver_dims *dims = (sim_algebraic_solver_dims *) dims_;
    dims->nu = *nu;
}

static void sim_algebraic_solver_set_nx(void *config_, void *dims_, const int *nx)
{
    sim_algebraic_solver_dims *dims = (sim_algebraic_solver_dims *) dims_;
    dims->nx = *nx;
}

static void sim_algebraic_solver_set_nz(void *config_, void *dims_, const int *nz)
{
    sim_algebraic_solver_dims *dims = (sim_algebraic_solver_dims *) dims_;
    dims->nz = *nz;
}


void sim_algebraic_solver_dims_set(void *config_, void *dims_, const char *field, const int* value)
{
    if (!strcmp(field, "nx"))
    {
        sim_algebraic_solver_set_nx(config_, dims_, value);
    }
    else if (!strcmp(field, "nu"))
    {
        sim_algebraic_solver_set_nu(config_, dims_, value);
    }
    else if (!strcmp(field, "nz"))
    {
        sim_algebraic_solver_set_nz(config_, dims_, value);
    }
    else
    {
        printf("\nerror: dimension type not available in module\n");
        exit(1);
    }
}

void sim_algebraic_solver_get_nx(void *dims_, int *nx)
{
    sim_algebraic_solver_dims *dims = (sim_algebraic_solver_dims *) dims_;
    *nx = dims->nx;
}

void sim_algebraic_solver_get_nu(void *dims_, int *nu)
{
    sim_algebraic_solver_dims *dims = (sim_algebraic_solver_dims *) dims_;
    *nu = dims->nu;
}

void sim_algebraic_solver_get_nz(void *dims_, int *nz)
{
    sim_algebraic_solver_dims *dims = (sim_algebraic_solver_dims *) dims_;
    *nz = dims->nz;
}

/************************************************
 * model
 ************************************************/

int sim_algebraic_solver_model_calculate_size(void *config, void *dims)
{
    int size = 0;

    size += sizeof(irk_model);

    return size;
}

void *sim_algebraic_solver_model_assign(void *config, void *dims, void *raw_memory)
{
    char *c_ptr = (char *) raw_memory;

    irk_model *data = (irk_model *) c_ptr;
    c_ptr += sizeof(irk_model);

    assert((char *) raw_memory + sim_algebraic_solver_model_calculate_size(config, dims) >= c_ptr);

    return data;
}

int sim_algebraic_solver_model_set_function(void *model_, sim_function_t fun_type, void *fun)
{
    irk_model *model = model_;

    switch (fun_type)
    {
        case IMPL_ODE_FUN:
            model->impl_ode_fun = (external_function_generic *) fun;
            break;
        case IMPL_ODE_FUN_JAC_X_XDOT:
            model->impl_ode_fun_jac_x_xdot_z = (external_function_generic *) fun;
            break;
        case IMPL_ODE_JAC_X_XDOT_U:
            model->impl_ode_jac_x_xdot_u_z = (external_function_generic *) fun;
            break;
        case IMPL_ODE_HESS:
            model->impl_ode_hess = (external_function_generic *) fun;
            break;
        default:
            return ACADOS_FAILURE;
    }
    return ACADOS_SUCCESS;
}

/************************************************
 * opts
 ************************************************/

int sim_algebraic_solver_opts_calculate_size(void *config_, void *dims)
{
    int ns_max = NS_MAX;

    int size = 0;

    size += sizeof(sim_rk_opts);

    size += ns_max * ns_max * sizeof(double);  // A_mat
    size += ns_max * sizeof(double);           // b_vec
    size += ns_max * sizeof(double);           // c_vec

    int tmp0 = gauss_nodes_work_calculate_size(ns_max);
    int tmp1 = butcher_table_work_calculate_size(ns_max);
    int work_size = tmp0 > tmp1 ? tmp0 : tmp1;
    size += work_size;  // work

    make_int_multiple_of(8, &size);
    size += 1 * 8;

    return size;
}

void *sim_algebraic_solver_opts_assign(void *config_, void *dims, void *raw_memory)
{
    int ns_max = NS_MAX;

    char *c_ptr = (char *) raw_memory;

    sim_rk_opts *opts = (sim_rk_opts *) c_ptr;
    c_ptr += sizeof(sim_rk_opts);

    align_char_to(8, &c_ptr);

    assign_and_advance_double(ns_max * ns_max, &opts->A_mat, &c_ptr);
    assign_and_advance_double(ns_max, &opts->b_vec, &c_ptr);
    assign_and_advance_double(ns_max, &opts->c_vec, &c_ptr);

    // work
    int tmp0 = gauss_nodes_work_calculate_size(ns_max);
    int tmp1 = butcher_table_work_calculate_size(ns_max);
    int work_size = tmp0 > tmp1 ? tmp0 : tmp1;
    opts->work = c_ptr;
    c_ptr += work_size;

    assert((char *) raw_memory + sim_algebraic_solver_opts_calculate_size(config_, dims) >= c_ptr);

    return (void *) opts;
}

void sim_algebraic_solver_opts_initialize_default(void *config_, void *dims_, void *opts_)
{
    sim_algebraic_solver_dims *dims = (sim_algebraic_solver_dims *) dims_;
    sim_rk_opts *opts = opts_;

    opts->ns = 3;  // GL 3
    int ns = opts->ns;

    assert(ns <= NS_MAX && "ns > NS_MAX!");

    // set tableau size
    opts->tableau_size = opts->ns;

    // gauss collocation nodes
    gauss_nodes(ns, opts->c_vec, opts->work);

    // butcher tableau
    butcher_table(ns, opts->c_vec, opts->b_vec, opts->A_mat, opts->work);

    // default options
    opts->newton_iter = 3;
    opts->scheme = NULL;
    opts->num_steps = 1;
    opts->num_forw_sens = dims->nx + dims->nu;
    opts->sens_forw = true;
    opts->sens_adj = false;
    opts->sens_hess = false;
    opts->jac_reuse = true;

    if (dims->nz > 0) {
        opts->output_z = true;
        opts->sens_algebraic = true;
    } else {
        opts->output_z = false;
        opts->sens_algebraic = false;
    }

    return;
}

void sim_algebraic_solver_opts_update(void *config_, void *dims, void *opts_)
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

    return;
}


int sim_algebraic_solver_opts_set(void *config_, void *opts_, const char *field, void *value)
{
    sim_rk_opts *opts = (sim_rk_opts *) opts_;
    return sim_rk_opts_set(opts, field, value);
}

/************************************************
 * memory
 ************************************************/

int sim_algebraic_solver_memory_calculate_size(void *config, void *dims, void *opts_) { return 0; }
void *sim_algebraic_solver_memory_assign(void *config, void *dims, void *opts_, void *raw_memory)
{
    return NULL;
}

/************************************************
 * workspace
 ************************************************/

int sim_algebraic_solver_workspace_calculate_size(void *config_, void *dims_, void *opts_)
{
    sim_algebraic_solver_dims *dims = (sim_algebraic_solver_dims *) dims_;
    sim_rk_opts *opts = opts_;

    int ns = opts->ns;

    int nx = dims->nx;
    int nu = dims->nu;
    int nz = dims->nz;

    int nK = (nx + nz) * ns;

    int steps = opts->num_steps;

    int size = sizeof(sim_algebraic_solver_workspace);

    if (opts->sens_algebraic || opts->output_z){
        size += (nx + nz) * sizeof(int);    // ipiv_one_stage
        size += ns * sizeof(double);        // Z_work
    }

    /* blasfeo structs */
    size += 4 * sizeof(struct blasfeo_dvec);          // rG, K, xt, xn

    if (opts->sens_adj || opts->sens_hess){
        size += 2 * steps * sizeof(struct blasfeo_dvec);  // **xn_traj, **K_traj
    }

    size += 2 * sizeof(struct blasfeo_dvec);  // lambda, lambdaK
    if (!opts->sens_hess){
        size += 4 * sizeof(struct blasfeo_dmat);  // dG_dxu, dG_dK, dK_dxu, S_forw
    } else {
        size += (4 * steps + 1) * sizeof(struct blasfeo_dmat);  // dG_dxu, dG_dK, dK_dxu, S_forw
    }

    /* blasfeo mem */
    size += blasfeo_memsize_dvec(nK);   // K
    size += blasfeo_memsize_dvec(nK);   // rG
    size += 3 * blasfeo_memsize_dvec(nx);           // xt, xn, xtdot
    size += 1 * blasfeo_memsize_dvec(nx + nu);      // lambda
    size += 1 * blasfeo_memsize_dvec(nK);           // lambdaK

    if (!opts->sens_hess){
        size += 1 * blasfeo_memsize_dmat(nK, nx + nu);  // dG_dxu
        size += 1 * blasfeo_memsize_dmat(nK, nK);       // dG_dK
        size += 1 * blasfeo_memsize_dmat(nK, nx + nu);  // dK_dxu
        size += 1 * blasfeo_memsize_dmat(nx, nx + nu);  // S_forw
        size += nK * sizeof(int);  // ipiv

    } else {
        size += steps * blasfeo_memsize_dmat(nK, nx + nu);      // dG_dxu
        size += steps * blasfeo_memsize_dmat(nK, nK);           // dG_dK
        size += steps * blasfeo_memsize_dmat(nK, nx + nu);      // dK_dxu
        size += (steps + 1) * blasfeo_memsize_dmat(nx, nx + nu);      // S_forw
        size += steps * nK * sizeof(int);  // ipiv

        size += 1 * blasfeo_memsize_dmat(nx + nu, nx + nu);  // f_hess
        size += 1 * blasfeo_memsize_dmat(2 * nx + nu + nz, nx + nu);  // dxkzu_dw0
        size += 1 * blasfeo_memsize_dmat(nx + nu, nx + nu);  // Hess
    }

    if ( opts->sens_adj || opts->sens_hess ){
        size += steps * blasfeo_memsize_dvec(nx);       // for xn_traj
        size += steps * blasfeo_memsize_dvec(nK);       // for K_traj
    }

    size += 2 * blasfeo_memsize_dmat(nx + nz, nx);  // df_dx, df_dxdot
    size += blasfeo_memsize_dmat(nx + nz, nu);      // df_du
    size += blasfeo_memsize_dmat(nx + nz, nz);      // df_dz

    if (opts->sens_algebraic){
        size += blasfeo_memsize_dmat(nx + nz, nx + nz);  // df_dxdotz
        size += blasfeo_memsize_dmat(nx + nz, nx + nu);  // dk0_dxu
    }

    make_int_multiple_of(64, &size);
    size += 1 * 64;

    return size;
}

static void *sim_algebraic_solver_workspace_cast(void *config_, void *dims_, void *opts_, void *raw_memory)
{
    sim_rk_opts *opts = opts_;
    sim_algebraic_solver_dims *dims = (sim_algebraic_solver_dims *) dims_;

    int ns = opts->ns;

    int nx = dims->nx;
    int nu = dims->nu;
    int nz = dims->nz;
    int nK = (nx + nz) * ns;

    int steps = opts->num_steps;

    char *c_ptr = (char *) raw_memory;

    sim_algebraic_solver_workspace *workspace = (sim_algebraic_solver_workspace *) c_ptr;
    c_ptr += sizeof(sim_algebraic_solver_workspace);

    if ( opts->sens_adj || opts->sens_hess ){
        assign_and_advance_blasfeo_dvec_structs(steps, &workspace->xn_traj, &c_ptr);
        assign_and_advance_blasfeo_dvec_structs(steps, &workspace->K_traj, &c_ptr);
    }

    assign_and_advance_blasfeo_dvec_structs(1, &workspace->rG, &c_ptr);
    assign_and_advance_blasfeo_dvec_structs(1, &workspace->K, &c_ptr);
    assign_and_advance_blasfeo_dvec_structs(1, &workspace->lambda, &c_ptr);
    assign_and_advance_blasfeo_dvec_structs(1, &workspace->lambdaK, &c_ptr);

    // dG_dxu, dG_dK, dK_dxu, S_forw
    if (!opts->sens_hess){
        assign_and_advance_blasfeo_dmat_structs(1, &workspace->dG_dxu, &c_ptr);
        assign_and_advance_blasfeo_dmat_structs(1, &workspace->dG_dK, &c_ptr);
        assign_and_advance_blasfeo_dmat_structs(1, &workspace->dK_dxu, &c_ptr);
        assign_and_advance_blasfeo_dmat_structs(1, &workspace->S_forw, &c_ptr);
    } else {
        assign_and_advance_blasfeo_dmat_structs(steps, &workspace->dG_dxu, &c_ptr);
        assign_and_advance_blasfeo_dmat_structs(steps, &workspace->dG_dK, &c_ptr);
        assign_and_advance_blasfeo_dmat_structs(steps, &workspace->dK_dxu, &c_ptr);
        assign_and_advance_blasfeo_dmat_structs(steps + 1, &workspace->S_forw, &c_ptr);
    }
    assign_and_advance_blasfeo_dvec_structs(1, &workspace->xt, &c_ptr);
    assign_and_advance_blasfeo_dvec_structs(1, &workspace->xn, &c_ptr);

    /* algin c_ptr to 64 blasfeo_dmat_mem has to be assigned directly after that  */
    align_char_to(64, &c_ptr);

    if (!opts->sens_hess){
        assign_and_advance_blasfeo_dmat_mem(nK, nx + nu, workspace->dG_dxu, &c_ptr);
        assign_and_advance_blasfeo_dmat_mem(nK, nK,      workspace->dG_dK, &c_ptr);
        assign_and_advance_blasfeo_dmat_mem(nK, nx + nu, workspace->dK_dxu, &c_ptr);
        assign_and_advance_blasfeo_dmat_mem(nx, nx + nu, workspace->S_forw, &c_ptr);
    }
    else
    {
        for (int ii = 0; ii < steps; ii++) {
            assign_and_advance_blasfeo_dmat_mem(nK, nx + nu, &workspace->dG_dxu[ii], &c_ptr);
            assign_and_advance_blasfeo_dmat_mem(nK, nK, &workspace->dG_dK[ii], &c_ptr);
            assign_and_advance_blasfeo_dmat_mem(nK, nx + nu, &workspace->dK_dxu[ii], &c_ptr);
            assign_and_advance_blasfeo_dmat_mem(nx, nx + nu, &workspace->S_forw[ii], &c_ptr);
        }
        assign_and_advance_blasfeo_dmat_mem(nx, nx + nu, &workspace->S_forw[steps], &c_ptr);
        assign_and_advance_blasfeo_dmat_mem(nx + nu, nx + nu, &workspace->f_hess, &c_ptr);
        assign_and_advance_blasfeo_dmat_mem(2*nx + nu + nz, nx + nu, &workspace->dxkzu_dw0, &c_ptr);
        assign_and_advance_blasfeo_dmat_mem(nx + nu, nx + nu, &workspace->Hess, &c_ptr);
    }

    assign_and_advance_blasfeo_dmat_mem(nx + nz, nx, &workspace->df_dx, &c_ptr);
    assign_and_advance_blasfeo_dmat_mem(nx + nz, nx, &workspace->df_dxdot, &c_ptr);
    assign_and_advance_blasfeo_dmat_mem(nx + nz, nu, &workspace->df_du, &c_ptr);
    assign_and_advance_blasfeo_dmat_mem(nx + nz, nz, &workspace->df_dz, &c_ptr);

    if (opts->sens_algebraic){
        assign_and_advance_blasfeo_dmat_mem(nx + nz, nx + nz, &workspace->df_dxdotz, &c_ptr);
        assign_and_advance_blasfeo_dmat_mem(nx + nz, nx + nu, &workspace->dk0_dxu, &c_ptr);
    }

    assign_and_advance_blasfeo_dvec_mem(nK, workspace->rG, &c_ptr);
    assign_and_advance_blasfeo_dvec_mem(nK, workspace->K, &c_ptr);
    assign_and_advance_blasfeo_dvec_mem(nx, workspace->xt, &c_ptr);
    assign_and_advance_blasfeo_dvec_mem(nx, workspace->xn, &c_ptr);
    assign_and_advance_blasfeo_dvec_mem(nx, &workspace->xtdot, &c_ptr);
    assign_and_advance_blasfeo_dvec_mem(nx + nu, workspace->lambda, &c_ptr);
    assign_and_advance_blasfeo_dvec_mem(nK, workspace->lambdaK, &c_ptr);


    if ( opts->sens_adj || opts->sens_hess ){
        for (int i = 0; i < steps; i++)
        {
            assign_and_advance_blasfeo_dvec_mem(nx, &workspace->xn_traj[i], &c_ptr);
            assign_and_advance_blasfeo_dvec_mem(nK, &workspace->K_traj[i], &c_ptr);
        }
    }

    if (opts->sens_algebraic || opts->output_z){
        assign_and_advance_double(ns, &workspace->Z_work, &c_ptr);
        assign_and_advance_int((nx + nz), &workspace->ipiv_one_stage, &c_ptr);
    }

    if (!opts->sens_hess){
        assign_and_advance_int(nK, &workspace->ipiv, &c_ptr);
    } else {
        assign_and_advance_int(steps * nK, &workspace->ipiv, &c_ptr);
    }

    // printf("\npointer moved - size calculated = %d bytes\n", c_ptr- (char*)raw_memory -
    // sim_algebraic_solver_calculate_workspace_size(dims, opts_));

    assert((char *) raw_memory + sim_algebraic_solver_workspace_calculate_size(config_, dims, opts_) >= c_ptr);

    return (void *) workspace;
}

int sim_algebraic_solver_precompute(void *config_, sim_in *in, sim_out *out, void *opts_, void *mem_,
                       void *work_)
{
    return ACADOS_SUCCESS;
}

/************************************************
 * algebraic_solver
 ************************************************/

int sim_algebraic_solver(void *config_, sim_in *in, sim_out *out, void *opts_, void *mem_, void *work_)
{
/* Get variables from workspace, etc; -- 
    cast pointers */
    sim_solver_config *config = config_;
    sim_rk_opts *opts = opts_;

    // if ( opts->ns != opts->tableau_size )
    // {
    //     printf("Error in sim_algebraic_solver: the Butcher tableau size does not match ns");
    //     return ACADOS_FAILURE;
    // }
    // int ns = opts->ns;

    void *dims_ = in->dims;
    sim_algebraic_solver_dims *dims = (sim_algebraic_solver_dims *) dims_;
    sim_algebraic_solver_workspace *workspace =
        (sim_algebraic_solver_workspace *) sim_algebraic_solver_workspace_cast(config, dims, opts, work_);

    irk_model *model = in->model;

    int nx = dims->nx;
    int nu = dims->nu;
    int nz = dims->nz;

    int nK = (nx + nz);

    double *u = in->u;

    int newton_iter = opts->newton_iter;
    double *A_mat = opts->A_mat;
    double *b_vec = opts->b_vec;
    int num_steps = opts->num_steps;
    double step = in->T / num_steps;

    int *ipiv = workspace->ipiv;
    int *ipiv_one_stage = workspace->ipiv_one_stage;
    double *Z_work = workspace->Z_work;

    struct blasfeo_dmat *dG_dK = workspace->dG_dK;
    struct blasfeo_dvec *rG = workspace->rG;
    struct blasfeo_dvec *K = workspace->K;
    struct blasfeo_dmat *dG_dxu = workspace->dG_dxu;
    struct blasfeo_dmat *dK_dxu = workspace->dK_dxu;
    struct blasfeo_dvec *xt = workspace->xt;
    struct blasfeo_dvec *xtdot = &workspace->xtdot;

    struct blasfeo_dvec *xn = workspace->xn;
    struct blasfeo_dmat *S_forw = workspace->S_forw;

    struct blasfeo_dmat df_dx = workspace->df_dx;
    struct blasfeo_dmat df_dxdot = workspace->df_dxdot;
    struct blasfeo_dmat df_du = workspace->df_du;
    struct blasfeo_dmat df_dz = workspace->df_dz;
    struct blasfeo_dmat f_hess = workspace->f_hess;
    struct blasfeo_dmat dxkzu_dw0 = workspace->dxkzu_dw0;

    struct blasfeo_dmat df_dxdotz = workspace->df_dxdotz;
    struct blasfeo_dmat dk0_dxu = workspace->dk0_dxu;

    // for adjoint
    struct blasfeo_dvec *lambda = workspace->lambda;
    struct blasfeo_dvec *lambdaK = workspace->lambdaK;
    struct blasfeo_dvec *xn_traj = workspace->xn_traj;
    struct blasfeo_dvec *K_traj = workspace->K_traj;

    // for hessians only
    struct blasfeo_dmat Hess = workspace->Hess;

    double *x_out = out->xn;
    double *S_forw_out = out->S_forw;
    double *S_adj_out = out->S_adj;
    double *S_algebraic = out->S_algebraic;

/* declare */
    acados_timer timer, timer_ad, timer_la;

    double a;
    // struct blasfeo_dmat *dG_dK;
    // struct blasfeo_dmat *dG_dxu;
    // struct blasfeo_dmat *dK_dxu;
    // struct blasfeo_dmat *S_forw;
    int *ipiv_ss;


/* SET FUNCTION IN- & OUTPUT TYPES */

    // INPUT: impl_ode
    ext_fun_arg_t impl_ode_type_in[4];
    void *impl_ode_in[4];

    impl_ode_type_in[0] = BLASFEO_DVEC;       // xt
    impl_ode_type_in[1] = BLASFEO_DVEC_ARGS;  // k_i
    impl_ode_type_in[2] = COLMAJ;             // u
    impl_ode_type_in[3] = BLASFEO_DVEC_ARGS;  // z_i

    struct blasfeo_dvec_args impl_ode_xdot_in;
    struct blasfeo_dvec_args impl_ode_z_in;

    impl_ode_in[0] = xt;  // 1st input is always xt
    impl_ode_in[1] =
        &impl_ode_xdot_in;  // 2nd input is part of K[ss], always update impl_ode_xdot_in
    impl_ode_in[2] = u;     // 3rd input is u (always)
    impl_ode_in[3] = &impl_ode_z_in;     // 4th input is part of Z[ss]

    // OUTPUT:
    // impl_ode_fun
    ext_fun_arg_t impl_ode_fun_type_out[1];
    void *impl_ode_fun_out[1];
    impl_ode_fun_type_out[0] = BLASFEO_DVEC_ARGS;

    struct blasfeo_dvec_args impl_ode_res_out;
    impl_ode_res_out.x = rG;

    impl_ode_fun_out[0] = &impl_ode_res_out;

    // impl_ode_fun_jac_x_xdot_z
    ext_fun_arg_t impl_ode_fun_jac_x_xdot_z_type_out[4];
    void *impl_ode_fun_jac_x_xdot_z_out[4];
    impl_ode_fun_jac_x_xdot_z_type_out[0] = BLASFEO_DVEC_ARGS;
    impl_ode_fun_jac_x_xdot_z_out[0] = &impl_ode_res_out;
    impl_ode_fun_jac_x_xdot_z_type_out[1] = BLASFEO_DMAT;
    impl_ode_fun_jac_x_xdot_z_out[1] = &df_dx;
    impl_ode_fun_jac_x_xdot_z_type_out[2] = BLASFEO_DMAT;
    impl_ode_fun_jac_x_xdot_z_out[2] = &df_dxdot;
    impl_ode_fun_jac_x_xdot_z_type_out[3] = BLASFEO_DMAT;
    impl_ode_fun_jac_x_xdot_z_out[3] = &df_dz;

    // impl_ode_jac_x_xdot_u_z
    ext_fun_arg_t impl_ode_jac_x_xdot_u_z_type_out[4];
    void *impl_ode_jac_x_xdot_u_z_out[4];
    impl_ode_jac_x_xdot_u_z_type_out[0] = BLASFEO_DMAT;
    impl_ode_jac_x_xdot_u_z_out[0] = &df_dx;
    impl_ode_jac_x_xdot_u_z_type_out[1] = BLASFEO_DMAT;
    impl_ode_jac_x_xdot_u_z_out[1] = &df_dxdot;
    impl_ode_jac_x_xdot_u_z_type_out[2] = BLASFEO_DMAT;
    impl_ode_jac_x_xdot_u_z_out[2] = &df_du;
    impl_ode_jac_x_xdot_u_z_type_out[3] = BLASFEO_DMAT;
    impl_ode_jac_x_xdot_u_z_out[3] = &df_dz;

    // impl_ode_hess
    // INPUT: impl_ode_hess
    ext_fun_arg_t impl_ode_hess_type_in[6];
    void *impl_ode_hess_in[6];

    struct blasfeo_dvec_args impl_ode_hess_lambda_in;

    impl_ode_hess_type_in[0] = BLASFEO_DVEC;       // xt
    impl_ode_hess_in[0] = xt;                      // 1st input is always xt
    impl_ode_hess_type_in[1] = BLASFEO_DVEC_ARGS;  // k_i
    impl_ode_hess_in[1] =  &impl_ode_xdot_in;      // 2nd input is part of K[ss]
    impl_ode_hess_type_in[2] = COLMAJ;             // u
    impl_ode_hess_in[2] = u;                       // 3rd input is u (always)
    impl_ode_hess_type_in[3] = BLASFEO_DVEC_ARGS;  // z_i
    impl_ode_hess_in[3] = &impl_ode_z_in;          // 4th input is part of Z[ss]
    impl_ode_hess_type_in[4] = BLASFEO_DVEC_ARGS;  // lambdaK component, direction
    impl_ode_hess_in[4] = &impl_ode_hess_lambda_in;     // 5th input is part of lambdaK[ss]
    impl_ode_hess_type_in[5] = BLASFEO_DMAT;       // dxkzu_dw0
    impl_ode_hess_in[5] = &dxkzu_dw0;              // 6th input is dxkzu_w0

    // OUTPUT
    ext_fun_arg_t impl_ode_hess_type_out[1];
    void *impl_ode_hess_out[1];
    impl_ode_hess_type_out[0] = BLASFEO_DMAT;
    impl_ode_hess_out[0] = &f_hess;


/* Initialize & Pack */
    // initialize
    double timing_ad = 0.0;
    double timing_la = 0.0;
    blasfeo_dvecse(nK, 0.0, lambdaK, 0);
    if (opts->sens_hess){
        blasfeo_dgese(nx + nu, nx + nu, 0.0, &Hess, 0, 0);
    }

    // pack
    blasfeo_pack_dvec(nx, in->x, xn, 0);
    blasfeo_pack_dmat(nx, nx + nu, in->S_forw, nx, S_forw, 0, 0);
    blasfeo_pack_dvec(nx + nu, in->S_adj, lambda, 0);

    // initialize integration variables
    //  state derivatives
    blasfeo_pack_dvec(nx, in->xdot, K, 0);
    //  algebraic variables
    blasfeo_pack_dvec(nz, in->z, K, nx);

    // TODO(dimitris, FreyJo): implement NF (number of forward sensis) properly, instead of nx+nu?

/************************************************
* Forward Sweep 
*       - (simulation & forward sensitivities)
************************************************/
    // set input for forward sweep
    impl_ode_xdot_in.x = K;
    impl_ode_z_in.x = K;

    // start the loop
    acados_tic(&timer);
    /* decide whether results from forward sensitivity propagation are stored,
        or if memory has to be reused --> set pointers accordingly */

    if ( opts->sens_adj || opts->sens_hess )  // store current xn
        assert(1 && "opts->sens_adj =1 and opts->sens_hess = 1  not implemented yet!");

    for (int iter = 0; iter < newton_iter; iter++)
    {
        if ((opts->jac_reuse && (iter == 0)) || (!opts->jac_reuse))
        {
            // if new jacobian gets computed, initialize dG_dK with zeros
            blasfeo_dgese(nK, nK, 0.0, dG_dK, 0, 0);
        }

        // take x(n); copy a strvec into a strvec
        blasfeo_dveccp(nx, xn, 0, xt, 0);

        impl_ode_xdot_in.xi = 0;  // use k_i of K = (k_1 ,z_1)
        impl_ode_z_in.xi    = nx;
                                      // use z_i of K = (k_1, z_1)
        impl_ode_res_out.xi = 0;  // store output in this position of rG

        // compute the residual of implicit ode at time t_ii
        if ((opts->jac_reuse && (iter == 0)) || (!opts->jac_reuse))
        {   // evaluate the ode function & jacobian w.r.t. x, xdot;
            //    &  compute jacobian dG_dK;
            acados_tic(&timer_ad);
            model->impl_ode_fun_jac_x_xdot_z->evaluate(
                model->impl_ode_fun_jac_x_xdot_z, impl_ode_type_in, impl_ode_in,
                impl_ode_fun_jac_x_xdot_z_type_out, impl_ode_fun_jac_x_xdot_z_out);
            timing_ad += acados_toc(&timer_ad);

            // compute the blocks of dG_dK
            blasfeo_dgead(nx + nz, nx, 1, &df_dxdot, 0, 0,
                          dG_dK, 0, 0);
            blasfeo_dgead(nx + nz, nz, 1, &df_dz,    0, 0,
                          dG_dK, 0, 0);
        }
        else // only eval function (without jacobian)
        {
            acados_tic(&timer_ad);
            model->impl_ode_fun->evaluate(model->impl_ode_fun, impl_ode_type_in,
                                          impl_ode_in, impl_ode_fun_type_out,
                                          impl_ode_fun_out);
            timing_ad += acados_toc(&timer_ad);
        }

        acados_tic(&timer_la);
        // DGETRF computes an LU factorization of a general M-by-N matrix A
        // using partial pivoting with row interchanges.
        // printf("dG_dK = (IRK) \n");
        // blasfeo_print_exp_dmat((nz+nx) *ns, (nz+nx) *ns, dG_dK, 0, 0);
        if ((opts->jac_reuse && (iter == 0)) || (!opts->jac_reuse))
        {
            blasfeo_dgetrf_rowpivot(nK, nK, dG_dK, 0, 0, dG_dK, 0, 0, ipiv);
        }

        // permute also the r.h.s
        blasfeo_dvecpe(nK, ipiv, rG, 0);

        // solve dG_dK * y = rG, dG_dK on the (l)eft, (l)ower-trian, (n)o-trans
        //                    (u)nit trian
        blasfeo_dtrsv_lnu(nK, dG_dK, 0, 0, rG, 0, rG, 0);

        // solve dG_dK * x = rG, dG_dK on the (l)eft, (u)pper-trian, (n)o-trans
        //                    (n)o unit trian , and store x in rG
        blasfeo_dtrsv_unn(nK, dG_dK, 0, 0, rG, 0, rG, 0);

        timing_la += acados_toc(&timer_la);

        // scale and add a generic strmat into a generic strmat // K = K - rG, where rG is
        // [DeltaK, DeltaZ]
        blasfeo_daxpy(nK, -1.0, rG, 0, K, 0, K, 0);
    }

    if ( opts->sens_adj || opts->sens_hess )
    {
        assert(1 && "opts->sens_adj =1 and opts->sens_hess = 1  not implemented yet!");
    }
    

    // TODO(andrea): these sensitivities are not necessary, right?
    /* evaluate forward sensitivities */
    // blasfeo_dgese(nK, nK, 0.0, dG_dK, 0, 0);
    //              // initialize dG_dK with zeros
    // // evaluate dG_dK(xn,Kn)

    // impl_ode_xdot_in.xi = 0;  // use k_i of K = (k_1, z_1)
    // impl_ode_z_in.xi    = nx;
    //                                 // use z_i of K = (k_1, z_1)
    // blasfeo_dveccp(nx, xn, 0, xt, 0);

    // acados_tic(&timer_ad);
    // model->impl_ode_jac_x_xdot_u_z->evaluate(
    //     model->impl_ode_jac_x_xdot_u_z, impl_ode_type_in, impl_ode_in,
    //     impl_ode_jac_x_xdot_u_z_type_out, impl_ode_jac_x_xdot_u_z_out);
    // timing_ad += acados_toc(&timer_ad);

    // blasfeo_dgecp(nx + nz, nx, &df_dx, 0, 0, dG_dxu, 0, 0);
    // blasfeo_dgecp(nx + nz, nu, &df_du, 0, 0, dG_dxu, 0, nx);

    // // compute the blocks of dG_dK
    // blasfeo_dgead(nx + nz, nx, 1, &df_dxdot, 0, 0,
    //                 dG_dK, 0, 0);
    // blasfeo_dgead(nx + nz, nz, 1, &df_dz,    0, 0,
    //                 dG_dK, 0, 0);

    // // factorize dG_dK
    // acados_tic(&timer_la);
    // blasfeo_dgetrf_rowpivot(nK, nK, dG_dK, 0, 0, dG_dK, 0, 0, ipiv_ss);
    // timing_la += acados_toc(&timer_la);
    // /* obtain dK_dxu */
    // // set up right hand side
    // // dK_dw = 0 * dK_dw + 1 * dG_dx * S_forw_old
    // blasfeo_dgemm_nn(nK, nx + nu, nx, 1.0, dG_dxu, 0, 0, S_forw, 0,
    //                      0, 0.0, dK_dxu, 0, 0, dK_dxu, 0, 0);
    // // printf("dG_dxu = \n");
    // // blasfeo_print_exp_dmat(nx + nz, nx+nu, dG_dxu, 0, 0);
    // // dK_du = dK_du + 1 * dG_du
    // blasfeo_dgead(nK, nu, 1.0, dG_dxu, 0, nx, dK_dxu, 0, nx);
    // // solve linear system
    // acados_tic(&timer_la);
    // blasfeo_drowpe(nK, ipiv_ss, dK_dxu);
    // blasfeo_dtrsm_llnu(nK, nx + nu, 1.0, dG_dK, 0, 0, dK_dxu, 0, 0, dK_dxu, 0, 0);
    // blasfeo_dtrsm_lunn(nK, nx + nu, 1.0, dG_dK, 0, 0, dK_dxu, 0, 0, dK_dxu, 0, 0);
    // timing_la += acados_toc(&timer_la);

    // // printf("dK_dxu (solved) = (IRK, ss = %d) \n", ss);
    // // blasfeo_print_exp_dmat(nK, nx + nu, dK_dxu, 0, 0);

    // // update forward sensitivity
    // for (int jj = 0; jj < ns; jj++)
    //     blasfeo_dgead(nx, nx + nu, 1.0, dK_dxu, 0, 0,
    //             S_forw, 0, 0);

    /* algebraic variables output and corresponding sensitivity propagation */
    // generate z output
    
    // copy corresponding values to out->zn
    blasfeo_unpack_dvec(nz, K, nx, out->zn);

    // set input for impl_ode
    impl_ode_type_in[0] = COLMAJ;
    impl_ode_type_in[1] = BLASFEO_DVEC;
    impl_ode_type_in[3] = COLMAJ;
    impl_ode_in[0] = in->x;  // 1st input is always xn
    impl_ode_in[1] = K;
    impl_ode_in[3] = &out->zn[0];

    // eval jacobians at interpolated values
    acados_tic(&timer_ad);
    model->impl_ode_jac_x_xdot_u_z->evaluate(
            model->impl_ode_jac_x_xdot_u_z, impl_ode_type_in, impl_ode_in,
            impl_ode_jac_x_xdot_u_z_type_out, impl_ode_jac_x_xdot_u_z_out);
    timing_ad += acados_toc(&timer_ad);

    // set up df_dxdotz
    blasfeo_dgecp(nx + nz, nx, &df_dxdot, 0, 0, &df_dxdotz, 0, 0);
    blasfeo_dgecp(nx + nz, nz, &df_dz,    0, 0, &df_dxdotz, 0, nx);
    // set up right hand side dk0_dxu
    blasfeo_dgecp(nx + nz, nx, &df_dx, 0, 0, &dk0_dxu, 0, 0);
    blasfeo_dgecp(nx + nz, nu, &df_du, 0, 0, &dk0_dxu, 0, nx);

    // solve linear system
    acados_tic(&timer_la);
    blasfeo_dgetrf_rowpivot(nx + nz, nx + nz, &df_dxdotz, 0, 0, &df_dxdotz, 0, 0,
                                                                ipiv_one_stage);
    blasfeo_drowpe(nx + nz, ipiv_one_stage, &dk0_dxu);
    blasfeo_dtrsm_llnu(nx + nz, nx + nu, 1.0, &df_dxdotz, 0, 0,
                       &dk0_dxu, 0, 0, &dk0_dxu, 0, 0);
    blasfeo_dtrsm_lunn(nx + nz, nx + nu, 1.0, &df_dxdotz, 0, 0,
                       &dk0_dxu, 0, 0, &dk0_dxu, 0, 0);
    timing_la += acados_toc(&timer_la);

    // solution has different sign
    blasfeo_dgesc(nx + nz, nx + nu, -1.0, &dk0_dxu, 0, 0);

    // extract output
    blasfeo_unpack_dmat(nz, nx + nu, &dk0_dxu, nx, 0, S_algebraic, nz);

    // Reset impl_ode inputs
    impl_ode_type_in[0] = BLASFEO_DVEC;       // xt
    impl_ode_type_in[1] = BLASFEO_DVEC_ARGS;  // k_i
    impl_ode_type_in[3] = BLASFEO_DVEC_ARGS;  // z_i
    impl_ode_in[0] = xt;  // 1st input is always xt
    impl_ode_in[1] = &impl_ode_xdot_in;
    impl_ode_in[3] = &impl_ode_z_in;     // 4th input is part of Z[ss]

    // extract results from forward sweep to output
    blasfeo_unpack_dvec(nx, xn, 0, x_out);
    if  ( opts->sens_forw || opts->sens_hess )
        blasfeo_unpack_dmat(nx, nx + nu, S_forw, 0, 0, S_forw_out, nx);

    if ( opts->sens_adj  || opts->sens_hess )
    {
        assert(1 && "opts->sens_adj =1 and opts->sens_hess = 1  not implemented yet!");
    }  // end if ( opts->sens_adj  || opts->sens_hess )

    out->info->CPUtime = acados_toc(&timer);
    out->info->LAtime =
        timing_la;  // note: this is the time for factorization and solving the linear systems
    out->info->ADtime = timing_ad;

    return ACADOS_SUCCESS;
}

void sim_algebraic_solver_config_initialize_default(void *config_)
{
    sim_solver_config *config = config_;

    config->evaluate = &sim_algebraic_solver;
    config->precompute = &sim_algebraic_solver_precompute;
    config->opts_calculate_size = &sim_algebraic_solver_opts_calculate_size;
    config->opts_assign = &sim_algebraic_solver_opts_assign;
    config->opts_initialize_default = &sim_algebraic_solver_opts_initialize_default;
    config->opts_update = &sim_algebraic_solver_opts_update;
    config->opts_set = &sim_algebraic_solver_opts_set;
    config->memory_calculate_size = &sim_algebraic_solver_memory_calculate_size;
    config->memory_assign = &sim_algebraic_solver_memory_assign;
    config->workspace_calculate_size = &sim_algebraic_solver_workspace_calculate_size;
    config->model_calculate_size = &sim_algebraic_solver_model_calculate_size;
    config->model_assign = &sim_algebraic_solver_model_assign;
    config->model_set_function = &sim_algebraic_solver_model_set_function;
    config->dims_calculate_size = &sim_algebraic_solver_dims_calculate_size;
    config->dims_assign = &sim_algebraic_solver_dims_assign;
    config->dims_set = &sim_algebraic_solver_dims_set;
    config->get_nx = &sim_algebraic_solver_get_nx;
    config->get_nu = &sim_algebraic_solver_get_nu;
    config->get_nz = &sim_algebraic_solver_get_nz;
    return;
}
