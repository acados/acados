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

#include "acados/sim/sim_irk_integrator.h"

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

int sim_irk_dims_calculate_size()
{
    int size = sizeof(sim_irk_dims);

    return size;
}

void *sim_irk_dims_assign(void *config_, void *raw_memory)
{
    char *c_ptr = raw_memory;

    sim_irk_dims *dims = (sim_irk_dims *) c_ptr;
    c_ptr += sizeof(sim_irk_dims);

    dims->nx = 0;
    dims->nu = 0;
    dims->nz = 0;

    assert((char *) raw_memory + sim_irk_dims_calculate_size() >= c_ptr);

    return dims;
}

void sim_irk_set_nx(void *dims_, int nx)
{
    sim_irk_dims *dims = (sim_irk_dims *) dims_;
    dims->nx = nx;
}

void sim_irk_set_nu(void *dims_, int nu)
{
    sim_irk_dims *dims = (sim_irk_dims *) dims_;
    dims->nu = nu;
}

void sim_irk_set_nz(void *dims_, int nz)
{
    sim_irk_dims *dims = (sim_irk_dims *) dims_;
    dims->nz = nz;
}

void sim_irk_get_nx(void *dims_, int *nx)
{
    sim_irk_dims *dims = (sim_irk_dims *) dims_;
    *nx = dims->nx;
}

void sim_irk_get_nu(void *dims_, int *nu)
{
    sim_irk_dims *dims = (sim_irk_dims *) dims_;
    *nu = dims->nu;
}

void sim_irk_get_nz(void *dims_, int *nz)
{
    sim_irk_dims *dims = (sim_irk_dims *) dims_;
    *nz = dims->nz;
}

/************************************************
 * model
 ************************************************/

int sim_irk_model_calculate_size(void *config, void *dims)
{
    int size = 0;

    size += sizeof(irk_model);

    return size;
}

void *sim_irk_model_assign(void *config, void *dims, void *raw_memory)
{
    char *c_ptr = (char *) raw_memory;

    irk_model *data = (irk_model *) c_ptr;
    c_ptr += sizeof(irk_model);

    assert((char *) raw_memory + sim_irk_model_calculate_size(config, dims) >= c_ptr);

    return data;
}

int sim_irk_model_set_function(void *model_, sim_function_t fun_type, void *fun)
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
        default:
            return ACADOS_FAILURE;
    }
    return ACADOS_SUCCESS;
}

/************************************************
 * opts
 ************************************************/

int sim_irk_opts_calculate_size(void *config_, void *dims)
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

void *sim_irk_opts_assign(void *config_, void *dims, void *raw_memory)
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

    assert((char *) raw_memory + sim_irk_opts_calculate_size(config_, dims) >= c_ptr);

    return (void *) opts;
}

void sim_irk_opts_initialize_default(void *config_, void *dims_, void *opts_)
{
    sim_irk_dims *dims = (sim_irk_dims *) dims_;
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
    opts->num_steps = 2;
    opts->num_forw_sens = dims->nx + dims->nu;
    opts->sens_forw = true;
    opts->sens_adj = false;
    opts->sens_hess = false;
    opts->jac_reuse = true;

    opts->output_z = false;
    opts->sens_algebraic = false;

    return;
}

void sim_irk_opts_update(void *config_, void *dims, void *opts_)
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

/************************************************
 * memory
 ************************************************/

int sim_irk_memory_calculate_size(void *config, void *dims, void *opts_) { return 0; }
void *sim_irk_memory_assign(void *config, void *dims, void *opts_, void *raw_memory)
{
    return NULL;
}

/************************************************
 * workspace
 ************************************************/

int sim_irk_workspace_calculate_size(void *config_, void *dims_, void *opts_)
{
    sim_irk_dims *dims = (sim_irk_dims *) dims_;
    sim_rk_opts *opts = opts_;

    int ns = opts->ns;

    int nx = dims->nx;
    int nu = dims->nu;
    int nz = dims->nz;

    int steps = opts->num_steps;

    int size = sizeof(sim_irk_workspace);

    size += 4 * sizeof(struct blasfeo_dmat);  // dG_dK, dG_dxu, dK_dxu, S_forw
    // size += steps*sizeof(struct blasfeo_dmat); // JG_traj

    size += 4 * sizeof(struct blasfeo_dvec);          // rG, K, xt, xn
    size += 2 * sizeof(struct blasfeo_dvec);          // lambda,lambdaK

    if (opts->sens_adj){
        size += 3 * steps * sizeof(struct blasfeo_dvec);  // **xn_traj, **K_traj, Ztraj
    }

    size += blasfeo_memsize_dmat((nx+nz) * ns, (nx+nz) * ns);      // dG_dK
    size += 2 * blasfeo_memsize_dmat((nx+nz) * ns, nx + nu);  // dG_dxu, dK_dxu
    size += blasfeo_memsize_dmat(nx, nx + nu);           // S_forw
    // size += steps * blasfeo_memsize_dmat(nx*ns, nx*ns); // for JG_traj

    size += blasfeo_memsize_dvec((nx + nz) * ns);   // K
    size += blasfeo_memsize_dvec((nx + nz) * ns);   // rG
    size += 3 * blasfeo_memsize_dvec(nx);           // xt, xn, xtdot
    size += blasfeo_memsize_dvec(nx + nu);          // lambda
    size += blasfeo_memsize_dvec((nx + nz) * ns);   // lambdaK

    if (opts->sens_adj){
        size += steps * blasfeo_memsize_dvec(nx);       // for xn_traj
        size += steps * blasfeo_memsize_dvec((nx + nz) * ns);  // for K_traj
    }

    size += (nx + nz) * ns * sizeof(int);  // ipiv
    size += (nx + nz) * sizeof(int);  // ipiv_one_stage
    size += ns * sizeof(double);  // Z_work

    size += 2 * blasfeo_memsize_dmat(nx + nz, nx);  // df_dx, df_dxdot
    size += blasfeo_memsize_dmat(nx + nz, nu);      // df_du
    size += blasfeo_memsize_dmat(nx + nz, nz);      // df_dz

    if (opts->sens_algebraic){
        size += blasfeo_memsize_dmat(nx + nz, nx + nz);  // df_dxdotz
    }

    make_int_multiple_of(64, &size);
    size += 1 * 64;

    return size;
}

static void *sim_irk_workspace_cast(void *config_, void *dims_, void *opts_, void *raw_memory)
{
    sim_rk_opts *opts = opts_;
    sim_irk_dims *dims = (sim_irk_dims *) dims_;

    int ns = opts->ns;

    int nx = dims->nx;
    int nu = dims->nu;
    int nz = dims->nz;

    int steps = opts->num_steps;

    char *c_ptr = (char *) raw_memory;

    sim_irk_workspace *workspace = (sim_irk_workspace *) c_ptr;
    c_ptr += sizeof(sim_irk_workspace);

    // assign_and_advance_blasfeo_dmat_structs(steps, &workspace->JG_traj, &c_ptr);

    workspace->dG_dK = (struct blasfeo_dmat *) c_ptr;
    c_ptr += sizeof(struct blasfeo_dmat);

    workspace->dG_dxu = (struct blasfeo_dmat *) c_ptr;
    c_ptr += sizeof(struct blasfeo_dmat);

    workspace->dK_dxu = (struct blasfeo_dmat *) c_ptr;
    c_ptr += sizeof(struct blasfeo_dmat);

    workspace->S_forw = (struct blasfeo_dmat *) c_ptr;
    c_ptr += sizeof(struct blasfeo_dmat);

    if (opts->sens_adj){
        assign_and_advance_blasfeo_dvec_structs(steps, &workspace->xn_traj, &c_ptr);
        assign_and_advance_blasfeo_dvec_structs(steps, &workspace->K_traj, &c_ptr);
    }

    workspace->rG = (struct blasfeo_dvec *) c_ptr;
    c_ptr += sizeof(struct blasfeo_dvec);

    assign_and_advance_blasfeo_dvec_structs(1, &workspace->K, &c_ptr);

    workspace->xt = (struct blasfeo_dvec *) c_ptr;
    c_ptr += sizeof(struct blasfeo_dvec);

    workspace->xn = (struct blasfeo_dvec *) c_ptr;
    c_ptr += sizeof(struct blasfeo_dvec);

    workspace->lambda = (struct blasfeo_dvec *) c_ptr;
    c_ptr += sizeof(struct blasfeo_dvec);

    workspace->lambdaK = (struct blasfeo_dvec *) c_ptr;
    c_ptr += sizeof(struct blasfeo_dvec);

    align_char_to(64, &c_ptr);

    assign_and_advance_blasfeo_dmat_mem((nx + nz) * ns, (nx+nz) * ns, workspace->dG_dK, &c_ptr);
    assign_and_advance_blasfeo_dmat_mem((nx + nz) * ns, nx + nu, workspace->dG_dxu, &c_ptr);
    assign_and_advance_blasfeo_dmat_mem((nx + nz) * ns, nx + nu, workspace->dK_dxu, &c_ptr);
    assign_and_advance_blasfeo_dmat_mem(nx, nx + nu, workspace->S_forw, &c_ptr);
    // for (int i=0;i<steps;i++){
    //     assign_and_advance_blasfeo_dmat_mem(nx*ns, nx*ns, &workspace->JG_traj[i], &c_ptr);
    // }

    assign_and_advance_blasfeo_dmat_mem(nx + nz, nx, &workspace->df_dx, &c_ptr);
    assign_and_advance_blasfeo_dmat_mem(nx + nz, nx, &workspace->df_dxdot, &c_ptr);
    assign_and_advance_blasfeo_dmat_mem(nx + nz, nu, &workspace->df_du, &c_ptr);
    assign_and_advance_blasfeo_dmat_mem(nx + nz, nz, &workspace->df_dz, &c_ptr);

    if (opts->sens_algebraic){
        assign_and_advance_blasfeo_dmat_mem(nx + nz, nx + nz, &workspace->df_dxdotz, &c_ptr);
    }

    assign_and_advance_blasfeo_dvec_mem((nx + nz) * ns, workspace->rG, &c_ptr);
    assign_and_advance_blasfeo_dvec_mem((nx + nz) * ns, workspace->K, &c_ptr);
    assign_and_advance_blasfeo_dvec_mem(nx, workspace->xt, &c_ptr);
    assign_and_advance_blasfeo_dvec_mem(nx, workspace->xn, &c_ptr);
    assign_and_advance_blasfeo_dvec_mem(nx, &workspace->xtdot, &c_ptr);
    assign_and_advance_blasfeo_dvec_mem(nx + nu, workspace->lambda, &c_ptr);
    assign_and_advance_blasfeo_dvec_mem((nx + nz) * ns, workspace->lambdaK, &c_ptr);

    if (opts->sens_adj){
        for (int i = 0; i < steps; i++)
        {
            assign_and_advance_blasfeo_dvec_mem(nx, &workspace->xn_traj[i], &c_ptr);
            assign_and_advance_blasfeo_dvec_mem((nx + nz) * ns, &workspace->K_traj[i], &c_ptr);
        }
    }
    assign_and_advance_double(ns, &workspace->Z_work, &c_ptr);
    assign_and_advance_int((nx + nz) * ns, &workspace->ipiv, &c_ptr);
    assign_and_advance_int((nx + nz), &workspace->ipiv_one_stage, &c_ptr);

    // printf("\npointer moved - size calculated = %d bytes\n", c_ptr- (char*)raw_memory -
    // sim_irk_calculate_workspace_size(dims, opts_));

    assert((char *) raw_memory + sim_irk_workspace_calculate_size(config_, dims, opts_) >= c_ptr);

    return (void *) workspace;
}

/************************************************
 * integrator
 ************************************************/

int sim_irk(void *config_, sim_in *in, sim_out *out, void *opts_, void *mem_, void *work_)
{
    /* INITIALIZE variables/pointers */
    sim_solver_config *config = config_;
    sim_rk_opts *opts = opts_;

    assert(opts->ns == opts->tableau_size && "the Butcher tableau size does not match ns");

    int ns = opts->ns;

    void *dims_ = in->dims;
    sim_irk_dims *dims = (sim_irk_dims *) dims_;
    sim_irk_workspace *workspace =
        (sim_irk_workspace *) sim_irk_workspace_cast(config, dims, opts, work_);

    int ii, jj, iter, ss;
    double a;

    int nx = dims->nx;
    int nu = dims->nu;
    int nz = dims->nz;

    int nK = (nx + nz) * ns;

    double *x = in->x;
    double *u = in->u;
    double *S_forw_in = in->S_forw;

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

    struct blasfeo_dmat df_dxdotz = workspace->df_dxdotz;

    // for adjoint
    struct blasfeo_dvec *lambda = workspace->lambda;
    struct blasfeo_dvec *lambdaK = workspace->lambdaK;
    struct blasfeo_dvec *xn_traj = workspace->xn_traj;
    struct blasfeo_dvec *K_traj = workspace->K_traj;
    // struct blasfeo_dmat *JG_traj = workspace->JG_traj;

    double *x_out = out->xn;
    double *S_forw_out = out->S_forw;
    double *S_adj_out = out->S_adj;
    double *S_algebraic = out->S_algebraic;

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

    irk_model *model = in->model;

    acados_timer timer, timer_ad, timer_la;
    double timing_ad = 0.0;
    double timing_la = 0.0;

    // initialize
    blasfeo_dgese(nK, nK, 0.0, dG_dK, 0, 0);
    blasfeo_dgese(nK, nx + nu, 0.0, dG_dxu, 0, 0);
    blasfeo_dgese(nK, nx + nu, 0.0, dK_dxu, 0, 0);

    blasfeo_dvecse(nK, 0.0, rG, 0);
    blasfeo_dvecse(nK, 0.0, K, 0);
    blasfeo_dvecse(nK, 0.0, lambdaK, 0);

    // pack
    blasfeo_pack_dvec(nx, x, xn, 0);
    blasfeo_pack_dmat(nx, nx + nu, S_forw_in, nx, S_forw, 0, 0);
    blasfeo_pack_dvec(nx + nu, in->S_adj, lambda, 0);

    // TODO(dimitris, FreyJo): implement NF (number of forward sensis) properly, instead of nx+nu?

    // start the loop
    acados_tic(&timer);
    for (ss = 0; ss < num_steps; ss++)
    {
        // obtain Kn
        // TODO(giaf): add exit condition on residuals ???
        // inf_norm_K = 1.0;
        // for(iter=0; inf_norm_K>tol_inf_norm_K & iter<newton_iter; iter++)
        impl_ode_xdot_in.x = K;
        impl_ode_z_in.x = K;
        for (iter = 0; iter < newton_iter; iter++)
        {
            if ((opts->jac_reuse && (ss == 0) && (iter == 0)) || (!opts->jac_reuse))
            {
                blasfeo_dgese(nK, nK, 0.0, dG_dK, 0,
                              0);  // if new jacobian gets computed, initialize dG_dK with zeros
            }

            if (opts->sens_adj)
            {
                blasfeo_dveccp(nx, xn, 0, &xn_traj[ss], 0);
            }

            for (ii = 0; ii < ns; ii++)
            {  // ii-th row of tableau
                // take x(n); copy a strvec into a strvec
                blasfeo_dveccp(nx, xn, 0, xt, 0);

                for (jj = 0; jj < ns; jj++)
                {  // jj-th col of tableau
                    a = A_mat[ii + ns * jj];
                    if (a != 0)
                    {  // xt = xt + T_int * a[i,j]*K_j
                        a *= step;
                        blasfeo_daxpy(nx, a, K, jj * nx, xt, 0, xt, 0);
                    }
                }
                impl_ode_xdot_in.xi = ii * nx;  // use k_i of K = (k_1,..., k_{ns},z_1,..., z_{ns})
                impl_ode_z_in.xi    = ns * nx + ii * nz;
                                              // use z_i of K = (k_1,..., k_{ns},z_1,..., z_{ns})
                impl_ode_res_out.xi = ii * (nx + nz);  // store output in this posistion of rG

                // compute the residual of implicit ode at time t_ii
                if ((opts->jac_reuse && (ss == 0) && (iter == 0)) || (!opts->jac_reuse))
                {  // evaluate the ode function & jacobian w.r.t. x, xdot; compute jacobian dG_dK;
                    acados_tic(&timer_ad);
                    model->impl_ode_fun_jac_x_xdot_z->evaluate(
                        model->impl_ode_fun_jac_x_xdot_z, impl_ode_type_in, impl_ode_in,
                        impl_ode_fun_jac_x_xdot_z_type_out, impl_ode_fun_jac_x_xdot_z_out);
                    timing_ad += acados_toc(&timer_ad);

                    // compute the blocks of dG_dK
                    for (jj = 0; jj < ns; jj++)
                    {  // compute the block (ii,jj)th block of dG_dK
                        a = A_mat[ii + ns * jj];
                        if (a != 0)
                        {
                            a *= step;
                            blasfeo_dgead(nx + nz, nx, a, &df_dx, 0, 0,
                                              dG_dK, ii * (nx + nz), jj * nx);
                        }
                        if (jj == ii)
                        {
                            blasfeo_dgead(nx + nz, nx, 1, &df_dxdot, 0, 0,
                                          dG_dK, ii * (nx + nz), jj * nx);
                            blasfeo_dgead(nx + nz, nz, 1, &df_dz,    0, 0,
                                          dG_dK, ii * (nx + nz), (nx * ns) + jj * nz);
                        }
                    }  // end jj
                }
                else // only eval function (without jacobian)
                {
                    acados_tic(&timer_ad);
                    model->impl_ode_fun->evaluate(model->impl_ode_fun, impl_ode_type_in,
                                                  impl_ode_in, impl_ode_fun_type_out,
                                                  impl_ode_fun_out);
                    timing_ad += acados_toc(&timer_ad);
                }
            }  // end ii

            acados_tic(&timer_la);
            // DGETRF computes an LU factorization of a general M-by-N matrix A
            // using partial pivoting with row interchanges.
            // printf("dG_dK = (IRK) \n");
            // blasfeo_print_exp_dmat((nz+nx) *ns, (nz+nx) *ns, dG_dK, 0, 0);
            if ((opts->jac_reuse && (ss == 0) && (iter == 0)) || (!opts->jac_reuse))
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
            // inf norm of K
            //       blasfeo_dvecnrm_inf(nx*ns, K, 0, &inf_norm_K);
        }

        if (opts->sens_adj)
        {
            blasfeo_dveccp(nK, K, 0, &K_traj[ss], 0);
        }

        // evaluate forward sensitivities
        if (opts->sens_forw)
        {
            blasfeo_dgese(nK, nK, 0.0, dG_dK, 0, 0);
                         // initialize dG_dK with zeros
            // evaluate dG_dK(xn,Kn)
            for (ii = 0; ii < ns; ii++)
            {
                impl_ode_xdot_in.xi = ii * nx;  // use k_i of K = (k_1,..., k_{ns},z_1,..., z_{ns})
                impl_ode_z_in.xi    = ns * nx + ii * nz;
                                                // use z_i of K = (k_1,..., k_{ns},z_1,..., z_{ns})
                blasfeo_dveccp(nx, xn, 0, xt, 0);

                for (jj = 0; jj < ns; jj++)
                {
                    a = A_mat[ii + ns * jj];
                    if (a != 0)
                    {  // xt = xt + T_int * a[i,j]*K_j
                        a *= step;
                        blasfeo_daxpy(nx, a, K, jj * nx, xt, 0, xt, 0);
                    }
                }

                acados_tic(&timer_ad);
                model->impl_ode_jac_x_xdot_u_z->evaluate(
                    model->impl_ode_jac_x_xdot_u_z, impl_ode_type_in, impl_ode_in,
                    impl_ode_jac_x_xdot_u_z_type_out, impl_ode_jac_x_xdot_u_z_out);
                timing_ad += acados_toc(&timer_ad);

                blasfeo_dgecp(nx + nz, nx, &df_dx, 0, 0, dG_dxu, ii * (nx + nz), 0);
                blasfeo_dgecp(nx + nz, nu, &df_du, 0, 0, dG_dxu, ii * (nx + nz), nx);

                // compute the blocks of dG_dK
                for (jj = 0; jj < ns; jj++)
                {  // compute the block (ii,jj)th block of dG_dK
                    a = A_mat[ii + ns * jj];
                    if (a != 0)
                    {
                        a *= step;
                        blasfeo_dgead(nx + nz, nx, a, &df_dx, 0, 0,
                                            dG_dK, ii * (nx + nz), jj * nx);
                    }
                    if (jj == ii)
                    {
                        blasfeo_dgead(nx + nz, nx, 1, &df_dxdot, 0, 0,
                                        dG_dK, ii * (nx + nz), jj * nx);
                        blasfeo_dgead(nx + nz, nz, 1, &df_dz,    0, 0,
                                        dG_dK, ii * (nx + nz), (nx * ns) + jj * nz);
                    }
                }  // end jj
            }  // end ii

            // factorize dG_dK
            acados_tic(&timer_la);
            blasfeo_dgetrf_rowpivot(nK, nK, dG_dK, 0, 0, dG_dK, 0, 0, ipiv);
            timing_la += acados_toc(&timer_la);
            // NOTE: it is possible to store the factorization and permutation of dG_dK and reuse it
            // in the adjoint propagation, but as in common MPC schemes only one of those is needed,
            // this is omited for now - JG_traj was used for this;

            // obtain dK_dxu
            // TODO(all): add the option to use VDE instead of dgemm ???
            // dK_dw = 0 * dK_dw + 1 * dG_dx * S_forw
            blasfeo_dgemm_nn(nK, nx + nu, nx, 1.0, dG_dxu, 0, 0, S_forw, 0,
                                 0, 0.0, dK_dxu, 0, 0, dK_dxu, 0, 0);
            // dK_du = dK_du + 1 * dG_du
            blasfeo_dgead(nK, nu, 1.0, dG_dxu, 0, nx, dK_dxu, 0, nx);

            acados_tic(&timer_la);
            blasfeo_drowpe(nK, ipiv, dK_dxu);
            blasfeo_dtrsm_llnu(nK, nx + nu, 1.0, dG_dK, 0, 0,
                                 dK_dxu, 0, 0, dK_dxu, 0, 0);
            blasfeo_dtrsm_lunn(nK, nx + nu, 1.0, dG_dK, 0, 0,
                                 dK_dxu, 0, 0, dK_dxu, 0, 0);
            timing_la += acados_toc(&timer_la);

            // update forward sensitivity
            for (jj = 0; jj < ns; jj++)
                blasfeo_dgead(nx, nx + nu, -step * b_vec[jj], dK_dxu, jj * nx, 0, S_forw, 0, 0);
        }  // end if sens_forw

        // obtain x(n+1)
        for (ii = 0; ii < ns; ii++){
            blasfeo_daxpy(nx, step * b_vec[ii], K, ii * nx, xn, 0, xn, 0);
        }

        if (ss == 0)
        {
            // generate z output
            if ((opts->output_z || opts->sens_algebraic))
            {
                for (int ii = 0; ii < nz; ii++)
                {
                    for (int jj = 0; jj < ns; jj++)
                    {
                        Z_work[jj] = blasfeo_dvecex1(K, nx * ns + nz * jj + ii);
                                // copy values of z_ii in first step, into Z_work
                    }
                    neville_algorithm(0.0, ns - 1, opts->c_vec, Z_work, &out->zn[ii]);
                                // eval polynomial through (c_jj, z_jj) at 0.
                }
            }

            if (opts->sens_algebraic)  // generate S_algebraic
            {
                double interpolated_value;
                for (int ii = 0; ii < nx; ii++)
                {
                    for (int jj = 0; jj < ns; jj++)
                    {
                        Z_work[jj] = blasfeo_dvecex1(K, nx * jj + ii);
                            // copy values of k_ii in first step, into Z_work
                    }
                    neville_algorithm(0.0, ns - 1, opts->c_vec, Z_work, &interpolated_value);
                                // eval polynomial through (c_jj, k_jj) at 0.
                    blasfeo_pack_dvec(1, &interpolated_value, xtdot, ii);
                }
                // xtdot now contains x_dot(0)

                // set input for impl_ode
                impl_ode_type_in[0] = COLMAJ;
                impl_ode_type_in[1] = BLASFEO_DVEC;
                impl_ode_type_in[3] = COLMAJ;
                impl_ode_in[0] = in->x;  // 1st input is always xn
                impl_ode_in[1] = xtdot;
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
                // set up right hand side dK_dxu (only first nx+nz columns used here)
                blasfeo_dgecp(nx + nz, nx, &df_dx, 0, 0, dK_dxu, 0, 0);
                blasfeo_dgecp(nx + nz, nu, &df_du, 0, 0, dK_dxu, 0, nx);

                // solve linear system
                acados_tic(&timer_la);
                blasfeo_dgetrf_rowpivot(nx + nz, nx + nz, &df_dxdotz, 0, 0, &df_dxdotz, 0, 0,
                                                                            ipiv_one_stage);
                blasfeo_drowpe(nx + nz, ipiv_one_stage, dK_dxu);
                blasfeo_dtrsm_llnu(nx + nz, nx + nu, 1.0, &df_dxdotz, 0, 0,
                                   dK_dxu, 0, 0, dK_dxu, 0, 0);
                blasfeo_dtrsm_lunn(nx + nz, nx + nu, 1.0, &df_dxdotz, 0, 0,
                                   dK_dxu, 0, 0, dK_dxu, 0, 0);
                timing_la += acados_toc(&timer_la);

                // solution has different sign
                blasfeo_dgesc(nx + nz, nx + nu, -1.0, dK_dxu, 0, 0);

                // extract output
                blasfeo_unpack_dmat(nz, nx + nu, dK_dxu, nx, 0, S_algebraic, nz);

                // Reset impl_ode inputs
                impl_ode_type_in[0] = BLASFEO_DVEC;       // xt
                impl_ode_type_in[1] = BLASFEO_DVEC_ARGS;  // k_i
                impl_ode_type_in[3] = BLASFEO_DVEC_ARGS;  // z_i
                impl_ode_in[0] = xt;  // 1st input is always xt
                impl_ode_in[1] = &impl_ode_xdot_in;
                impl_ode_in[3] = &impl_ode_z_in;     // 4th input is part of Z[ss]
            }  // end if sens_algebraic
        }  //  end if (ss == 0)
    }  // end step loop (ss)



    // evaluate backwards
    if (opts->sens_adj)
    {
        for (ss = num_steps - 1; ss > -1; ss--)
        {
            impl_ode_xdot_in.x = &K_traj[ss];              // use K values of step ss
            impl_ode_z_in.x = &K_traj[ss];                 // use Z values of step ss

            blasfeo_dgese(nK, nK, 0.0, dG_dK, 0, 0);
                                                    // initialize dG_dK with zeros

            for (ii = 0; ii < ns; ii++)
            {
                impl_ode_xdot_in.xi = ii * nx;  // use k_i of K = (k_1,..., k_{ns},z_1,..., z_{ns})
                impl_ode_z_in.xi    = ns * nx + ii * nz;
                            // use z_i of K = (k_1,..., k_{ns},z_1,..., z_{ns})

                // build stage value
                blasfeo_dveccp(nx, &xn_traj[ss], 0, xt, 0);
                for (jj = 0; jj < ns; jj++)
                {
                    a = A_mat[ii + ns * jj];
                    if (a != 0)
                    {
                        a *= step;
                        blasfeo_daxpy(nx, a, &K_traj[ss], jj * nx, xt, 0, xt, 0);
                    }
                }

                acados_tic(&timer_ad);
                model->impl_ode_jac_x_xdot_u_z->evaluate(
                    model->impl_ode_jac_x_xdot_u_z, impl_ode_type_in, impl_ode_in,
                    impl_ode_jac_x_xdot_u_z_type_out, impl_ode_jac_x_xdot_u_z_out);
                timing_ad += acados_toc(&timer_ad);

                blasfeo_dgecp(nx + nz, nx, &df_dx, 0, 0, dG_dxu, ii * (nx + nz), 0);
                blasfeo_dgecp(nx + nz, nu, &df_du, 0, 0, dG_dxu, ii * (nx + nz), nx);

                // compute the blocks of dG_dK
                for (jj = 0; jj < ns; jj++)
                {  // compute the block (ii,jj)th block of dG_dK
                    a = A_mat[ii + ns * jj];
                    if (a != 0)
                    {
                        a *= step;
                        blasfeo_dgead(nx + nz, nx, a, &df_dx, 0, 0,
                                            dG_dK, ii * (nx + nz), jj * nx);
                    }
                    if (jj == ii)
                    {
                        blasfeo_dgead(nx + nz, nx, 1, &df_dxdot, 0, 0,
                                        dG_dK, ii * (nx + nz), jj * nx);
                        blasfeo_dgead(nx + nz, nz, 1, &df_dz,    0, 0,
                                        dG_dK, ii * (nx + nz), (nx * ns) + jj * nz);
                    }
                }  // end jj
            }  // end ii

            // factorize dG_dK
            acados_tic(&timer_la);
            blasfeo_dgetrf_rowpivot(nK, nK, dG_dK, 0, 0, dG_dK, 0, 0, ipiv);
            timing_la += acados_toc(&timer_la);

            blasfeo_dvecse(nK, 0.0, lambdaK, 0);
            for (jj = 0; jj < ns; jj++)
                blasfeo_dveccpsc(nx, -step * b_vec[jj], lambda, 0, lambdaK, jj * nx);
                // lambdaK_jj = -step b_jj * lambda_x

            // obtain lambda_K
            acados_tic(&timer_la);
            blasfeo_dtrsv_utn(nK, dG_dK, 0, 0, lambdaK, 0, lambdaK, 0);
            blasfeo_dtrsv_ltu(nK, dG_dK, 0, 0, lambdaK, 0, lambdaK, 0);
            blasfeo_dvecpei(nK, ipiv, lambdaK, 0);
            timing_la += acados_toc(&timer_la);

            // lambda = 1 * lambda + 1 * dG_dxu' * lambdaK
            blasfeo_dgemv_t(nK, nx + nu, 1.0, dG_dxu, 0, 0, lambdaK, 0, 1.0, lambda,
                             0, lambda, 0);
        }
    }
    // extract output

    blasfeo_unpack_dvec(nx, xn, 0, x_out);

    blasfeo_unpack_dmat(nx, nx + nu, S_forw, 0, 0, S_forw_out, nx);

    blasfeo_unpack_dvec(nx + nu, lambda, 0, S_adj_out);

    out->info->CPUtime = acados_toc(&timer);
    out->info->LAtime =
        timing_la;  // note: this is the time for factorization and solving the linear systems
    out->info->ADtime = timing_ad;

    return 0;
}

void sim_irk_config_initialize_default(void *config_)
{
    sim_solver_config *config = config_;

    config->evaluate = &sim_irk;
    config->opts_calculate_size = &sim_irk_opts_calculate_size;
    config->opts_assign = &sim_irk_opts_assign;
    config->opts_initialize_default = &sim_irk_opts_initialize_default;
    config->opts_update = &sim_irk_opts_update;
    config->memory_calculate_size = &sim_irk_memory_calculate_size;
    config->memory_assign = &sim_irk_memory_assign;
    config->workspace_calculate_size = &sim_irk_workspace_calculate_size;
    config->model_calculate_size = &sim_irk_model_calculate_size;
    config->model_assign = &sim_irk_model_assign;
    config->model_set_function = &sim_irk_model_set_function;
    config->dims_calculate_size = &sim_irk_dims_calculate_size;
    config->dims_assign = &sim_irk_dims_assign;
    config->set_nx = &sim_irk_set_nx;
    config->set_nu = &sim_irk_set_nu;
    config->set_nz = &sim_irk_set_nz;
    config->get_nx = &sim_irk_get_nx;
    config->get_nu = &sim_irk_get_nu;
    config->get_nz = &sim_irk_get_nz;
    return;
}
