/*
 * Copyright (c) The acados authors.
 *
 * This file is part of acados.
 *
 * The 2-Clause BSD License
 */

#include "acados/ocp_nlp/ocp_nlp_dynamics_cont_with_cost.h"
#include "acados/ocp_nlp/ocp_nlp_cost_ls.h"
#include "acados/utils/mem.h"

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "blasfeo_d_aux.h"
#include "blasfeo_d_blas.h"

static void ocp_nlp_dynamics_cont_with_cost_dims_set(void *config_, void *dims_,
        const char *field, int *value)
{
    ocp_nlp_dynamics_cont_dims *dims = dims_;
    ocp_nlp_dynamics_config *config = config_;

    if (!strcmp(field, "nx"))
    {
        dims->nx = *value;
        int nx_rk = *value + 1;
        config->sim_solver->dims_set(config->sim_solver, dims->sim, "nx", &nx_rk);
    }
    else if (!strcmp(field, "nx1"))
    {
        dims->nx1 = *value;
    }
    else if (!strcmp(field, "nu1"))
    {
        dims->nu1 = *value;
    }
    else if (!strcmp(field, "nz"))
    {
        dims->nz = *value;
        config->sim_solver->dims_set(config->sim_solver, dims->sim, field, value);
    }
    else if (!strcmp(field, "nu"))
    {
        dims->nu = *value;
        config->sim_solver->dims_set(config->sim_solver, dims->sim, field, value);
    }
    else if (!strcmp(field, "np"))
    {
        dims->np = *value;
        config->sim_solver->dims_set(config->sim_solver, dims->sim, field, value);
    }
    else
    {
        config->sim_solver->dims_set(config->sim_solver, dims->sim, field, value);
    }
}

static void ocp_nlp_dynamics_cont_with_cost_dims_get(void *config_, void *dims_,
        const char *field, int *value)
{
    ocp_nlp_dynamics_cont_dims *dims = dims_;

    if (!strcmp(field, "nx"))
        *value = dims->nx;
    else if (!strcmp(field, "nx1"))
        *value = dims->nx1;
    else if (!strcmp(field, "nu"))
        *value = dims->nu;
    else if (!strcmp(field, "nu1"))
        *value = dims->nu1;
    else if (!strcmp(field, "np"))
        *value = dims->np;
    else
    {
        ocp_nlp_dynamics_config *config = config_;
        config->sim_solver->dims_get(config->sim_solver, dims->sim, field, value);
    }
}

static void ocp_nlp_dynamics_cont_with_cost_opts_initialize_default(void *config_,
        void *dims_, void *opts_)
{
    ocp_nlp_dynamics_cont_opts_initialize_default(config_, dims_, opts_);

    ocp_nlp_dynamics_config *config = config_;
    ocp_nlp_dynamics_cont_opts *opts = opts_;
    bool cost_computation = true;
    config->sim_solver->opts_set(config->sim_solver, opts->sim_solver,
            "cost_computation", &cost_computation);
}

static void ocp_nlp_dynamics_cont_with_cost_opts_set(void *config_, void *opts_,
        const char *field, void *value)
{
    if (!strcmp(field, "cost_computation"))
        return;

    ocp_nlp_dynamics_cont_opts_set(config_, opts_, field, value);
}

static void ocp_nlp_dynamics_cont_with_cost_opts_get(void *config_, void *opts_,
        const char *field, void *value)
{
    ocp_nlp_dynamics_config *config = config_;
    ocp_nlp_dynamics_cont_opts *opts = opts_;

    if (!strcmp(field, "cost_computation"))
    {
        *(int *) value = 1;
        return;
    }
    if (!strcmp(field, "compute_adj"))
        *(int *) value = opts->compute_adj;
    else if (!strcmp(field, "compute_hess"))
        *(int *) value = opts->compute_hess;
    else
        config->sim_solver->opts_get(config->sim_solver, opts->sim_solver, field, value);
}

static void ocp_nlp_dynamics_cont_with_cost_memory_set(void *config_, void *dims_,
        void *mem_, const char *field, void *value)
{
    ocp_nlp_dynamics_cont_memory *mem = mem_;

    if (!strcmp(field, "cost_capsule_ptr"))
    {
        mem->cost_capsule = value;
        return;
    }
    if (!strcmp(field, "ux_ptr"))
        mem->ux = value;
    else if (!strcmp(field, "ux1_ptr"))
        mem->ux1 = value;
    else if (!strcmp(field, "pi_ptr"))
        mem->pi = value;
    else if (!strcmp(field, "BAbt_ptr"))
        mem->BAbt = value;
    else if (!strcmp(field, "RSQrq_ptr"))
        mem->RSQrq = value;
    else if (!strcmp(field, "dzduxt_ptr"))
        mem->dzduxt = value;
    else if (!strcmp(field, "sim_guess"))
        mem->sim_guess = value;
    else if (!strcmp(field, "set_sim_guess"))
        mem->set_sim_guess = value;
    else if (!strcmp(field, "z_alg_ptr"))
        mem->z_alg = value;
    else if (!strcmp(field, "dyn_jac_p_global_ptr") ||
            !strcmp(field, "jac_lag_stat_p_global_ptr") ||
            !strcmp(field, "adj_lag_p_global_ptr") || !strcmp(field, "seed_ux_ptr") ||
            !strcmp(field, "seed_pi_ptr"))
        return;
    else
    {
        printf("\nerror: ocp_nlp_dynamics_cont_with_cost_memory_set: field %s not available\n", field);
        exit(1);
    }
}

static void ocp_nlp_dynamics_cont_with_cost_memory_get(void *config_, void *dims_,
        void *mem_, const char *field, void *value)
{
    ocp_nlp_dynamics_cont_memory *mem = mem_;
    ocp_nlp_dynamics_config *config = config_;
    ocp_nlp_dynamics_cont_dims *dims = dims_;

    if (!strcmp(field, "time_sim") || !strcmp(field, "time_sim_ad") ||
            !strcmp(field, "time_sim_la") || !strcmp(field, "S_p"))
    {
        config->sim_solver->memory_get(config->sim_solver, dims->sim, mem->sim_solver,
                field, value);
        return;
    }

    printf("\nerror: ocp_nlp_dynamics_cont_with_cost_memory_get: field %s not available\n", field);
    exit(1);
}

static void ocp_nlp_dynamics_cont_with_cost_set_sim_input(
        ocp_nlp_dynamics_cont_dims *dims, ocp_nlp_dynamics_cont_model *model,
        ocp_nlp_dynamics_cont_memory *mem, ocp_nlp_dynamics_cont_workspace *work)
{
    int nx = dims->nx;
    int nu = dims->nu;
    int nx_rk = nx + 1;

    work->sim_in->model = model->sim_model;
    work->sim_in->T = model->T;
    memset(work->sim_in->x, 0, nx_rk * sizeof(double));
    blasfeo_unpack_dvec(nu, mem->ux, 0, work->sim_in->u, 1);
    blasfeo_unpack_dvec(nx, mem->ux, nu, work->sim_in->x, 1);

    if (mem->set_sim_guess != NULL && mem->set_sim_guess[0])
    {
        /* The augmented cost state has the known initial value zero. */
        mem->set_sim_guess[0] = false;
    }
}

static void ocp_nlp_dynamics_cont_with_cost_set_identity_seed(
        ocp_nlp_dynamics_cont_dims *dims, ocp_nlp_dynamics_cont_workspace *work)
{
    int nx_rk = dims->nx + 1;
    int nu = dims->nu;

    for (int ii = 0; ii < nx_rk * (nx_rk + nu); ii++)
        work->sim_in->S_forw[ii] = 0.0;
    for (int ii = 0; ii < nx_rk; ii++)
        work->sim_in->S_forw[ii * (nx_rk + 1)] = 1.0;
    work->sim_in->identity_seed = true;
}

static void ocp_nlp_dynamics_cont_with_cost_cast_workspace(void *config_, void *dims_,
        void *opts_, void *work_, void *mem_)
{
    ocp_nlp_dynamics_config *config = config_;
    ocp_nlp_dynamics_cont_dims *dims = dims_;
    ocp_nlp_dynamics_cont_memory *mem = mem_;
    ocp_nlp_dynamics_cont_workspace *work = work_;
    char *c_ptr = (char *) work_;
    int nx = dims->nx;
    int nu = dims->nu;

    c_ptr += sizeof(ocp_nlp_dynamics_cont_workspace);
    sim_in_assign_and_advance(config->sim_solver, dims->sim, &work->sim_in, &c_ptr);
    work->sim_out = sim_out_assign(config->sim_solver, dims->sim, c_ptr);
    c_ptr += sim_out_calculate_size(config->sim_solver, dims->sim);
    work->sim_solver = c_ptr;
    c_ptr += mem->sim_workspace_size;
    align_char_to(64, &c_ptr);
    assign_and_advance_blasfeo_dmat_mem(nu + nx, nu + nx, &work->hess, &c_ptr);

    assert((char *) work + mem->workspace_size >= c_ptr);
}

static ocp_nlp_cost_ls_memory *ocp_nlp_dynamics_cont_with_cost_get_cost_memory(
        ocp_nlp_dynamics_cont_memory *mem)
{
    assert(mem->cost_capsule != NULL);
    ocp_nlp_cost_capsule *cost_capsule = mem->cost_capsule;
    return cost_capsule->memory;
}

void ocp_nlp_dynamics_cont_with_cost_update_qp_matrices(void *config_, void *dims_,
        void *model_, void *opts_, void *mem_, void *work_)
{
    ocp_nlp_dynamics_cont_with_cost_cast_workspace(config_, dims_, opts_, work_, mem_);

    ocp_nlp_dynamics_config *config = config_;
    ocp_nlp_dynamics_cont_dims *dims = dims_;
    ocp_nlp_dynamics_cont_model *model = model_;
    ocp_nlp_dynamics_cont_opts *opts = opts_;
    ocp_nlp_dynamics_cont_memory *mem = mem_;
    ocp_nlp_dynamics_cont_workspace *work = work_;
    ocp_nlp_cost_ls_memory *cost_memory = ocp_nlp_dynamics_cont_with_cost_get_cost_memory(mem);

    int nx = dims->nx;
    int nu = dims->nu;
    int nz = dims->nz;
    int nx1 = dims->nx1;
    int nu1 = dims->nu1;
    int nx_rk = nx + 1;
    int nf_rk = nx_rk + nu;

    ocp_nlp_dynamics_cont_with_cost_set_sim_input(dims, model, mem, work);
    ocp_nlp_dynamics_cont_with_cost_set_identity_seed(dims, work);

    for (int ii = 0; ii < nf_rk; ii++)
        work->sim_in->S_adj[ii] = 0.0;
    blasfeo_unpack_dvec(nx1, mem->pi, 0, work->sim_in->S_adj, 1);
    double cost_scaling;
    ocp_nlp_cost_capsule *cost_capsule = mem->cost_capsule;
    ocp_nlp_cost_config *cost_config = cost_capsule->config;
    cost_config->model_get(cost_config, cost_capsule->dims, cost_capsule->model,
            "scaling", &cost_scaling);
    work->sim_in->S_adj[nx] = cost_scaling;

    config->sim_solver->evaluate(config->sim_solver, work->sim_in, work->sim_out,
            opts->sim_solver, mem->sim_solver, work->sim_solver);

    blasfeo_pack_tran_dmat(nx1, nu, work->sim_out->S_forw + nx_rk * nx_rk,
            nx_rk, mem->BAbt, 0, 0);
    blasfeo_pack_tran_dmat(nx1, nx, work->sim_out->S_forw, nx_rk,
            mem->BAbt, nu, 0);

    blasfeo_pack_tran_dmat(nz, nu, work->sim_out->S_algebraic + nz * nx_rk,
            nz, mem->dzduxt, 0, 0);
    blasfeo_pack_tran_dmat(nz, nx, work->sim_out->S_algebraic, nz,
            mem->dzduxt, nu, 0);
    blasfeo_pack_dvec(nz, work->sim_out->zn, 1, mem->z_alg, 0);

    blasfeo_pack_dvec(nx1, work->sim_out->xn, 1, &mem->fun, 0);
    blasfeo_daxpy(nx1, -1.0, mem->ux1, nu1, &mem->fun, 0, &mem->fun, 0);

    if (opts->compute_adj)
    {
        blasfeo_dgemv_n(nu + nx, nx1, -1.0, mem->BAbt, 0, 0, mem->pi, 0,
                0.0, &mem->adj, 0, &mem->adj, 0);
        blasfeo_dveccp(nx1, mem->pi, 0, &mem->adj, nu + nx);
    }

    cost_memory->common->fun = work->sim_out->xn[nx];
    blasfeo_pack_dvec(nu, work->sim_out->S_forw + nx + nx_rk * nx_rk, nx_rk,
            &cost_memory->common->grad, 0);
    blasfeo_pack_dvec(nx, work->sim_out->S_forw + nx, nx_rk,
            &cost_memory->common->grad, nu);

    if (opts->compute_hess)
    {
        blasfeo_pack_dmat(nu, nu, work->sim_out->S_hess + nx_rk + nf_rk * nx_rk,
                nf_rk, &work->hess, 0, 0);
        blasfeo_pack_dmat(nx, nu, work->sim_out->S_hess + nf_rk * nx_rk,
                nf_rk, &work->hess, nu, 0);
        blasfeo_pack_dmat(nx, nx, work->sim_out->S_hess, nf_rk,
                &work->hess, nu, nu);
        blasfeo_dtrcp_l(nu + nx, &work->hess, 0, 0, mem->RSQrq, 0, 0);
    }
}

void ocp_nlp_dynamics_cont_with_cost_compute_fun(void *config_, void *dims_,
        void *model_, void *opts_, void *mem_, void *work_)
{
    ocp_nlp_dynamics_cont_with_cost_cast_workspace(config_, dims_, opts_, work_, mem_);

    ocp_nlp_dynamics_config *config = config_;
    ocp_nlp_dynamics_cont_dims *dims = dims_;
    ocp_nlp_dynamics_cont_model *model = model_;
    ocp_nlp_dynamics_cont_opts *opts = opts_;
    ocp_nlp_dynamics_cont_memory *mem = mem_;
    ocp_nlp_dynamics_cont_workspace *work = work_;
    ocp_nlp_cost_ls_memory *cost_memory = ocp_nlp_dynamics_cont_with_cost_get_cost_memory(mem);

    int nx1 = dims->nx1;
    int nu1 = dims->nu1;
    ocp_nlp_dynamics_cont_with_cost_set_sim_input(dims, model, mem, work);

    bool sens_forw, sens_adj, sens_hess, disabled = false;
    config->sim_solver->opts_get(config->sim_solver, opts->sim_solver, "sens_forw", &sens_forw);
    config->sim_solver->opts_get(config->sim_solver, opts->sim_solver, "sens_adj", &sens_adj);
    config->sim_solver->opts_get(config->sim_solver, opts->sim_solver, "sens_hess", &sens_hess);
    config->sim_solver->opts_set(config->sim_solver, opts->sim_solver, "sens_forw", &disabled);
    config->sim_solver->opts_set(config->sim_solver, opts->sim_solver, "sens_adj", &disabled);
    config->sim_solver->opts_set(config->sim_solver, opts->sim_solver, "sens_hess", &disabled);

    config->sim_solver->evaluate(config->sim_solver, work->sim_in, work->sim_out,
            opts->sim_solver, mem->sim_solver, work->sim_solver);

    config->sim_solver->opts_set(config->sim_solver, opts->sim_solver, "sens_forw", &sens_forw);
    config->sim_solver->opts_set(config->sim_solver, opts->sim_solver, "sens_adj", &sens_adj);
    config->sim_solver->opts_set(config->sim_solver, opts->sim_solver, "sens_hess", &sens_hess);

    blasfeo_pack_dvec(nx1, work->sim_out->xn, 1, &mem->fun, 0);
    blasfeo_daxpy(nx1, -1.0, mem->ux1, nu1, &mem->fun, 0, &mem->fun, 0);
    cost_memory->common->fun = work->sim_out->xn[dims->nx];
}

void ocp_nlp_dynamics_cont_with_cost_compute_fun_and_adj(void *config_, void *dims_,
        void *model_, void *opts_, void *mem_, void *work_)
{
    printf("\nerror: ocp_nlp_dynamics_cont_with_cost_compute_fun_and_adj not implemented\n");
    exit(1);
}

void ocp_nlp_dynamics_cont_with_cost_compute_adj_sol_sens_pdiff(void *config_, void *dims_,
        void *model_, void *opts_, void *mem_, void *work_)
{
    printf("\nerror: ocp_nlp_dynamics_cont_with_cost_compute_adj_sol_sens_pdiff not implemented\n");
    exit(1);
}

void ocp_nlp_dynamics_cont_with_cost_compute_adj_p(void *config_, void *dims_,
        void *model_, void *opts_, void *mem_, struct blasfeo_dvec *out)
{
    printf("\nerror: ocp_nlp_dynamics_cont_with_cost_compute_adj_p not implemented\n");
    exit(1);
}

void ocp_nlp_dynamics_cont_with_cost_config_initialize_default(void *config_, int stage)
{
    ocp_nlp_dynamics_config *config = config_;

    ocp_nlp_dynamics_cont_config_initialize_default(config_, stage);
    config->dims_set = ocp_nlp_dynamics_cont_with_cost_dims_set;
    config->dims_get = ocp_nlp_dynamics_cont_with_cost_dims_get;
    config->opts_initialize_default = ocp_nlp_dynamics_cont_with_cost_opts_initialize_default;
    config->opts_set = ocp_nlp_dynamics_cont_with_cost_opts_set;
    config->opts_get = ocp_nlp_dynamics_cont_with_cost_opts_get;
    config->memory_set = ocp_nlp_dynamics_cont_with_cost_memory_set;
    config->memory_get = ocp_nlp_dynamics_cont_with_cost_memory_get;
    config->update_qp_matrices = ocp_nlp_dynamics_cont_with_cost_update_qp_matrices;
    config->compute_fun = ocp_nlp_dynamics_cont_with_cost_compute_fun;
    config->compute_fun_and_adj = ocp_nlp_dynamics_cont_with_cost_compute_fun_and_adj;
    config->compute_adj_sol_sens_pdiff = ocp_nlp_dynamics_cont_with_cost_compute_adj_sol_sens_pdiff;
    config->compute_adj_p = ocp_nlp_dynamics_cont_with_cost_compute_adj_p;
}