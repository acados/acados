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

#include "acados/ocp_nlp/ocp_nlp_constraints_bgh.h"

#include <assert.h>
#include <stdlib.h>
#include <string.h>

// blasfeo
#include "blasfeo/include/blasfeo_d_aux.h"
#include "blasfeo/include/blasfeo_d_blas.h"
// acados
#include "acados/ocp_qp/ocp_qp_common.h"
#include "acados/utils/mem.h"



/************************************************
 * dims
 ************************************************/

int ocp_nlp_constraints_bgh_dims_calculate_size(void *config_)
{
    int size = sizeof(ocp_nlp_constraints_bgh_dims);

    return size;
}



void *ocp_nlp_constraints_bgh_dims_assign(void *config_, void *raw_memory)
{
    char *c_ptr = (char *) raw_memory;

    ocp_nlp_constraints_bgh_dims *dims = (ocp_nlp_constraints_bgh_dims *) c_ptr;
    c_ptr += sizeof(ocp_nlp_constraints_bgh_dims);

    assert((char *) raw_memory + ocp_nlp_constraints_bgh_dims_calculate_size(config_) >= c_ptr);

    return dims;
}



void ocp_nlp_constraints_bgh_dims_initialize(void *config_, void *dims_, int nx, int nu, int nbx,
                                             int nbu, int ng, int nh, int dummy0, int ns)
{
    ocp_nlp_constraints_bgh_dims *dims = dims_;

    dims->nx = nx;
    dims->nu = nu;
    dims->nbx = nbx;
    dims->nbu = nbu;
    dims->nb = nbx + nbu;
    dims->ng = ng;
    dims->nh = nh;
    dims->ns = ns;

    return;
}



/************************************************
 * model
 ************************************************/

int ocp_nlp_constraints_bgh_model_calculate_size(void *config, void *dims_)
{
    ocp_nlp_constraints_bgh_dims *dims = dims_;

    // extract dims
    int nx = dims->nx;
    int nu = dims->nu;
    int nb = dims->nb;
    int ng = dims->ng;
    int nh = dims->nh;
    int ns = dims->ns;

    int size = 0;

    size += sizeof(ocp_nlp_constraints_bgh_model);

    size += sizeof(int) * nb;                                         // idxb
    size += sizeof(int) * ns;                                         // idxs
    size += blasfeo_memsize_dvec(2 * nb + 2 * ng + 2 * nh + 2 * ns);  // d
    size += blasfeo_memsize_dmat(nu + nx, ng);                        // DCt

    size += 64;  // blasfeo_mem align

    return size;
}



void *ocp_nlp_constraints_bgh_model_assign(void *config, void *dims_, void *raw_memory)
{
    ocp_nlp_constraints_bgh_dims *dims = dims_;

    char *c_ptr = (char *) raw_memory;

    // extract sizes
    int nx = dims->nx;
    int nu = dims->nu;
    int nb = dims->nb;
    int ng = dims->ng;
    int nh = dims->nh;
    int ns = dims->ns;

    // struct
    ocp_nlp_constraints_bgh_model *model = (ocp_nlp_constraints_bgh_model *) c_ptr;
    c_ptr += sizeof(ocp_nlp_constraints_bgh_model);

    // dims
    //  model->dims = dims;

    // blasfeo_mem align
    align_char_to(64, &c_ptr);

    // blasfeo_dmat
    // DCt
    assign_and_advance_blasfeo_dmat_mem(nu + nx, ng, &model->DCt, &c_ptr);

    // blasfeo_dvec
    // d
    assign_and_advance_blasfeo_dvec_mem(2 * nb + 2 * ng + 2 * nh + 2 * ns, &model->d, &c_ptr);

    // int
    // idxb
    assign_and_advance_int(nb, &model->idxb, &c_ptr);
    // idxs
    assign_and_advance_int(ns, &model->idxs, &c_ptr);

    // h
    //  model->h = NULL;

    // assert
    assert((char *) raw_memory + ocp_nlp_constraints_bgh_model_calculate_size(config, dims) >=
           c_ptr);

    return model;
}



/************************************************
 * options
 ************************************************/

int ocp_nlp_constraints_bgh_opts_calculate_size(void *config_, void *dims_)
{
    int size = 0;

    size += sizeof(ocp_nlp_constraints_bgh_opts);

    return size;
}



void *ocp_nlp_constraints_bgh_opts_assign(void *config_, void *dims_, void *raw_memory)
{
    char *c_ptr = (char *) raw_memory;

    ocp_nlp_constraints_bgh_opts *opts = (ocp_nlp_constraints_bgh_opts *) c_ptr;
    c_ptr += sizeof(ocp_nlp_constraints_bgh_opts);

    assert((char *) raw_memory + ocp_nlp_constraints_bgh_opts_calculate_size(config_, dims_) >=
           c_ptr);

    return opts;
}



void ocp_nlp_constraints_bgh_opts_initialize_default(void *config_, void *dims_, void *opts_)
{
    ocp_nlp_constraints_bgh_opts *opts = opts_;

    opts->compute_adj = 1;

    return;
}



void ocp_nlp_constraints_bgh_opts_update(void *config_, void *dims_, void *opts_)
{
    //  ocp_nlp_constraints_bgh_opts *opts = opts_;

    return;
}



void ocp_nlp_constraints_bgh_opts_set(void *config_, void *dims_, void *opts_,
    enum acados_opts name, void *ptr_value)
{

    ocp_nlp_constraints_bgh_opts *opts = opts_;

    if (name == COMPUTE_ADJ)
    {
        int *compute_adj = ptr_value;
        opts->compute_adj = *compute_adj;
    }
    else
    {
        // TODO(fuck_you_lint): something better tha this print-and-exit
        printf("\nocp_nlp_constraints_bgh_opts_set: unknown opts name !\n");
        exit(1);
    }

    return;

}



/************************************************
 * memory
 ************************************************/

int ocp_nlp_constraints_bgh_memory_calculate_size(void *config_, void *dims_, void *opts_)
{
    ocp_nlp_constraints_bgh_dims *dims = dims_;

    // extract dims
    int nx = dims->nx;
    int nu = dims->nu;
    int nb = dims->nb;
    int ng = dims->ng;
    int nh = dims->nh;
    int ns = dims->ns;

    int size = 0;

    size += sizeof(ocp_nlp_constraints_bgh_memory);

    size += 1 * blasfeo_memsize_dvec(2 * nb + 2 * ng + 2 * nh + 2 * ns);  // fun
    size += 1 * blasfeo_memsize_dvec(nu + nx + 2 * ns);                   // adj

    size += 1 * 64;  // blasfeo_mem align

    return size;
}



void *ocp_nlp_constraints_bgh_memory_assign(void *config_, void *dims_, void *opts_,
                                            void *raw_memory)
{
    ocp_nlp_constraints_bgh_dims *dims = dims_;

    char *c_ptr = (char *) raw_memory;

    // extract dims
    int nx = dims->nx;
    int nu = dims->nu;
    int nb = dims->nb;
    int ng = dims->ng;
    int nh = dims->nh;
    int ns = dims->ns;

    // struct
    ocp_nlp_constraints_bgh_memory *memory = (ocp_nlp_constraints_bgh_memory *) c_ptr;
    c_ptr += sizeof(ocp_nlp_constraints_bgh_memory);

    // blasfeo_mem align
    align_char_to(64, &c_ptr);

    // fun
    assign_and_advance_blasfeo_dvec_mem(2 * nb + 2 * ng + 2 * nh + 2 * ns, &memory->fun, &c_ptr);
    // adj
    assign_and_advance_blasfeo_dvec_mem(nu + nx + 2 * ns, &memory->adj, &c_ptr);

    assert((char *) raw_memory +
               ocp_nlp_constraints_bgh_memory_calculate_size(config_, dims, opts_) >=
           c_ptr);

    return memory;
}



struct blasfeo_dvec *ocp_nlp_constraints_bgh_memory_get_fun_ptr(void *memory_)
{
    ocp_nlp_constraints_bgh_memory *memory = memory_;

    return &memory->fun;
}



struct blasfeo_dvec *ocp_nlp_constraints_bgh_memory_get_adj_ptr(void *memory_)
{
    ocp_nlp_constraints_bgh_memory *memory = memory_;

    return &memory->adj;
}



void ocp_nlp_constraints_bgh_memory_set_ux_ptr(struct blasfeo_dvec *ux, void *memory_)
{
    ocp_nlp_constraints_bgh_memory *memory = memory_;

    memory->ux = ux;
}



void ocp_nlp_constraints_bgh_memory_set_lam_ptr(struct blasfeo_dvec *lam, void *memory_)
{
    ocp_nlp_constraints_bgh_memory *memory = memory_;

    memory->lam = lam;
}



void ocp_nlp_constraints_bgh_memory_set_DCt_ptr(struct blasfeo_dmat *DCt, void *memory_)
{
    ocp_nlp_constraints_bgh_memory *memory = memory_;

    memory->DCt = DCt;
}



void ocp_nlp_constraints_bgh_memory_set_RSQrq_ptr(struct blasfeo_dmat *RSQrq, void *memory_)
{
    ocp_nlp_constraints_bgh_memory *memory = memory_;

    memory->RSQrq = RSQrq;
}



void ocp_nlp_constraints_bgh_memory_set_idxb_ptr(int *idxb, void *memory_)
{
    ocp_nlp_constraints_bgh_memory *memory = memory_;

    memory->idxb = idxb;
}



void ocp_nlp_constraints_bgh_memory_set_idxs_ptr(int *idxs, void *memory_)
{
    ocp_nlp_constraints_bgh_memory *memory = memory_;

    memory->idxs = idxs;
}



/************************************************
 * workspace
 ************************************************/

int ocp_nlp_constraints_bgh_workspace_calculate_size(void *config_, void *dims_, void *opts_)
{
    ocp_nlp_constraints_bgh_dims *dims = dims_;

    // extract dims
    int nx = dims->nx;
    int nu = dims->nu;
    int nb = dims->nb;
    int ng = dims->ng;
    int nh = dims->nh;
    int ns = dims->ns;

    int size = 0;

    size += sizeof(ocp_nlp_constraints_bgh_workspace);

    size += 1 * blasfeo_memsize_dvec(nb + ng + nh + ns);  // tmp_ni

    size += 1 * 64;  // blasfeo_mem align

    return size;
}



static void ocp_nlp_constraints_bgh_cast_workspace(void *config_, void *dims_, void *opts_,
                                                   void *work_)
{
    ocp_nlp_constraints_bgh_dims *dims = dims_;
    ocp_nlp_constraints_bgh_workspace *work = work_;

    // extract dims
    int nx = dims->nx;
    int nu = dims->nu;
    int nb = dims->nb;
    int ng = dims->ng;
    int nh = dims->nh;
    int ns = dims->ns;

    char *c_ptr = (char *) work_;
    c_ptr += sizeof(ocp_nlp_constraints_bgh_workspace);

    // blasfeo_mem align
    align_char_to(64, &c_ptr);

    // tmp_ni
    assign_and_advance_blasfeo_dvec_mem(nb + ng + nh + ns, &work->tmp_ni, &c_ptr);

    assert((char *) work + ocp_nlp_constraints_bgh_workspace_calculate_size(config_, dims, opts_) >=
           c_ptr);

    return;
}



/************************************************
 * functions
 ************************************************/

void ocp_nlp_constraints_bgh_initialize(void *config_, void *dims_, void *model_, void *opts,
                                        void *memory_, void *work_)
{
    ocp_nlp_constraints_bgh_dims *dims = dims_;
    ocp_nlp_constraints_bgh_model *model = model_;
    ocp_nlp_constraints_bgh_memory *memory = memory_;

    // loop index
    int j;

    // extract dims
    int nx = dims->nx;
    int nu = dims->nu;
    int nb = dims->nb;
    int ng = dims->ng;
    int ns = dims->ns;

    // initialize idxb
    for (j = 0; j < nb; j++)
    {
        memory->idxb[j] = model->idxb[j];
    }

    // initialize idxs
    for (j = 0; j < ns; j++)
    {
        memory->idxs[j] = model->idxs[j];
    }

    // initialize general constraints matrix
    blasfeo_dgecp(nu + nx, ng, &model->DCt, 0, 0, memory->DCt, 0, 0);

    return;
}



void ocp_nlp_constraints_bgh_update_qp_matrices(void *config_, void *dims_, void *model_,
                                                void *opts_, void *memory_, void *work_)
{
    ocp_nlp_constraints_bgh_dims *dims = dims_;
    ocp_nlp_constraints_bgh_model *model = model_;
    ocp_nlp_constraints_bgh_opts *opts = opts_;
    ocp_nlp_constraints_bgh_memory *memory = memory_;
    ocp_nlp_constraints_bgh_workspace *work = work_;

    ocp_nlp_constraints_bgh_cast_workspace(config_, dims, opts_, work_);

    // extract dims
    int nx = dims->nx;
    int nu = dims->nu;
    int nb = dims->nb;
    int ng = dims->ng;
    int nh = dims->nh;
    int ns = dims->ns;

    // XXX large enough ?
    ext_fun_arg_t ext_fun_type_in[2];
    void *ext_fun_in[2];
    ext_fun_arg_t ext_fun_type_out[2];
    void *ext_fun_out[2];

    // box
    blasfeo_dvecex_sp(nb, 1.0, model->idxb, memory->ux, 0, &work->tmp_ni, 0);

    // general linear
    blasfeo_dgemv_t(nu + nx, ng, 1.0, memory->DCt, 0, 0, memory->ux, 0, 0.0, &work->tmp_ni, nb,
                    &work->tmp_ni, nb);

    // nonlinear
    if (nh > 0)
    {
        //
        ext_fun_type_in[0] = BLASFEO_DVEC;
        ext_fun_in[0] = memory->ux;  // ux: nu+nx

        //
        ext_fun_type_out[0] = BLASFEO_DVEC_ARGS;
        struct blasfeo_dvec_args h_args;
        h_args.x = &work->tmp_ni;
        h_args.xi = nb + ng;
        ext_fun_out[0] = &h_args;  // fun: nh
        //
        ext_fun_type_out[1] = BLASFEO_DMAT_ARGS;
        struct blasfeo_dmat_args Jht_args;
        Jht_args.A = memory->DCt;
        Jht_args.ai = 0;
        Jht_args.aj = ng;
        ext_fun_out[1] = &Jht_args;  // jac': (nu+nx) * nh

        model->h->evaluate(model->h, ext_fun_type_in, ext_fun_in, ext_fun_type_out, ext_fun_out);
    }

    blasfeo_daxpy(nb + ng + nh, -1.0, &work->tmp_ni, 0, &model->d, 0, &memory->fun, 0);
    blasfeo_daxpy(nb + ng + nh, -1.0, &model->d, nb + ng + nh, &work->tmp_ni, 0, &memory->fun,
                  nb + ng + nh);

    // soft
    blasfeo_dvecad_sp(ns, -1.0, memory->ux, nu + nx, model->idxs, &memory->fun, 0);
    blasfeo_dvecad_sp(ns, -1.0, memory->ux, nu + nx + ns, model->idxs, &memory->fun, nb + ng + nh);

    blasfeo_daxpy(2 * ns, -1.0, memory->ux, nu + nx, &model->d, 2 * nb + 2 * ng + 2 * nh,
                  &memory->fun, 2 * nb + 2 * ng + 2 * nh);

    // nlp_mem: ineq_adj
    if (opts->compute_adj)
    {
        blasfeo_dvecse(nu + nx + 2 * ns, 0.0, &memory->adj, 0);
        blasfeo_daxpy(nb + ng + nh, -1.0, memory->lam, nb+ng+nh, memory->lam, 0, &work->tmp_ni, 0);
        blasfeo_dvecad_sp(nb, 1.0, &work->tmp_ni, 0, model->idxb, &memory->adj, 0);
        blasfeo_dgemv_n(nu+nx, ng+nh, 1.0, memory->DCt, 0, 0, &work->tmp_ni, nb, 1.0, &memory->adj,
                        0, &memory->adj, 0);
        // soft
        blasfeo_dvecex_sp(ns, 1.0, model->idxs, memory->lam, 0, &memory->adj, nu + nx);
        blasfeo_dvecex_sp(ns, 1.0, model->idxs, memory->lam, nb + ng + nh, &memory->adj, nu+nx+ns);
        blasfeo_daxpy(2 * ns, 1.0, memory->lam, 2 * nb + 2 * ng + 2 * nh, &memory->adj, nu + nx,
                      &memory->adj, nu + nx);
    }

    return;
}



void ocp_nlp_constraints_bgh_config_initialize_default(void *config_)
{
    ocp_nlp_constraints_config *config = config_;

    config->dims_calculate_size = &ocp_nlp_constraints_bgh_dims_calculate_size;
    config->dims_assign = &ocp_nlp_constraints_bgh_dims_assign;
    config->dims_initialize = &ocp_nlp_constraints_bgh_dims_initialize;
    config->model_calculate_size = &ocp_nlp_constraints_bgh_model_calculate_size;
    config->model_assign = &ocp_nlp_constraints_bgh_model_assign;
    config->opts_calculate_size = &ocp_nlp_constraints_bgh_opts_calculate_size;
    config->opts_assign = &ocp_nlp_constraints_bgh_opts_assign;
    config->opts_initialize_default = &ocp_nlp_constraints_bgh_opts_initialize_default;
    config->opts_update = &ocp_nlp_constraints_bgh_opts_update;
    config->opts_set = &ocp_nlp_constraints_bgh_opts_set;
    config->memory_calculate_size = &ocp_nlp_constraints_bgh_memory_calculate_size;
    config->memory_assign = &ocp_nlp_constraints_bgh_memory_assign;
    config->memory_get_fun_ptr = &ocp_nlp_constraints_bgh_memory_get_fun_ptr;
    config->memory_get_adj_ptr = &ocp_nlp_constraints_bgh_memory_get_adj_ptr;
    config->memory_set_ux_ptr = &ocp_nlp_constraints_bgh_memory_set_ux_ptr;
    config->memory_set_lam_ptr = &ocp_nlp_constraints_bgh_memory_set_lam_ptr;
    config->memory_set_DCt_ptr = &ocp_nlp_constraints_bgh_memory_set_DCt_ptr;
    config->memory_set_RSQrq_ptr = &ocp_nlp_constraints_bgh_memory_set_RSQrq_ptr;
    config->memory_set_idxb_ptr = &ocp_nlp_constraints_bgh_memory_set_idxb_ptr;
    config->memory_set_idxs_ptr = &ocp_nlp_constraints_bgh_memory_set_idxs_ptr;
    config->workspace_calculate_size = &ocp_nlp_constraints_bgh_workspace_calculate_size;
    config->initialize = &ocp_nlp_constraints_bgh_initialize;
    config->update_qp_matrices = &ocp_nlp_constraints_bgh_update_qp_matrices;
    config->config_initialize_default = &ocp_nlp_constraints_bgh_config_initialize_default;

    return;
}
