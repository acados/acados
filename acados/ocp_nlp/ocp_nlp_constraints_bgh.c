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

    // initialize to zero
    dims->nx = 0;
    dims->nu = 0;
    dims->nb = 0;
    dims->nbx = 0;
    dims->nbu = 0;
    dims->ng = 0;
    dims->nh = 0;
    dims->ns = 0;
    dims->nsbu = 0;
    dims->nsbx = 0;
    dims->nsg = 0;
    dims->nsh = 0;

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


/* dimension setters */
static void ocp_nlp_constraints_bgh_set_nx(void *config_, void *dims_, const int *nx)
{
    ocp_nlp_constraints_bgh_dims *dims = (ocp_nlp_constraints_bgh_dims *) dims_;
    dims->nx = *nx;
}



static void ocp_nlp_constraints_bgh_set_nu(void *config_, void *dims_, const int *nu)
{
    ocp_nlp_constraints_bgh_dims *dims = (ocp_nlp_constraints_bgh_dims *) dims_;
    dims->nu = *nu;
}



static void ocp_nlp_constraints_bgh_set_nbx(void *config_, void *dims_, const int *nbx)
{
    ocp_nlp_constraints_bgh_dims *dims = (ocp_nlp_constraints_bgh_dims *) dims_;
    dims->nbx = *nbx;
    dims->nb = dims->nbu + dims->nbx;
}



static void ocp_nlp_constraints_bgh_set_nbu(void *config_, void *dims_, const int *nbu)
{
    ocp_nlp_constraints_bgh_dims *dims = (ocp_nlp_constraints_bgh_dims *) dims_;
    dims->nbu = *nbu;
    dims->nb = dims->nbu + dims->nbx;
}



static void ocp_nlp_constraints_bgh_set_ng(void *config_, void *dims_, const int *ng)
{
    ocp_nlp_constraints_bgh_dims *dims = (ocp_nlp_constraints_bgh_dims *) dims_;
    dims->ng = *ng;
}



static void ocp_nlp_constraints_bgh_set_nh(void *config_, void *dims_, const int *nh)
{
    ocp_nlp_constraints_bgh_dims *dims = (ocp_nlp_constraints_bgh_dims *) dims_;
    dims->nh = *nh;
}



static void ocp_nlp_constraints_bgh_set_nsbu(void *config_, void *dims_, const int *nsbu)
{
    ocp_nlp_constraints_bgh_dims *dims = (ocp_nlp_constraints_bgh_dims *) dims_;
    dims->nsbu = *nsbu;
    dims->ns = dims->nsbu + dims->nsbx + dims->nsg + dims->nsh;
}



static void ocp_nlp_constraints_bgh_set_nsbx(void *config_, void *dims_, const int *nsbx)
{
    ocp_nlp_constraints_bgh_dims *dims = (ocp_nlp_constraints_bgh_dims *) dims_;
    dims->nsbx = *nsbx;
    dims->ns = dims->nsbu + dims->nsbx + dims->nsg + dims->nsh;
}



static void ocp_nlp_constraints_bgh_set_nsg(void *config_, void *dims_, const int *nsg)
{
    ocp_nlp_constraints_bgh_dims *dims = (ocp_nlp_constraints_bgh_dims *) dims_;
    dims->nsg = *nsg;
    dims->ns = dims->nsbu + dims->nsbx + dims->nsg + dims->nsh;
}



static void ocp_nlp_constraints_bgh_set_nsh(void *config_, void *dims_, const int *nsh)
{
    ocp_nlp_constraints_bgh_dims *dims = (ocp_nlp_constraints_bgh_dims *) dims_;
    dims->nsh = *nsh;
    dims->ns = dims->nsbu + dims->nsbx + dims->nsg + dims->nsh;
}



void ocp_nlp_constraints_bgh_dims_set(void *config_, void *dims_, const char *field,
                                             const int* value)
{
    if (!strcmp(field, "nx"))
    {
        ocp_nlp_constraints_bgh_set_nx(config_, dims_, value);
    }
    else if (!strcmp(field, "nu"))
    {
        ocp_nlp_constraints_bgh_set_nu(config_, dims_, value);
    }
    else if (!strcmp(field, "nz"))
    {
        // do nothing
        // TODO(all): implement constraints with daes
    }
    else if (!strcmp(field, "nbx"))
    {
        ocp_nlp_constraints_bgh_set_nbx(config_, dims_, value);
    }
    else if (!strcmp(field, "nbu"))
    {
        ocp_nlp_constraints_bgh_set_nbu(config_, dims_, value);
    }
    else if (!strcmp(field, "ng"))
    {
        ocp_nlp_constraints_bgh_set_ng(config_, dims_, value);
    }
    else if (!strcmp(field, "nh"))
    {
        ocp_nlp_constraints_bgh_set_nh(config_, dims_, value);
    }
    else if (!strcmp(field, "nsbu"))
    {
        ocp_nlp_constraints_bgh_set_nsbu(config_, dims_, value);
    }
    else if (!strcmp(field, "nsbx"))
    {
        ocp_nlp_constraints_bgh_set_nsbx(config_, dims_, value);
    }
    else if (!strcmp(field, "nsg"))
    {
        ocp_nlp_constraints_bgh_set_nsg(config_, dims_, value);
    }
    else if (!strcmp(field, "nsh"))
    {
        ocp_nlp_constraints_bgh_set_nsh(config_, dims_, value);
    }
    else
    {
        printf("\nerror: dims type not available in module ocp_nlp_constraints_bgh: %s\n", field);
        exit(1);
    }
}



/* dimension getters */
static void ocp_nlp_constraints_bgh_get_ni(void *config_, void *dims_, int* value)
{
    ocp_nlp_constraints_bgh_dims *dims = (ocp_nlp_constraints_bgh_dims *) dims_;
    *value = dims->nbx + dims->nbu + dims->ng + dims->nh + dims->ns;
}



static void ocp_nlp_constraints_bgh_get_nb(void *config_, void *dims_, int* value)
{
    ocp_nlp_constraints_bgh_dims *dims = (ocp_nlp_constraints_bgh_dims *) dims_;
    *value = dims->nb;
}



static void ocp_nlp_constraints_bgh_get_ng(void *config_, void *dims_, int* value)
{
    ocp_nlp_constraints_bgh_dims *dims = (ocp_nlp_constraints_bgh_dims *) dims_;
    *value = dims->ng;
}



static void ocp_nlp_constraints_bgh_get_nh(void *config_, void *dims_, int* value)
{
    ocp_nlp_constraints_bgh_dims *dims = (ocp_nlp_constraints_bgh_dims *) dims_;
    *value = dims->nh;
}



static void ocp_nlp_constraints_bgh_get_ns(void *config_, void *dims_, int* value)
{
    ocp_nlp_constraints_bgh_dims *dims = (ocp_nlp_constraints_bgh_dims *) dims_;
    *value = dims->ns;
}



void ocp_nlp_constraints_bgh_dims_get(void *config_, void *dims_, const char *field, int* value)
{
    if (!strcmp(field, "ni"))
    {
        ocp_nlp_constraints_bgh_get_ni(config_, dims_, value);
    }
    else if (!strcmp(field, "nb"))
    {
        ocp_nlp_constraints_bgh_get_nb(config_, dims_, value);
    }
    else if (!strcmp(field, "ng"))
    {
        ocp_nlp_constraints_bgh_get_ng(config_, dims_, value);
    }
    else if (!strcmp(field, "nh"))
    {
        ocp_nlp_constraints_bgh_get_nh(config_, dims_, value);
    }
    else if (!strcmp(field, "ns"))
    {
        ocp_nlp_constraints_bgh_get_ns(config_, dims_, value);
    }
    else
    {
        printf("error: attempt to get dimension from constraint model, that is not there");
        exit(1);
    }
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
    // default initialization to zero
    blasfeo_dvecse(2*nb+2*ng+2*nh+2*ns, 0.0, &model->d, 0);

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



int ocp_nlp_constraints_bgh_model_set(void *config_, void *dims_,
                         void *model_, const char *field, void *value)
{
    ocp_nlp_constraints_bgh_dims *dims = (ocp_nlp_constraints_bgh_dims *) dims_;
    ocp_nlp_constraints_bgh_model *model = (ocp_nlp_constraints_bgh_model *) model_;

    int status = ACADOS_FAILURE;

    int ii;
    int *ptr_i;

    if (!dims || !model || !field || !value) return ACADOS_FAILURE;

    int nu = dims->nu;
    int nx = dims->nx;
    int nb = dims->nb;
    int ng = dims->ng;
    int nh = dims->nh;
    int ns = dims->ns;
    int nsbu = dims->nsbu;
    int nsbx = dims->nsbx;
    int nsg = dims->nsg;
    int nsh = dims->nsh;
    int nbx = dims->nbx;
    int nbu = dims->nbu;
    // TODO(oj): document which strings mean what! - adapted from prev implementation..
    if (!strcmp(field, "lb")) // TODO(fuck_lint) remove !!!
    {
        blasfeo_pack_dvec(nb, value, &model->d, 0);
        status = ACADOS_SUCCESS;
    }
    else if (!strcmp(field, "ub")) // TODO(fuck_lint) remove !!!
    {
        blasfeo_pack_dvec(nb, value, &model->d, nb+ng+nh);
        status = ACADOS_SUCCESS;
    }
    else if (!strcmp(field, "idxbx"))
    {
        ptr_i = (int *) value;
        for (ii=0; ii < nbx; ii++)
            model->idxb[nbu+ii] = nu+ptr_i[ii];
        status = ACADOS_SUCCESS;
    }
    else if (!strcmp(field, "lbx"))
    {
        blasfeo_pack_dvec(nbx, value, &model->d, nbu);
        status = ACADOS_SUCCESS;
    }
    else if (!strcmp(field, "ubx"))
    {
        blasfeo_pack_dvec(nbx, value, &model->d, nb + ng + nh + nbu);
        status = ACADOS_SUCCESS;
    }
    else if (!strcmp(field, "idxbu"))
    {
        ptr_i = (int *) value;
        for (ii=0; ii < nbu; ii++)
            model->idxb[ii] = ptr_i[ii];
        status = ACADOS_SUCCESS;
    }
    else if (!strcmp(field, "lbu"))
    {
        blasfeo_pack_dvec(nbu, value, &model->d, 0);
        status = ACADOS_SUCCESS;
    }
    else if (!strcmp(field, "ubu"))
    {
        blasfeo_pack_dvec(nbu, value, &model->d, nb + ng + nh);
        status = ACADOS_SUCCESS;
    }
    else if (!strcmp(field, "C"))
    {
        blasfeo_pack_tran_dmat(ng, nx, value, ng, &model->DCt, nu, 0);
        status = ACADOS_SUCCESS;
    }
    else if (!strcmp(field, "D"))
    {
        blasfeo_pack_tran_dmat(ng, nu, value, ng, &model->DCt, 0, 0);
        status = ACADOS_SUCCESS;
    }
    else if (!strcmp(field, "lg"))
    {
        blasfeo_pack_dvec(ng, value, &model->d, nb);
        status = ACADOS_SUCCESS;
    }
    else if (!strcmp(field, "ug"))
    {
        blasfeo_pack_dvec(ng, value, &model->d, 2*nb+ng+nh);
        status = ACADOS_SUCCESS;
    }
    else if (!strcmp(field, "h"))
    {
        model->h = value;
        status = ACADOS_SUCCESS;
    }
    else if (!strcmp(field, "lh")) // TODO(fuck_lint) remove
    {
        blasfeo_pack_dvec(nh, value, &model->d, nb+ng);
        status = ACADOS_SUCCESS;
    }
    else if (!strcmp(field, "uh")) // TODO(fuck_lint) remove
    {
        blasfeo_pack_dvec(nh, value, &model->d, 2*nb+2*ng+nh);
        status = ACADOS_SUCCESS;
    }
    else if (!strcmp(field, "idxsbu"))
    {
        ptr_i = (int *) value;
        for (ii=0; ii < nsbu; ii++)
            model->idxs[ii] = ptr_i[ii];
        status = ACADOS_SUCCESS;
    }
    else if (!strcmp(field, "lsbu"))
    {
        blasfeo_pack_dvec(nsbu, value, &model->d, 2*nb+2*ng+2*nh);
        status = ACADOS_SUCCESS;
    }
    else if (!strcmp(field, "usbu"))
    {
        blasfeo_pack_dvec(nsbu, value, &model->d, 2*nb+2*ng+2*nh+ns);
        status = ACADOS_SUCCESS;
    }
    else if (!strcmp(field, "idxsbx"))
    {
        ptr_i = (int *) value;
        for (ii=0; ii < nsbx; ii++)
            model->idxs[nsbu+ii] = nbu+ptr_i[ii];
        status = ACADOS_SUCCESS;
    }
    else if (!strcmp(field, "lsbx"))
    {
        blasfeo_pack_dvec(nsbx, value, &model->d, 2*nb+2*ng+2*nh+nsbu);
        status = ACADOS_SUCCESS;
    }
    else if (!strcmp(field, "usbx"))
    {
        blasfeo_pack_dvec(nsbx, value, &model->d, 2*nb+2*ng+2*nh+ns+nsbu);
        status = ACADOS_SUCCESS;
    }
    else if (!strcmp(field, "idxsg"))
    {
        ptr_i = (int *) value;
        for (ii=0; ii < nsg; ii++)
            model->idxs[nsbu+nsbx+ii] = nbu+nbx+ptr_i[ii];
        status = ACADOS_SUCCESS;
    }
    else if (!strcmp(field, "lsg"))
    {
        blasfeo_pack_dvec(nsg, value, &model->d, 2*nb+2*ng+2*nh+nsbu+nsbx);
        status = ACADOS_SUCCESS;
    }
    else if (!strcmp(field, "usg"))
    {
        blasfeo_pack_dvec(nsg, value, &model->d, 2*nb+2*ng+2*nh+ns+nsbu+nsbx);
        status = ACADOS_SUCCESS;
    }
    else if (!strcmp(field, "idxsh"))
    {
        ptr_i = (int *) value;
        for (ii=0; ii < nsh; ii++)
            model->idxs[nsbu+nsbx+nsg+ii] = nbu+nbx+ng+ptr_i[ii];
        status = ACADOS_SUCCESS;
    }
    else if (!strcmp(field, "lsh"))
    {
        blasfeo_pack_dvec(nsh, value, &model->d, 2*nb+2*ng+2*nh+nsbu+nsbx+nsg);
        status = ACADOS_SUCCESS;
    }
    else if (!strcmp(field, "ush"))
    {
        blasfeo_pack_dvec(nsh, value, &model->d, 2*nb+2*ng+2*nh+ns+nsbu+nsbx+nsg);
        status = ACADOS_SUCCESS;
    }
    else
    {
        printf("\nerror: model field not available in module ocp_nlp_constraints_bgh: %s\n", field);
        exit(1);
    }

    return status;
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
    // int nx = dims->nx;
    // int nu = dims->nu;
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
    // int nx = dims->nx;
    // int nu = dims->nu;
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
    config->dims_set = &ocp_nlp_constraints_bgh_dims_set;
    config->get_dims = &ocp_nlp_constraints_bgh_dims_get;
    config->model_calculate_size = &ocp_nlp_constraints_bgh_model_calculate_size;
    config->model_assign = &ocp_nlp_constraints_bgh_model_assign;
    config->model_set = &ocp_nlp_constraints_bgh_model_set;
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
