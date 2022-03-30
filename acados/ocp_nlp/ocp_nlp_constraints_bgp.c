/*
 * Copyright 2019 Gianluca Frison, Dimitris Kouzoupis, Robin Verschueren,
 * Andrea Zanelli, Niels van Duijkeren, Jonathan Frey, Tommaso Sartor,
 * Branimir Novoselnik, Rien Quirynen, Rezart Qelibari, Dang Doan,
 * Jonas Koenemann, Yutao Chen, Tobias Schöls, Jonas Schlagenhauf, Moritz Diehl
 *
 * This file is part of acados.
 *
 * The 2-Clause BSD License
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 * this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 * this list of conditions and the following disclaimer in the documentation
 * and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.;
 */


#include "acados/ocp_nlp/ocp_nlp_constraints_bgp.h"

#include <assert.h>
#include <stdlib.h>
#include <string.h>

// blasfeo
#include "blasfeo/include/blasfeo_d_aux.h"
#include "blasfeo/include/blasfeo_d_blas.h"
// acados
#include "acados/ocp_qp/ocp_qp_common.h"
#include "acados/utils/mem.h"



/* dims */

acados_size_t ocp_nlp_constraints_bgp_dims_calculate_size(void *config_)
{
    acados_size_t size = sizeof(ocp_nlp_constraints_bgp_dims);

    return size;
}



void *ocp_nlp_constraints_bgp_dims_assign(void *config_, void *raw_memory)
{
    char *c_ptr = (char *) raw_memory;

    ocp_nlp_constraints_bgp_dims *dims = (ocp_nlp_constraints_bgp_dims *) c_ptr;
    c_ptr += sizeof(ocp_nlp_constraints_bgp_dims);

    assert((char *) raw_memory + ocp_nlp_constraints_bgp_dims_calculate_size(config_) >= c_ptr);

    // initialize to zero
    dims->nx = 0;
    dims->nu = 0;
    dims->nz = 0;
    dims->nb = 0;
    dims->nbx = 0;
    dims->nbu = 0;
    dims->ng = 0;
    dims->nphi = 0;
    dims->ns = 0;
    dims->nsbu = 0;
    dims->nsbx = 0;
    dims->nsg = 0;
    dims->nsphi = 0;
    dims->nr = 0;
    dims->nbxe = 0;
    dims->nbue = 0;
    dims->nge = 0;
    dims->nphie = 0;

    return dims;
}



// TODO outdated? to remove?
void ocp_nlp_constraints_bgp_dims_initialize(void *config_, void *dims_, int nx, int nu, int nz,
        int nbx, int nbu, int ng, int nphi, int nr, int ns)
{
    ocp_nlp_constraints_bgp_dims *dims = dims_;

    dims->nx = nx;
    dims->nu = nu;
    dims->nz = nz;
    dims->nbx = nbx;
    dims->nbu = nbu;
    dims->nb = nbx + nbu;
    dims->ng = ng;
    dims->nphi = nphi;
    dims->nr = nr;
    dims->ns = ns;

    return;
}



/* dimension setters */
static void ocp_nlp_constraints_bgp_set_nx(void *config_, void *dims_, const int *nx)
{
    ocp_nlp_constraints_bgp_dims *dims = (ocp_nlp_constraints_bgp_dims *) dims_;
    dims->nx = *nx;
}



static void ocp_nlp_constraints_bgp_set_nu(void *config_, void *dims_, const int *nu)
{
    ocp_nlp_constraints_bgp_dims *dims = (ocp_nlp_constraints_bgp_dims *) dims_;
    dims->nu = *nu;
}



static void ocp_nlp_constraints_bgp_set_nz(void *config_, void *dims_, const int *nz)
{
    ocp_nlp_constraints_bgp_dims *dims = (ocp_nlp_constraints_bgp_dims *) dims_;
    dims->nz = *nz;
}



static void ocp_nlp_constraints_bgp_set_nbx(void *config_, void *dims_, const int *nbx)
{
    ocp_nlp_constraints_bgp_dims *dims = (ocp_nlp_constraints_bgp_dims *) dims_;
    dims->nbx = *nbx;
    dims->nb = *nbx + dims->nbu;
}



static void ocp_nlp_constraints_bgp_set_nbu(void *config_, void *dims_, const int *nbu)
{
    ocp_nlp_constraints_bgp_dims *dims = (ocp_nlp_constraints_bgp_dims *) dims_;
    dims->nbu = *nbu;
    dims->nb = *nbu + dims->nbx;
}



static void ocp_nlp_constraints_bgp_set_ng(void *config_, void *dims_, const int *ng)
{
    ocp_nlp_constraints_bgp_dims *dims = (ocp_nlp_constraints_bgp_dims *) dims_;
    dims->ng = *ng;
}



static void ocp_nlp_constraints_bgp_set_nphi(void *config_, void *dims_, const int *nphi)
{
    ocp_nlp_constraints_bgp_dims *dims = (ocp_nlp_constraints_bgp_dims *) dims_;
    dims->nphi = *nphi;
}



static void ocp_nlp_constraints_bgp_set_nsbu(void *config_, void *dims_, const int *nsbu)
{
    ocp_nlp_constraints_bgp_dims *dims = (ocp_nlp_constraints_bgp_dims *) dims_;
    dims->nsbu = *nsbu;
    dims->ns = dims->nsbu + dims->nsbx + dims->nsg + dims->nsphi;
}



static void ocp_nlp_constraints_bgp_set_nsbx(void *config_, void *dims_, const int *nsbx)
{
    ocp_nlp_constraints_bgp_dims *dims = (ocp_nlp_constraints_bgp_dims *) dims_;
    dims->nsbx = *nsbx;
    dims->ns = dims->nsbu + dims->nsbx + dims->nsg + dims->nsphi;
}



static void ocp_nlp_constraints_bgp_set_nsg(void *config_, void *dims_, const int *nsg)
{
    ocp_nlp_constraints_bgp_dims *dims = (ocp_nlp_constraints_bgp_dims *) dims_;
    dims->nsg = *nsg;
    dims->ns = dims->nsbu + dims->nsbx + dims->nsg + dims->nsphi;
}



static void ocp_nlp_constraints_bgp_set_nsphi(void *config_, void *dims_, const int *nsphi)
{
    ocp_nlp_constraints_bgp_dims *dims = (ocp_nlp_constraints_bgp_dims *) dims_;
    dims->nsphi = *nsphi;
    dims->ns = dims->nsbu + dims->nsbx + dims->nsg + dims->nsphi;
}



static void ocp_nlp_constraints_bgp_set_nr(void *config_, void *dims_, const int *nr)
{
    ocp_nlp_constraints_bgp_dims *dims = (ocp_nlp_constraints_bgp_dims *) dims_;
    dims->nr = *nr;
}



static void ocp_nlp_constraints_bgp_set_nbxe(void *config_, void *dims_, const int *nbxe)
{
    ocp_nlp_constraints_bgp_dims *dims = (ocp_nlp_constraints_bgp_dims *) dims_;
    dims->nbxe = *nbxe;
}



static void ocp_nlp_constraints_bgp_set_nbue(void *config_, void *dims_, const int *nbue)
{
    ocp_nlp_constraints_bgp_dims *dims = (ocp_nlp_constraints_bgp_dims *) dims_;
    dims->nbue = *nbue;
}



static void ocp_nlp_constraints_bgp_set_nge(void *config_, void *dims_, const int *nge)
{
    ocp_nlp_constraints_bgp_dims *dims = (ocp_nlp_constraints_bgp_dims *) dims_;
    dims->nge = *nge;
}



static void ocp_nlp_constraints_bgp_set_nphie(void *config_, void *dims_, const int *nphie)
{
    ocp_nlp_constraints_bgp_dims *dims = (ocp_nlp_constraints_bgp_dims *) dims_;
    dims->nphie = *nphie;
}



void ocp_nlp_constraints_bgp_dims_set(void *config_, void *dims_,
                                       const char *field, const int* value)
{
    if (!strcmp(field, "nx"))
    {
        ocp_nlp_constraints_bgp_set_nx(config_, dims_, value);
    }
    else if (!strcmp(field, "nz"))
    {
        ocp_nlp_constraints_bgp_set_nz(config_, dims_, value);
    }
    else if (!strcmp(field, "nu"))
    {
        ocp_nlp_constraints_bgp_set_nu(config_, dims_, value);
    }
    else if (!strcmp(field, "nbx"))
    {
        ocp_nlp_constraints_bgp_set_nbx(config_, dims_, value);
    }
    else if (!strcmp(field, "nbu"))
    {
        ocp_nlp_constraints_bgp_set_nbu(config_, dims_, value);
    }
    else if (!strcmp(field, "ng"))
    {
        ocp_nlp_constraints_bgp_set_ng(config_, dims_, value);
    }
    else if (!strcmp(field, "nphi"))
    {
        ocp_nlp_constraints_bgp_set_nphi(config_, dims_, value);
    }
    else if (!strcmp(field, "nsbu"))
    {
        ocp_nlp_constraints_bgp_set_nsbu(config_, dims_, value);
    }
    else if (!strcmp(field, "nsbx"))
    {
        ocp_nlp_constraints_bgp_set_nsbx(config_, dims_, value);
    }
    else if (!strcmp(field, "nsg"))
    {
        ocp_nlp_constraints_bgp_set_nsg(config_, dims_, value);
    }
    else if (!strcmp(field, "nsphi"))
    {
        ocp_nlp_constraints_bgp_set_nsphi(config_, dims_, value);
    }
    else if (!strcmp(field, "nr"))
    {
        ocp_nlp_constraints_bgp_set_nr(config_, dims_, value);
    }
    else if (!strcmp(field, "nbxe"))
    {
        ocp_nlp_constraints_bgp_set_nbxe(config_, dims_, value);
    }
    else if (!strcmp(field, "nbue"))
    {
        ocp_nlp_constraints_bgp_set_nbue(config_, dims_, value);
    }
    else if (!strcmp(field, "nge"))
    {
        ocp_nlp_constraints_bgp_set_nge(config_, dims_, value);
    }
    else if (!strcmp(field, "nphie"))
    {
        ocp_nlp_constraints_bgp_set_nphie(config_, dims_, value);
    }
    else
    {
        printf("\nerror: dims type not available in module ocp_nlp_constraints_bgp: %s\n", field);
        exit(1);
    }
}



/* dimension getters */
static void ocp_nlp_constraints_bgp_get_ni(void *config_, void *dims_, int* value)
{
    ocp_nlp_constraints_bgp_dims *dims = (ocp_nlp_constraints_bgp_dims *) dims_;
    *value = dims->nbx + dims->nbu + dims->ng + dims->nphi + dims->ns;
}



static void ocp_nlp_constraints_bgp_get_nb(void *config_, void *dims_, int* value)
{
    ocp_nlp_constraints_bgp_dims *dims = (ocp_nlp_constraints_bgp_dims *) dims_;
    *value = dims->nb;
}



static void ocp_nlp_constraints_bgp_get_nbx(void *config_, void *dims_, int* value)
{
    ocp_nlp_constraints_bgp_dims *dims = (ocp_nlp_constraints_bgp_dims *) dims_;
    *value = dims->nbx;
}



static void ocp_nlp_constraints_bgp_get_nbu(void *config_, void *dims_, int* value)
{
    ocp_nlp_constraints_bgp_dims *dims = (ocp_nlp_constraints_bgp_dims *) dims_;
    *value = dims->nbu;
}



static void ocp_nlp_constraints_bgp_get_ng(void *config_, void *dims_, int* value)
{
    ocp_nlp_constraints_bgp_dims *dims = (ocp_nlp_constraints_bgp_dims *) dims_;
    *value = dims->ng;
}



static void ocp_nlp_constraints_bgp_get_nphi(void *config_, void *dims_, int* value)
{
    ocp_nlp_constraints_bgp_dims *dims = (ocp_nlp_constraints_bgp_dims *) dims_;
    *value = dims->nphi;
}



static void ocp_nlp_constraints_bgp_get_ns(void *config_, void *dims_, int* value)
{
    ocp_nlp_constraints_bgp_dims *dims = (ocp_nlp_constraints_bgp_dims *) dims_;
    *value = dims->ns;
}


static void ocp_nlp_constraints_bgp_get_nsg(void *config_, void *dims_, int* value)
{
    ocp_nlp_constraints_bgp_dims *dims = (ocp_nlp_constraints_bgp_dims *) dims_;
    *value = dims->nsg;
}


static void ocp_nlp_constraints_bgp_get_nr(void *config_, void *dims_, int* value)
{
    ocp_nlp_constraints_bgp_dims *dims = (ocp_nlp_constraints_bgp_dims *) dims_;
    *value = dims->nr;
}



static void ocp_nlp_constraints_bgp_get_nsphi(void *config_, void *dims_, int* value)
{
    ocp_nlp_constraints_bgp_dims *dims = (ocp_nlp_constraints_bgp_dims *) dims_;
    *value = dims->nsphi;
}



static void ocp_nlp_constraints_bgp_get_nbxe(void *config_, void *dims_, int* value)
{
    ocp_nlp_constraints_bgp_dims *dims = (ocp_nlp_constraints_bgp_dims *) dims_;
    *value = dims->nbxe;
}



static void ocp_nlp_constraints_bgp_get_nbue(void *config_, void *dims_, int* value)
{
    ocp_nlp_constraints_bgp_dims *dims = (ocp_nlp_constraints_bgp_dims *) dims_;
    *value = dims->nbue;
}



static void ocp_nlp_constraints_bgp_get_nge(void *config_, void *dims_, int* value)
{
    ocp_nlp_constraints_bgp_dims *dims = (ocp_nlp_constraints_bgp_dims *) dims_;
    *value = dims->nge;
}



static void ocp_nlp_constraints_bgp_get_nphie(void *config_, void *dims_, int* value)
{
    ocp_nlp_constraints_bgp_dims *dims = (ocp_nlp_constraints_bgp_dims *) dims_;
    *value = dims->nphie;
}



void ocp_nlp_constraints_bgp_dims_get(void *config_, void *dims_, const char *field, int* value)
{

    if (!strcmp(field, "ni"))
    {
        ocp_nlp_constraints_bgp_get_ni(config_, dims_, value);
    }
    else if (!strcmp(field, "nb"))
    {
        ocp_nlp_constraints_bgp_get_nb(config_, dims_, value);
    }
    else if (!strcmp(field, "nbx"))
    {
        ocp_nlp_constraints_bgp_get_nbx(config_, dims_, value);
    }
    else if (!strcmp(field, "nbu"))
    {
        ocp_nlp_constraints_bgp_get_nbu(config_, dims_, value);
    }
    else if (!strcmp(field, "ng"))
    {
        ocp_nlp_constraints_bgp_get_ng(config_, dims_, value);
    }
    else if (!strcmp(field, "nphi"))
    {
        ocp_nlp_constraints_bgp_get_nphi(config_, dims_, value);
    }
    else if (!strcmp(field, "ns"))
    {
        ocp_nlp_constraints_bgp_get_ns(config_, dims_, value);
    }
    else if (!strcmp(field, "nsphi"))
    {
        ocp_nlp_constraints_bgp_get_nsphi(config_, dims_, value);
    }
    else if (!strcmp(field, "nsg"))
    {
        ocp_nlp_constraints_bgp_get_nsg(config_, dims_, value);
    }
    else if (!strcmp(field, "nr"))
    {
        ocp_nlp_constraints_bgp_get_nr(config_, dims_, value);
    }
    else if (!strcmp(field, "ng_qp_solver"))
    {
        int ng, nphi;
        ocp_nlp_constraints_bgp_get_ng(config_, dims_, &ng);
        ocp_nlp_constraints_bgp_get_nphi(config_, dims_, &nphi);
        *value = ng + nphi;
    }
    else if (!strcmp(field, "nsg_qp_solver"))
    {
        int nsg, nsphi;
        ocp_nlp_constraints_bgp_get_nsg(config_, dims_, &nsg);
        ocp_nlp_constraints_bgp_get_nsphi(config_, dims_, &nsphi);
        *value = nsg + nsphi;
    }
    else if (!strcmp(field, "nbxe"))
    {
        ocp_nlp_constraints_bgp_get_nbxe(config_, dims_, value);
    }
    else if (!strcmp(field, "nbue"))
    {
        ocp_nlp_constraints_bgp_get_nbue(config_, dims_, value);
    }
    else if (!strcmp(field, "nge"))
    {
        ocp_nlp_constraints_bgp_get_nge(config_, dims_, value);
    }
    else if (!strcmp(field, "nphie"))
    {
        ocp_nlp_constraints_bgp_get_nphie(config_, dims_, value);
    }
    else if (!strcmp(field, "nge_qp_solver"))
    {
        int nge, nphie;
        ocp_nlp_constraints_bgp_get_nge(config_, dims_, &nge);
        ocp_nlp_constraints_bgp_get_nphie(config_, dims_, &nphie);
        *value = nge + nphie;
    }
    else
    {
        printf("error: attempt to get dimension %s from constraint model, that is not there", field);
        exit(1);
    }
}


/* model */

acados_size_t ocp_nlp_constraints_bgp_model_calculate_size(void *config, void *dims_)
{
    ocp_nlp_constraints_bgp_dims *dims = dims_;

    // extract dims
    int nx = dims->nx;
    int nu = dims->nu;
    int nb = dims->nb;
    int ng = dims->ng;
    int nphi = dims->nphi;
    int ns = dims->ns;
    int nbue = dims->nbue;
    int nbxe = dims->nbxe;
    int nge = dims->nge;
    int nphie = dims->nphie;

    acados_size_t size = 0;

    size += sizeof(ocp_nlp_constraints_bgp_model);

    size += sizeof(int) * nb;                                         // idxb
    size += sizeof(int) * ns;                                         // idxs
    size += sizeof(int)*(nbue+nbxe+nge+nphie);                        // idxe
    size += blasfeo_memsize_dvec(2 * nb + 2 * ng + 2 * nphi + 2 * ns);  // d
    size += blasfeo_memsize_dmat(nu + nx, ng);                        // DCt

    size += 64;  // blasfeo_mem align

    return size;
}



void *ocp_nlp_constraints_bgp_model_assign(void *config, void *dims_, void *raw_memory)
{
    ocp_nlp_constraints_bgp_dims *dims = dims_;

    char *c_ptr = (char *) raw_memory;

    // extract sizes
    int nx = dims->nx;
    int nu = dims->nu;
    int nb = dims->nb;
    int ng = dims->ng;
    int nphi = dims->nphi;
    int ns = dims->ns;
    int nbue = dims->nbue;
    int nbxe = dims->nbxe;
    int nge = dims->nge;
    int nphie = dims->nphie;

	int ii;

    // struct
    ocp_nlp_constraints_bgp_model *model = (ocp_nlp_constraints_bgp_model *) c_ptr;
    c_ptr += sizeof(ocp_nlp_constraints_bgp_model);

    // dims
    //  model->dims = dims;

    // blasfeo_mem align
    align_char_to(64, &c_ptr);

    // blasfeo_dmat
    // DCt
    assign_and_advance_blasfeo_dmat_mem(nu + nx, ng, &model->DCt, &c_ptr);

    // blasfeo_dvec
    // d
    assign_and_advance_blasfeo_dvec_mem(2 * nb + 2 * ng + 2 * nphi + 2 * ns, &model->d, &c_ptr);
    // default initialization to zero
    blasfeo_dvecse(2*nb+2*ng+2*nphi+2*ns, 0.0, &model->d, 0);

    // int
    // idxb
    assign_and_advance_int(nb, &model->idxb, &c_ptr);
    // idxs
    assign_and_advance_int(ns, &model->idxs, &c_ptr);
    // idxe
    assign_and_advance_int(nbue+nbxe+nge+nphie, &model->idxe, &c_ptr);

    // h
    //  model->nl_constr_phi_o_r_fun_phi_jac_ux_z_phi_hess_r_jac_ux = NULL;

	// default initialization
	for(ii=0; ii<nbue+nbxe+nge+nphie; ii++)
		model->idxe[ii] = 0;

    // assert
    assert((char *) raw_memory + ocp_nlp_constraints_bgp_model_calculate_size(config, dims) >=
           c_ptr);

    return model;
}


int ocp_nlp_constraints_bgp_model_set(void *config_, void *dims_,
                         void *model_, const char *field, void *value)
{
    // NOTE(oj): this is adapted from the bgh module, maybe something has to be changed here.
    ocp_nlp_constraints_bgp_dims *dims = (ocp_nlp_constraints_bgp_dims *) dims_;
    ocp_nlp_constraints_bgp_model *model = (ocp_nlp_constraints_bgp_model *) model_;

    int ii;
    int *ptr_i;

    if (!dims || !model || !field || !value)
    {
        printf("ocp_nlp_constraints_bgp_model_set: got Null pointer \n");
        exit(1);
    }

    int nu = dims->nu;
    int nx = dims->nx;
    int nb = dims->nb;
    int ng = dims->ng;
    int nphi = dims->nphi;
    int ns = dims->ns;
    int nsbu = dims->nsbu;
    int nsbx = dims->nsbx;
    int nsg = dims->nsg;
    int nsphi = dims->nsphi;
    int nbx = dims->nbx;
    int nbu = dims->nbu;
    int nbue = dims->nbue;
    int nbxe = dims->nbxe;
    int nge = dims->nge;
    int nphie = dims->nphie;

    if (!strcmp(field, "lb")) // TODO(fuck_lint) remove !!!
    {
        blasfeo_pack_dvec(nb, value, 1, &model->d, 0);
    }
    else if (!strcmp(field, "ub")) // TODO(fuck_lint) remove !!!
    {
        blasfeo_pack_dvec(nb, value, 1, &model->d, nb+ng+nphi);
    }
    else if (!strcmp(field, "idxbx"))
    {
        ptr_i = (int *) value;
        for (ii=0; ii < nbx; ii++)
            model->idxb[nbu+ii] = nu+ptr_i[ii];
    }
    else if (!strcmp(field, "lbx"))
    {
        blasfeo_pack_dvec(nbx, value, 1, &model->d, nbu);
    }
    else if (!strcmp(field, "ubx"))
    {
        blasfeo_pack_dvec(nbx, value, 1, &model->d, nb + ng + nphi + nbu);
    }
    else if (!strcmp(field, "idxbu"))
    {
        ptr_i = (int *) value;
        for (ii=0; ii < nbu; ii++)
            model->idxb[ii] = ptr_i[ii];
    }
    else if (!strcmp(field, "lbu"))
    {
        blasfeo_pack_dvec(nbu, value, 1, &model->d, 0);
    }
    else if (!strcmp(field, "ubu"))
    {
        blasfeo_pack_dvec(nbu, value, 1, &model->d, nb + ng + nphi);
    }
    else if (!strcmp(field, "C"))
    {
        blasfeo_pack_tran_dmat(ng, nx, value, ng, &model->DCt, nu, 0);
    }
    else if (!strcmp(field, "D"))
    {
        blasfeo_pack_tran_dmat(ng, nu, value, ng, &model->DCt, 0, 0);
    }
    else if (!strcmp(field, "lg"))
    {
        blasfeo_pack_dvec(ng, value, 1, &model->d, nb);
    }
    else if (!strcmp(field, "ug"))
    {
        blasfeo_pack_dvec(ng, value, 1, &model->d, 2*nb+ng+nphi);
    }
    else if (!strcmp(field, "nl_constr_phi_o_r_fun_phi_jac_ux_z_phi_hess_r_jac_ux"))
    {
        model->nl_constr_phi_o_r_fun_phi_jac_ux_z_phi_hess_r_jac_ux = value;
    }
    else if (!strcmp(field, "nl_constr_phi_o_r_fun"))
    {
        model->nl_constr_phi_o_r_fun = value;
    }
    else if (!strcmp(field, "lphi")) // TODO(fuck_lint) remove
    {
        blasfeo_pack_dvec(nphi, value, 1, &model->d, nb+ng);
    }
    else if (!strcmp(field, "uphi"))
    {
        blasfeo_pack_dvec(nphi, value, 1, &model->d, 2*nb+2*ng+nphi);
    }
    else if (!strcmp(field, "idxsbu"))
    {
        ptr_i = (int *) value;
        for (ii=0; ii < nsbu; ii++)
            model->idxs[ii] = ptr_i[ii];
    }
    else if (!strcmp(field, "lsbu"))
    {
        blasfeo_pack_dvec(nsbu, value, 1, &model->d, 2*nb+2*ng+2*nphi);
    }
    else if (!strcmp(field, "usbu"))
    {
        blasfeo_pack_dvec(nsbu, value, 1, &model->d, 2*nb+2*ng+2*nphi+ns);
    }
    else if (!strcmp(field, "idxsbx"))
    {
        ptr_i = (int *) value;
        for (ii=0; ii < nsbx; ii++)
            model->idxs[nsbu+ii] = nbu+ptr_i[ii];
    }
    else if (!strcmp(field, "lsbx"))
    {
        blasfeo_pack_dvec(nsbx, value, 1, &model->d, 2*nb+2*ng+2*nphi+nsbu);
    }
    else if (!strcmp(field, "usbx"))
    {
        blasfeo_pack_dvec(nsbx, value, 1, &model->d, 2*nb+2*ng+2*nphi+ns+nsbu);
    }
    else if (!strcmp(field, "idxsg"))
    {
        ptr_i = (int *) value;
        for (ii=0; ii < nsg; ii++)
            model->idxs[nsbu+nsbx+ii] = nbu+nbx+ptr_i[ii];
    }
    else if (!strcmp(field, "lsg"))
    {
        blasfeo_pack_dvec(nsg, value, 1, &model->d, 2*nb+2*ng+2*nphi+nsbu+nsbx);
    }
    else if (!strcmp(field, "usg"))
    {
        blasfeo_pack_dvec(nsg, value, 1, &model->d, 2*nb+2*ng+2*nphi+ns+nsbu+nsbx);
    }
    else if (!strcmp(field, "idxsphi"))
    {
        ptr_i = (int *) value;
        for (ii=0; ii < nsphi; ii++)
            model->idxs[nsbu+nsbx+nsg+ii] = nbu+nbx+ng+ptr_i[ii];
    }
    else if (!strcmp(field, "lsphi"))
    {
        blasfeo_pack_dvec(nsphi, value, 1, &model->d, 2*nb+2*ng+2*nphi+nsbu+nsbx+nsg);
    }
    else if (!strcmp(field, "usphi"))
    {
        blasfeo_pack_dvec(nsphi, value, 1, &model->d, 2*nb+2*ng+2*nphi+ns+nsbu+nsbx+nsg);
    }
    else if (!strcmp(field, "idxbue"))
    {
        ptr_i = (int *) value;
        for (ii=0; ii < nbue; ii++)
            model->idxe[ii] = ptr_i[ii];
    }
    else if (!strcmp(field, "idxbxe"))
    {
        ptr_i = (int *) value;
        for (ii=0; ii < nbxe; ii++)
            model->idxe[nbue+ii] = nbu+ptr_i[ii];
    }
    else if (!strcmp(field, "idxge"))
    {
        ptr_i = (int *) value;
        for (ii=0; ii < nge; ii++)
            model->idxe[nbue+nbxe+ii] = nbu+nbx+ptr_i[ii];
    }
    else if (!strcmp(field, "idxphie"))
    {
        ptr_i = (int *) value;
        for (ii=0; ii < nphie; ii++)
            model->idxe[nbue+nbxe+nge+ii] = nbu+nbx+ng+ptr_i[ii];
    }
    else
    {
        printf("\nerror: model field not available in module ocp_nlp_constraints_bgp: %s\n",
            field);
        exit(1);
    }

    return ACADOS_SUCCESS;
}



/* options */

acados_size_t ocp_nlp_constraints_bgp_opts_calculate_size(void *config_, void *dims_)
{
    acados_size_t size = 0;

    size += sizeof(ocp_nlp_constraints_bgp_opts);

    return size;
}



void *ocp_nlp_constraints_bgp_opts_assign(void *config_, void *dims_, void *raw_memory)
{
    char *c_ptr = (char *) raw_memory;

    ocp_nlp_constraints_bgp_opts *opts = (ocp_nlp_constraints_bgp_opts *) c_ptr;
    c_ptr += sizeof(ocp_nlp_constraints_bgp_opts);

    assert((char *) raw_memory + ocp_nlp_constraints_bgp_opts_calculate_size(config_, dims_) >=
           c_ptr);

    return opts;
}



void ocp_nlp_constraints_bgp_opts_initialize_default(void *config_, void *dims_, void *opts_)
{
    ocp_nlp_constraints_bgp_opts *opts = opts_;

    opts->compute_adj = 1;
    opts->compute_hess = 0;

    return;
}



void ocp_nlp_constraints_bgp_opts_update(void *config_, void *dims_, void *opts_)
{
    //  ocp_nlp_constraints_bgp_opts *opts = opts_;

    return;
}



void ocp_nlp_constraints_bgp_opts_set(void *config_, void *opts_, char *field, void *value)
{

    ocp_nlp_constraints_bgp_opts *opts = opts_;

    if(!strcmp(field, "compute_adj"))
    {
        int *compute_adj = value;
        opts->compute_adj = *compute_adj;
    }
    else if(!strcmp(field, "compute_hess"))
    {
        int *compute_hess = value;
        opts->compute_hess = *compute_hess;
    }
    else
    {
        printf("\nerror: field %s not available in ocp_nlp_constraints_bgp_opts_set\n", field);
        exit(1);
    }

    return;

}



/* memory */

acados_size_t ocp_nlp_constraints_bgp_memory_calculate_size(void *config_, void *dims_, void *opts_)
{
    ocp_nlp_constraints_bgp_dims *dims = dims_;

    // extract dims
    int nx = dims->nx;
    int nu = dims->nu;
    int nb = dims->nb;
    int ng = dims->ng;
    int nphi = dims->nphi;
    int ns = dims->ns;

    acados_size_t size = 0;

    size += sizeof(ocp_nlp_constraints_bgp_memory);

    size += 1 * blasfeo_memsize_dvec(2 * nb + 2 * ng + 2 * nphi + 2 * ns);  // fun
    size += 1 * blasfeo_memsize_dvec(nu + nx + 2 * ns);                   // adj

    size += 1 * 64;  // blasfeo_mem align

    return size;
}



void *ocp_nlp_constraints_bgp_memory_assign(void *config_, void *dims_, void *opts_,
                                            void *raw_memory)
{
    ocp_nlp_constraints_bgp_dims *dims = dims_;

    char *c_ptr = (char *) raw_memory;

    // extract dims
    int nx = dims->nx;
    int nu = dims->nu;
    int nb = dims->nb;
    int ng = dims->ng;
    int nphi = dims->nphi;
    int ns = dims->ns;

    // struct
    ocp_nlp_constraints_bgp_memory *memory = (ocp_nlp_constraints_bgp_memory *) c_ptr;
    c_ptr += sizeof(ocp_nlp_constraints_bgp_memory);

    // blasfeo_mem align
    align_char_to(64, &c_ptr);

    // fun
    assign_and_advance_blasfeo_dvec_mem(2 * nb + 2 * ng + 2 * nphi + 2 * ns, &memory->fun, &c_ptr);
    // adj
    assign_and_advance_blasfeo_dvec_mem(nu + nx + 2 * ns, &memory->adj, &c_ptr);

    assert((char *) raw_memory +
               ocp_nlp_constraints_bgp_memory_calculate_size(config_, dims, opts_) >=
           c_ptr);

    return memory;
}



struct blasfeo_dvec *ocp_nlp_constraints_bgp_memory_get_fun_ptr(void *memory_)
{
    ocp_nlp_constraints_bgp_memory *memory = memory_;

    return &memory->fun;
}



struct blasfeo_dvec *ocp_nlp_constraints_bgp_memory_get_adj_ptr(void *memory_)
{
    ocp_nlp_constraints_bgp_memory *memory = memory_;

    return &memory->adj;
}



void ocp_nlp_constraints_bgp_memory_set_ux_ptr(struct blasfeo_dvec *ux, void *memory_)
{
    ocp_nlp_constraints_bgp_memory *memory = memory_;

    memory->ux = ux;
}



void ocp_nlp_constraints_bgp_memory_set_tmp_ux_ptr(struct blasfeo_dvec *tmp_ux, void *memory_)
{
    ocp_nlp_constraints_bgp_memory *memory = memory_;

    memory->tmp_ux = tmp_ux;
}



void ocp_nlp_constraints_bgp_memory_set_lam_ptr(struct blasfeo_dvec *lam, void *memory_)
{
    ocp_nlp_constraints_bgp_memory *memory = memory_;

    memory->lam = lam;
}



void ocp_nlp_constraints_bgp_memory_set_tmp_lam_ptr(struct blasfeo_dvec *tmp_lam, void *memory_)
{
    ocp_nlp_constraints_bgp_memory *memory = memory_;

    memory->tmp_lam = tmp_lam;
}



void ocp_nlp_constraints_bgp_memory_set_DCt_ptr(struct blasfeo_dmat *DCt, void *memory_)
{
    ocp_nlp_constraints_bgp_memory *memory = memory_;

    memory->DCt = DCt;
}



void ocp_nlp_constraints_bgp_memory_set_RSQrq_ptr(struct blasfeo_dmat *RSQrq, void *memory_)
{
    ocp_nlp_constraints_bgp_memory *memory = memory_;

    memory->RSQrq = RSQrq;
}



void ocp_nlp_constraints_bgp_memory_set_z_alg_ptr(struct blasfeo_dvec *z_alg, void *memory_)
{
    ocp_nlp_constraints_bgp_memory *memory = memory_;

    memory->z_alg = z_alg;
}



void ocp_nlp_constraints_bgp_memory_set_dzduxt_ptr(struct blasfeo_dmat *dzduxt, void *memory_)
{
    ocp_nlp_constraints_bgp_memory *memory = memory_;

    memory->dzduxt = dzduxt;
}




void ocp_nlp_constraints_bgp_memory_set_idxb_ptr(int *idxb, void *memory_)
{
    ocp_nlp_constraints_bgp_memory *memory = memory_;

    memory->idxb = idxb;
}



void ocp_nlp_constraints_bgp_memory_set_idxs_rev_ptr(int *idxs_rev, void *memory_)
{
    ocp_nlp_constraints_bgp_memory *memory = memory_;

    memory->idxs_rev = idxs_rev;
}



void ocp_nlp_constraints_bgp_memory_set_idxe_ptr(int *idxe, void *memory_)
{
    ocp_nlp_constraints_bgp_memory *memory = memory_;

    memory->idxe = idxe;
}



/* workspace */

acados_size_t ocp_nlp_constraints_bgp_workspace_calculate_size(void *config_, void *dims_, void *opts_)
{
    ocp_nlp_constraints_bgp_dims *dims = dims_;

    // extract dims
    int nx = dims->nx;
    int nu = dims->nu;
    int nz = dims->nz;
    int nb = dims->nb;
    int ng = dims->ng;
    int nphi = dims->nphi;
    int ns = dims->ns;
    int nr = dims->nr;

    int nv = nx + nu;

    acados_size_t size = 0;

    size += sizeof(ocp_nlp_constraints_bgp_workspace);

    size += 1 * blasfeo_memsize_dvec(nb + ng + nphi + ns); // tmp_ni
    size += nr * (nx + nu) * sizeof(double);
    size += 1 * blasfeo_memsize_dmat(nx + nu, nr);
    size += 1 * blasfeo_memsize_dmat(nr * nphi, nr);       // tmp_nr_nphi_nr
    size += 1 * blasfeo_memsize_dmat(nv, nr);              // tmp_nv_nr
    size += 1 * blasfeo_memsize_dmat(nz, nphi);            // tmp_nz_nphi
    size += 1 * blasfeo_memsize_dmat(nv, nphi);            // tmp_nv_nphi

    size += 2 * 64;  // blasfeo_mem align

    return size;
}



static void ocp_nlp_constraints_bgp_cast_workspace(void *config_, void *dims_, void *opts_,
                                                   void *work_)
{
    ocp_nlp_constraints_bgp_dims *dims = dims_;
    ocp_nlp_constraints_bgp_workspace *work = work_;

    // extract dims
    int nx = dims->nx;
    int nu = dims->nu;
    int nz = dims->nz;
    int nb = dims->nb;
    int ng = dims->ng;
    int nphi = dims->nphi;
    int ns = dims->ns;
    int nr = dims->nr;

    int nv = nu + nx;

    char *c_ptr = (char *) work_;
    c_ptr += sizeof(ocp_nlp_constraints_bgp_workspace);

    // blasfeo_mem align
    align_char_to(64, &c_ptr);

    // tmp_ni
    assign_and_advance_blasfeo_dvec_mem(nb + ng + nphi + ns, &work->tmp_ni, &c_ptr);
    c_ptr += nr * (nx + nu) * sizeof(double);
    align_char_to(64, &c_ptr);
    assign_and_advance_blasfeo_dmat_mem(nx + nu, nr, &work->jac_r_ux_tran, &c_ptr);
    assign_and_advance_blasfeo_dmat_mem(nr * nphi, nr, &work->tmp_nr_nphi_nr, &c_ptr);
    assign_and_advance_blasfeo_dmat_mem(nv, nr, &work->tmp_nv_nr, &c_ptr);
    assign_and_advance_blasfeo_dmat_mem(nz, nphi, &work->tmp_nz_nphi, &c_ptr);
    assign_and_advance_blasfeo_dmat_mem(nv, nphi, &work->tmp_nv_nphi, &c_ptr);
    assert((char *) work + ocp_nlp_constraints_bgp_workspace_calculate_size(config_, dims, opts_)
           >= c_ptr);

    return;
}



/* functions */

void ocp_nlp_constraints_bgp_initialize(void *config_, void *dims_, void *model_, void *opts,
                                        void *memory_, void *work_)
{
    ocp_nlp_constraints_bgp_dims *dims = dims_;
    ocp_nlp_constraints_bgp_model *model = model_;
    ocp_nlp_constraints_bgp_memory *memory = memory_;

    // loop index
    int j;

    // extract dims
    int nx = dims->nx;
    int nu = dims->nu;
    int nb = dims->nb;
    int ng = dims->ng;
    int ns = dims->ns;
    int nbue = dims->nbue;
    int nbxe = dims->nbxe;
    int nge = dims->nge;
    int nphie = dims->nphie;

    // initialize idxb
    for (j = 0; j < nb; j++)
    {
        memory->idxb[j] = model->idxb[j];
    }

    // initialize idxs_rev
    for (j = 0; j < ns; j++)
    {
//        memory->idxs[j] = model->idxs[j];
        memory->idxs_rev[model->idxs[j]] = j;
    }

    // initialize idxe
    for (j = 0; j < nbue+nbxe+nge+nphie; j++)
    {
        memory->idxe[j] = model->idxe[j];
    }

    // initialize general constraints matrix
    blasfeo_dgecp(nu + nx, ng, &model->DCt, 0, 0, memory->DCt, 0, 0);

    return;
}



void ocp_nlp_constraints_bgp_update_qp_matrices(void *config_, void *dims_, void *model_,
                                                void *opts_, void *memory_, void *work_)
{
    ocp_nlp_constraints_bgp_dims *dims = dims_;
    ocp_nlp_constraints_bgp_model *model = model_;
    ocp_nlp_constraints_bgp_opts *opts = opts_;
    ocp_nlp_constraints_bgp_memory *memory = memory_;
    ocp_nlp_constraints_bgp_workspace *work = work_;

    ocp_nlp_constraints_bgp_cast_workspace(config_, dims, opts_, work_);

    // extract dims
    int nx = dims->nx;
    int nu = dims->nu;
    int nz = dims->nz;
    int nb = dims->nb;
    int ng = dims->ng;
    int nphi = dims->nphi;
    int ns = dims->ns;
    int nr = dims->nr;

    int nv = nx + nu;

    ext_fun_arg_t ext_fun_type_in[3];
    void *ext_fun_in[3];
    ext_fun_arg_t ext_fun_type_out[5];
    void *ext_fun_out[5];

    // box
    blasfeo_dvecex_sp(nb, 1.0, model->idxb, memory->ux, 0, &work->tmp_ni, 0);

    // general linear
    blasfeo_dgemv_t(nu + nx, ng, 1.0, memory->DCt, 0, 0, memory->ux, 0, 0.0, &work->tmp_ni, nb,
                    &work->tmp_ni, nb);

    // TODO(andrea): nz > 0 supported, but Hessian contribution associated with algebraic variables is neglected. 

    // nonlinear
    if (nphi > 0)
    {
        struct blasfeo_dvec_args x_in;  // input x of external fun;
        x_in.x = memory->ux;
        x_in.xi = nu;

        struct blasfeo_dvec_args u_in;  // input u of external fun;
        u_in.x = memory->ux;
        u_in.xi = 0;

        struct blasfeo_dvec_args z_in;  // input z of external fun;
        z_in.x = memory->z_alg;
        z_in.xi = 0;

        ext_fun_type_in[0] = BLASFEO_DVEC_ARGS;
        ext_fun_in[0] = &x_in;
        ext_fun_type_in[1] = BLASFEO_DVEC_ARGS;
        ext_fun_in[1] = &u_in;
        ext_fun_type_in[2] = BLASFEO_DVEC_ARGS;
        ext_fun_in[2] = &z_in;

        struct blasfeo_dvec_args fun_out;
        fun_out.x = &work->tmp_ni;
        fun_out.xi = nb + ng;
        ext_fun_type_out[0] = BLASFEO_DVEC_ARGS;
        ext_fun_out[0] = &fun_out;  // fun: nphi

        struct blasfeo_dmat_args jac_phi_tran_out;
        jac_phi_tran_out.A = memory->DCt;
        jac_phi_tran_out.ai = 0;
        jac_phi_tran_out.aj = ng;
        ext_fun_type_out[1] = BLASFEO_DMAT_ARGS;
        ext_fun_out[1] = &jac_phi_tran_out;  // jac': (nu+nx) * nphi

        struct blasfeo_dmat_args jac_phi_z_tran_out; // Jacobian dphidz treated separately
        jac_phi_z_tran_out.A = &work->tmp_nz_nphi;
        jac_phi_z_tran_out.ai = 0;
        jac_phi_z_tran_out.aj = 0;
        ext_fun_type_out[2] = BLASFEO_DMAT_ARGS;
        ext_fun_out[2] = &jac_phi_z_tran_out;  // jac': nz * nphi

        struct blasfeo_dmat_args hess_out;
        hess_out.A = &work->tmp_nr_nphi_nr;
        hess_out.ai = 0;
        hess_out.aj = 0;
        ext_fun_type_out[3] = BLASFEO_DMAT_ARGS;
        ext_fun_out[3] = &hess_out;  // hess: nphi * nr * nr

        struct blasfeo_dmat_args jac_r_tran_out;
        jac_r_tran_out.A = &work->jac_r_ux_tran;
        jac_r_tran_out.ai = 0;
        jac_r_tran_out.aj = 0;
        ext_fun_type_out[4] = BLASFEO_DMAT_ARGS;
        ext_fun_out[4] = &jac_r_tran_out;  // jac': (nu+nx) * nr

        model->nl_constr_phi_o_r_fun_phi_jac_ux_z_phi_hess_r_jac_ux->evaluate(
                model->nl_constr_phi_o_r_fun_phi_jac_ux_z_phi_hess_r_jac_ux, 
                ext_fun_type_in, ext_fun_in, ext_fun_type_out, ext_fun_out);
        
        // expand phi:
        // phi(x, u, z) ~
        // phi(\bar{x}, \bar{u}, \bar{z}) + 
        // dphidx*(x - \bar{x}) + 
        // dphidu*(u - \bar{u}) + 
        // dphidz*(z - \bar{z}) =
        //
        // phi(\bar{x}, \bar{u}, \bar{z}) - dphidz*dzdx*\bar{x} - dphidz*dzdu*\bar{u} + 
        // (dphidx + dphidz*dzdx)*(x - \bar{x}) +  
        // (dphidu + dphidz*dzdu)*(u - \bar{u})  
        
        // update DCt
        blasfeo_dgemm_nn(nu+nx, nphi, nz, 1.0, memory->dzduxt, 0, 0, &work->tmp_nz_nphi, 0, 0, 0.0, &work->tmp_nv_nphi, 0, 0, &work->tmp_nv_nphi, 0, 0);
        blasfeo_dgead(nu+nx, nphi, 1.0, &work->tmp_nv_nphi, 0, 0, memory->DCt, ng, 0);	
        // update memory->fun	
        blasfeo_dgemv_t(nu+nx, nphi, -1.0, &work->tmp_nv_nphi, 0, 0, memory->ux, 0, 1.0, &memory->fun, 0, &memory->fun, 0);
    }

    // add SCQP Hessian contribution
    for (int i = 0; i < nphi; i++) { 
        double lam_i = blasfeo_dvecex1(memory->lam, 
                2 * (nb + ng) + nphi + i);

        blasfeo_dgemm_nt(nv, nr, nr, lam_i, &work->jac_r_ux_tran, 
                0, 0, &work->tmp_nr_nphi_nr, nr * i, 0, 0.0, 
                &work->tmp_nv_nr, 0, 0, &work->tmp_nv_nr, 0, 0);

        blasfeo_dsyrk_ln(nv, nr, 1.0, &work->tmp_nv_nr, 0, 0,
                &work->jac_r_ux_tran, 0, 0, 1.0, memory->RSQrq,
                0, 0, memory->RSQrq, 0, 0);
    }

    blasfeo_daxpy(nb + ng + nphi, -1.0, &work->tmp_ni, 0, &model->d, 0, 
            &memory->fun, 0);
    blasfeo_daxpy(nb + ng + nphi, -1.0, &model->d, nb + ng + nphi, 
            &work->tmp_ni, 0, &memory->fun, nb + ng + nphi);

    // soft
    blasfeo_dvecad_sp(ns, -1.0, memory->ux, nu + nx, model->idxs, &memory->fun, 0);
    blasfeo_dvecad_sp(ns, -1.0, memory->ux, nu + nx + ns, model->idxs, &memory->fun, nb + ng + nphi);

    blasfeo_daxpy(2 * ns, -1.0, memory->ux, nu + nx, &model->d, 2 * nb + 2 * ng + 2 * nphi,
                  &memory->fun, 2 * nb + 2 * ng + 2 * nphi);

    // nlp_mem: ineq_adj
    if (opts->compute_adj)
    {
        blasfeo_dvecse(nu + nx + 2 * ns, 0.0, &memory->adj, 0);
        blasfeo_daxpy(nb+ng+nphi, -1.0, memory->lam, nb + ng + nphi, memory->lam, 0, &work->tmp_ni, 0);
        blasfeo_dvecad_sp(nb, 1.0, &work->tmp_ni, 0, model->idxb, &memory->adj, 0);
        blasfeo_dgemv_n(nu+nx, ng+nphi, 1.0, memory->DCt, 0, 0, &work->tmp_ni, nb, 1.0, &memory->adj,
                        0, &memory->adj, 0);
        // soft
        blasfeo_dvecex_sp(ns, 1.0, model->idxs, memory->lam, 0, &memory->adj, nu + nx);
        blasfeo_dvecex_sp(ns, 1.0, model->idxs, memory->lam, nb+ng+nphi, &memory->adj, nu+nx+ns);
        blasfeo_daxpy(2 * ns, 1.0, memory->lam, 2 * nb + 2 * ng + 2 * nphi, &memory->adj, nu + nx,
                      &memory->adj, nu + nx);
    }

    if (opts->compute_hess)
    {
        printf("\nerror: compute_hess!=0 not supported (yet) in ocp_nlp_constraints_bgp\n");
        exit(1);
    }

    return;
}



void ocp_nlp_constraints_bgp_compute_fun(void *config_, void *dims_, void *model_, void *opts_, void *memory_, void *work_)
{
    ocp_nlp_constraints_bgp_dims *dims = dims_;
    ocp_nlp_constraints_bgp_model *model = model_;
    // ocp_nlp_constraints_bgp_opts *opts = opts_;
    ocp_nlp_constraints_bgp_memory *memory = memory_;
    ocp_nlp_constraints_bgp_workspace *work = work_;

    ocp_nlp_constraints_bgp_cast_workspace(config_, dims, opts_, work_);

    // extract dims
    int nx = dims->nx;
    int nu = dims->nu;
    int nz = dims->nz;
    int nb = dims->nb;
    int ng = dims->ng;
    int nphi = dims->nphi;
    int ns = dims->ns;

    ext_fun_arg_t ext_fun_type_in[3];
    void *ext_fun_in[3];
    ext_fun_arg_t ext_fun_type_out[3];
    void *ext_fun_out[3];

    // box
    blasfeo_dvecex_sp(nb, 1.0, model->idxb, memory->tmp_ux, 0, &work->tmp_ni, 0);

    // general linear
    blasfeo_dgemv_t(nu+nx, ng, 1.0, memory->DCt, 0, 0, memory->tmp_ux, 0, 0.0, &work->tmp_ni, nb, &work->tmp_ni, nb);

    // nonlinear
    if (nphi > 0)
    {
		if(nz > 0)
		{
			// TODO
			printf("\nerror: ocp_nlp_constraints_bgp_compute_fun: not implemented yet for nz>0\n");
			exit(1);
		}

        struct blasfeo_dvec_args x_in;  // input x of external fun;
        x_in.x = memory->tmp_ux;
        x_in.xi = nu;

        struct blasfeo_dvec_args u_in;  // input u of external fun;
        u_in.x = memory->tmp_ux;
        u_in.xi = 0;

		// TODO tmp_z_alg !!!
        struct blasfeo_dvec_args z_in;  // input z of external fun;
        z_in.x = memory->z_alg;
        z_in.xi = 0;

        struct blasfeo_dvec_args fun_out;
        fun_out.x = &work->tmp_ni;
        fun_out.xi = nb + ng;

		ext_fun_type_in[0] = BLASFEO_DVEC_ARGS;
		ext_fun_in[0] = &x_in;
		ext_fun_type_in[1] = BLASFEO_DVEC_ARGS;
		ext_fun_in[1] = &u_in;
		ext_fun_type_in[2] = BLASFEO_DVEC_ARGS;
		ext_fun_in[2] = &z_in;

		ext_fun_type_out[0] = BLASFEO_DVEC_ARGS;
		ext_fun_out[0] = &fun_out;  // fun: nphi

        if (model->nl_constr_phi_o_r_fun == 0)
        {
            printf("ocp_nlp_constraints_bgp_compute_fun: nl_constr_phi_o_r_fun is not provided. Exiting.\n");
            exit(1);
        }
		model->nl_constr_phi_o_r_fun->evaluate(model->nl_constr_phi_o_r_fun, ext_fun_type_in, ext_fun_in, ext_fun_type_out, ext_fun_out);
	}

    blasfeo_daxpy(nb+ng+nphi, -1.0, &work->tmp_ni, 0, &model->d, 0, &memory->fun, 0);
    blasfeo_daxpy(nb+ng+nphi, -1.0, &model->d, nb+ng+nphi, &work->tmp_ni, 0, &memory->fun, nb+ng+nphi);

    // soft
    blasfeo_dvecad_sp(ns, -1.0, memory->ux, nu+nx, model->idxs, &memory->fun, 0);
    blasfeo_dvecad_sp(ns, -1.0, memory->ux, nu+nx+ns, model->idxs, &memory->fun, nb+ng+nphi);

    blasfeo_daxpy(2*ns, -1.0, memory->ux, nu+nx, &model->d, 2*nb+2*ng+2*nphi, &memory->fun, 2*nb+2*ng+2*nphi);

    return;

}



void ocp_nlp_constraints_bgp_bounds_update(void *config_, void *dims_, void *model_,
                                            void *opts_, void *memory_, void *work_)
{
    
    ocp_nlp_constraints_bgp_dims *dims = dims_;
    ocp_nlp_constraints_bgp_model *model = model_;
    // ocp_nlp_constraints_bgp_opts *opts = opts_;
    ocp_nlp_constraints_bgp_memory *memory = memory_;
    ocp_nlp_constraints_bgp_workspace *work = work_;

    ocp_nlp_constraints_bgp_cast_workspace(config_, dims, opts_, work_);

    // extract dims
    // int nx = dims->nx;
    // int nu = dims->nu;
    int nb = dims->nb;
    int ng = dims->ng;
    int nphi = dims->nphi;

    // box
    blasfeo_dvecex_sp(nb, 1.0, model->idxb, memory->ux, 0, &work->tmp_ni, 0);
    
    blasfeo_daxpy(nb, -1.0, &work->tmp_ni, 0, &model->d, 0, &memory->fun, 0);
    blasfeo_daxpy(nb, -1.0, &model->d, nb+ng+nphi, &work->tmp_ni, 0, &memory->fun, nb+ng+nphi);

    return;
}



void ocp_nlp_constraints_bgp_config_initialize_default(void *config_)
{
    ocp_nlp_constraints_config *config = config_;

    config->dims_calculate_size = &ocp_nlp_constraints_bgp_dims_calculate_size;
    config->dims_assign = &ocp_nlp_constraints_bgp_dims_assign;
    config->dims_initialize = &ocp_nlp_constraints_bgp_dims_initialize;
    config->dims_set = &ocp_nlp_constraints_bgp_dims_set;
    config->dims_get = &ocp_nlp_constraints_bgp_dims_get;
    config->model_calculate_size = &ocp_nlp_constraints_bgp_model_calculate_size;
    config->model_assign = &ocp_nlp_constraints_bgp_model_assign;
    config->model_set = &ocp_nlp_constraints_bgp_model_set;
    config->opts_calculate_size = &ocp_nlp_constraints_bgp_opts_calculate_size;
    config->opts_assign = &ocp_nlp_constraints_bgp_opts_assign;
    config->opts_initialize_default = &ocp_nlp_constraints_bgp_opts_initialize_default;
    config->opts_update = &ocp_nlp_constraints_bgp_opts_update;
    config->opts_set = &ocp_nlp_constraints_bgp_opts_set;
    config->memory_calculate_size = &ocp_nlp_constraints_bgp_memory_calculate_size;
    config->memory_assign = &ocp_nlp_constraints_bgp_memory_assign;
    config->memory_get_fun_ptr = &ocp_nlp_constraints_bgp_memory_get_fun_ptr;
    config->memory_get_adj_ptr = &ocp_nlp_constraints_bgp_memory_get_adj_ptr;
    config->memory_set_ux_ptr = &ocp_nlp_constraints_bgp_memory_set_ux_ptr;
    config->memory_set_tmp_ux_ptr = &ocp_nlp_constraints_bgp_memory_set_tmp_ux_ptr;
    config->memory_set_lam_ptr = &ocp_nlp_constraints_bgp_memory_set_lam_ptr;
    config->memory_set_tmp_lam_ptr = &ocp_nlp_constraints_bgp_memory_set_tmp_lam_ptr;
    config->memory_set_DCt_ptr = &ocp_nlp_constraints_bgp_memory_set_DCt_ptr;
    config->memory_set_RSQrq_ptr = &ocp_nlp_constraints_bgp_memory_set_RSQrq_ptr;
    config->memory_set_z_alg_ptr = &ocp_nlp_constraints_bgp_memory_set_z_alg_ptr;
    config->memory_set_dzdux_tran_ptr = &ocp_nlp_constraints_bgp_memory_set_dzduxt_ptr;
    config->memory_set_idxb_ptr = &ocp_nlp_constraints_bgp_memory_set_idxb_ptr;
    config->memory_set_idxs_rev_ptr = &ocp_nlp_constraints_bgp_memory_set_idxs_rev_ptr;
    config->memory_set_idxe_ptr = &ocp_nlp_constraints_bgp_memory_set_idxe_ptr;
    config->workspace_calculate_size = &ocp_nlp_constraints_bgp_workspace_calculate_size;
    config->initialize = &ocp_nlp_constraints_bgp_initialize;
    config->update_qp_matrices = &ocp_nlp_constraints_bgp_update_qp_matrices;
    config->compute_fun = &ocp_nlp_constraints_bgp_compute_fun;
    config->bounds_update = &ocp_nlp_constraints_bgp_bounds_update;
    config->config_initialize_default = &ocp_nlp_constraints_bgp_config_initialize_default;


    return;
}
