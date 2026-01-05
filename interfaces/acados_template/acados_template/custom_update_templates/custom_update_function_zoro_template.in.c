/*
 * Copyright (c) The acados authors.
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

// This is a template based custom_update function
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>

#include "custom_update_function.h"
#include "acados_solver_{{ model.name }}.h"
#include "acados_c/ocp_nlp_interface.h"
#include "acados/utils/mem.h"

#include "blasfeo_d_aux_ext_dep.h"
#include "blasfeo_d_blasfeo_api.h"


// number of tightened path constraints
#define NCT {{ zoro_description.nlbu_t + zoro_description.nlbx_t + zoro_description.nlg_t + zoro_description.nlh_t + zoro_description.nubu_t + zoro_description.nubx_t + zoro_description.nug_t + zoro_description.nuh_t }}
// number of tightened terminal constraints
#define NCT_E {{zoro_description.nlbx_e_t + zoro_description.nlg_e_t + zoro_description.nlh_e_t + zoro_description.nubx_e_t + zoro_description.nug_e_t + zoro_description.nuh_e_t}}

typedef struct custom_memory
{
    // covariance matrics
    struct blasfeo_dmat *uncertainty_matrix_buffer;      // shape = (N+1) * (nx, nx)
    // covariance matrix of the additive disturbance
    struct blasfeo_dmat W_mat;                           // shape = (nw, nw)
    struct blasfeo_dmat unc_jac_G_mat;                   // shape = (nx, nw)
    struct blasfeo_dmat temp_GW_mat;                     // shape = (nx, nw)
    struct blasfeo_dmat GWG_mat;                         // shape = (nx, nx)
    // NOTE: Covariance matrix of the additive noise (used if input_W_add_diag)
    struct blasfeo_dmat W_stage_mat;  // shape = (nw, nw)
    // sensitivity matrices
    struct blasfeo_dmat A_mat;                           // shape = (nx, nx)
    struct blasfeo_dmat B_mat;                           // shape = (nx, nu)
    // matrix in linear constraints
    struct blasfeo_dmat Cg_mat;                          // shape = (ng, nx)
    struct blasfeo_dmat Dg_mat;                          // shape = (ng, nu)
    struct blasfeo_dmat Cg_e_mat;                        // shape = (ng_e, nx)
    struct blasfeo_dmat dummy_Dgh_e_mat;                 // shape = (ngh_e_max, nu)
    // matrix in nonlinear constraints
    struct blasfeo_dmat Ch_mat;                          // shape = (nh, nx)
    struct blasfeo_dmat Dh_mat;                          // shape = (nh, nu)
    struct blasfeo_dmat Ch_e_mat;                        // shape = (nh_e, nx)
    // feedback gain matrix
    struct blasfeo_dmat K_mat;                           // shape = (nu, nx)

    struct blasfeo_dmat dct_dux;             // shape = (nlbu_t + nlbx_t + nlg_t + nlh_t + nubu_t + nubx_t + nug_t + nuh_t, nx + nu)
    // Jacobian of tightened constraints wrt [u,x].
    struct blasfeo_dmat scaled_dct_dux;      // shape = (nlbu_t + nlbx_t + nlg_t + nlh_t + nubu_t + nubx_t + nug_t + nuh_t, nx + nu)

    struct blasfeo_dmat dcet_dx;           // shape = (nlbx_e_t + nlg_e_t + nlh_e_t + nubx_e_t + nug_e_t + nuh_e_t,       nx)
    // Jacobian of tightened constraints wrt x at terminal node

    struct blasfeo_dmat scaled_dcet_dx;    // shape = (nlbx_e_t + nlg_e_t + nlh_e_t + nubx_e_t + nug_e_t + nuh_e_t,       nx)
    struct blasfeo_dvec *ineq_backoff_sq_buffer;         // shape = N * (nbu + nbx + ng + nh,) + (nbx_e + ng_e + nh_e, )
    struct blasfeo_dvec ricc_ones;         // shape = (nbu + nbx + ng + nh,) + (nbx_e + ng_e + nh_e, )

    // AK = A - B@K
    struct blasfeo_dmat AK_mat;                          // shape = (nx, nx)
    // A@P_k
    struct blasfeo_dmat temp_AP_mat;                     // shape = (nx, nx)
    // K@P_k, K@P_k@K^T
    struct blasfeo_dmat temp_KP_mat;                     // shape = (nu, nx)
    struct blasfeo_dmat temp_KPK_mat;                    // shape = (nu, nu)
    // C + D @ K, (C + D @ K) @ P_k
    struct blasfeo_dmat temp_CaDK_mat;                   // shape = (ngh_me_max, nx)
    struct blasfeo_dmat temp_CaDKmP_mat;                 // shape = (ngh_me_max, nx)
    struct blasfeo_dmat temp_beta_mat;                   // shape = (ngh_me_max, ngh_me_max)


    // matrices for Riccati recursion
    struct blasfeo_dmat *riccati_K_buffer;               // shape = N * (nu, nx)
    struct blasfeo_dmat riccati_Q_mat;                   // shape = (nx, nx)
    struct blasfeo_dmat riccati_Q_const;              // shape = (nx, nx)
    struct blasfeo_dmat riccati_Q_const_e;            // shape = (nx, nx)
    struct blasfeo_dmat riccati_R_mat;                   // shape = (nu, nu)
    struct blasfeo_dmat riccati_R_const;              // shape = (nu, nu)
    struct blasfeo_dmat riccati_S_mat;                   // shape = (nu, nx)
    struct blasfeo_dmat riccati_S_const;              // shape = (nu, nx)

    struct blasfeo_dmat temp_riccati_P_mat;              // shape = (nx, nx)
    struct blasfeo_dmat temp_riccati_P_plus_mat;         // shape = (nx, nx)
    // B.T @ P_k
    struct blasfeo_dmat temp_riccati_BP_mat;             // shape = (nu, nx)
    // R + B.T @ P_k @ B
    struct blasfeo_dmat temp_riccati_RaBPB_mat;          // shape = (nu, nu)
    // S + B.T @ P_k @ A
    struct blasfeo_dmat temp_riccati_SaBPA_mat;          // shape = (nu, nx)
    struct blasfeo_dmat temp_riccati_chol_mat;           // shape = (nu, nu)
    struct blasfeo_dmat temp_riccati_cholinvSaBPA_mat;   // shape = (nu, nx)
    struct blasfeo_dmat temp_riccati_AP_mat;             // shape = (nx, nx)

    double *d_A_mat;                                     // shape = (nx, nx)
    double *d_B_mat;                                     // shape = (nx, nu)
    double *d_Cg_mat;                                    // shape = (ng, nx)
    double *d_Dg_mat;                                    // shape = (ng, nu)
    double *d_Cg_e_mat;                                  // shape = (ng_e, nx)
    double *d_Cgh_mat;                                   // shape = (ng+nh, nx)
    double *d_Dgh_mat;                                   // shape = (ng+nh, nu)
    double *d_Cgh_e_mat;                                 // shape = (ng_e+nh_e, nx)
    double *d_state_vec;
    // upper and lower bounds on state variables
    double *d_lbx;                                       // shape = (nbx,)
    double *d_ubx;                                       // shape = (nbx,)
    double *d_lbx_e;                                     // shape = (nbx_e,)
    double *d_ubx_e;                                     // shape = (nbx_e,)
    // tightened upper and lower bounds on state variables
    double *d_lbx_tightened;                             // shape = (nbx,)
    double *d_ubx_tightened;                             // shape = (nbx,)
    double *d_lbx_e_tightened;                           // shape = (nbx_e,)
    double *d_ubx_e_tightened;                           // shape = (nbx_e,)
    // upper and lower bounds on control inputs
    double *d_lbu;                                       // shape = (nbu,)
    double *d_ubu;                                       // shape = (nbu,)
    // tightened upper and lower bounds on control inputs
    double *d_lbu_tightened;                             // shape = (nbu,)
    double *d_ubu_tightened;                             // shape = (nbu,)
    // upper and lower bounds on polytopic constraints
    double *d_lg;                                        // shape = (ng,)
    double *d_ug;                                        // shape = (ng,)
    double *d_lg_e;                                      // shape = (ng_e,)
    double *d_ug_e;                                      // shape = (ng_e,)
    // tightened lower bounds on polytopic constraints
    double *d_lg_tightened;                              // shape = (ng,)
    double *d_ug_tightened;                              // shape = (ng,)
    double *d_lg_e_tightened;                            // shape = (ng_e,)
    double *d_ug_e_tightened;                            // shape = (ng_e,)
    // upper and lower bounds on nonlinear constraints
    double *d_lh;                                        // shape = (nh,)
    double *d_uh;                                        // shape = (nh,)
    double *d_lh_e;                                      // shape = (nh_e,)
    double *d_uh_e;                                      // shape = (nh_e,)
    // tightened upper and lower bounds on nonlinear constraints
    double *d_lh_tightened;                              // shape = (nh,)
    double *d_uh_tightened;                              // shape = (nh,)
    double *d_lh_e_tightened;                            // shape = (nh_e,)
    double *d_uh_e_tightened;                            // shape = (nh_e,)
    // ineq constraint values, used for riccati
    double *d_ineq_val;                                  // shape = (nbu + nbx + ng + nh, ), actually only the max of (nbu, nbx, ng + nh) is needed
    double *d_ineq_e_val;                                // shape = (nbx_e + ng_e + nh_e, )

    int *idxbx;                                          // shape = (nbx,)
    int *idxbu;                                          // shape = (nbu,)
    int *idxbx_e;                                        // shape = (nbx_e,)

    int offset_W_diag;
    int offset_W_add_diag;
    int offset_P_out;

    double tau;

    void *raw_memory; // Pointer to allocated memory, to be used for freeing
} custom_memory;

static int int_max(int num1, int num2)
{
    return (num1 > num2 ) ? num1 : num2;
}


static int custom_memory_calculate_size(ocp_nlp_config *nlp_config, ocp_nlp_dims *nlp_dims)
{
    int N = nlp_dims->N;
    int nx = {{ dims.nx }};
    int nu = {{ dims.nu }};
    int nw = {{ zoro_description.nw }};

    int ng = {{ dims.ng }};
    int nh = {{ dims.nh }};
    int nbx = {{ dims.nbx }};
    int nbu = {{ dims.nbu }};

    int ng_e = {{ dims.ng_e }};
    int nh_e = {{ dims.nh_e }};
    int ngh_e_max = int_max(ng_e, nh_e);
    int ngh_me_max = int_max(ngh_e_max, int_max(ng, nh));
    int nbx_e = {{ dims.nbx_e }};

    assert({{zoro_description.nlbx_t}} <= nbx);
    assert({{zoro_description.nubx_t}} <= nbx);
    assert({{zoro_description.nlbu_t}} <= nbu);
    assert({{zoro_description.nubu_t}} <= nbu);
    assert({{zoro_description.nlg_t}} <= ng);
    assert({{zoro_description.nug_t}} <= ng);
    assert({{zoro_description.nlh_t}} <= nh);
    assert({{zoro_description.nuh_t}} <= nh);
    assert({{zoro_description.nlbx_e_t}} <= nbx_e);
    assert({{zoro_description.nubx_e_t}} <= nbx_e);
    assert({{zoro_description.nlg_e_t}} <= ng_e);
    assert({{zoro_description.nug_e_t}} <= ng_e);
    assert({{zoro_description.nlh_e_t}} <= nh_e);
    assert({{zoro_description.nuh_e_t}} <= nh_e);

    acados_size_t size = sizeof(custom_memory);
    size += nbx * sizeof(int);
    /* blasfeo structs */
    size += (2 * N + 1) * sizeof(struct blasfeo_dmat); // uncertainty_matrix_buffer, riccati_K_buffer
    /* blasfeo mem: mat */
    size += (N + 1) * blasfeo_memsize_dmat(nx, nx); // uncertainty_matrix_buffer
    size += blasfeo_memsize_dmat(nw, nw);           // W_mat
    size += 2 * blasfeo_memsize_dmat(nx, nw);       // unc_jac_G_mat, temp_GW_mat
    size += 4 * blasfeo_memsize_dmat(nx, nx);       // GWG_mat, A_mat, AK_mat, temp_AP_mat
    size += blasfeo_memsize_dmat(nx, nu);           // B_mat
    size += 2 * blasfeo_memsize_dmat(nu, nx); // K_mat, temp_KP_mat
    size += blasfeo_memsize_dmat(nu, nu);       // temp_KPK_mat
    size += blasfeo_memsize_dmat(ng, nx);           // Cg_mat
    size += blasfeo_memsize_dmat(ng, nu);           // Dg_mat
    size += blasfeo_memsize_dmat(ng_e, nx);         // Cg_e_mat
    size += blasfeo_memsize_dmat(ngh_e_max, nu);    // dummy_Dgh_e_mat
    size += blasfeo_memsize_dmat(nh, nx);           // Ch_mat
    size += blasfeo_memsize_dmat(nh, nu);           // Dh_mat
    size += blasfeo_memsize_dmat(nh_e, nx);         // Ch_e_mat
    size += 2 * blasfeo_memsize_dmat(ngh_me_max, nx);           // temp_CaDK_mat, temp_CaDKmP_mat
    size += blasfeo_memsize_dmat(ngh_me_max, ngh_me_max);       // temp_beta_mat
    // NOTE: Covariance matrix of the additive noise (used if input_W_add_diag)
    size += blasfeo_memsize_dmat(nw, nw);  // W_stage_mat

{%- if zoro_description.feedback_optimization_mode != "CONSTANT_FEEDBACK" %}
    // Riccati specific memory
    size += 6 * blasfeo_memsize_dmat(nx, nx);       // riccati_Q_mat, riccati_Q_const, riccati_Q_const_e, temp_riccati_P_mat, temp_riccati_P_plus_mat, temp_riccati_AP_mat
    size += (N + 5) * blasfeo_memsize_dmat(nu, nx); // riccati_K_buffer, riccati_S_mat, riccati_S_const, temp_riccati_BP_mat, temp_riccati_SaBPA_mat, temp_riccati_cholinvSaBPA_mat
    size += 4 * blasfeo_memsize_dmat(nu, nu);       // riccati_R_mat, riccati_R_const, temp_riccati_RaBPB_mat, temp_riccati_chol_mat
    size += 2 * blasfeo_memsize_dmat(NCT, nx + nu);  // dct_dux, scaled_dct_dux
    size += 2 * blasfeo_memsize_dmat(NCT_E, nx);  // dcet_dx, scaled_dcet_dx
{% endif %}

    size += (N+1) * sizeof(struct blasfeo_dvec);                // ineq_backoff_sq_buffer
    size += N * blasfeo_memsize_dvec(nbu + nbx + ng + nh);      // ineq_backoff_sq_buffer--stage
    size += blasfeo_memsize_dvec(nbx_e + ng_e + nh_e);          // ineq_backoff_sq_buffer--terminal
    size += blasfeo_memsize_dvec(nbu + nbx + ng + nh + nbx_e + ng_e + nh_e);          // ricc_ones

    /* blasfeo mem: vec */
    /* Arrays */
    size += nx*nx *sizeof(double);                  // d_A_mat
    size += nx*nu *sizeof(double);                  // d_B_mat
    size += (ng + ng_e) * nx * sizeof(double);      // d_Cg_mat, d_Cg_e_mat
    size += (ng) * nu * sizeof(double);             // d_Dg_mat
    size += (nh + nh_e + ng + ng_e) * nx * sizeof(double);      // d_Cgh_mat, d_Cgh_e_mat
    size += (nh + ng) * nu * sizeof(double);        // d_Dgh_mat
    // d_state_vec
    size += nx *sizeof(double);
    // constraints and tightened constraints
    size += 5 * (nbx + nbu + ng + nh)*sizeof(double);
    size += 5 * (nbx_e + ng_e + nh_e)*sizeof(double);
    size += (nbx + nbu + nbx_e)*sizeof(int);        // idxbx, idxbu, idxbx_e

    size += 1 * 8; // initial alignment
    make_int_multiple_of(64, &size);
    size += 1 * 64;

    return size;
}


static custom_memory *custom_memory_assign(ocp_nlp_config *nlp_config, ocp_nlp_dims *nlp_dims, void *raw_memory)
{
    int N = nlp_dims->N;
    int nx = {{ dims.nx }};
    int nu = {{ dims.nu }};
    int nw = {{ zoro_description.nw }};

    int ng = {{ dims.ng }};
    int nh = {{ dims.nh }};
    int nbx = {{ dims.nbx }};
    int nbu = {{ dims.nbu }};

    int ng_e = {{ dims.ng_e }};
    int nh_e = {{ dims.nh_e }};
    int ngh_e_max = int_max(ng_e, nh_e);
    int ngh_me_max = int_max(ngh_e_max, int_max(ng, nh));
    int nbx_e = {{ dims.nbx_e }};

    char *c_ptr = (char *) raw_memory;
    custom_memory *mem = (custom_memory *) c_ptr;
    c_ptr += sizeof(custom_memory);

    align_char_to(8, &c_ptr);
    assign_and_advance_blasfeo_dmat_structs(N+1, &mem->uncertainty_matrix_buffer, &c_ptr);
    assign_and_advance_blasfeo_dmat_structs(N, &mem->riccati_K_buffer, &c_ptr);
    assign_and_advance_blasfeo_dvec_structs(N+1, &mem->ineq_backoff_sq_buffer, &c_ptr);

    align_char_to(64, &c_ptr);

    for (int ii = 0; ii <= N; ii++)
    {
        assign_and_advance_blasfeo_dmat_mem(nx, nx, &mem->uncertainty_matrix_buffer[ii], &c_ptr);
    }
    // Disturbance Dynamics
    assign_and_advance_blasfeo_dmat_mem(nw, nw, &mem->W_mat, &c_ptr);
    assign_and_advance_blasfeo_dmat_mem(nx, nw, &mem->unc_jac_G_mat, &c_ptr);
    assign_and_advance_blasfeo_dmat_mem(nx, nw, &mem->temp_GW_mat, &c_ptr);
    assign_and_advance_blasfeo_dmat_mem(nx, nx, &mem->GWG_mat, &c_ptr);
    // NOTE: Covariance matrix of the additive noise disturbance (used if input_W_add_diag)
    assign_and_advance_blasfeo_dmat_mem(nw, nw, &mem->W_stage_mat, &c_ptr);
    // System Dynamics
    assign_and_advance_blasfeo_dmat_mem(nx, nx, &mem->A_mat, &c_ptr);
    assign_and_advance_blasfeo_dmat_mem(nx, nu, &mem->B_mat, &c_ptr);
    assign_and_advance_blasfeo_dmat_mem(ng, nx, &mem->Cg_mat, &c_ptr);
    assign_and_advance_blasfeo_dmat_mem(ng, nu, &mem->Dg_mat, &c_ptr);
    assign_and_advance_blasfeo_dmat_mem(ng_e, nx, &mem->Cg_e_mat, &c_ptr);
    assign_and_advance_blasfeo_dmat_mem(ngh_e_max, nu, &mem->dummy_Dgh_e_mat, &c_ptr);
    assign_and_advance_blasfeo_dmat_mem(nh, nx, &mem->Ch_mat, &c_ptr);
    assign_and_advance_blasfeo_dmat_mem(nh, nu, &mem->Dh_mat, &c_ptr);
    assign_and_advance_blasfeo_dmat_mem(nh_e, nx, &mem->Ch_e_mat, &c_ptr);
    assign_and_advance_blasfeo_dmat_mem(nu, nx, &mem->K_mat, &c_ptr);

    assign_and_advance_blasfeo_dmat_mem(nx, nx, &mem->AK_mat, &c_ptr);
    assign_and_advance_blasfeo_dmat_mem(nx, nx, &mem->temp_AP_mat, &c_ptr);
    assign_and_advance_blasfeo_dmat_mem(nu, nx, &mem->temp_KP_mat, &c_ptr);
    assign_and_advance_blasfeo_dmat_mem(nu, nu, &mem->temp_KPK_mat, &c_ptr);
    assign_and_advance_blasfeo_dmat_mem(ngh_me_max, nx, &mem->temp_CaDK_mat, &c_ptr);
    assign_and_advance_blasfeo_dmat_mem(ngh_me_max, nx, &mem->temp_CaDKmP_mat, &c_ptr);
    assign_and_advance_blasfeo_dmat_mem(ngh_me_max, ngh_me_max, &mem->temp_beta_mat, &c_ptr);

{%- if zoro_description.feedback_optimization_mode != "CONSTANT_FEEDBACK" %}
    for (int ii = 0; ii < N; ii++)
    {
        assign_and_advance_blasfeo_dmat_mem(nu, nx, &mem->riccati_K_buffer[ii], &c_ptr);
    }
    assign_and_advance_blasfeo_dmat_mem(nx, nx, &mem->riccati_Q_mat, &c_ptr);
    assign_and_advance_blasfeo_dmat_mem(nx, nx, &mem->riccati_Q_const, &c_ptr);
    assign_and_advance_blasfeo_dmat_mem(nx, nx, &mem->riccati_Q_const_e, &c_ptr);
    assign_and_advance_blasfeo_dmat_mem(nu, nu, &mem->riccati_R_mat, &c_ptr);
    assign_and_advance_blasfeo_dmat_mem(nu, nu, &mem->riccati_R_const, &c_ptr);
    assign_and_advance_blasfeo_dmat_mem(nu, nx, &mem->riccati_S_mat, &c_ptr);
    assign_and_advance_blasfeo_dmat_mem(nu, nx, &mem->riccati_S_const, &c_ptr);

    assign_and_advance_blasfeo_dmat_mem(nx, nx, &mem->temp_riccati_P_mat, &c_ptr);
    assign_and_advance_blasfeo_dmat_mem(nx, nx, &mem->temp_riccati_P_plus_mat, &c_ptr);
    assign_and_advance_blasfeo_dmat_mem(nx, nx, &mem->temp_riccati_AP_mat, &c_ptr);
    assign_and_advance_blasfeo_dmat_mem(nu, nx, &mem->temp_riccati_BP_mat, &c_ptr);
    assign_and_advance_blasfeo_dmat_mem(nu, nu, &mem->temp_riccati_RaBPB_mat, &c_ptr);
    assign_and_advance_blasfeo_dmat_mem(nu, nx, &mem->temp_riccati_SaBPA_mat, &c_ptr);
    assign_and_advance_blasfeo_dmat_mem(nu, nu, &mem->temp_riccati_chol_mat, &c_ptr);
    assign_and_advance_blasfeo_dmat_mem(nu, nx, &mem->temp_riccati_cholinvSaBPA_mat, &c_ptr);
    assign_and_advance_blasfeo_dmat_mem(NCT, nx + nu, &mem->dct_dux, &c_ptr);
    assign_and_advance_blasfeo_dmat_mem(NCT, nx + nu, &mem->scaled_dct_dux, &c_ptr);
    assign_and_advance_blasfeo_dmat_mem(NCT_E, nx, &mem->dcet_dx, &c_ptr);
    assign_and_advance_blasfeo_dmat_mem(NCT_E, nx, &mem->scaled_dcet_dx, &c_ptr);
{% endif %}

    for (int ii = 0; ii < N; ii++)
    {
        assign_and_advance_blasfeo_dvec_mem(nbu + nbx + ng + nh, &mem->ineq_backoff_sq_buffer[ii], &c_ptr);
    }
    assign_and_advance_blasfeo_dvec_mem(nbx_e + ng_e + nh_e, &mem->ineq_backoff_sq_buffer[N], &c_ptr);

    assign_and_advance_blasfeo_dvec_mem(nbu + nbx + ng + nh + nbx_e + ng_e + nh_e, &mem->ricc_ones, &c_ptr);
    blasfeo_dvecse(nbu + nbx + ng + nh + nbx_e + ng_e + nh_e, 1.0, &mem->ricc_ones, 0);

    assign_and_advance_double(nx*nx, &mem->d_A_mat, &c_ptr);
    assign_and_advance_double(nx*nu, &mem->d_B_mat, &c_ptr);
    assign_and_advance_double(ng*nx, &mem->d_Cg_mat, &c_ptr);
    assign_and_advance_double(ng*nu, &mem->d_Dg_mat, &c_ptr);
    assign_and_advance_double(ng_e*nx, &mem->d_Cg_e_mat, &c_ptr);
    assign_and_advance_double((ng + nh)*nx, &mem->d_Cgh_mat, &c_ptr);
    assign_and_advance_double((ng + nh)*nu, &mem->d_Dgh_mat, &c_ptr);
    assign_and_advance_double((ng_e + nh_e)*nx, &mem->d_Cgh_e_mat, &c_ptr);
    assign_and_advance_double(nx, &mem->d_state_vec, &c_ptr);
    assign_and_advance_double(nbx, &mem->d_lbx, &c_ptr);
    assign_and_advance_double(nbx, &mem->d_ubx, &c_ptr);
    assign_and_advance_double(nbx_e, &mem->d_lbx_e, &c_ptr);
    assign_and_advance_double(nbx_e, &mem->d_ubx_e, &c_ptr);
    assign_and_advance_double(nbx, &mem->d_lbx_tightened, &c_ptr);
    assign_and_advance_double(nbx, &mem->d_ubx_tightened, &c_ptr);
    assign_and_advance_double(nbx_e, &mem->d_lbx_e_tightened, &c_ptr);
    assign_and_advance_double(nbx_e, &mem->d_ubx_e_tightened, &c_ptr);
    assign_and_advance_double(nbu, &mem->d_lbu, &c_ptr);
    assign_and_advance_double(nbu, &mem->d_ubu, &c_ptr);
    assign_and_advance_double(nbu, &mem->d_lbu_tightened, &c_ptr);
    assign_and_advance_double(nbu, &mem->d_ubu_tightened, &c_ptr);
    assign_and_advance_double(ng, &mem->d_lg, &c_ptr);
    assign_and_advance_double(ng, &mem->d_ug, &c_ptr);
    assign_and_advance_double(ng_e, &mem->d_lg_e, &c_ptr);
    assign_and_advance_double(ng_e, &mem->d_ug_e, &c_ptr);
    assign_and_advance_double(ng, &mem->d_lg_tightened, &c_ptr);
    assign_and_advance_double(ng, &mem->d_ug_tightened, &c_ptr);
    assign_and_advance_double(ng_e, &mem->d_lg_e_tightened, &c_ptr);
    assign_and_advance_double(ng_e, &mem->d_ug_e_tightened, &c_ptr);
    assign_and_advance_double(nh, &mem->d_lh, &c_ptr);
    assign_and_advance_double(nh, &mem->d_uh, &c_ptr);
    assign_and_advance_double(nh_e, &mem->d_lh_e, &c_ptr);
    assign_and_advance_double(nh_e, &mem->d_uh_e, &c_ptr);
    assign_and_advance_double(nh, &mem->d_lh_tightened, &c_ptr);
    assign_and_advance_double(nh, &mem->d_uh_tightened, &c_ptr);
    assign_and_advance_double(nh_e, &mem->d_lh_e_tightened, &c_ptr);
    assign_and_advance_double(nh_e, &mem->d_uh_e_tightened, &c_ptr);
    assign_and_advance_double(nbu + nbx + ng + nh, &mem->d_ineq_val, &c_ptr);
    assign_and_advance_double(nbx_e + ng_e + nh_e, &mem->d_ineq_e_val, &c_ptr);

    assign_and_advance_int(nbx, &mem->idxbx, &c_ptr);
    assign_and_advance_int(nbu, &mem->idxbu, &c_ptr);
    assign_and_advance_int(nbx_e, &mem->idxbx_e, &c_ptr);

    // compute input and output data offsets:
    mem->offset_W_diag = 0;
{%- if zoro_description.input_P0_diag %}
    mem->offset_W_diag += nx;  // P0_diag
{%- elif zoro_description.input_P0 %}
    mem->offset_W_diag += nx*nx;  // P0_diag
{% endif %}

    mem->offset_W_add_diag = mem->offset_W_diag;
{%- if zoro_description.input_W_diag %}
    mem->offset_W_add_diag += nw;  // W_diag
{% endif %}

    mem->offset_P_out = mem->offset_W_add_diag;
{%- if zoro_description.input_W_add_diag %}
    mem->offset_P_out += N * nw;
{% endif %}

    mem->tau = 1.0;

    assert((char *) raw_memory + custom_memory_calculate_size(nlp_config, nlp_dims) >= c_ptr);
    mem->raw_memory = raw_memory;

    return mem;
}



static void *custom_memory_create({{ model.name }}_solver_capsule* capsule)
{
    // printf("\nin custom_memory_create_function\n");

    ocp_nlp_dims *nlp_dims = {{ model.name }}_acados_get_nlp_dims(capsule);
    ocp_nlp_config *nlp_config = {{ model.name }}_acados_get_nlp_config(capsule);
    acados_size_t bytes = custom_memory_calculate_size(nlp_config, nlp_dims);

    void *ptr = acados_calloc(1, bytes);

    custom_memory *custom_mem = custom_memory_assign(nlp_config, nlp_dims, ptr);
    custom_mem->raw_memory = ptr;

    return custom_mem;
}


static void custom_val_init_function(ocp_nlp_dims *nlp_dims, ocp_nlp_in *nlp_in, ocp_nlp_solver *nlp_solver, custom_memory *custom_mem)
{
    int N = nlp_dims->N;
    int nx = {{ dims.nx }};
    int nu = {{ dims.nu }};
    int nw = {{ zoro_description.nw }};

    int ng = {{ dims.ng }};
    int nh = {{ dims.nh }};
    int nbx = {{ dims.nbx }};
    int nbu = {{ dims.nbu }};

    int ng_e = {{ dims.ng_e }};
    int nh_e = {{ dims.nh_e }};
    int ngh_e_max = int_max(ng_e, nh_e);
    int nbx_e = {{ dims.nbx_e }};

    /* Get the state constraint bounds */
    ocp_nlp_constraints_model_get(nlp_solver->config, nlp_dims, nlp_in, 1, "idxbx", custom_mem->idxbx);
    ocp_nlp_constraints_model_get(nlp_solver->config, nlp_dims, nlp_in, N, "idxbx", custom_mem->idxbx_e);
    ocp_nlp_constraints_model_get(nlp_solver->config, nlp_dims, nlp_in, 1, "lbx", custom_mem->d_lbx);
    ocp_nlp_constraints_model_get(nlp_solver->config, nlp_dims, nlp_in, 1, "ubx", custom_mem->d_ubx);
    ocp_nlp_constraints_model_get(nlp_solver->config, nlp_dims, nlp_in, N, "lbx", custom_mem->d_lbx_e);
    ocp_nlp_constraints_model_get(nlp_solver->config, nlp_dims, nlp_in, N, "ubx", custom_mem->d_ubx_e);
    ocp_nlp_constraints_model_get(nlp_solver->config, nlp_dims, nlp_in, 1, "idxbu", custom_mem->idxbu);
    ocp_nlp_constraints_model_get(nlp_solver->config, nlp_dims, nlp_in, 1, "lbu", custom_mem->d_lbu);
    ocp_nlp_constraints_model_get(nlp_solver->config, nlp_dims, nlp_in, 1, "ubu", custom_mem->d_ubu);
    // Get the Jacobians and the bounds of the linear constraints
    ocp_nlp_constraints_model_get(nlp_solver->config, nlp_dims, nlp_in, 1, "lg", custom_mem->d_lg);
    ocp_nlp_constraints_model_get(nlp_solver->config, nlp_dims, nlp_in, 1, "ug", custom_mem->d_ug);
    ocp_nlp_constraints_model_get(nlp_solver->config, nlp_dims, nlp_in, N, "lg", custom_mem->d_lg_e);
    ocp_nlp_constraints_model_get(nlp_solver->config, nlp_dims, nlp_in, N, "ug", custom_mem->d_ug_e);
    ocp_nlp_constraints_model_get(nlp_solver->config, nlp_dims, nlp_in, 1, "C", custom_mem->d_Cg_mat);
    ocp_nlp_constraints_model_get(nlp_solver->config, nlp_dims, nlp_in, 1, "D", custom_mem->d_Dg_mat);
    ocp_nlp_constraints_model_get(nlp_solver->config, nlp_dims, nlp_in, N, "C", custom_mem->d_Cg_e_mat);
    blasfeo_pack_dmat(ng, nx, custom_mem->d_Cg_mat, ng, &custom_mem->Cg_mat, 0, 0);
    blasfeo_pack_dmat(ng, nu, custom_mem->d_Dg_mat, ng, &custom_mem->Dg_mat, 0, 0);
    blasfeo_pack_dmat(ng_e, nx, custom_mem->d_Cg_e_mat, ng_e, &custom_mem->Cg_e_mat, 0, 0);
    blasfeo_dgese(ngh_e_max, nu, 0., &custom_mem->dummy_Dgh_e_mat, 0, 0); //fill with zeros
    // NOTE: fixed lower and upper bounds of nonlinear constraints
    ocp_nlp_constraints_model_get(nlp_solver->config, nlp_dims, nlp_in, 1, "lh", custom_mem->d_lh);
    ocp_nlp_constraints_model_get(nlp_solver->config, nlp_dims, nlp_in, 1, "uh", custom_mem->d_uh);
    ocp_nlp_constraints_model_get(nlp_solver->config, nlp_dims, nlp_in, N, "lh", custom_mem->d_lh_e);
    ocp_nlp_constraints_model_get(nlp_solver->config, nlp_dims, nlp_in, N, "uh", custom_mem->d_uh_e);

    /* Initilize tightened constraints*/
    // NOTE: tightened constraints are only initialized once
    ocp_nlp_constraints_model_get(nlp_solver->config, nlp_dims, nlp_in, 1, "lbx", custom_mem->d_lbx_tightened);
    ocp_nlp_constraints_model_get(nlp_solver->config, nlp_dims, nlp_in, 1, "ubx", custom_mem->d_ubx_tightened);
    ocp_nlp_constraints_model_get(nlp_solver->config, nlp_dims, nlp_in, N, "lbx", custom_mem->d_lbx_e_tightened);
    ocp_nlp_constraints_model_get(nlp_solver->config, nlp_dims, nlp_in, N, "ubx", custom_mem->d_ubx_e_tightened);
    ocp_nlp_constraints_model_get(nlp_solver->config, nlp_dims, nlp_in, 1, "lbu", custom_mem->d_lbu_tightened);
    ocp_nlp_constraints_model_get(nlp_solver->config, nlp_dims, nlp_in, 1, "ubu", custom_mem->d_ubu_tightened);
    ocp_nlp_constraints_model_get(nlp_solver->config, nlp_dims, nlp_in, 1, "lg", custom_mem->d_lg_tightened);
    ocp_nlp_constraints_model_get(nlp_solver->config, nlp_dims, nlp_in, 1, "ug", custom_mem->d_ug_tightened);
    ocp_nlp_constraints_model_get(nlp_solver->config, nlp_dims, nlp_in, N, "lg", custom_mem->d_lg_e_tightened);
    ocp_nlp_constraints_model_get(nlp_solver->config, nlp_dims, nlp_in, N, "ug", custom_mem->d_ug_e_tightened);
    ocp_nlp_constraints_model_get(nlp_solver->config, nlp_dims, nlp_in, 1, "lh", custom_mem->d_lh_tightened);
    ocp_nlp_constraints_model_get(nlp_solver->config, nlp_dims, nlp_in, 1, "uh", custom_mem->d_uh_tightened);
    ocp_nlp_constraints_model_get(nlp_solver->config, nlp_dims, nlp_in, N, "lh", custom_mem->d_lh_e_tightened);
    ocp_nlp_constraints_model_get(nlp_solver->config, nlp_dims, nlp_in, N, "uh", custom_mem->d_uh_e_tightened);

    /* Initialize the W matrix */
    // blasfeo_dgese(nw, nw, 0., &custom_mem->W_mat, 0, 0);
{%- for ir in range(end=zoro_description.nw) %}
    {%- for ic in range(end=zoro_description.nw) %}
    blasfeo_dgein1({{zoro_description.W_mat[ir][ic]}}, &custom_mem->W_mat, {{ir}}, {{ic}});
    {%- endfor %}
{%- endfor %}

{%- for ir in range(end=dims.nx) %}
    {%- for ic in range(end=zoro_description.nw) %}
    blasfeo_dgein1({{zoro_description.unc_jac_G_mat[ir][ic]}}, &custom_mem->unc_jac_G_mat, {{ir}}, {{ic}});
    {%- endfor %}
{%- endfor %}

{%- if not zoro_description.input_W_diag %}
    // NOTE: G, W are not changing -> precompute GWG
    // temp_GW_mat = unc_jac_G_mat * W_mat
    blasfeo_dgemm_nn(nx, nw, nw, 1.0, &custom_mem->unc_jac_G_mat, 0, 0,
                        &custom_mem->W_mat, 0, 0, 0.0,
                        &custom_mem->temp_GW_mat, 0, 0, &custom_mem->temp_GW_mat, 0, 0);
    // GWG_mat = temp_GW_mat * unc_jac_G_mat^T
    blasfeo_dgemm_nt(nx, nx, nw, 1.0, &custom_mem->temp_GW_mat, 0, 0,
                        &custom_mem->unc_jac_G_mat, 0, 0, 0.0,
                        &custom_mem->GWG_mat, 0, 0, &custom_mem->GWG_mat, 0, 0);
{%- endif %}

{%- if zoro_description.input_W_add_diag %}
    // NOTE: Initialize additive covariance matrix to 0.
    blasfeo_dgese(nw, nw, 0.0, &custom_mem->W_stage_mat, 0, 0);
{%- endif %}

    /* Initialize the uncertainty_matrix_buffer[0] */
{%- for ir in range(end=dims.nx) %}
    {%- for ic in range(end=dims.nx) %}
    blasfeo_dgein1({{zoro_description.P0_mat[ir][ic]}}, &custom_mem->uncertainty_matrix_buffer[0], {{ir}}, {{ic}});
    {%- endfor %}
{%- endfor %}

    /* Initialize the feedback gain matrix */
{%- for ir in range(end=dims.nu) %}
    {%- for ic in range(end=dims.nx) %}
    blasfeo_dgein1({{zoro_description.fdbk_K_mat[ir][ic]}}, &custom_mem->K_mat, {{ir}}, {{ic}});
    {%- endfor %}
{%- endfor %}

{%- if not zoro_description.feedback_optimization_mode == "CONSTANT_FEEDBACK" %}
{%- for ir in range(end=dims.nx) %}
    {%- for ic in range(end=dims.nx) %}
    blasfeo_dgein1({{zoro_description.riccati_Q_const[ir][ic]}}, &custom_mem->riccati_Q_const, {{ir}}, {{ic}});
    {%- endfor %}
{%- endfor %}
{%- for ir in range(end=dims.nx) %}
    {%- for ic in range(end=dims.nx) %}
    blasfeo_dgein1({{zoro_description.riccati_Q_const_e[ir][ic]}}, &custom_mem->riccati_Q_const_e, {{ir}}, {{ic}});
    {%- endfor %}
{%- endfor %}
{%- for ir in range(end=dims.nu) %}
    {%- for ic in range(end=dims.nu) %}
    blasfeo_dgein1({{zoro_description.riccati_R_const[ir][ic]}}, &custom_mem->riccati_R_const, {{ir}}, {{ic}});
    {%- endfor %}
{%- endfor %}
{%- for ir in range(end=dims.nu) %}
    {%- for ic in range(end=dims.nx) %}
    blasfeo_dgein1({{zoro_description.riccati_S_const[ir][ic]}}, &custom_mem->riccati_S_const, {{ir}}, {{ic}});
    {%- endfor %}
{%- endfor %}
    // initialize riccati_Q_mat, riccati_R_mat, riccati_S_mat <- they dont change for RICCATI_CONSTANT_COST
    blasfeo_dgecp(nx, nx, &custom_mem->riccati_Q_const, 0, 0, &custom_mem->riccati_Q_mat, 0, 0);
    blasfeo_dgecp(nu, nu, &custom_mem->riccati_R_const, 0, 0, &custom_mem->riccati_R_mat, 0, 0);
    blasfeo_dgecp(nu, nx, &custom_mem->riccati_S_const, 0, 0, &custom_mem->riccati_S_mat, 0, 0);

    // Set constant values of dct_dux (nlbu_t + nlbx_t + nlg_t + nlh_t + nubu_t + nubx_t + nug_t + nuh_t, nu + nx)
    // the gradients of u is on left of x
    blasfeo_dgese({{zoro_description.nlbu_t}} + {{zoro_description.nlbx_t}} + {{zoro_description.nlg_t}} + {{zoro_description.nlh_t}}
            + {{zoro_description.nubu_t}} + {{zoro_description.nubx_t}} + {{zoro_description.nug_t}} + {{zoro_description.nuh_t}}, nx + nu, 0.0,
            &custom_mem->dct_dux, 0, 0);
    blasfeo_dgese({{zoro_description.nlbx_e_t}} + {{zoro_description.nlg_e_t}} + {{zoro_description.nlh_e_t}}
            + {{zoro_description.nubx_e_t}} + {{zoro_description.nug_e_t}} + {{zoro_description.nuh_e_t}}, nx, 0.0,
            &custom_mem->dcet_dx, 0, 0);
    int ir = 0;
    // lbu
    {%- for it in zoro_description.idx_lbu_t %}
    blasfeo_dgein1(1.0, &custom_mem->dct_dux, ir, custom_mem->idxbu[{{it}}]);
    ir += 1;
    {%- endfor %}
    // lbx
    {%- for it in zoro_description.idx_lbx_t %}
    blasfeo_dgein1(1.0, &custom_mem->dct_dux, ir, nu + custom_mem->idxbx[{{it}}]);
    ir += 1;
    {%- endfor %}
    // lg
    {%- for it in zoro_description.idx_lg_t %}
    blasfeo_dgecp(1, nu, &custom_mem->Dg_mat, it, 0, &custom_mem->dct_dux, ir, 0);
    blasfeo_dgecp(1, nx, &custom_mem->Cg_mat, it, 0, &custom_mem->dct_dux, ir, nu);
    ir += 1;
    {%- endfor %}
    // lh
    ir += {{zoro_description.nlh_t}};
    // ubu
    {%- for it in zoro_description.idx_ubu_t %}
    blasfeo_dgein1(1.0, &custom_mem->dct_dux, ir, custom_mem->idxbu[{{it}}]);
    ir += 1;
    {%- endfor %}
    // ubx
    {%- for it in zoro_description.idx_ubx_t %}
    blasfeo_dgein1(1.0, &custom_mem->dct_dux, ir, nu + custom_mem->idxbx[{{it}}]);
    ir += 1;
    {%- endfor %}
    // ug
    {%- for it in zoro_description.idx_ug_t %}
    blasfeo_dgecp(1, nu, &custom_mem->Dg_mat, it, 0, &custom_mem->dct_dux, ir, 0);
    blasfeo_dgecp(1, nx, &custom_mem->Cg_mat, it, 0, &custom_mem->dct_dux, ir, nu);
    ir += 1;
    {%- endfor %}

    ir = 0;
    // lbx_e
    {%- for it in zoro_description.idx_lbx_e_t %}
    blasfeo_dgein1(1.0, &custom_mem->dcet_dx, ir, custom_mem->idxbx[{{it}}]);
    ir += 1;
    {%- endfor %}
    // lg_e
    {%- for it in zoro_description.idx_lg_e_t %}
    blasfeo_dgecp(1, nx, &custom_mem->Cg_e_mat, it, 0, &custom_mem->dcet_dx, ir, 0);
    ir += 1;
    {%- endfor %}
    // lh_e (just add offset)
    ir += {{zoro_description.nlh_e_t}};
    // ubx_e
    {%- for it in zoro_description.idx_ubx_e_t %}
    blasfeo_dgein1(1.0, &custom_mem->dcet_dx, ir, custom_mem->idxbx[{{it}}]);
    ir += 1;
    {%- endfor %}
    // ug_e
    {%- for it in zoro_description.idx_ug_e_t %}
    blasfeo_dgecp(1, nx, &custom_mem->Cg_e_mat, it, 0, &custom_mem->dcet_dx, ir, 0);
    ir += 1;
    {%- endfor %}

    // fill with some values, but just for very first iteration, when back-offs are not available yet.
    for (int ii = 0; ii < N; ii++)
    {
        blasfeo_dvecse(nbu + nbx + ng + nh, 1e-4, &custom_mem->ineq_backoff_sq_buffer[ii], 0);
    }
    blasfeo_dvecse(nbx_e + ng_e + nh_e, 1e-4, &custom_mem->ineq_backoff_sq_buffer[N], 0);
    {%- endif %}
}


int custom_update_init_function({{ model.name }}_solver_capsule* capsule)
{
    capsule->custom_update_memory = custom_memory_create(capsule);
    ocp_nlp_in *nlp_in = {{ model.name }}_acados_get_nlp_in(capsule);

    ocp_nlp_dims *nlp_dims = {{ model.name }}_acados_get_nlp_dims(capsule);
    ocp_nlp_solver *nlp_solver = {{ model.name }}_acados_get_nlp_solver(capsule);
    custom_val_init_function(nlp_dims, nlp_in, nlp_solver, capsule->custom_update_memory);
    return 1;
}

static void compute_gh_beta(struct blasfeo_dmat* K_mat, struct blasfeo_dmat* C_mat,
                         struct blasfeo_dmat* D_mat, struct blasfeo_dmat* CaDK_mat,
                         struct blasfeo_dmat* CaDKmP_mat, struct blasfeo_dmat* beta_mat,
                         struct blasfeo_dmat* P_mat,
                         int n_cstr, int nx, int nu)
{
    // (C+DK)@P@(C^T+K^TD^T)
    // CaDK_mat = C_mat + D_mat @ K_mat
    blasfeo_dgemm_nn(n_cstr, nx, nu, 1.0, D_mat, 0, 0,
                        K_mat, 0, 0, 1.0,
                        C_mat, 0, 0, CaDK_mat, 0, 0);
    // CaDKmP_mat = CaDK_mat @ P_mat
    blasfeo_dgemm_nn(n_cstr, nx, nx, 1.0, CaDK_mat, 0, 0,
                        P_mat, 0, 0, 0.0,
                        CaDKmP_mat, 0, 0, CaDKmP_mat, 0, 0);
    // NOTE: here we also compute cross-terms which are not needed.
    //       Only diag(beta_mat) is used later
    // beta_mat = CaDKmP_mat @ CaDK_mat^T
    // blasfeo_dgemm_nt(n_cstr, n_cstr, nx, 1.0, CaDKmP_mat, 0, 0,
    //                     CaDK_mat, 0, 0, 0.0,
    //                     beta_mat, 0, 0, beta_mat, 0, 0);
    for (int ii = 0; ii < n_cstr; ii++)
    {
        blasfeo_dgemm_nt(1, 1, nx, 1.0, CaDKmP_mat, ii, 0,
                CaDK_mat, ii, 0, 0.0,
                beta_mat, ii, ii, beta_mat, ii, ii);
    }
}

static void compute_KPK(struct blasfeo_dmat* K_mat, struct blasfeo_dmat* temp_KP_mat,
                        struct blasfeo_dmat* temp_KPK_mat, struct blasfeo_dmat* P_mat,
                        int nx, int nu)
{
    // K @ P_k @ K^T
    // temp_KP_mat = K_mat @ P_mat
    blasfeo_dgemm_nn(nu, nx, nx, 1.0, K_mat, 0, 0,
                     P_mat, 0, 0, 0.0,
                     temp_KP_mat, 0, 0, temp_KP_mat, 0, 0);
    // temp_KPK_mat = temp_KP_mat @ K_mat^T
    blasfeo_dgemm_nt(nu, nu, nx, 1.0, temp_KP_mat, 0, 0,
                     K_mat, 0, 0, 0.0,
                     temp_KPK_mat, 0, 0, temp_KPK_mat, 0, 0);
}

static void compute_next_P_matrix(struct blasfeo_dmat* P_mat, struct blasfeo_dmat* P_next_mat,
                                  struct blasfeo_dmat* A_mat, struct blasfeo_dmat* B_mat,
                                  struct blasfeo_dmat* K_mat, struct blasfeo_dmat* W_mat,
                                  struct blasfeo_dmat* AK_mat, struct blasfeo_dmat* temp_AP_mat, int nx, int nu)
{
    // TODO: exploit symmetry of P, however, only blasfeo_dtrmm_rlnn is implemented in high-performance BLAFEO variant.
    // AK_mat = -B@K + A
    blasfeo_dgemm_nn(nx, nx, nu, -1.0, B_mat, 0, 0, K_mat, 0, 0,
                        1.0, A_mat, 0, 0, AK_mat, 0, 0);
    // temp_AP_mat = AK_mat @ P_k
    blasfeo_dgemm_nn(nx, nx, nx, 1.0, AK_mat, 0, 0,
                        P_mat, 0, 0, 0.0,
                        temp_AP_mat, 0, 0, temp_AP_mat, 0, 0);
    // P_{k+1} = temp_AP_mat @ AK_mat^T + GWG_mat
    blasfeo_dgemm_nt(nx, nx, nx, 1.0, temp_AP_mat, 0, 0,
                        AK_mat, 0, 0, 1.0,
                        W_mat, 0, 0, P_next_mat, 0, 0);
}

/**
 * @brief Sets the initial uncertainty P_0.
 *
 * We initialize the initial uncertainty matrix P_0 using a diagonal covariance matrix where
 * the elements are read from the first nx elements in the data buffer.
*/
static void reset_P0_matrix(ocp_nlp_dims *nlp_dims, struct blasfeo_dmat* P_mat, double* data)
{
    const int nx = nlp_dims->nx[0];
{%- if zoro_description.input_P0_diag %}
    for (int i = 0; i < nx; ++i)
    {
        if (data[i] < 0)
        {
            printf("skipping P0_diag[%d] = %f\n", i, data[i]);
            continue;
        }

        blasfeo_dgein1(data[i], P_mat, i, i);
    }
{%- elif zoro_description.input_P0 %}
    blasfeo_pack_dmat(nx, nx, data, nx, P_mat, 0, 0);
{% endif %}
}


{%- if zoro_description.input_W_diag %}
/**
 * @brief Update the process noise
 *
 * We update the process noise matrix W using the elements from the data buffer.
*/
static void reset_process_noise_matrix(custom_memory* custom_mem, ocp_nlp_dims *nlp_dims, struct blasfeo_dmat* W_mat, double* data)
{
    const int nx = nlp_dims->nx[0];
    const int nw = {{ zoro_description.nw }};

    for (int i = 0; i < nw; ++i)
    {
        if (data[custom_mem->offset_W_diag + i] < 0)
        {
            printf("skipping W_diag[%d] = %f\n", i, data[i]);
            continue;
        }
        blasfeo_dgein1(data[custom_mem->offset_W_diag + i], W_mat, i, i);
    }
}
{% endif %}



{%- if zoro_description.input_W_add_diag %}

/**
 * @brief Computes the adjusted GWG based on the uncertainties given W + stagewise varying part.
 */
static void compute_GWG_stagewise_varying(ocp_nlp_solver* solver, custom_memory* custom_mem, double* data, const int current_stage)
{
    ocp_nlp_dims* nlp_dims = solver->dims;

    // int N = nlp_dims->N;
    const int nx = nlp_dims->nx[0];
    const int nw = {{ zoro_description.nw }};

    // NOTE: update noise covariance terms for current stage
    for (int i = 0; i < nw; ++i)
    {
        blasfeo_dgein1(data[custom_mem->offset_W_add_diag + nw * current_stage + i], &custom_mem->W_stage_mat, i, i);
    }
    //   blasfeo_print_exp_dmat(nw, nw, &custom_mem->W_stage_mat, 0, 0);

    // NOTE: Combine constant W part with additive one
    blasfeo_dgead(nw, nw, 1.0, &custom_mem->W_mat, 0, 0, &custom_mem->W_stage_mat, 0, 0);
    //   blasfeo_print_exp_dmat(nw, nw, &custom_mem->W_stage_mat, 0, 0);

    // NOTE: Compute G@W@G^T term with W_stage_mat
    // temp_GW_mat = unc_jac_G_mat * W_stage_mat
    blasfeo_dgemm_nn(nx, nw, nw, 1.0, &custom_mem->unc_jac_G_mat, 0, 0, &custom_mem->W_stage_mat, 0, 0, 0.0,
                    &custom_mem->temp_GW_mat, 0, 0, &custom_mem->temp_GW_mat, 0, 0);
    // GWG_mat = temp_GW_stage_mat * unc_jac_G_mat^T
    blasfeo_dgemm_nt(nx, nx, nw, 1.0, &custom_mem->temp_GW_mat, 0, 0, &custom_mem->unc_jac_G_mat, 0, 0, 0.0,
                    &custom_mem->GWG_mat, 0, 0, &custom_mem->GWG_mat, 0, 0);
}
{% endif %}


static void update_riccati_quad_matrices(ocp_nlp_solver *solver, ocp_nlp_memory *nlp_mem, custom_memory* custom_mem, int ii)
{
    // extract dim
    ocp_nlp_dims* nlp_dims = solver->dims;
    int N = nlp_dims->N;
	int nbu = {{ dims.nbu}};
	int nbx = {{ dims.nbx}};
    int ng = {{ dims.ng }};
    int nh = {{ dims.nh }};
    int nx = {{ dims.nx }};
    int nu = {{ dims.nu }};

    blasfeo_dgecp(NCT, nu+nx, &custom_mem->dct_dux, 0, 0, &custom_mem->scaled_dct_dux, 0, 0);

    double temp_nominal_val;
    int ir = 0;
    // NOTE: this could be enhanced by using blasfeo_dgecpsc instead of blasfeo_dgesc.
// lbu
{%- if zoro_description.nlbu_t > 0 %}
    d_ocp_qp_get_lbu(ii, nlp_mem->qp_in, custom_mem->d_ineq_val);
{%- for it in zoro_description.idx_lbu_t %}
    temp_nominal_val = custom_mem->d_ineq_val[{{it}}] - sqrt(blasfeo_dvecex1(custom_mem->ineq_backoff_sq_buffer+ii, {{it}}));
{%- if zoro_description.feedback_optimization_mode == "RICCATI_BARRIER_1" %}
    blasfeo_dgesc(1, nu+nx, custom_mem->tau / (temp_nominal_val * temp_nominal_val), &custom_mem->scaled_dct_dux, ir, 0);
{%- else %}
    blasfeo_dgesc(1, nu+nx, -1.0 * custom_mem->tau / temp_nominal_val / (2 * sqrt(blasfeo_dvecex1(custom_mem->ineq_backoff_sq_buffer+ii, {{it}}))), &custom_mem->scaled_dct_dux, ir, 0);
{%- endif %}
    ir += 1;
{%- endfor %}
{%- endif %}
// lbx
{%- if zoro_description.nlbx_t > 0 %}
    d_ocp_qp_get_lbx(ii, nlp_mem->qp_in, custom_mem->d_ineq_val);
{%- for it in zoro_description.idx_lbx_t %}
    temp_nominal_val = custom_mem->d_ineq_val[{{it}}] - sqrt(blasfeo_dvecex1(custom_mem->ineq_backoff_sq_buffer+ii, nbu + {{it}}));
{%- if zoro_description.feedback_optimization_mode == "RICCATI_BARRIER_1" %}
    blasfeo_dgesc(1, nu+nx, custom_mem->tau / (temp_nominal_val * temp_nominal_val), &custom_mem->scaled_dct_dux, ir, 0);
{%- else %}
    blasfeo_dgesc(1, nu+nx, -1.0 * custom_mem->tau / temp_nominal_val / (2 * sqrt(blasfeo_dvecex1(custom_mem->ineq_backoff_sq_buffer+ii, nbu + {{it}}))), &custom_mem->scaled_dct_dux, ir, 0);
{%- endif %}
    ir += 1;
{%- endfor %}
{%- endif %}

{%- if zoro_description.nlh_t + zoro_description.nuh_t > 0 %}
    // Get C_k and D_k
    ocp_nlp_get_at_stage(solver, ii, "D", custom_mem->d_Dgh_mat);
    ocp_nlp_get_at_stage(solver, ii, "C", custom_mem->d_Cgh_mat);
{%- endif %}

{%- if zoro_description.nlg_t + zoro_description.nlh_t > 0 %}
    // lg
    d_ocp_qp_get_lg(ii, nlp_mem->qp_in, custom_mem->d_ineq_val);
{%- for it in zoro_description.idx_lg_t %}
    temp_nominal_val = custom_mem->d_ineq_val[{{it}}] - sqrt(blasfeo_dvecex1(custom_mem->ineq_backoff_sq_buffer+ii, nbu + nbx + {{it}}));
{%- if zoro_description.feedback_optimization_mode == "RICCATI_BARRIER_1" %}
    blasfeo_dgesc(1, nu+nx, custom_mem->tau / (temp_nominal_val * temp_nominal_val), &custom_mem->scaled_dct_dux, ir, 0);
{%- else %}
    blasfeo_dgesc(1, nu+nx, -1.0 * custom_mem->tau / temp_nominal_val / (2 * sqrt(blasfeo_dvecex1(custom_mem->ineq_backoff_sq_buffer+ii, nbu + nbx + {{it}}))),  &custom_mem->scaled_dct_dux, ir, 0);
{%- endif %}
    ir += 1;
{%- endfor %}
    // lh
{%- for it in zoro_description.idx_lh_t %}
    // NOTE: the d_Cgh_mat is column-major, the first ng rows are the Jacobians of the linear constraints
    blasfeo_pack_dmat(1, nu, custom_mem->d_Dgh_mat + ng + {{it}}, ng+nh, &custom_mem->dct_dux, 0, 0);
    blasfeo_pack_dmat(1, nu, custom_mem->d_Dgh_mat + ng + {{it}}, ng+nh, &custom_mem->scaled_dct_dux, 0, 0);
    blasfeo_pack_dmat(1, nx, custom_mem->d_Cgh_mat + ng + {{it}}, ng+nh, &custom_mem->dct_dux, ir, nu);
    blasfeo_pack_dmat(1, nx, custom_mem->d_Cgh_mat + ng + {{it}}, ng+nh, &custom_mem->scaled_dct_dux, ir, nu);
    temp_nominal_val = custom_mem->d_ineq_val[ng + {{it}}] - sqrt(blasfeo_dvecex1(custom_mem->ineq_backoff_sq_buffer+ii, nbu + nbx + ng + {{it}}));
{%- if zoro_description.feedback_optimization_mode == "RICCATI_BARRIER_1" %}
    blasfeo_dgesc(1, nu+nx, custom_mem->tau / (temp_nominal_val * temp_nominal_val), &custom_mem->scaled_dct_dux, ir, 0);
{%- else %}
    blasfeo_dgesc(1, nu+nx, -1.0 * custom_mem->tau / temp_nominal_val / (2 * sqrt(blasfeo_dvecex1(custom_mem->ineq_backoff_sq_buffer+ii, nbu + nbx + ng + {{it}}))), &custom_mem->scaled_dct_dux, ir, 0);
{%- endif %}
    ir += 1;
{%- endfor %}
{%- endif %}

// ubu
{%- if zoro_description.nubu_t > 0 %}
    d_ocp_qp_get_ubu(ii, nlp_mem->qp_in, custom_mem->d_ineq_val);
{%- for it in zoro_description.idx_ubu_t %}
    temp_nominal_val = -custom_mem->d_ineq_val[{{it}}] - sqrt(blasfeo_dvecex1(custom_mem->ineq_backoff_sq_buffer+ii, {{it}}));
{%- if zoro_description.feedback_optimization_mode == "RICCATI_BARRIER_1" %}
    blasfeo_dgesc(1, nu+nx, custom_mem->tau / (temp_nominal_val * temp_nominal_val), &custom_mem->scaled_dct_dux, ir, 0);
{%- else %}
    blasfeo_dgesc(1, nu+nx, -1.0 * custom_mem->tau / temp_nominal_val / (2 * sqrt(blasfeo_dvecex1(custom_mem->ineq_backoff_sq_buffer+ii, {{it}}))), &custom_mem->scaled_dct_dux, ir, 0);
{%- endif %}
    ir += 1;
{%- endfor %}
{%- endif %}
// ubx
{%- if zoro_description.nubx_t > 0 %}
    d_ocp_qp_get_ubx(ii, nlp_mem->qp_in, custom_mem->d_ineq_val);
{%- for it in zoro_description.idx_ubx_t %}
    temp_nominal_val = -custom_mem->d_ineq_val[{{it}}] - sqrt(blasfeo_dvecex1(custom_mem->ineq_backoff_sq_buffer+ii, nbu + {{it}}));
{%- if zoro_description.feedback_optimization_mode == "RICCATI_BARRIER_1" %}
    blasfeo_dgesc(1, nu+nx, custom_mem->tau / (temp_nominal_val * temp_nominal_val), &custom_mem->scaled_dct_dux, ir, 0);
{%- else %}
    blasfeo_dgesc(1, nu+nx, -1.0 * custom_mem->tau / temp_nominal_val / (2 * sqrt(blasfeo_dvecex1(custom_mem->ineq_backoff_sq_buffer+ii, nbu + {{it}}))), &custom_mem->scaled_dct_dux, ir, 0);
{%- endif %}
    ir += 1;
{%- endfor %}
{%- endif %}
// ug
{%- if zoro_description.nug_t + zoro_description.nuh_t > 0 %}
    d_ocp_qp_get_ug(ii, nlp_mem->qp_in, custom_mem->d_ineq_val);
{%- for it in zoro_description.idx_ug_t %}
    temp_nominal_val = - custom_mem->d_ineq_val[{{it}}] - sqrt(blasfeo_dvecex1(custom_mem->ineq_backoff_sq_buffer+ii, nbu + nbx + {{it}})) ;
{%- if zoro_description.feedback_optimization_mode == "RICCATI_BARRIER_1" %}
    blasfeo_dgesc(1, nu+nx, custom_mem->tau / (temp_nominal_val * temp_nominal_val), &custom_mem->scaled_dct_dux, ir, 0);
{%- else %}
    blasfeo_dgesc(1, nu+nx, -1.0 * custom_mem->tau / temp_nominal_val / (2 * sqrt(blasfeo_dvecex1(custom_mem->ineq_backoff_sq_buffer+ii, nbu + nbx + {{it}}))), &custom_mem->scaled_dct_dux, ir, 0);
{%- endif %}
    ir += 1;
{%- endfor %}
// uh
{%- for it in zoro_description.idx_uh_t %}
    // NOTE: the d_Cgh_mat is column-major, the first ng rows are the Jacobians of the linear constraints
    blasfeo_pack_dmat(1, nu, custom_mem->d_Dgh_mat + ng + {{it}}, ng+nh, &custom_mem->dct_dux, 0, 0);
    blasfeo_pack_dmat(1, nu, custom_mem->d_Dgh_mat + ng + {{it}}, ng+nh, &custom_mem->scaled_dct_dux, 0, 0);
    blasfeo_pack_dmat(1, nx, custom_mem->d_Cgh_mat + ng + {{it}}, ng+nh, &custom_mem->dct_dux, ir, nu);
    blasfeo_pack_dmat(1, nx, custom_mem->d_Cgh_mat + ng + {{it}}, ng+nh, &custom_mem->scaled_dct_dux, ir, nu);
    temp_nominal_val = -custom_mem->d_ineq_val[ng + {{it}}] - sqrt(blasfeo_dvecex1(custom_mem->ineq_backoff_sq_buffer+ii, nbu + nbx + ng + {{it}}));
{%- if zoro_description.feedback_optimization_mode == "RICCATI_BARRIER_1" %}
    blasfeo_dgesc(1, nu+nx, custom_mem->tau / (temp_nominal_val * temp_nominal_val), &custom_mem->scaled_dct_dux, ir, 0);
{%- else %}
    blasfeo_dgesc(1, nu+nx, -1.0 * custom_mem->tau / temp_nominal_val / (2 * sqrt(blasfeo_dvecex1(custom_mem->ineq_backoff_sq_buffer+ii, nbu + nbx + ng + {{it}}))), &custom_mem->scaled_dct_dux, ir, 0);
{%- endif %}
    ir += 1;
{%- endfor %}
{%- endif %}

    blasfeo_dgemm_tn(nx, nx, NCT, 1.0, &custom_mem->dct_dux, 0, nu,
        &custom_mem->scaled_dct_dux, 0, nu, 1.0, &custom_mem->riccati_Q_const, 0, 0, &custom_mem->riccati_Q_mat, 0, 0);
    blasfeo_dgemm_tn(nu, nu, NCT, 1.0, &custom_mem->dct_dux, 0, 0,
        &custom_mem->scaled_dct_dux, 0, 0, 1.0, &custom_mem->riccati_R_const, 0, 0, &custom_mem->riccati_R_mat, 0, 0);
    blasfeo_dgemm_tn(nu, nx, NCT, 1.0, &custom_mem->dct_dux, 0, 0,
        &custom_mem->scaled_dct_dux, 0, nu, 1.0, &custom_mem->riccati_S_const, 0, 0, &custom_mem->riccati_S_mat, 0, 0);

}


static void update_riccati_quad_matrices_terminal(ocp_nlp_solver *solver, ocp_nlp_memory *nlp_mem, custom_memory* custom_mem)
{
    // extract dim
    ocp_nlp_dims* nlp_dims = solver->dims;
    int N = nlp_dims->N;
	int nbx_e = {{ dims.nbx_e}};
    int ng_e = {{ dims.ng_e }};
    int nh_e = {{ dims.nh_e }};
    int nx = {{ dims.nx }};
    int nu = {{ dims.nu }};

    blasfeo_dgecp(NCT_E, nx, &custom_mem->dcet_dx, 0, 0, &custom_mem->scaled_dcet_dx, 0, 0);

    double temp_nominal_val;
    int ir = 0;
// lbx_e
{%- if zoro_description.nlbx_e_t > 0 %}
    d_ocp_qp_get_lbx(N, nlp_mem->qp_in, custom_mem->d_ineq_e_val);
{%- for it in zoro_description.idx_lbx_e_t %}
    temp_nominal_val = custom_mem->d_ineq_e_val[{{it}}] - sqrt(blasfeo_dvecex1(custom_mem->ineq_backoff_sq_buffer+N, {{it}}));
{%- if zoro_description.feedback_optimization_mode == "RICCATI_BARRIER_1" %}
    blasfeo_dgesc(1, nx, custom_mem->tau / (temp_nominal_val * temp_nominal_val), &custom_mem->scaled_dcet_dx, ir, 0);
{%- else %}
    blasfeo_dgesc(1, nx, -1.0 * custom_mem->tau / temp_nominal_val / (2 * sqrt(blasfeo_dvecex1(custom_mem->ineq_backoff_sq_buffer+N, {{it}}))), &custom_mem->scaled_dcet_dx, ir, 0);
{%- endif %}
    ir += 1;
{%- endfor %}
{%- endif %}

{%- if zoro_description.nlh_e_t + zoro_description.nuh_e_t > 0 %}
    // Get C_N
    ocp_nlp_get_at_stage(solver, N, "C", custom_mem->d_Cgh_mat);
{%- endif %}

{%- if zoro_description.nlg_e_t + zoro_description.nlh_e_t > 0 %}
    d_ocp_qp_get_lg(N, nlp_mem->qp_in, custom_mem->d_ineq_e_val);
// lg_e
{%- for it in zoro_description.idx_lg_e_t %}
    temp_nominal_val = custom_mem->d_ineq_e_val[{{it}}] - sqrt(blasfeo_dvecex1(custom_mem->ineq_backoff_sq_buffer+N, nbx_e + {{it}}));
{%- if zoro_description.feedback_optimization_mode == "RICCATI_BARRIER_1" %}
    blasfeo_dgesc(1, nx, custom_mem->tau / (temp_nominal_val * temp_nominal_val), &custom_mem->scaled_dcet_dx, ir, 0);
{%- else %}
    blasfeo_dgesc(1, nx, -1.0 * custom_mem->tau / temp_nominal_val / (2 * sqrt(blasfeo_dvecex1(custom_mem->ineq_backoff_sq_buffer+N, nbx_e + {{it}}))), &custom_mem->scaled_dcet_dx, ir, 0);
{%- endif %}
    ir += 1;
{%- endfor %}
// lh_e
{%- for it in zoro_description.idx_lh_e_t %}
    // NOTE: the d_Cgh_mat is column-major, the first ng_e rows are the Jacobians of the linear constraints
    blasfeo_pack_dmat(1, nx, custom_mem->d_Cgh_mat + ng_e + {{it}}, ng_e+nh_e, &custom_mem->dcet_dx, ir, 0);
    blasfeo_pack_dmat(1, nx, custom_mem->d_Cgh_mat + ng_e + {{it}}, ng_e+nh_e, &custom_mem->scaled_dcet_dx, ir, 0);
    temp_nominal_val = custom_mem->d_ineq_e_val[ng_e + {{it}}] - sqrt(blasfeo_dvecex1(custom_mem->ineq_backoff_sq_buffer+N, nbx_e + ng_e + {{it}}));
{%- if zoro_description.feedback_optimization_mode == "RICCATI_BARRIER_1" %}
    blasfeo_dgesc(1, nx, custom_mem->tau / (temp_nominal_val * temp_nominal_val), &custom_mem->scaled_dcet_dx, ir, 0);
{%- else %}
    blasfeo_dgesc(1, nx, -1.0 * custom_mem->tau / temp_nominal_val / (2 * sqrt(blasfeo_dvecex1(custom_mem->ineq_backoff_sq_buffer+N, nbx_e + ng_e + {{it}}))), &custom_mem->scaled_dcet_dx, ir, 0);
{%- endif %}
    ir += 1;
{%- endfor %}
{%- endif %}

// ubx_e
{%- if zoro_description.nubx_e_t > 0 %}
    d_ocp_qp_get_ubx(N, nlp_mem->qp_in, custom_mem->d_ineq_e_val);
{%- for it in zoro_description.idx_ubx_e_t %}
    temp_nominal_val = -custom_mem->d_ineq_e_val[{{it}}] - sqrt(blasfeo_dvecex1(custom_mem->ineq_backoff_sq_buffer+N, {{it}}));
{%- if zoro_description.feedback_optimization_mode == "RICCATI_BARRIER_1" %}
    blasfeo_dgesc(1, nx, custom_mem->tau / (temp_nominal_val * temp_nominal_val), &custom_mem->scaled_dcet_dx, ir, 0);
{%- else %}
    blasfeo_dgesc(1, nx, -1.0 * custom_mem->tau / temp_nominal_val / (2 * sqrt(blasfeo_dvecex1(custom_mem->ineq_backoff_sq_buffer+N, {{it}}))), &custom_mem->scaled_dcet_dx, ir, 0);
{%- endif %}
    ir += 1;
{%- endfor %}
{%- endif %}

{%- if zoro_description.nug_e_t + zoro_description.nuh_e_t > 0 %}
    d_ocp_qp_get_ug(N, nlp_mem->qp_in, custom_mem->d_ineq_e_val);
// ug_e
{%- for it in zoro_description.idx_ug_e_t %}
    temp_nominal_val = -custom_mem->d_ineq_e_val[{{it}}] - sqrt(blasfeo_dvecex1(custom_mem->ineq_backoff_sq_buffer+N, nbx_e + {{it}}));
{%- if zoro_description.feedback_optimization_mode == "RICCATI_BARRIER_1" %}
    blasfeo_dgesc(1, nx, custom_mem->tau / (temp_nominal_val * temp_nominal_val), &custom_mem->scaled_dcet_dx, ir, 0);
{%- else %}
    blasfeo_dgesc(1, nx, -1.0 * custom_mem->tau / temp_nominal_val / (2 * sqrt(blasfeo_dvecex1(custom_mem->ineq_backoff_sq_buffer+N, nbx_e + {{it}}))), &custom_mem->scaled_dcet_dx, ir, 0);
{%- endif %}
    ir += 1;
{%- endfor %}
// uh_e
{%- for it in zoro_description.idx_uh_e_t %}
    // NOTE: the d_Cgh_mat is column-major, the first ng_e rows are the Jacobians of the linear constraints
    blasfeo_pack_dmat(1, nx, custom_mem->d_Cgh_mat + ng_e + {{it}}, ng_e+nh_e, &custom_mem->dcet_dx, ir, 0);
    blasfeo_pack_dmat(1, nx, custom_mem->d_Cgh_mat + ng_e + {{it}}, ng_e+nh_e, &custom_mem->scaled_dcet_dx, ir, 0);
    temp_nominal_val = - custom_mem->d_ineq_e_val[ng_e + {{it}}] - sqrt(blasfeo_dvecex1(custom_mem->ineq_backoff_sq_buffer+N, nbx_e + ng_e + {{it}}));
{%- if zoro_description.feedback_optimization_mode == "RICCATI_BARRIER_1" %}
    blasfeo_dgesc(1, nx, custom_mem->tau / (temp_nominal_val * temp_nominal_val), &custom_mem->scaled_dcet_dx, ir, 0);
{%- else %}
    blasfeo_dgesc(1, nx, -1.0 * custom_mem->tau / temp_nominal_val / (2 * sqrt(blasfeo_dvecex1(custom_mem->ineq_backoff_sq_buffer+N, nbx_e + ng_e + {{it}}))), &custom_mem->scaled_dcet_dx, ir, 0);
{%- endif %}
    ir += 1;
{%- endfor %}
{%- endif %}

    blasfeo_dgemm_tn(nx, nx, NCT_E, 1.0, &custom_mem->dcet_dx, 0, 0, &custom_mem->scaled_dcet_dx, 0, 0, 1.0, &custom_mem->riccati_Q_const_e, 0, 0, &custom_mem->riccati_Q_mat, 0, 0);

}


static void riccati_recursion(ocp_nlp_solver* solver, ocp_nlp_memory *nlp_mem, custom_memory* custom_mem)
{
    ocp_nlp_dims* nlp_dims = solver->dims;
    int N = nlp_dims->N;
    int nx = nlp_dims->nx[0];
    int nu = nlp_dims->nu[0];

{%- if zoro_description.feedback_optimization_mode is containing("BARRIER") %}
    update_riccati_quad_matrices_terminal(solver, nlp_mem, custom_mem);
    blasfeo_dgecp(nx, nx, &custom_mem->riccati_Q_mat, 0, 0, &custom_mem->temp_riccati_P_plus_mat, 0, 0);
{%- else %}
    blasfeo_dgecp(nx, nx, &custom_mem->riccati_Q_const_e, 0, 0, &custom_mem->temp_riccati_P_plus_mat, 0, 0);
{%- endif %}

    for (int ii = N-1; ii >= 1; ii--)
    {
        // NOTE: It would be most efficient to get pointers to A, B matrices in blasfeo form directly to avoid unpacking AND packing.
        ocp_nlp_get_at_stage(solver, ii, "A", custom_mem->d_A_mat);
        blasfeo_pack_dmat(nx, nx, custom_mem->d_A_mat, nx, &custom_mem->A_mat, 0, 0);
        ocp_nlp_get_at_stage(solver, ii, "B", custom_mem->d_B_mat);
        blasfeo_pack_dmat(nx, nu, custom_mem->d_B_mat, nx, &custom_mem->B_mat, 0, 0);

    {%- if zoro_description.feedback_optimization_mode is containing("BARRIER") %}
        update_riccati_quad_matrices(solver, nlp_mem, custom_mem, ii);
    {%- endif %}

        // temp_riccati_BP_mat = B^T @ P
        blasfeo_dgemm_tn(nu, nx, nx, 1.0, &custom_mem->B_mat, 0, 0, &custom_mem->temp_riccati_P_plus_mat, 0, 0,
                            0.0, &custom_mem->temp_riccati_BP_mat, 0, 0, &custom_mem->temp_riccati_BP_mat, 0, 0);
        // temp_riccati_RaBPB_mat = R + temp_riccati_BP_mat @ B
        blasfeo_dgemm_nn(nu, nu, nx, 1.0, &custom_mem->temp_riccati_BP_mat, 0, 0, &custom_mem->B_mat, 0, 0,
                            1.0, &custom_mem->riccati_R_mat, 0, 0, &custom_mem->temp_riccati_RaBPB_mat, 0, 0);
        // temp_riccati_SaBPA_mat = S + temp_riccati_BP_mat @ A
        blasfeo_dgemm_nn(nu, nx, nx, 1.0, &custom_mem->temp_riccati_BP_mat, 0, 0, &custom_mem->A_mat, 0, 0,
                            1.0, &custom_mem->riccati_S_mat, 0, 0, &custom_mem->temp_riccati_SaBPA_mat, 0, 0);

        // temp_riccati_chol_mat = chol( temp_riccati_RaBPB_mat )
        blasfeo_dpotrf_l(nu, &custom_mem->temp_riccati_RaBPB_mat, 0, 0, &custom_mem->temp_riccati_chol_mat, 0, 0);
        // temp_riccati_cholinvSaBPA_mat = inv( temp_riccati_chol_mat ) @ temp_riccati_SaBPA_mat
        blasfeo_dtrsm_llnn(nu, nx, 1.0, &custom_mem->temp_riccati_chol_mat, 0, 0,
                            &custom_mem->temp_riccati_SaBPA_mat, 0, 0, &custom_mem->temp_riccati_cholinvSaBPA_mat, 0, 0);
        // custom_mem->riccati_K_buffer[ii] = inv( temp_riccati_chol_mat ).T @ temp_riccati_cholinvSaBPA_mat
        blasfeo_dtrsm_lltn(nu, nx, 1.0, &custom_mem->temp_riccati_chol_mat, 0, 0,
                            &custom_mem->temp_riccati_cholinvSaBPA_mat, 0, 0, &custom_mem->riccati_K_buffer[ii], 0, 0);

        // temp_riccati_P_mat = Q - temp_riccati_cholinvSaBPA_mat^T @ temp_riccati_cholinvSaBPA_mat
        blasfeo_dgemm_tn(nx, nx, nu, -1.0, &custom_mem->temp_riccati_cholinvSaBPA_mat, 0, 0, &custom_mem->temp_riccati_cholinvSaBPA_mat, 0, 0,
            1.0, &custom_mem->riccati_Q_mat, 0, 0, &custom_mem->temp_riccati_P_mat, 0, 0);
        // temp_riccati_AP_mat = A^T @ P
        blasfeo_dgemm_tn(nx, nx, nx, 1.0, &custom_mem->A_mat, 0, 0, &custom_mem->temp_riccati_P_plus_mat, 0, 0,
                            0.0, &custom_mem->temp_riccati_AP_mat, 0, 0, &custom_mem->temp_riccati_AP_mat, 0, 0);
        // temp_riccati_P_mat += temp_riccati_AP_mat @ A
        blasfeo_dgemm_nn(nx, nx, nx, 1.0, &custom_mem->temp_riccati_AP_mat, 0, 0, &custom_mem->A_mat, 0, 0,
                            1.0, &custom_mem->temp_riccati_P_mat, 0, 0, &custom_mem->temp_riccati_P_mat, 0, 0);
        // temp_riccati_P_plus_mat = temp_riccati_P_mat
        blasfeo_dgecp(nx, nx, &custom_mem->temp_riccati_P_mat, 0, 0, &custom_mem->temp_riccati_P_plus_mat, 0, 0);
    }

    // NOTE: K[0] is forced to be zero.
    blasfeo_dgese(nu, nx, 0., &custom_mem->riccati_K_buffer[0], 0, 0);
}


static void uncertainty_propagate_and_update(ocp_nlp_solver *solver, ocp_nlp_in *nlp_in, ocp_nlp_out *nlp_out, custom_memory *custom_mem, double* data, int data_len)
{
    ocp_nlp_config *nlp_config = solver->config;
    ocp_nlp_dims *nlp_dims = solver->dims;

    int N = nlp_dims->N;
    int nx = nlp_dims->nx[0];
    int nu = nlp_dims->nu[0];
    int nx_sqr = nx*nx;
    int nbx = {{ dims.nbx }};
    int nbu = {{ dims.nbu }};
    int ng = {{ dims.ng }};
    int nh = {{ dims.nh }};
    int ng_e = {{ dims.ng_e }};
    int nh_e = {{ dims.nh_e }};
    int nbx_e = {{ dims.nbx_e }};
    double backoff_scaling_gamma = {{ zoro_description.backoff_scaling_gamma }};
    double backoff_eps = 1e-8;
    struct blasfeo_dmat *K_mat = &custom_mem->K_mat;

{%- if zoro_description.feedback_optimization_mode != "CONSTANT_FEEDBACK" %}
    K_mat = &custom_mem->uncertainty_matrix_buffer[0];
{%- endif %}

    /* First Stage */
    // NOTE: lbx_0 and ubx_0 should not be tightened.
    // NOTE: lg_0 and ug_0 are not tightened.
    // NOTE: lh_0 and uh_0 are not tightened.
{%- if zoro_description.nlbu_t + zoro_description.nubu_t > 0 %}
    compute_KPK(K_mat, &custom_mem->temp_KP_mat,
                &custom_mem->temp_KPK_mat, &custom_mem->uncertainty_matrix_buffer[0], nx, nu);
    blasfeo_ddiaex_sp(nbu, backoff_scaling_gamma*backoff_scaling_gamma, custom_mem->idxbu, &custom_mem->temp_KPK_mat, 0, 0, &custom_mem->ineq_backoff_sq_buffer[0], 0);
{%- if zoro_description.feedback_optimization_mode is containing("BARRIER") %}
    blasfeo_dvecad(nbu, backoff_eps, &custom_mem->ricc_ones, 0, &custom_mem->ineq_backoff_sq_buffer[0], 0);
{%- endif %}

{%- if zoro_description.nlbu_t > 0 %}
    // backoff lbu
    {%- for it in zoro_description.idx_lbu_t %}
    custom_mem->d_lbu_tightened[{{it}}]
        = custom_mem->d_lbu[{{it}}] + sqrt(blasfeo_dvecex1(&custom_mem->ineq_backoff_sq_buffer[0], {{it}}));
    {%- endfor %}
    ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, nlp_out, 0, "lbu", custom_mem->d_lbu_tightened);
{%- endif %}
{%- if zoro_description.nubu_t > 0 %}
    // backoff ubu
    {%- for it in zoro_description.idx_ubu_t %}
    custom_mem->d_ubu_tightened[{{it}}]
        = custom_mem->d_ubu[{{it}}] - sqrt(blasfeo_dvecex1(&custom_mem->ineq_backoff_sq_buffer[0], {{it}}));
    {%- endfor %}
    ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, nlp_out, 0, "ubu", custom_mem->d_ubu_tightened);
{%- endif %}
{%- endif %}


    /* Middle Stages */
    // constraint tightening: for next stage based on dynamics of ii stage
    // P[ii+1] = (A-B@K) @ P[ii] @ (A-B@K).T + G@W@G.T
    for (int ii = 0; ii < N-1; ii++)
    {
{%- if zoro_description.feedback_optimization_mode != "CONSTANT_FEEDBACK" %}
        K_mat = &custom_mem->riccati_K_buffer[ii];
{%- endif %}

        // get and pack: A, B
        ocp_nlp_get_at_stage(solver, ii, "A", custom_mem->d_A_mat);
        blasfeo_pack_dmat(nx, nx, custom_mem->d_A_mat, nx, &custom_mem->A_mat, 0, 0);
        ocp_nlp_get_at_stage(solver, ii, "B", custom_mem->d_B_mat);
        blasfeo_pack_dmat(nx, nu, custom_mem->d_B_mat, nx, &custom_mem->B_mat, 0, 0);

{% if zoro_description.input_W_add_diag %}
        compute_GWG_stagewise_varying(solver, custom_mem, data, ii);
{% endif %}

        compute_next_P_matrix(&custom_mem->uncertainty_matrix_buffer[ii],
                              &custom_mem->uncertainty_matrix_buffer[ii+1],
                              &custom_mem->A_mat, &custom_mem->B_mat,
                              K_mat, &custom_mem->GWG_mat,
                              &custom_mem->AK_mat, &custom_mem->temp_AP_mat, nx, nu);

        // state constraints
{%- if zoro_description.nlbx_t + zoro_description.nubx_t > 0 %}
    blasfeo_ddiaex_sp(nbx, backoff_scaling_gamma*backoff_scaling_gamma, custom_mem->idxbx, &custom_mem->uncertainty_matrix_buffer[ii+1], 0, 0, &custom_mem->ineq_backoff_sq_buffer[ii+1], nbu);
{%- if zoro_description.feedback_optimization_mode is containing("BARRIER") %}
    blasfeo_dvecad(nbx, backoff_eps, &custom_mem->ricc_ones, 0, &custom_mem->ineq_backoff_sq_buffer[ii+1], nbu);
{%- endif %}
    {%- if zoro_description.nlbx_t > 0 %}
        // lbx
        {%- for it in zoro_description.idx_lbx_t %}
        custom_mem->d_lbx_tightened[{{it}}]
            = custom_mem->d_lbx[{{it}}] + sqrt(blasfeo_dvecex1(&custom_mem->ineq_backoff_sq_buffer[ii+1], nbu + {{it}}));
        {%- endfor %}
        ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, nlp_out, ii+1, "lbx", custom_mem->d_lbx_tightened);
    {%- endif %}
    {% if zoro_description.nubx_t > 0 %}
        // ubx
        {%- for it in zoro_description.idx_ubx_t %}
        custom_mem->d_ubx_tightened[{{it}}]
            = custom_mem->d_ubx[{{it}}] - sqrt(blasfeo_dvecex1(&custom_mem->ineq_backoff_sq_buffer[ii+1], nbu + {{it}}));
        {%- endfor %}
        ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, nlp_out, ii+1, "ubx", custom_mem->d_ubx_tightened);
    {%- endif %}
{%- endif %}

{%- if zoro_description.nlbu_t + zoro_description.nubu_t > 0 %}
        // input constraints
        compute_KPK(K_mat, &custom_mem->temp_KP_mat,
            &custom_mem->temp_KPK_mat, &custom_mem->uncertainty_matrix_buffer[ii+1], nx, nu);
        blasfeo_ddiaex_sp(nbu, backoff_scaling_gamma*backoff_scaling_gamma, custom_mem->idxbu, &custom_mem->temp_KPK_mat, 0, 0, &custom_mem->ineq_backoff_sq_buffer[ii+1], 0);
{%- if zoro_description.feedback_optimization_mode is containing("BARRIER") %}
        blasfeo_dvecad(nbu, backoff_eps, &custom_mem->ricc_ones, 0, &custom_mem->ineq_backoff_sq_buffer[ii+1], 0);
{%- endif %}
    {%- if zoro_description.nlbu_t > 0 %}
        {%- for it in zoro_description.idx_lbu_t %}
        custom_mem->d_lbu_tightened[{{it}}] = custom_mem->d_lbu[{{it}}] + sqrt(blasfeo_dvecex1(&custom_mem->ineq_backoff_sq_buffer[ii+1], {{it}}));
        {%- endfor %}

        ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, nlp_out, ii+1, "lbu", custom_mem->d_lbu_tightened);
    {%- endif %}
    {%- if zoro_description.nubu_t > 0 %}
        {%- for it in zoro_description.idx_ubu_t %}
        custom_mem->d_ubu_tightened[{{it}}] = custom_mem->d_ubu[{{it}}] - sqrt(blasfeo_dvecex1(&custom_mem->ineq_backoff_sq_buffer[ii+1], {{it}}));
        {%- endfor %}
        ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, nlp_out, ii+1, "ubu", custom_mem->d_ubu_tightened);
    {%- endif %}
{%- endif %}

{%- if zoro_description.nlg_t + zoro_description.nug_t > 0 %}
        // Linear constraints: g
        compute_gh_beta(K_mat, &custom_mem->Cg_mat,
                     &custom_mem->Dg_mat, &custom_mem->temp_CaDK_mat,
                     &custom_mem->temp_CaDKmP_mat, &custom_mem->temp_beta_mat,
                     &custom_mem->uncertainty_matrix_buffer[ii+1], ng, nx, nu);
        blasfeo_ddiaex(ng, backoff_scaling_gamma*backoff_scaling_gamma, &custom_mem->temp_beta_mat, 0, 0, &custom_mem->ineq_backoff_sq_buffer[ii+1], nbu + nbx);
{%- if zoro_description.feedback_optimization_mode is containing("BARRIER") %}
        blasfeo_dvecad(ng, backoff_eps, &custom_mem->ricc_ones, 0, &custom_mem->ineq_backoff_sq_buffer[ii+1], nbu + nbx);
{%- endif %}
    {%- if zoro_description.nlg_t > 0 %}
        {%- for it in zoro_description.idx_lg_t %}
        custom_mem->d_lg_tightened[{{it}}]
            = custom_mem->d_lg[{{it}}] + sqrt(blasfeo_dvecex1(&custom_mem->ineq_backoff_sq_buffer[ii+1], nbu + nbx + {{it}}));
        {%- endfor %}
        ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, nlp_out, ii+1, "lg", custom_mem->d_lg_tightened);
    {%- endif %}
    {%- if zoro_description.nug_t > 0 %}
        {%- for it in zoro_description.idx_ug_t %}
        custom_mem->d_ug_tightened[{{it}}]
            = custom_mem->d_ug[{{it}}] - sqrt(blasfeo_dvecex1(&custom_mem->ineq_backoff_sq_buffer[ii+1], nbu + nbx + {{it}}));
        {%- endfor %}
        ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, nlp_out, ii+1, "ug", custom_mem->d_ug_tightened);
    {%- endif %}
{%- endif %}


{%- if zoro_description.nlh_t + zoro_description.nuh_t > 0 %}
        // nonlinear constraints: h
        // Get C_{k+1} and D_{k+1}
        ocp_nlp_get_at_stage(solver, ii+1, "C", custom_mem->d_Cgh_mat);
        ocp_nlp_get_at_stage(solver, ii+1, "D", custom_mem->d_Dgh_mat);
        // NOTE: the d_Cgh_mat is column-major, the first ng rows are the Jacobians of the linear constraints
        blasfeo_pack_dmat(nh, nx, custom_mem->d_Cgh_mat+ng, ng+nh, &custom_mem->Ch_mat, 0, 0);
        blasfeo_pack_dmat(nh, nu, custom_mem->d_Dgh_mat+ng, ng+nh, &custom_mem->Dh_mat, 0, 0);

        compute_gh_beta(K_mat, &custom_mem->Ch_mat,
                     &custom_mem->Dh_mat, &custom_mem->temp_CaDK_mat,
                     &custom_mem->temp_CaDKmP_mat, &custom_mem->temp_beta_mat,
                     &custom_mem->uncertainty_matrix_buffer[ii+1], nh, nx, nu);
        blasfeo_ddiaex(nh, backoff_scaling_gamma*backoff_scaling_gamma, &custom_mem->temp_beta_mat, 0, 0, &custom_mem->ineq_backoff_sq_buffer[ii+1], nbu + nbx + ng);
{%- if zoro_description.feedback_optimization_mode is containing("BARRIER") %}
        blasfeo_dvecad(nh, backoff_eps, &custom_mem->ricc_ones, 0, &custom_mem->ineq_backoff_sq_buffer[ii+1], nbu + nbx + ng);
{%- endif %}
        // TODO: eval hessian(h) -> H_hess (nh*(nx+nu)**2)
        // temp_Kt_hhess = h_i_hess[:nx, :] + K^T * h_i_hess[nx:nx+nu, :]
        // tempCD = temp_CaDKmP_mat * temp_Kt_hhess
        // acados_CD += tempCD
        // += for upper or -= for lower bound

        // only works if 1 bound is trivially satisfied.

    {%- if zoro_description.nlh_t > 0 %}
        {%- for it in zoro_description.idx_lh_t %}
        custom_mem->d_lh_tightened[{{it}}]
            = custom_mem->d_lh[{{it}}] + sqrt(blasfeo_dvecex1(&custom_mem->ineq_backoff_sq_buffer[ii+1], nbu + nbx + ng + {{it}}));
        {%- endfor %}
        ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, nlp_out, ii+1, "lh", custom_mem->d_lh_tightened);
    {%- endif %}
    {%- if zoro_description.nuh_t > 0 %}
        {%- for it in zoro_description.idx_uh_t %}
        custom_mem->d_uh_tightened[{{it}}]
            = custom_mem->d_uh[{{it}}] - sqrt(blasfeo_dvecex1(&custom_mem->ineq_backoff_sq_buffer[ii+1], nbu + nbx + ng + {{it}}));
        {%- endfor %}
        ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, nlp_out, ii+1, "uh", custom_mem->d_uh_tightened);
    {%- endif %}
{%- endif %}
    }

    // Last stage
    // get and pack: A, B
    ocp_nlp_get_at_stage(solver, N-1, "A", custom_mem->d_A_mat);
    blasfeo_pack_dmat(nx, nx, custom_mem->d_A_mat, nx, &custom_mem->A_mat, 0, 0);
    ocp_nlp_get_at_stage(solver, N-1, "B", custom_mem->d_B_mat);
    blasfeo_pack_dmat(nx, nu, custom_mem->d_B_mat, nx, &custom_mem->B_mat, 0, 0);

{%- if zoro_description.feedback_optimization_mode != "CONSTANT_FEEDBACK" %}
    K_mat = &custom_mem->riccati_K_buffer[N-1];
{%- endif %}

{% if zoro_description.input_W_add_diag %}
    compute_GWG_stagewise_varying(solver, custom_mem, data, N - 1);
{%- endif %}

    // AK_mat = -B*K + A
    compute_next_P_matrix(&custom_mem->uncertainty_matrix_buffer[N-1],
                        &custom_mem->uncertainty_matrix_buffer[N],
                        &custom_mem->A_mat, &custom_mem->B_mat,
                        K_mat, &custom_mem->GWG_mat,
                        &custom_mem->AK_mat, &custom_mem->temp_AP_mat, nx, nu);

    // state constraints nlbx_e_t
{%- if zoro_description.nlbx_e_t + zoro_description.nubx_e_t > 0 %}
    blasfeo_ddiaex_sp(nbx_e, backoff_scaling_gamma*backoff_scaling_gamma, custom_mem->idxbx_e, &custom_mem->uncertainty_matrix_buffer[N], 0, 0, &custom_mem->ineq_backoff_sq_buffer[N], 0);
{%- if zoro_description.feedback_optimization_mode is containing("BARRIER") %}
    blasfeo_dvecad(nbx_e, backoff_eps, &custom_mem->ricc_ones, 0, &custom_mem->ineq_backoff_sq_buffer[N], 0);
{%- endif %}
{%- if zoro_description.nlbx_e_t > 0 %}
    // lbx_e
    {%- for it in zoro_description.idx_lbx_e_t %}
    custom_mem->d_lbx_e_tightened[{{it}}]
        = custom_mem->d_lbx_e[{{it}}] + sqrt(blasfeo_dvecex1(&custom_mem->ineq_backoff_sq_buffer[N], {{it}}));
    {%- endfor %}
    ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, nlp_out, N, "lbx", custom_mem->d_lbx_e_tightened);
{%- endif %}
{% if zoro_description.nubx_e_t > 0 %}
    // ubx_e
    {%- for it in zoro_description.idx_ubx_e_t %}
    custom_mem->d_ubx_e_tightened[{{it}}]
        = custom_mem->d_ubx_e[{{it}}] - sqrt(blasfeo_dvecex1(&custom_mem->ineq_backoff_sq_buffer[N], {{it}}));
    {%- endfor %}
    ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, nlp_out, N, "ubx", custom_mem->d_ubx_e_tightened);
{%- endif %}
{%- endif %}

{%- if zoro_description.nlg_e_t + zoro_description.nug_e_t > 0 %}
    // Linear constraints: g
    compute_gh_beta(K_mat, &custom_mem->Cg_mat,
                    &custom_mem->dummy_Dgh_e_mat, &custom_mem->temp_CaDK_mat,
                    &custom_mem->temp_CaDKmP_mat, &custom_mem->temp_beta_mat,
                    &custom_mem->uncertainty_matrix_buffer[N], ng_e, nx, nu);
    blasfeo_ddiaex(ng_e, backoff_scaling_gamma*backoff_scaling_gamma, &custom_mem->temp_beta_mat, 0, 0, &custom_mem->ineq_backoff_sq_buffer[N], nbx_e);
{%- if zoro_description.feedback_optimization_mode is containing("BARRIER") %}
    blasfeo_dvecad(ng_e, backoff_eps, &custom_mem->ricc_ones, 0, &custom_mem->ineq_backoff_sq_buffer[N], nbx_e);
{%- endif %}
{%- if zoro_description.nlg_e_t > 0 %}
    {%- for it in zoro_description.idx_lg_e_t %}
    custom_mem->d_lg_e_tightened[{{it}}]
        = custom_mem->d_lg_e[{{it}}] + sqrt(blasfeo_dvecex1(&custom_mem->ineq_backoff_sq_buffer[N], nbx_e + {{it}}));
    {%- endfor %}
    ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, nlp_out, N, "lg", custom_mem->d_lg_e_tightened);
{%- endif %}
{%- if zoro_description.nug_e_t > 0 %}
    {%- for it in zoro_description.idx_ug_e_t %}
    custom_mem->d_ug_e_tightened[{{it}}]
        = custom_mem->d_ug_e[{{it}}] - sqrt(blasfeo_dvecex1(&custom_mem->ineq_backoff_sq_buffer[N], nbx_e + {{it}}));
    {%- endfor %}
    ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, nlp_out, N, "ug", custom_mem->d_ug_e_tightened);
{%- endif %}
{%- endif %}


{%- if zoro_description.nlh_e_t + zoro_description.nuh_e_t > 0 %}
    // nonlinear constraints: h
    // Get C_{k+1} and D_{k+1}
    ocp_nlp_get_at_stage(solver, N, "C", custom_mem->d_Cgh_e_mat);
    // NOTE: the d_Cgh_e_mat is column-major, the first ng_e rows are the Jacobians of the linear constraints
    blasfeo_pack_dmat(nh_e, nx, custom_mem->d_Cgh_e_mat+ng_e, ng_e+nh_e, &custom_mem->Ch_mat, 0, 0);

    compute_gh_beta(K_mat, &custom_mem->Ch_mat,
                    &custom_mem->dummy_Dgh_e_mat, &custom_mem->temp_CaDK_mat,
                    &custom_mem->temp_CaDKmP_mat, &custom_mem->temp_beta_mat,
                    &custom_mem->uncertainty_matrix_buffer[N], nh_e, nx, nu);
    blasfeo_ddiaex(nh_e, backoff_scaling_gamma*backoff_scaling_gamma, &custom_mem->temp_beta_mat, 0, 0, &custom_mem->ineq_backoff_sq_buffer[N], nbx_e + ng_e);
{%- if zoro_description.feedback_optimization_mode is containing("BARRIER") %}
    blasfeo_dvecad(nh_e, backoff_eps, &custom_mem->ricc_ones, 0, &custom_mem->ineq_backoff_sq_buffer[N], nbx_e + ng_e);
{%- endif %}
    {%- if zoro_description.nlh_e_t > 0 %}
        {%- for it in zoro_description.idx_lh_e_t %}
        custom_mem->d_lh_e_tightened[{{it}}]
            = custom_mem->d_lh_e[{{it}}] + sqrt(blasfeo_dvecex1(&custom_mem->ineq_backoff_sq_buffer[N], nbx_e + ng_e + {{it}}));
        {%- endfor %}
        ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, nlp_out, N, "lh", custom_mem->d_lh_e_tightened);
    {%- endif %}
    {%- if zoro_description.nuh_e_t > 0 %}
        {%- for it in zoro_description.idx_uh_e_t %}
        custom_mem->d_uh_e_tightened[{{it}}]
            = custom_mem->d_uh_e[{{it}}] - sqrt(blasfeo_dvecex1(&custom_mem->ineq_backoff_sq_buffer[N], nbx_e + ng_e + {{it}}));
        {%- endfor %}
        ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, nlp_out, N, "uh", custom_mem->d_uh_e_tightened);
    {%- endif %}
{%- endif %}

}


int custom_update_function({{ model.name }}_solver_capsule* capsule, double* data, int data_len)
{
    custom_memory *custom_mem = (custom_memory *) capsule->custom_update_memory;
    ocp_nlp_config *nlp_config = {{ model.name }}_acados_get_nlp_config(capsule);
    ocp_nlp_dims *nlp_dims = {{ model.name }}_acados_get_nlp_dims(capsule);
    ocp_nlp_in *nlp_in = {{ model.name }}_acados_get_nlp_in(capsule);
    ocp_nlp_out *nlp_out = {{ model.name }}_acados_get_nlp_out(capsule);
    ocp_nlp_solver *nlp_solver = {{ model.name }}_acados_get_nlp_solver(capsule);
    void *nlp_opts = {{ model.name }}_acados_get_nlp_opts(capsule);
    ocp_nlp_memory *nlp_mem;
    nlp_config->get(nlp_config, nlp_dims, nlp_solver->mem, "nlp_mem", &nlp_mem);

    int N = nlp_dims->N;
    int nx = {{ dims.nx }};
    int nw = {{ zoro_description.nw }};

{%- if zoro_description.output_P_matrices or zoro_description.output_riccati_t or zoro_description.input_W_diag and not zoro_description.input_W_add_diag -%}
    if (data_len != {{ zoro_description.data_size }})
    {
        printf("custom_update_zoro: data_length does not match expected one. Got %d, expected {{ zoro_description.data_size }}\n", data_len);
        exit(1);
    }
{%- endif %}

{%- if zoro_description.input_P0_diag or zoro_description.input_P0 %}
    if (data_len > 0)
    {
        reset_P0_matrix(nlp_dims, &custom_mem->uncertainty_matrix_buffer[0], data);
    }
{%- endif %}

{%- if zoro_description.input_W_diag %}
    reset_process_noise_matrix(custom_mem, nlp_dims, &custom_mem->W_mat, data);
{%- endif %}

{%- if zoro_description.input_W_diag and not zoro_description.input_W_add_diag %}
    // compute GWG with updated W
    blasfeo_dgemm_nn(nx, nw, nw, 1.0, &custom_mem->unc_jac_G_mat, 0, 0,
                        &custom_mem->W_mat, 0, 0, 0.0,
                        &custom_mem->temp_GW_mat, 0, 0, &custom_mem->temp_GW_mat, 0, 0);
    // GWG_mat = temp_GW_mat * unc_jac_G_mat^T
    blasfeo_dgemm_nt(nx, nx, nw, 1.0, &custom_mem->temp_GW_mat, 0, 0,
                        &custom_mem->unc_jac_G_mat, 0, 0, 0.0,
                        &custom_mem->GWG_mat, 0, 0, &custom_mem->GWG_mat, 0, 0);
{%- endif %}

{%- if zoro_description.feedback_optimization_mode != "CONSTANT_FEEDBACK" %}
    acados_timer timer0;
    acados_tic(&timer0);
    riccati_recursion(nlp_solver, nlp_mem, custom_mem);
    double time_riccati = acados_toc(&timer0);
{%- endif %}
    uncertainty_propagate_and_update(nlp_solver, nlp_in, nlp_out, custom_mem, data, data_len);


{%- if zoro_description.output_P_matrices %}
    for (int i = 0; i < N+1; ++i)
    {
        blasfeo_unpack_dmat(nx, nx, &custom_mem->uncertainty_matrix_buffer[i], 0, 0,
                    &data[custom_mem->offset_P_out + i * nx * nx], nx);
    }
    {%- if zoro_description.output_riccati_t %}
        data[custom_mem->offset_P_out + (N+1) * nx * nx] = time_riccati;
    {%- endif %}
{%- else %}
    {%- if zoro_description.output_riccati_t %}
        data[custom_mem->offset_P_out] = time_riccati;
    {%- endif %}
{%- endif %}

    return 1;
}


int custom_update_terminate_function({{ model.name }}_solver_capsule* capsule)
{
    custom_memory *mem = capsule->custom_update_memory;

    free(mem->raw_memory);
    return 1;
}

// useful prints for debugging

/*
printf("A_mat:\n");
blasfeo_print_exp_dmat(nx, nx, &custom_mem->A_mat, 0, 0);
printf("B_mat:\n");
blasfeo_print_exp_dmat(nx, nu, &custom_mem->B_mat, 0, 0);
printf("K_mat:\n");
blasfeo_print_exp_dmat(nu, nx, &custom_mem->K_mat, 0, 0);
printf("AK_mat:\n");
blasfeo_print_exp_dmat(nx, nx, &custom_mem->AK_mat, 0, 0);
printf("temp_AP_mat:\n");
blasfeo_print_exp_dmat(nx, nx, &custom_mem->temp_AP_mat, 0, 0);
printf("W_mat:\n");
blasfeo_print_exp_dmat(nx, nx, &custom_mem->W_mat, 0, 0);
printf("P_k+1:\n");
blasfeo_print_exp_dmat(nx, nx, &custom_mem->uncertainty_matrix_buffer[ii+1], 0, 0);*/