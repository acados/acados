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

#ifndef ACADOS_SIM_SIM_COMMON_H_
#define ACADOS_SIM_SIM_COMMON_H_

#include <stdbool.h>

#include "acados/utils/timing.h"
#include "acados/utils/types.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    int_t nx;   // NX
    int_t nu;   // NU
    real_t *x;  // x[NX]
    real_t *u;  // u[NU]

    real_t *S_forw;  // forward seed
    real_t *S_adj;   // backward seed

    bool sens_forw;
    bool sens_adj;
    bool sens_hess;
    int_t num_forw_sens;

    casadi_function_t vde;
    void (*forward_vde_wrapper)(const int_t, const int_t, const real_t *, real_t *,
                                casadi_function_t);
    casadi_function_t jac;
    void (*jacobian_wrapper)(const int_t, const real_t *, real_t *, casadi_function_t);
    casadi_function_t vde_adj;
    void (*adjoint_vde_wrapper)(const int_t, const int_t, const real_t *, real_t *,
                                casadi_function_t);

    real_t step;
    uint num_steps;

    real_t *grad_K;  // gradient correction
} sim_in;

typedef struct {
    real_t CPUtime;
    real_t LAtime;
    real_t ADtime;
} sim_info;

typedef struct {
    real_t *xn;      // xn[NX]
    real_t *S_forw;  // S_forw[NX*(NX+NU)]
    real_t *S_adj;   //
    real_t *S_hess;  //

    real_t *grad;  // gradient correction

    sim_info *info;
} sim_out;

typedef struct {
    int_t (*fun)(const sim_in *, sim_out *, void *, void *, void *);
    sim_in *in;
    sim_out *out;
    void *args;
    void *mem;
    void *work;
} sim_solver;

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // ACADOS_SIM_SIM_COMMON_H_
