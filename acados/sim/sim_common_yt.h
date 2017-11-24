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

#ifndef ACADOS_SIM_SIM_COMMON_YT_H_
#define ACADOS_SIM_SIM_COMMON_YT_H_

#include <stdbool.h>

#include "acados/utils/timing.h"
#include "acados/utils/types.h"

typedef struct {
    int_t nx;   // NX
    int_t nu;   // NU
    int_t NF;   // NO. of forward sens

    // int_t nz;   // ALGEBRAIC VARIABLES: currently only internal, similar to ACADO code generation
    real_t *x;  // x[NX]
    real_t *u;  // u[NU]

    double *S_forw;  // forward seed
    double *S_adj;   // backward seed

    bool sens_forw;
    bool sens_adj;
    bool sens_hess;

    casadi_function_t vde;
    void (*VDE_forw)(const int_t, const int_t, const real_t *, real_t *, casadi_function_t);

    casadi_function_t adj;
    void (*VDE_adj)(const int_t, const int_t, const real_t *, real_t *, casadi_function_t);
    
    // casadi_function_t jac;
    // void (*jac_fun)(int, double *, double *, casadi_function_t);

    casadi_function_t hess;
    void (*Hess_fun)(const int_t, const int_t, const real_t *, real_t *, casadi_function_t);

    casadi_function_t impl_ode;
    void (*eval_impl_res)(const int, const int, const double *, double *, casadi_function_t); // function pointer to residuals of implicit ode

    casadi_function_t impl_jac_x;
    void (*eval_impl_jac_x)(const int, const int, const double *, double *, casadi_function_t); // function pointer to jacobian of implicit ode

    casadi_function_t impl_jac_xdot;
    void (*eval_impl_jac_xdot)(const int, const int, const double *, double *, casadi_function_t); // function pointer to jacobian of implicit ode

    casadi_function_t impl_jac_u;
    void (*eval_impl_jac_u)(const int, const int, const double *, double *, casadi_function_t); // function pointer to jacobian of implicit ode

    real_t step;
    int_t num_steps;

} sim_in;

typedef struct {
    double CPUtime;
    double LAtime;
    double ADtime;
} sim_info;

typedef struct {
    double *xn;      // xn[NX]
    double *S_forw;  // S_forw[NX*(NX+NU)]
    double *S_adj;   //
    double *S_hess;  //

    sim_info *info;
} sim_out;

int_t sim_in_calculate_size(int_t nx, int_t nu, int_t NF);
char *assign_sim_in(int_t nx, int_t nu, int_t NF, sim_in **in, void *ptr);
sim_in *create_sim_in(int_t nx, int_t nu, int_t NF);

int_t sim_out_calculate_size(int_t nx, int_t nu, int_t NF);
char *assign_sim_out(int_t nx, int_t nu, int_t NF, sim_out **out, void *ptr);
sim_out *create_sim_out(int_t nx, int_t nu, int_t NF);
#endif  // ACADOS_SIM_SIM_YT_COMMON_H_
