/*    /home/dang/acados/acados/sim/sim_common.h
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

#include "timing.h"
#include "types.h"

typedef struct sim_in_ {
    int_t nx;   // NX
    int_t nu;   // NU
    real_t *x;  // x[NX]
    real_t *u;  // u[NU]

    real_t *S_forw;     // forward seed
    real_t *S_adj;      // backward seed

    bool sens_forw;
    bool sens_adj;
    bool sens_hess;
    int_t nsens_forw;

    void (*VDE_forw)(const real_t*, real_t*);
    void (*VDE_adj)(const real_t*, real_t*);
    void (*jac_fun)(const real_t*, real_t*);

    real_t step;
    uint nSteps;

    void *opts;
} sim_in;

typedef struct sim_info_ {
    real_t CPUtime;
    real_t LAtime;
    real_t ADtime;
} sim_info;

typedef struct sim_out_ {
    real_t *xn;         // xn[NX]
    real_t *S_forw;     // S_forw[NX*(NX+NU)]
    real_t *S_adj;      //
    real_t *S_hess;     //

    sim_info *info;
} sim_out;

typedef struct sim_solver_ {
    void (*fun)(const sim_in*, sim_out*, void*, void*);
    sim_in *in;
    sim_out *out;
    void *mem;
    void *work;
} sim_solver;

#endif  // ACADOS_SIM_SIM_COMMON_H_
