/*
 *    This file is part of ACADOS.
 *
 *    ACADOS is free software; you can redistribute it and/or
 *    modify it under the terms of the GNU Lesser General Public
 *    License as published by the Free Software Foundation; either
 *    version 3 of the License, or (at your option) any later version.
 *
 *    ACADOS is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *    Lesser General Public License for more details.
 *
 *    You should have received a copy of the GNU Lesser General Public
 *    License along with ACADOS; if not, write to the Free Software Foundation,
 *    Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 *
 */

#ifndef ACADOS_SIM_SIM_COMMON_H_
#define ACADOS_SIM_SIM_COMMON_H_

#include "acados/utils/types.h"
#include "acados/utils/timing.h"

#define FIXED_STEP_SIZE 0

typedef struct sim_in_ {
    int_t nx;   // NX
    int_t nu;   // NU
    real_t *x;  // x[NX]
    real_t *u;  // u[NU]

    void (*VDE_fun)(const real_t*, real_t*);
    void (*jac_fun)(const real_t*, real_t*);

#if FIXED_STEP_SIZE == 0
    real_t step;
    unsigned int nSteps;
#endif
} sim_in;

typedef struct sim_info_ {
    real_t CPUtime;
    real_t LAtime;
    real_t ADtime;
} sim_info;

typedef struct sim_out_ {
    real_t *xn;     // xn[NX]
    real_t *Sx;     // Sx[NX*NX]
    real_t *Su;     // Su[NX*NU]

    sim_info *info;
} sim_out;


#endif  // ACADOS_SIM_SIM_COMMON_H_
