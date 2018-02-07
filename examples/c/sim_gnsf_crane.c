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

// external
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
// acados
// #include <acados_c/sim.h>
// #include <acados_c/options.h>

#include <acados/sim/sim_gnsf.h>
//#include <acados/sim/sim_erk_integrator.h>
//#include <acados/sim/sim_casadi_wrapper.h>

#include "examples/c/gnsf_crane_model/gnsf_crane_model.h"

// // blasfeo
// #include <blasfeo/include/blasfeo_target.h>
// #include <blasfeo/include/blasfeo_common.h>
// #include <blasfeo/include/blasfeo_d_aux.h>
// #include <blasfeo/include/blasfeo_d_aux_ext_dep.h>
// #include <blasfeo/include/blasfeo_v_aux_ext_dep.h>
// #include <blasfeo/include/blasfeo_d_blas.h>
int main() {

    int ii;
    int jj;

    int nx = 6;
    int nu = 2;
    int NF = nx + nu; // columns of forward seed

    double T = 0.1;
    int num_stages = 4;
/*
    sim_dims dims;
    dims.num_stages = num_stages;
    dims.nx = nx;
    dims.nu = nu; */

    gnsf_dims dims;
    dims.nu = 2;
    dims.nx = 9;
    printf("dims.nx = %d \n", dims.nx);
    print_gnsf_dims(dims);
    //sim_rk_opts *erk_opts = (sim_rk_opts *) args;
    //erk_opts->num_steps = 4;

    return nx;
}