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
#include <acados/sim/sim_common.h>
#include <acados/sim/sim_gnsf_casadi_wrapper.h>

#include "examples/c/gnsf_crane_model/gnsf_crane_model.h"

// blasfeo
// #include <blasfeo/include/blasfeo_target.h>
// #include <blasfeo/include/blasfeo_common.h>
// #include <blasfeo/include/blasfeo_d_aux.h>
// #include <blasfeo/include/blasfeo_d_aux_ext_dep.h>
// #include <blasfeo/include/blasfeo_v_aux_ext_dep.h>
// #include <blasfeo/include/blasfeo_d_blas.h>
int main() {

    gnsf_dims dims;
    gnsf_get_dims(&dims, get_ints_fun);
    sim_dims simdim;
    simdim.nx = 9;
    simdim.nu = 2;
    simdim.num_stages = 4;

    gnsf_in in;
    gnsf_opts opts;
    opts.sens_forw = (bool) calloc(1, sizeof(bool));
    opts.sens_adj = (bool) calloc(1, sizeof(bool));
    opts.sens_forw = 1;
    opts.sens_adj = 1;

    int gnsf_fixed_size = gnsf_fixed_calculate_size(&dims, &opts);
    void *raw_memory_ptr = malloc(gnsf_fixed_size);

    gnsf_fixed fix = *gnsf_fixed_assign(&dims, raw_memory_ptr);
    gnsf_import(&dims, &fix, But_KK_ZZ_LO_fun);
    printf("left_import \n");

    // opts.A_mat = (double*) calloc(dims.num_stages * dims.num_stages, sizeof(double));
    // *opts.A_mat = 1;
    
    // gnsf_allocate_fixed(&dims,&fix);
    // dims.num_steps = 2;
    

    in.res_inc_Jff = res_inc_Jff_fun;
    in.f_LO_inc_J_x1k1uz = f_LO_inc_J_x1k1uz_fun;
    in.jac_res_ffx1u = jac_res_ffx1u_fun;
    sim_out* out;

    int sim_out_size = sim_out_calculate_size(&simdim);
    void* sim_out_ptr = (void*) calloc(sim_out_size, sizeof(double));
    out = assign_sim_out(&simdim, sim_out_ptr);
    
    in.u = (double*) calloc(dims.nu, sizeof(double));
    in.x = (double*) calloc(dims.nx, sizeof(double));
    in.S_forw = (double*) calloc(dims.nx * (dims.nx + dims.nu), sizeof(double));

    for (int ii = 0; ii < dims.nx; ii++) {
        in.S_forw[ii+ ii*dims.nx] = 1.0;
    }
    in.x[2] = 0.8;
    in.u[0] = 40.108149413030752;
    in.u[1] = -50.446662212534974;

    int gnsf_workspace_size = gnsf_calculate_workspace_size(&dims, &opts);
    void *work_ = malloc(gnsf_workspace_size);
    // void* work_ = (void*) calloc(8, sizeof(double));

    printf("test\n");
    int num_executions = 100;
    double gnsf_time = 0;
    double casadi_time = 0;
    for (int i = 0; i < num_executions; i++) {
        gnsf_simulate( &dims, &fix, &in, out, &opts, work_);
        gnsf_time += out->info->CPUtime;
        casadi_time += out->info->ADtime;
    }
    gnsf_time = gnsf_time/num_executions;
    casadi_time = casadi_time/num_executions;
    printf("\n gnsf _time =  %f \n", gnsf_time);
    printf("casadi_time =  %f \n", casadi_time);

    free(work_);
    // free(dims);
    free(raw_memory_ptr);
    free(sim_out_ptr);

    return 0;
}