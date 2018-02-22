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
#include <assert.h>
// acados
// #include <acados_c/sim.h>
// #include <acados_c/options.h>

#include <acados/sim/sim_gnsf.h>
#include <acados/sim/sim_common.h>
#include <acados/sim/sim_gnsf_casadi_wrapper.h>

#include "examples/c/gnsf_crane_model/gnsf_crane_model.h"

int main() {
    // set up gnsf_dims
    int dims_size = gnsf_dims_calculate_size();
    void *dims_memory = malloc(dims_size);
    gnsf_dims *dims = assign_gnsf_dims(dims_memory);
    gnsf_get_dims(dims, get_ints_fun);

    // set up sim_dims
    int sim_dims_size = sim_dims_calculate_size();
    void *sim_dims_mem = malloc(sim_dims_size);
    sim_dims *simdim = assign_sim_dims(sim_dims_mem);
    simdim->nx = 9;
    simdim->nu = 2;
    simdim->num_stages = 4;

    // set up gnsf_in
    int gnsf_in_size = gnsf_in_calculate_size(dims);
    void *gnsf_in_mem = malloc(gnsf_in_size);
    gnsf_in *in = assign_gnsf_in(dims, gnsf_in_mem);
    in->res_inc_Jff = res_inc_Jff_fun;
    in->f_LO_inc_J_x1k1uz = f_LO_inc_J_x1k1uz_fun;
    in->jac_res_ffx1u = jac_res_ffx1u_fun;    
    for (int ii = 0; ii < dims->nx *(dims->nx +dims->nu); ii++) {
        in->S_forw[ii] = 0.0;
    }
    for (int ii = 0; ii < dims->nx; ii++) {
        in->S_forw[ii+ ii*dims->nx] = 1.0;
        in->x[ii] = 0.0;
    }
    in->x[2] = 0.8;
    in->u[0] = 40.108149413030752;
    in->u[1] = -50.446662212534974;

    // set up gnsf_opts
    int gnsf_opts_size = gnsf_opts_calculate_size(dims);
    void *gnsf_opts_mem = malloc(gnsf_opts_size);
    gnsf_opts *opts = assign_gnsf_opts(dims, gnsf_opts_mem);
    opts->sens_forw = 1;
    opts->sens_adj = 1;
    opts->newton_max = 10;

    // set up gnsf_fixed
    int gnsf_fixed_size = gnsf_fixed_calculate_size(dims, opts);
    void *gnsf_fixed_mem = malloc(gnsf_fixed_size);
    gnsf_fixed* fix = gnsf_fixed_assign(dims, gnsf_fixed_mem, gnsf_fixed_size);
    gnsf_import(dims, fix, But_KK_ZZ_LO_fun);

    // set up sm_out
    int sim_out_size = sim_out_calculate_size(simdim);
    void* sim_out_ptr = (void*) malloc(sim_out_size);
    sim_out* out = assign_sim_out(simdim, sim_out_ptr);

    // setup workspace
    int gnsf_workspace_size = gnsf_calculate_workspace_size(dims, opts);
    void *work_ = malloc(gnsf_workspace_size);

    int num_executions = 1000;
    double gnsf_time = 0;
    double casadi_time = 0;
    for (int i = 0; i < num_executions; i++) {
        gnsf_simulate( dims, fix, in, out, opts, work_);
        gnsf_time += out->info->CPUtime;
        casadi_time += out->info->ADtime;
    }
    printf("\nforw_Sensitivities = \n");
    d_print_e_mat(dims->nx, dims->nx + dims->nu, out->S_forw, dims->nx);
    gnsf_time = gnsf_time/num_executions;
    casadi_time = casadi_time/num_executions;
    printf("\n gnsf _time =  %f \n", gnsf_time);
    printf("casadi_time =  %f \n", casadi_time);

    free(dims_memory);
    free(sim_dims_mem);
    free(gnsf_in_mem);
    free(gnsf_opts_mem);
    free(gnsf_fixed_mem);
    free(sim_out_ptr);
    free(work_);

    return 0;
}