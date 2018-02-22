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

// blasfeo
// #include <blasfeo/include/blasfeo_target.h>
// #include <blasfeo/include/blasfeo_common.h>
// #include <blasfeo/include/blasfeo_d_aux.h>
// #include <blasfeo/include/blasfeo_d_aux_ext_dep.h>
// #include <blasfeo/include/blasfeo_v_aux_ext_dep.h>
// #include <blasfeo/include/blasfeo_d_blas.h>
int main() {
    int dims_size = gnsf_dims_calculate_size();
    void *dims_memory = malloc(dims_size);
    printf("%d\n", dims_size);
    gnsf_dims *dims = assign_gnsf_dims(dims_memory);
    gnsf_get_dims(dims, get_ints_fun);

    int sim_dims_size = sim_dims_calculate_size();
    void *sim_dims_mem = malloc(sim_dims_size);
    printf("%d\n", sim_dims_size);
    sim_dims *simdim =assign_sim_dims(sim_dims_mem);

    simdim->nx = 9;
    simdim->nu = 2;
    simdim->num_stages = 4;

    int gnsf_in_size = gnsf_in_calculate_size(dims);
    void *gnsf_in_mem = malloc(gnsf_in_size);
    printf("%d\n", gnsf_in_size);
    gnsf_in *in = assign_gnsf_in(dims, gnsf_in_mem);

    gnsf_opts opts;
    opts.sens_forw = (bool) calloc(1, sizeof(bool));
    opts.sens_adj = (bool) calloc(1, sizeof(bool));
    opts.sens_forw = 1;
    opts.sens_adj = 1;

    int gnsf_fixed_size = gnsf_fixed_calculate_size(dims, &opts);
    void *gnsf_fixed_mem = malloc(gnsf_fixed_size);
    gnsf_fixed* fix = gnsf_fixed_assign(dims, gnsf_fixed_mem, gnsf_fixed_size);
    gnsf_import(dims, fix, But_KK_ZZ_LO_fun);

    in->res_inc_Jff = res_inc_Jff_fun;
    in->f_LO_inc_J_x1k1uz = f_LO_inc_J_x1k1uz_fun;
    in->jac_res_ffx1u = jac_res_ffx1u_fun;    
    for (int ii = 0; ii < dims->nx *(dims->nx +dims->nu); ii++) {
        in->S_forw[ii] = 0.0;
    }

    d_print_mat(dims->nx, 1, in->x, 1);
    printf("\n loop :  \n");
    for (int ii = 0; ii < dims->nx; ii++) {
        // in->S_forw[ii+ ii*dims->nx] = 1.0;
        // printf("\n adress \n%p",(void*)&fix->KKf);
        // printf("ii = %d \n", ii);
        // printf("nx = %d \n", dims->nx);
        // d_print_mat(dims->nx, dims->nx + dims->nu, in->S_forw, dims->nx);
        // printf("\n adress %p\n",(void*)&fix->KKf);
        // blasfeo_print_dmat(dims->nx1 * dims->num_stages, dims->n_out * dims->num_stages, &fix->KKf, 0, 0);
        in->x[ii] = 0.0;
    }
    in->x[2] = 0.8;
    in->u[0] = 40.108149413030752;
    in->u[1] = -50.446662212534974;

    sim_out* out;

    int sim_out_size = sim_out_calculate_size(simdim);
    void* sim_out_ptr = (void*) calloc(sim_out_size, sizeof(double));
    out = assign_sim_out(simdim, sim_out_ptr);


    int gnsf_workspace_size = gnsf_calculate_workspace_size(dims, &opts);
    void *work_ = malloc(gnsf_workspace_size);

    int num_executions = 1;
    double gnsf_time = 0;
    double casadi_time = 0;
    for (int i = 0; i < num_executions; i++) {
        gnsf_simulate( dims, fix, in, out, &opts, work_);
        gnsf_time += out->info->CPUtime;
        casadi_time += out->info->ADtime;
    }
    gnsf_time = gnsf_time/num_executions;
    casadi_time = casadi_time/num_executions;
    printf("\n gnsf _time =  %f \n", gnsf_time);
    printf("casadi_time =  %f \n", casadi_time);

    free(work_);
    free(gnsf_fixed_mem);
    free(sim_out_ptr);
    free(dims_memory);
    free(sim_dims_mem);
    free(gnsf_in_mem);

    return 0;
}