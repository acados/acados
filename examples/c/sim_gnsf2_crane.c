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

#include <acados/sim/sim_gnsf2.h>
#include <acados/sim/sim_common.h>
#include <acados/sim/sim_gnsf_casadi_wrapper.h> // todo remove
#include "acados/utils/external_function_generic.h"

#include "examples/c/gnsf2_crane_model/gnsf2_crane_model.h"

int main() {

    // set up gnsf_dims
    int gnsf2_dims_size = gnsf2_dims_calculate_size();
    void *dims_memory = malloc(gnsf2_dims_size);
    gnsf2_dims *dims = gnsf2_dims_assign(dims_memory);
    gnsf2_get_dims(dims, get_ints_fun);

    // set up sim_dims
    sim_dims *simdim = (sim_dims *) dims; // typecasting works as gnsf_dims has entries of sim_dims at the beginning

    // set up sim config
    int config_size = sim_solver_config_calculate_size();
    void *config_mem = malloc(config_size);
    sim_solver_config *config = sim_solver_config_assign(config_mem);
    sim_gnsf2_config_initialize_default(config);


    // // set up gnsf_in
    int gnsf2_in_size = gnsf2_in_calculate_size(dims);
    void *gnsf2_in_mem = malloc(gnsf2_in_size);
    gnsf2_in *in = gnsf2_in_assign(dims, gnsf2_in_mem);


    for (int ii = 0; ii < dims->nx *(dims->nx +dims->nu); ii++) {
        in->S_forw[ii] = 0.0;
    }
    for (int ii = 0; ii < dims->nx; ii++) {
        in->S_forw[ii+ ii*dims->nx] = 1.0;
        in->x[ii] = 0.0;
    }
    for (int ii = 0; ii < dims->nx; ii++) {
        in->S_adj[ii] = 1.0;
    }
    in->x[2] = 0.8;
    in->u[0] = 40.108149413030752;
    in->u[1] = -50.446662212534974;

    // set up gnsf_opts
    int gnsf2_opts_size = gnsf2_opts_calculate_size(dims);
    void *gnsf2_opts_mem = malloc(gnsf2_opts_size);
    gnsf2_opts *opts = gnsf2_opts_assign(dims, gnsf2_opts_mem);
    opts->sens_forw = 1;
    opts->sens_adj = 1;
    opts->newton_max = 3;

    // set up gnsf_fixed
    int gnsf2_fixed_size = gnsf2_fixed_calculate_size(dims, opts);
    void *gnsf2_fixed_mem = malloc(gnsf2_fixed_size);
    gnsf2_fixed* fix = gnsf2_fixed_assign(dims, gnsf2_fixed_mem, gnsf2_fixed_size);
    gnsf2_import(dims, fix, But_KK_YY_ZZ_LO_fun);

    // NEW STYLE CASADI FCNZ:::
    // Phi_inc_dy
    external_function_casadi Phi_inc_dy;
    Phi_inc_dy.casadi_fun = &Phi_inc_dy_fun;
    Phi_inc_dy.casadi_work = &Phi_inc_dy_fun_work;
    Phi_inc_dy.casadi_sparsity_in  = &Phi_inc_dy_fun_sparsity_in;
    Phi_inc_dy.casadi_sparsity_out = &Phi_inc_dy_fun_sparsity_out;

    int Phi_inc_dy_size = external_function_casadi_calculate_size(&Phi_inc_dy);
    void *Phi_inc_dy_mem = malloc(Phi_inc_dy_size);
    external_function_casadi_assign(&Phi_inc_dy, Phi_inc_dy_mem);
    fix->Phi_inc_dy = (external_function_generic *) &Phi_inc_dy;

    // jac_Phi_y_fun
    external_function_casadi jac_Phi_y;
    jac_Phi_y.casadi_fun = &jac_Phi_y_fun;
    jac_Phi_y.casadi_work = &jac_Phi_y_fun_work;
    jac_Phi_y.casadi_sparsity_in  = &jac_Phi_y_fun_sparsity_in;
    jac_Phi_y.casadi_sparsity_out = &jac_Phi_y_fun_sparsity_out;

    int jac_Phi_y_size = external_function_casadi_calculate_size(&jac_Phi_y);
    void *jac_Phi_y_mem = malloc(jac_Phi_y_size);
    external_function_casadi_assign(&jac_Phi_y, jac_Phi_y_mem);
    fix->jac_Phi_y = (external_function_generic *) &jac_Phi_y;

    // f_LO_inc_J_x1k1uz
    external_function_casadi f_LO_inc_J_x1k1uz;
    f_LO_inc_J_x1k1uz.casadi_fun = &f_LO_inc_J_x1k1uz_fun;
    f_LO_inc_J_x1k1uz.casadi_work = &f_LO_inc_J_x1k1uz_fun_work;
    f_LO_inc_J_x1k1uz.casadi_sparsity_in  = &f_LO_inc_J_x1k1uz_fun_sparsity_in;
    f_LO_inc_J_x1k1uz.casadi_sparsity_out = &f_LO_inc_J_x1k1uz_fun_sparsity_out;

    int f_LO_inc_J_x1k1uz_size = external_function_casadi_calculate_size(&f_LO_inc_J_x1k1uz);
    void *f_LO_inc_J_x1k1uz_mem = malloc(f_LO_inc_J_x1k1uz_size);
    external_function_casadi_assign(&f_LO_inc_J_x1k1uz, f_LO_inc_J_x1k1uz_mem);
    fix->f_LO_inc_J_x1k1uz = (external_function_generic *) &f_LO_inc_J_x1k1uz;

    // // set up sim_out
    int sim_out_size = sim_out_calculate_size(config, simdim);
    void* sim_out_ptr = (void*) malloc(sim_out_size);
    sim_out* out = sim_out_assign(config, simdim, sim_out_ptr);

    // setup workspace
    int gnsf2_workspace_size = gnsf2_calculate_workspace_size(dims, opts);
    void *work_ = malloc(gnsf2_workspace_size);

    printf("Newton_iter = %d,\t num_steps = %d \n", opts->newton_max, dims->num_steps);
    int NREP = 10000;
    double casadi_times[NREP];
    double gnsf_times[NREP];

    for (int i = 0; i < NREP; i++) {
        gnsf2_simulate( dims, fix, in, out, opts, work_);

        casadi_times[i] = out->info->ADtime;
        gnsf_times[i] = out->info->CPUtime;
    }
    double casadi_time = minimum_of_doubles(casadi_times, NREP);
    double gnsf_time = minimum_of_doubles(gnsf_times, NREP);


    printf("xf =\n");
    d_print_e_mat(1, dims->nx, out->xn, 1);
    printf("forw_Sensitivities = \n");
    d_print_e_mat(dims->nx, dims->nx + dims->nu, out->S_forw, dims->nx);
    printf("adj Sensitivities =\n");
    d_print_e_mat(1, dims->nx + dims->nu, out->S_adj, 1);
    
    printf("gnsf2_time  =  %f [ms]\n", gnsf_time*1000);
    printf("casadi_time =  %f [ms]\n", casadi_time*1000);

    
    free(dims_memory);
    free(gnsf2_in_mem);
    free(gnsf2_fixed_mem);
    free(sim_out_ptr);
    free(gnsf2_opts_mem);
    free(work_);
    free(config_mem);

    free(jac_Phi_y_mem);
    free(Phi_inc_dy_mem);
    free(f_LO_inc_J_x1k1uz_mem);

    // return 0;
}