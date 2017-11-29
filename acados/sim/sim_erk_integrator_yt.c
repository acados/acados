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

// standard
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
// acados
#include "acados/utils/mem.h"
#include "acados/utils/print.h"
#include "acados/sim/sim_rk_common_yt.h"
#include "acados/sim/sim_common_yt.h"
#include "acados/sim/sim_erk_integrator_yt.h"

#include "acados/sim/sim_casadi_wrapper.h"

// TODO(dimitris): rename sim_erk_calculate_memory_size
int erk_calculate_memory_size(sim_dims *dims, void *opts_)
{
    sim_RK_opts *opts = (sim_RK_opts *) opts_;

    int nx = dims->nx;
    int nu = dims->nu;
    int NF = opts->num_forw_sens;

    int num_stages = opts->num_stages; // number of stages
    int nX = nx*(1+NF); // (nx) for ODE and (NF*nx) for VDE
    int nhess = (NF + 1) * NF / 2;
    uint num_steps = opts->num_steps;  // number of steps

    int size = sizeof(sim_erk_memory);

    size += (nX + nu) * sizeof(double); // rhs_forw_in

    if(opts->sens_adj){
        size += num_steps * num_stages * nX * sizeof(double); // K_traj
        size += (num_steps + 1) * nX *sizeof(double); // out_forw_traj
    }else{
        size += num_stages * nX * sizeof(double); // K_traj
        size += nX *sizeof(double); // out_forw_traj
    }

    if (opts->sens_hess && opts->sens_adj){
        size += (nx + nX + nu) * sizeof(double); //rhs_adj_in
        size += (nx + nu + nhess) * sizeof(double); //out_adj_tmp
        size += num_stages * (nx + nu + nhess) * sizeof(double); //adj_traj
    }else if (opts->sens_adj){
        size += (nx * 2 + nu) * sizeof(double); //rhs_adj_in
        size += (nx + nu)* sizeof(double); //out_adj_tmp
        size += num_stages * (nx + nu) * sizeof(double); //adj_traj
    }

    make_int_multiple_of(8, &size);
    size += 1 * 8;

    return size;
}


// TODO(dimitris): rename sim_erk_assign_memory
void *assign_erk_memory(sim_dims *dims, void *opts_, void *raw_memory)
{
    sim_RK_opts *opts = (sim_RK_opts *) opts_;

    int nx = dims->nx;
    int nu = dims->nu;
    int NF = opts->num_forw_sens;

    int num_stages = opts->num_stages; // number of stages
    int nX = nx*(1+NF); // (nx) for ODE and (NF*nx) for VDE
    int nhess = (NF + 1) * NF / 2;
    int num_steps = opts->num_steps;  // number of steps

    char *c_ptr = (char *)raw_memory;

    sim_erk_memory *mem = (sim_erk_memory *) c_ptr;
    c_ptr += sizeof(sim_erk_memory);

    // align memory to typical cache line size
    align_char_to(8, &c_ptr);

    assign_double(nX + nu, &mem->rhs_forw_in, &c_ptr);

    if(opts->sens_adj)
    {
        assign_double(num_steps*num_stages*nX, &mem->K_traj, &c_ptr);
        assign_double((num_steps + 1)*nX, &mem->out_forw_traj, &c_ptr);
    } else
    {
        assign_double(num_stages*nX, &mem->K_traj, &c_ptr);
        assign_double(nX, &mem->out_forw_traj, &c_ptr);
    }

    if (opts->sens_hess && opts->sens_adj)
    {
        assign_double(nx+nX+nu, &mem->rhs_adj_in, &c_ptr);
        assign_double(nx+nu+nhess, &mem->out_adj_tmp, &c_ptr);
        assign_double(num_stages*(nx+nu+nhess), &mem->adj_traj, &c_ptr);
    } else if (opts->sens_adj)
    {
        assign_double((nx*2+nu), &mem->rhs_adj_in, &c_ptr);
        assign_double(nx+nu, &mem->out_adj_tmp, &c_ptr);
        assign_double(num_stages*(nx+nu), &mem->adj_traj, &c_ptr);
    }

    assert((char*)raw_memory + erk_calculate_memory_size(dims, opts_) >= c_ptr);

    return (void *)mem;
}



void *sim_erk_create_memory(sim_dims *dims, void *opts_)
{
    int bytes = erk_calculate_memory_size(dims, opts_);
    void *ptr = malloc(bytes);
    sim_erk_memory *memory = assign_erk_memory(dims, opts_, ptr);

    return (void *)memory;
}



int sim_erk_calculate_workspace_size(sim_dims *dims, void *opts_)
{
    return 0;
}



int sim_erk_yt(sim_in *in, sim_out *out, void *opts_, void *mem_, void *work_)
{
    sim_RK_opts *opts = (sim_RK_opts *) opts_;
    sim_erk_memory *mem = (sim_erk_memory *) mem_;
    sim_erk_workspace *work = (sim_erk_workspace *) work_;

    int i, j, s, istep;
    double a = 0, b =0; // temp values of A_mat and b_vec
    int nx = in->nx;
    int nu = in->nu;

    int NF = opts->num_forw_sens;
    if (!opts->sens_forw)
        NF = 0;

    int nhess = (NF + 1) * NF / 2;
    int nX = nx * (1 + NF);

    double *x = in->x;
    double *u = in->u;
    double *S_forw_in = in->S_forw;
    int num_steps = opts->num_steps;
    double step = in->step;

    double *S_adj_in = in->S_adj;

    double *A_mat = opts->A_mat;
    double *b_vec = opts->b_vec;
    //    double *c_vec = opts->c_vec;
    int num_stages = opts->num_stages;

    double *K_traj = mem->K_traj;
    double *forw_traj = mem->out_forw_traj;
    double *rhs_forw_in = mem->rhs_forw_in;

    double *adj_tmp = mem->out_adj_tmp;
    double *adj_traj = mem->adj_traj;
    double *rhs_adj_in = mem->rhs_adj_in;

    double *xn = out->xn;
    double *S_forw_out = out->S_forw;
    double *S_adj_out = out->S_adj;
    double *S_hess_out = out->S_hess;

    acados_timer timer, timer_ad;
    double timing_ad = 0.0;

    acados_tic(&timer);
    for (i = 0; i < nx; i++)
        forw_traj[i] = x[i];  // x0
    if (opts->sens_forw) {
        for (i = 0; i < nx * NF; i++)
            forw_traj[nx + i] = S_forw_in[i];  // sensitivities
    }

    for (i = 0; i < nu; i++)
        rhs_forw_in[nX + i] = u[i]; // controls

    // FORWARD SWEEP:
    for (istep = 0; istep < num_steps; istep++) {
        if (opts->sens_adj) {
            K_traj = mem->K_traj + istep * num_stages * nX;
            forw_traj = mem->out_forw_traj + (istep + 1) * nX;
            for (i = 0; i < nX; i++)
                forw_traj[i] = forw_traj[i - nX];
        }

        for (s = 0; s < num_stages; s++) {
            for (i = 0; i < nX; i++)
                rhs_forw_in[i] = forw_traj[i];
            for (j = 0; j < s; j++){
                a = A_mat[j * num_stages + s];
                if (a!=0){
                    a *= step;
                    for (i = 0; i < nX; i++)
                        rhs_forw_in[i] += a * K_traj[j * nX + i];
                }
            }

            acados_tic(&timer_ad);
            in->forward_vde_wrapper(nx, nu, rhs_forw_in, K_traj+s*nX, in->vde);  // k evaluation
            timing_ad += acados_toc(&timer_ad)*1000;
        }
        for (s = 0; s < num_stages; s++){
            b = step * b_vec[s];
            for (i = 0; i < nX; i++)
                forw_traj[i] += b * K_traj[s * nX + i];  // ERK step
        }
    }

    for (i = 0; i < nx; i++)
        xn[i] = forw_traj[i];
    if (opts->sens_forw) {
        for (i = 0; i < nx * NF; i++)
            S_forw_out[i] = forw_traj[nx + i];
    }

    // ADJOINT SWEEP:
    if (opts->sens_adj) {
        for (i = 0; i < nx; i++)
            adj_tmp[i] = S_adj_in[i];
        for (i = 0; i < nu; i++)
            adj_tmp[nx+i] = 0.0;

        int nForw = nx;
        int nAdj = nx + nu;
        if (opts->sens_hess) {
            nForw = nX;
            nAdj = nx + nu + nhess;
            for (i = 0; i < nhess; i++)
                adj_tmp[nx + nu + i] = 0.0;
        }

        printf("\nnFOrw=%d nAdj=%d\n", nForw, nAdj);

        for (i = 0; i < nu; i++)
            rhs_adj_in[nForw + nx + i] = u[i];

        for (istep = num_steps - 1; istep > -1; istep--) {
            K_traj = mem->K_traj + istep * num_stages * nX;
            forw_traj = mem->out_forw_traj + istep * nX;
            for (s = num_stages - 1; s > -1; s--) {
                // forward variables:
                for (i = 0; i < nForw; i++)
                    rhs_adj_in[i] = forw_traj[i]; // extract x trajectory
                for (j = 0; j < s; j++)
                    a = A_mat[j * num_stages + s];
                    if (a!=0){
                        a*=step;
                        for (i = 0; i < nForw; i++)
                            rhs_adj_in[i] += a *K_traj[j * nX + i];
                    } // plus k traj
                // adjoint variables:
                b = step * b_vec[s];
                for (i = 0; i < nx; i++)
                    rhs_adj_in[nForw + i] = b * adj_tmp[i];
                for (j = s + 1; j < num_stages; j++){
                    a = A_mat[s * num_stages + j];
                    if (a!=0){
                        a *= step;
                        for (i = 0; i < nx; i++)
                            rhs_adj_in[nForw + i] += a * adj_traj[j * nAdj + i];
                    }
                }
                acados_tic(&timer_ad);
                if (opts->sens_hess){
                    in->Hess_fun(nx, nu, rhs_adj_in, adj_traj+s*nAdj, in->hess);
                }else{
                    in->adjoint_vde_wrapper(nx, nu, rhs_adj_in, adj_traj+s*nAdj, in->vde_adj); // adjoint VDE evaluation
                }
                timing_ad += acados_toc(&timer_ad);

                // printf("\nadj_traj:\n");
                // for (int ii=0;ii<num_stages*nAdj;ii++)
                //     printf("%3.1f ", adj_traj[ii]);
            }
            for (s = 0; s < num_stages; s++)
                for (i = 0; i < nAdj; i++)
                    adj_tmp[i] += adj_traj[s * nAdj + i];  // ERK step
        }
        for (i = 0; i < nx + nu; i++)
            S_adj_out[i] = adj_tmp[i];
        if (opts->sens_hess) {
            for (i = 0; i < nhess; i++)
                S_hess_out[i] = adj_tmp[nx + nu + i];
        }
    }
    out->info->CPUtime = acados_toc(&timer)*1000;
    out->info->LAtime = 0.0;
    out->info->ADtime = timing_ad;
    return 0;  // success
}
