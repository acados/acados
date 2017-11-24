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
#include "acados/utils/print.h"
#include "acados/sim/sim_rk_common_yt.h"
#include "acados/sim/sim_common_yt.h"
#include "acados/sim/sim_erk_integrator_yt.h"

#include "acados/sim/sim_casadi_wrapper.h"

int_t erk_calculate_memory_size(sim_RK_opts *opts, sim_in *in)
{

int_t nx = in->nx;
int_t nu = in->nu; 
int_t NF = in->NF;
 
int_t num_stages = opts->num_stages; // number of stages
int_t nX = nx*(1+NF); // (nx) for ODE and (NF*nx) for VDE
int_t nhess = (NF + 1) * NF / 2;
uint num_steps = in->num_steps;  // number of steps

int_t size = sizeof(sim_erk_memory);

size += (nX + nu) * sizeof(real_t); // rhs_forw_in

if(in->sens_adj){
    size += num_steps * num_stages * nX * sizeof(real_t); // K_traj
    size += (num_steps + 1) * nX *sizeof(real_t); // out_forw_traj
}else{
    size += num_stages * nX * sizeof(real_t); // K_traj
    size += nX *sizeof(real_t); // out_forw_traj
}

if (in->sens_hess && in->sens_adj){
    size += (nx + nX + nu) * sizeof(real_t); //rhs_adj_in
    size += (nx + nu + nhess) * sizeof(real_t); //out_adj_tmp
    size += num_stages * (nx + nu + nhess) * sizeof(real_t); //adj_traj
}else if (in->sens_adj){
    size += (nx * 2 + nu) * sizeof(real_t); //rhs_adj_in
    size += (nx + nu)* sizeof(real_t); //out_adj_tmp
    size += num_stages * (nx + nu) * sizeof(real_t); //adj_traj
}

size = (size + 63) / 64 * 64;
size += 1 * 64;

return size;
}

char *assign_erk_memory(sim_RK_opts *opts, sim_in *in, sim_erk_memory **memory, void *raw_memory)
{

int_t nx = in->nx;
int_t nu = in->nu; 
int_t NF = in->NF;

int_t num_stages = opts->num_stages; // number of stages
int_t nX = nx*(1+NF); // (nx) for ODE and (NF*nx) for VDE
int_t nhess = (NF + 1) * NF / 2;
int_t num_steps = in->num_steps;  // number of steps

char *c_ptr = (char *)raw_memory;

*memory = (sim_erk_memory *) c_ptr;
c_ptr += sizeof(sim_erk_memory);

// align memory to typical cache line size
size_t s_ptr = (size_t)c_ptr;
s_ptr = (s_ptr + 63) / 64 * 64;
c_ptr = (char *)s_ptr;

(*memory)->rhs_forw_in = (real_t *)c_ptr;
c_ptr += (nX + nu) * sizeof(real_t);

if(in->sens_adj){
    (*memory)->K_traj = (real_t *)c_ptr;
    c_ptr += num_steps * num_stages * nX * sizeof(real_t);
    (*memory)->out_forw_traj = (real_t *)c_ptr;
    c_ptr += (num_steps + 1) * nX * sizeof(real_t);
}else{
    (*memory)->K_traj = (real_t *)c_ptr;
    c_ptr += num_stages * nX * sizeof(real_t);
    (*memory)->out_forw_traj = (real_t *)c_ptr;
    c_ptr += nX * sizeof(real_t);
}

if (in->sens_hess && in->sens_adj){
    (*memory)->rhs_adj_in = (real_t *)c_ptr;
    c_ptr += (nx + nX + nu) * sizeof(real_t);
    (*memory)->out_adj_tmp = (real_t *)c_ptr;
    c_ptr += (nx + nu + nhess) * sizeof(real_t); 
    (*memory)->adj_traj = (real_t *)c_ptr;
    c_ptr += num_stages * (nx + nu + nhess) * sizeof(real_t); 
    }else if (in->sens_adj){
        (*memory)->rhs_adj_in = (real_t *)c_ptr;
        c_ptr += (nx * 2 + nu) * sizeof(real_t);
        (*memory)->out_adj_tmp = (real_t *)c_ptr;
        c_ptr += (nx + nu) * sizeof(real_t);
        (*memory)->adj_traj = (real_t *)c_ptr;
        c_ptr += num_stages * (nx + nu) * sizeof(real_t);
    }

    return c_ptr;
}

sim_erk_memory *sim_erk_create_memory(sim_RK_opts *opts, sim_in *in)
{
    sim_erk_memory *memory;

    int_t bytes = erk_calculate_memory_size(opts, in);
    void *ptr = malloc(bytes);
    char *ptr_end = assign_erk_memory(opts, in, &memory, ptr);
    assert((char*)ptr + bytes >= ptr_end); (void) ptr_end;

    return memory;
}

int_t sim_erk_yt(const sim_in *in, sim_out *out, void *opts_, void *mem_) {

    sim_RK_opts *opts = (sim_RK_opts *) opts_;
    sim_erk_memory *mem = (sim_erk_memory *) mem_;

    int_t i, j, s, istep;
    real_t a = 0, b =0; // temp values of A_mat and b_vec
    int_t nx = in->nx;
    int_t nu = in->nu;

    int_t NF = in->NF;
    if (!in->sens_forw)
        NF = 0;
    
    int_t nhess = (NF + 1) * NF / 2;
    int_t nX = nx * (1 + NF);

    real_t *x = in->x;
    real_t *u = in->u;
    real_t *S_forw_in = in->S_forw; 
    int_t num_steps = in->num_steps;
    real_t step = in->step;

    real_t *S_adj_in = in->S_adj;

    real_t *A_mat = opts->A_mat;
    real_t *b_vec = opts->b_vec;
    //    real_t *c_vec = opts->c_vec;
    int_t num_stages = opts->num_stages;

    real_t *K_traj = mem->K_traj;
    real_t *forw_traj = mem->out_forw_traj;
    real_t *rhs_forw_in = mem->rhs_forw_in;

    real_t *adj_tmp = mem->out_adj_tmp;
    real_t *adj_traj = mem->adj_traj;
    real_t *rhs_adj_in = mem->rhs_adj_in;

    real_t *xn = out->xn;
    real_t *S_forw_out = out->S_forw;
    real_t *S_adj_out = out->S_adj;
    real_t *S_hess_out = out->S_hess;

    acados_timer timer, timer_ad;
    double timing_ad = 0.0;

    acados_tic(&timer);
    for (i = 0; i < nx; i++)
        forw_traj[i] = x[i];  // x0
    if (in->sens_forw) {
        for (i = 0; i < nx * NF; i++)
            forw_traj[nx + i] = S_forw_in[i];  // sensitivities
    }

    for (i = 0; i < nu; i++)
        rhs_forw_in[nX + i] = u[i]; // controls

    // FORWARD SWEEP:
    for (istep = 0; istep < num_steps; istep++) {
        if (in->sens_adj) {
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
            in->VDE_forw(nx, nu, rhs_forw_in, K_traj+s*nX, in->vde);  // k evaluation
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
    if (in->sens_forw) {
        for (i = 0; i < nx * NF; i++)
            S_forw_out[i] = forw_traj[nx + i];
    }

    // ADJOINT SWEEP:
    if (in->sens_adj) {
        for (i = 0; i < nx; i++)
            adj_tmp[i] = S_adj_in[i];
        for (i = 0; i < nu; i++)
            adj_tmp[nx+i] = 0.0;

        int nForw = nx;
        int nAdj = nx + nu;
        if (in->sens_hess) {
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
                if (in->sens_hess){
                    in->Hess_fun(nx, nu, rhs_adj_in, adj_traj+s*nAdj, in->hess);                   
                }else{
                    in->VDE_adj(nx, nu, rhs_adj_in, adj_traj+s*nAdj, in->adj); // adjoint VDE evaluation
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
        if (in->sens_hess) {
            for (i = 0; i < nhess; i++)
                S_hess_out[i] = adj_tmp[nx + nu + i];
        }
    }
    out->info->CPUtime = acados_toc(&timer)*1000;
    out->info->LAtime = 0.0;
    out->info->ADtime = timing_ad;
    return 0;  // success
}
