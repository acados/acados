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

#include "acados/sim/sim_erk_integrator.h"

#include <stdlib.h>
#include <string.h>

#include "hpmpc/include/aux_d.h"

#include "acados/utils/print.h"


static void sim_erk_cast_workspace(sim_erk_workspace *work,
         const sim_in *in, void *args) {
     int_t nx = in->nx;
     int_t nu = in->nu;
     sim_RK_opts *opts = (sim_RK_opts*) args;
     int_t num_stages = opts->num_stages;
     int_t NF = in->nsens_forw;
     int_t nSteps = in->nSteps;
     if (!in->sens_forw) {
         NF = 0;
     }
     int_t nhess = (int_t)(NF+1)*(real_t)NF/2.0;

     char *ptr = (char *)work;
     ptr += sizeof(sim_erk_workspace);
     work->rhs_forw_in = (real_t*)ptr;
     ptr += (nx*(1+NF)+nu)*sizeof(real_t);  // rhs_forw_in

     if (!in->sens_adj) {
         work->K_traj = (real_t*)ptr;
         ptr += (num_stages*nx*(1+NF))*sizeof(real_t);  // K_traj
         work->out_forw_traj = (real_t*)ptr;
         ptr += (nx*(1+NF))*sizeof(real_t);  // out_forw_traj
     } else {
         work->K_traj = (real_t*)ptr;
         ptr += (nSteps*num_stages*nx*(1+NF))*sizeof(real_t);  // K_traj
         work->out_forw_traj = (real_t*)ptr;
         ptr += ((nSteps+1)*nx*(1+NF))*sizeof(real_t);  // out_forw_traj
     }

     if (in->sens_hess && in->sens_adj) {
         work->rhs_adj_in = (real_t*)ptr;
         ptr += (nx*(2+NF)+nu)*sizeof(real_t);  // rhs_adj_in
         work->out_adj_tmp = (real_t*)ptr;
         ptr += (nx+nu+nhess)*sizeof(real_t);  // out_adj_tmp
         work->adj_traj = (real_t*)ptr;
         ptr += (num_stages*(nx+nu+nhess))*sizeof(real_t);  // adj_traj
     } else if (in->sens_adj) {
         work->rhs_adj_in = (real_t*)ptr;
         ptr += (nx*2+nu)*sizeof(real_t);  // rhs_adj_in
         work->out_adj_tmp = (real_t*)ptr;
         ptr += (nx+nu)*sizeof(real_t);  // out_adj_tmp
         work->adj_traj = (real_t*)ptr;
         ptr += (num_stages*(nx+nu))*sizeof(real_t);  // adj_traj
     }
}

int_t sim_erk(const sim_in *in, sim_out *out, void *args, void *mem, void *work_) {
    int_t nx = in->nx;
    int_t nu = in->nu;
    sim_RK_opts *opts = (sim_RK_opts*) args;
    int_t num_stages = opts->num_stages;
    int_t i, s, j, istep;
    real_t H_INT = in->step;
    int_t NSTEPS = in->nSteps;
    int_t NF = in->nsens_forw;
    sim_erk_workspace *work = (sim_erk_workspace*) work_;
    sim_erk_cast_workspace(work, in, args);
    if (!in->sens_forw) {
        NF = 0;
    }
    int_t nhess = (int_t)(NF+1)*(real_t)NF/2.0;

    mem = 0; mem += 0;

    real_t *A_mat = opts->A_mat;
    real_t *b_vec = opts->b_vec;
//    real_t *c_vec = opts->c_vec;

//    print_matrix("stdout", A_mat, num_stages, num_stages);

    real_t *K_traj = work->K_traj;
    real_t *forw_traj = work->out_forw_traj;
    real_t *rhs_forw_in = work->rhs_forw_in;

    real_t *adj_tmp = work->out_adj_tmp;
    real_t *adj_traj = work->adj_traj;
    real_t *rhs_adj_in = work->rhs_adj_in;

#ifdef MEASURE_TIMINGS
    acado_timer timer, timer_ad;
    real_t timing_ad = 0.0;
    acado_tic(&timer);
#endif
    for (i = 0; i < nx; i++) forw_traj[i] = in->x[i];
    if (in->sens_forw) {
        for (i = 0; i < nx*NF; i++) forw_traj[nx+i] = in->S_forw[i];  // sensitivities
    }

    for (i = 0; i < nu; i++) rhs_forw_in[nx*(1+NF)+i] = in->u[i];

    // FORWARD SWEEP:
    for (istep = 0; istep < NSTEPS; istep++) {
        if (in->sens_adj) {
            K_traj = &work->K_traj[istep*num_stages*nx*(1+NF)];
            forw_traj = &work->out_forw_traj[(istep+1)*nx*(1+NF)];
            for (i = 0; i < nx*(1+NF); i++) {
                forw_traj[i] = forw_traj[i-nx*(1+NF)];
            }
        }

        for (s = 0; s < num_stages; s++) {
            for (i = 0; i < nx*(1+NF); i++) {
                rhs_forw_in[i] = forw_traj[i];
            }
            for (j = 0; j < s; j++) {
                if (A_mat[j*num_stages+s] != 0) {
                    for (i = 0; i < nx*(1+NF); i++) {
                        rhs_forw_in[i] += H_INT*A_mat[j*num_stages+s]*K_traj[j*nx*(1+NF)+i];
                    }
                }
            }
#ifdef MEASURE_TIMINGS
            acado_tic(&timer_ad);
#endif
            in->VDE_forw(rhs_forw_in, &(K_traj[s*nx*(1+NF)]), in->vde);  // k evaluation
#ifdef MEASURE_TIMINGS
            timing_ad += acado_toc(&timer_ad);
#endif
        }
        for (s = 0; s < num_stages; s++) {
            for (i = 0; i < nx*(1+NF); i++) {
                forw_traj[i] += H_INT*b_vec[s]*K_traj[s*nx*(1+NF)+i];  // ERK step
            }
        }
    }
    for (i = 0; i < nx; i++)    out->xn[i] = forw_traj[i];
    if (in->sens_forw) {
        for (i = 0; i < nx*NF; i++) out->S_forw[i] = forw_traj[nx+i];
    }

    // ADJOINT SWEEP:
    if (in->sens_adj) {
        for (i = 0; i < nx+nu; i++) adj_tmp[i] = in->S_adj[i];

        int_t nForw = nx;
        int_t nAdj = nx+nu;
        if (in->sens_hess) {
            nForw = nx*(1+NF);
            nAdj = nx+nu+nhess;
            for (i = 0; i < nhess; i++) adj_tmp[nx+nu+i] = 0.0;
        }
        for (i = 0; i < nu; i++) rhs_adj_in[nForw+nx+i] = in->u[i];

        for (istep = NSTEPS-1; istep > -1; istep--) {
            K_traj = &work->K_traj[istep*num_stages*nx*(1+NF)];
            forw_traj = &work->out_forw_traj[istep*nx*(1+NF)];

            for (s = num_stages-1; s > -1; s--) {
                // forward variables:
                for (i = 0; i < nForw; i++) {
                    rhs_adj_in[i] = forw_traj[i];
                }
                for (j = 0; j < s; j++) {
                    if (A_mat[j*num_stages+s] != 0) {
                        for (i = 0; i < nForw; i++) {
                            rhs_adj_in[i] += H_INT*A_mat[j*num_stages+s]*K_traj[j*nx*(1+NF)+i];
                        }
                    }
                }
                // adjoint variables:
                for (i = 0; i < nx; i++) {
                    rhs_adj_in[nForw+i] = H_INT*b_vec[s]*adj_tmp[i];
                }
                for (j = s+1; j < num_stages; j++) {
                    if (A_mat[s*num_stages+j] != 0) {
                        for (i = 0; i < nx; i++) {
                            rhs_adj_in[nForw+i] += H_INT*A_mat[s*num_stages+j]*adj_traj[j*nAdj+i];
                        }
                    }
                }
#ifdef MEASURE_TIMINGS
                acado_tic(&timer_ad);
#endif
                in->VDE_adj(rhs_adj_in, &(adj_traj[s*nAdj]));  // adjoint VDE evaluation
#ifdef MEASURE_TIMINGS
                timing_ad += acado_toc(&timer_ad);
#endif
            }
            for (s = 0; s < num_stages; s++) {
                for (i = 0; i < nAdj; i++) {
                    adj_tmp[i] += adj_traj[s*nAdj+i];  // ERK step
                }
            }
        }
        for (i = 0; i < nx+nu; i++) out->S_adj[i] = adj_tmp[i];
        if (in->sens_hess) {
            for (i = 0; i < nhess; i++) out->S_hess[i] = adj_tmp[nx+nu+i];
        }
    }
#ifdef MEASURE_TIMINGS
    out->info->CPUtime = acado_toc(&timer);
    out->info->LAtime = 0.0;
    out->info->ADtime = timing_ad;
#endif
    return 0;  // success
}

int_t sim_erk_calculate_workspace_size(const sim_in *in, void *args) {
    int_t nx = in->nx;
    int_t nu = in->nu;
    sim_RK_opts *opts = (sim_RK_opts*) args;
    int_t num_stages = opts->num_stages;
    int_t NF = in->nsens_forw;
    int_t nSteps = in->nSteps;
    if (!in->sens_forw) {
        NF = 0;
    }
    int_t nhess = (int_t)(NF+1)*(real_t)NF/2.0;

    int_t size = sizeof(sim_erk_workspace);
    size += (nx*(1+NF)+nu)*sizeof(real_t);  // rhs_forw_in

    if (!in->sens_adj) {
        size += (num_stages*nx*(1+NF))*sizeof(real_t);  // K_traj
        size += (nx*(1+NF))*sizeof(real_t);  // out_forw_traj
    } else {
        size += (nSteps*num_stages*nx*(1+NF))*sizeof(real_t);  // K_traj
        size += ((nSteps+1)*nx*(1+NF))*sizeof(real_t);  // out_forw_traj
    }

    if (in->sens_hess && in->sens_adj) {
        size += (nx*(2+NF)+nu)*sizeof(real_t);  // rhs_adj_in
        size += (nx+nu+nhess)*sizeof(real_t);  // out_adj_tmp
        size += (num_stages*(nx+nu+nhess))*sizeof(real_t);  // adj_traj
    } else if (in->sens_adj) {
        size += (nx*2+nu)*sizeof(real_t);  // rhs_adj_in
        size += (nx+nu)*sizeof(real_t);  // out_adj_tmp
        size += (num_stages*(nx+nu))*sizeof(real_t);  // adj_traj
    }
    return size;
}


void sim_erk_create_arguments(void *args, const int_t num_stages) {
    sim_RK_opts *opts = (sim_RK_opts*) args;
    opts->scheme.type = exact;
    if ( num_stages == 1 ) {
        opts->num_stages = 1;       // explicit Euler
        opts->A_mat = calloc(num_stages*num_stages, sizeof(*opts->A_mat));
        opts->b_vec = calloc(num_stages, sizeof(*opts->b_vec));
        opts->c_vec = calloc(num_stages, sizeof(*opts->c_vec));
        opts->A_mat[0] = 0;
        opts->b_vec[0] = 1.0;
        opts->c_vec[0] = 0;
    } else if ( num_stages == 4 ) {
        opts->num_stages = 4;       // 4-stage, 4th order ERK method
        opts->A_mat = calloc(num_stages*num_stages, sizeof(*opts->A_mat));
        opts->b_vec = calloc(num_stages, sizeof(*opts->b_vec));
        opts->c_vec = calloc(num_stages, sizeof(*opts->c_vec));

        memcpy(opts->A_mat,
                ((real_t[]) {0, 0.5, 0, 0, 0, 0, 0.5, 0, 0, 0, 0, 1, 0, 0, 0, 0}),
                sizeof(*opts->A_mat) * (num_stages*num_stages));
        memcpy(opts->b_vec,
                ((real_t[]) {1.0/6, 2.0/6, 2.0/6, 1.0/6}), sizeof(*opts->b_vec) * (num_stages));
        memcpy(opts->c_vec,
                ((real_t[]) {0.0, 0.5, 0.5, 1.0}), sizeof(*opts->c_vec) * (num_stages));
    } else {
        // throw error somehow?
    }
}


void sim_erk_initialize(const sim_in *in, void *args, void **work) {
    sim_RK_opts *opts = (sim_RK_opts*) args;

    // TODO(dimitris): opts should be an input to initialize
    if (opts->num_stages > 0) {
        sim_erk_create_arguments(args, opts->num_stages);
    } else {
        sim_erk_create_arguments(args, 4);
    }
    int_t work_space_size = sim_erk_calculate_workspace_size(in, args);
    *work = (void *) malloc(work_space_size);
}


void sim_erk_destroy(void *work) {
    free(work);
}
