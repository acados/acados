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

#include <stdlib.h>
#include <string.h>

#include "hpmpc/include/aux_d.h"
#include "acados/sim/sim_erk_integrator.h"
#include "acados/utils/print.h"

void sim_erk(const sim_in *in, sim_out *out, const sim_RK_opts *opts, sim_erk_workspace *work) {
    int_t nx = in->nx;
    int_t nu = in->nu;
    int_t num_stages = opts->num_stages;
    int_t i, s, j, istep;
    real_t H_INT = in->step;
    int_t NSTEPS = in->nSteps;
    int_t NF = in->nsens_forw;
    if (!in->sens_forw) {
        NF = 0;
    }
    int_t nhess = (int_t)(NF+1)*(real_t)NF/2.0;

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

    acado_timer timer, timer_ad;
    acado_tic(&timer);
    real_t timing_ad = 0.0;

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
            acado_tic(&timer_ad);
            in->VDE_forw(rhs_forw_in, &(K_traj[s*nx*(1+NF)]));  // k evaluation
            timing_ad += acado_toc(&timer_ad);
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

                acado_tic(&timer_ad);
                in->VDE_adj(rhs_adj_in, &(adj_traj[s*nAdj]));  // adjoint VDE evaluation
                timing_ad += acado_toc(&timer_ad);
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

    out->info->CPUtime = acado_toc(&timer);
    out->info->LAtime = 0.0;
    out->info->ADtime = timing_ad;
}


void sim_erk_create_workspace(const sim_in *in, sim_RK_opts *opts, sim_erk_workspace *work) {
    int_t nx = in->nx;
    int_t nu = in->nu;
    int_t num_stages = opts->num_stages;
    int_t NF = in->nsens_forw;
    int_t nSteps = in->nSteps;
    if (!in->sens_forw) {
        NF = 0;
    }
    int_t nhess = (int_t)(NF+1)*(real_t)NF/2.0;


    work->rhs_forw_in = malloc(sizeof(*work->rhs_forw_in) * (nx*(1+NF)+nu));
    if (!in->sens_adj) {
        work->K_traj = malloc(sizeof(*work->K_traj) * (num_stages*nx*(1+NF)));
        work->out_forw_traj = malloc(sizeof(*work->out_forw_traj) * (nx*(1+NF)));
    } else {
        work->K_traj = malloc(sizeof(*work->K_traj) * (nSteps*num_stages*nx*(1+NF)));
        work->out_forw_traj = malloc(sizeof(*work->out_forw_traj) * ((nSteps+1)*nx*(1+NF)));
    }

    if (in->sens_adj) {
    }
    if (in->sens_hess && in->sens_adj) {
        work->rhs_adj_in = malloc(sizeof(*work->rhs_adj_in) * (nx*(2+NF)+nu));
        work->out_adj_tmp = malloc(sizeof(*work->out_adj_tmp) * (nx+nu+nhess));
        work->adj_traj = malloc(sizeof(*work->adj_traj) * (num_stages*(nx+nu+nhess)));
    } else if (in->sens_adj) {
        work->rhs_adj_in = malloc(sizeof(*work->rhs_adj_in) * (nx*2+nu));
        work->out_adj_tmp = malloc(sizeof(*work->out_adj_tmp) * (nx+nu));
        work->adj_traj = malloc(sizeof(*work->adj_traj) * (num_stages*(nx+nu)));
    }
}


void sim_erk_create_opts(const int_t num_stages, sim_RK_opts *opts) {
    if ( num_stages == 1 ) {
        opts->num_stages = 1;       // explicit Euler
        opts->A_mat = malloc(sizeof(*opts->A_mat) * (num_stages*num_stages));
        opts->b_vec = malloc(sizeof(*opts->b_vec) * (num_stages));
        opts->c_vec = malloc(sizeof(*opts->c_vec) * (num_stages));
        opts->A_mat[0] = 0;
        opts->b_vec[0] = 1.0;
        opts->c_vec[0] = 0;
    } else if ( num_stages == 4 ) {
        opts->num_stages = 4;       // 4-stage, 4th order ERK method
        opts->A_mat = malloc(sizeof(*opts->A_mat) * (num_stages*num_stages));
        opts->b_vec = malloc(sizeof(*opts->b_vec) * (num_stages));
        opts->c_vec = malloc(sizeof(*opts->c_vec) * (num_stages));

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
