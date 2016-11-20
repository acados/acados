#include "test/sim/pendulum/casadi/casadi_pendulum.h"
#include "acados/utils/types.h"
#include "acados/utils/timing.h"
#include "acados/utils/print.h"
#include "acados/sim/sim_erk_integrator.h"
#include "test/test_utils/eigen.h"

void create_ERK_integrator(sim_in* sim_in, sim_out* sim_out,
        sim_info* info, sim_RK_opts* rk_opts, sim_erk_workspace* erk_work,
        int_t NX, int_t NU, real_t T, bool hessian) {
        sim_in->nSteps = 10;
        sim_in->step = T/sim_in->nSteps;
        sim_in->nx = NX;
        sim_in->nu = NU;

        int_t nhess = (int_t)(NX+NU+1)*(real_t)(NX+NU)/2.0;

        sim_in->sens_forw = true;
        sim_in->sens_adj = true;
        sim_in->sens_hess = hessian;
        sim_in->nsens_forw = NX+NU;

        sim_in->VDE_forw = &VDE_forw_pendulum;
        if ( hessian ) {
            sim_in->VDE_adj = &VDE_hess_pendulum;
        } else {
            sim_in->VDE_adj = &VDE_adj_pendulum;
        }
        sim_in->jac_fun = &jac_fun_pendulum;

        sim_in->x = (real_t*) malloc(sizeof(*sim_in->x) * (NX));
        sim_in->u = (real_t*) malloc(sizeof(*sim_in->u) * (NU));
        sim_in->S_forw = (real_t*) malloc(sizeof(*sim_in->S_forw) * (NX*(NX+NU)));
        sim_in->S_adj = (real_t*) malloc(sizeof(*sim_in->S_adj) * (NX+NU));
        for (int_t i = 0; i < NX*(NX+NU); i++) sim_in->S_forw[i] = 0.0;
        for (int_t i = 0; i < NX; i++) sim_in->S_forw[i*(NX+1)] = 1.0;

        sim_out->xn = (real_t*) malloc(sizeof(*sim_out->xn) * (NX));
        sim_out->S_forw = (real_t*) malloc(sizeof(*sim_out->S_forw) * (NX*(NX+NU)));
        sim_out->S_adj = (real_t*) malloc(sizeof(*sim_out->S_adj) * (NX+NU));
        if ( hessian ) {
            sim_out->S_hess = (real_t*) malloc(sizeof(*sim_out->S_hess) * (nhess));
        }
        sim_out->info = info;

        sim_erk_create_opts(4, rk_opts);
        sim_erk_create_workspace(sim_in, rk_opts, erk_work);
}
