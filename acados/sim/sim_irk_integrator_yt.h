#ifndef ACADOS_SIM_SIM_IRK_INTEGRATOR_YT_H_
#define ACADOS_SIM_SIM_IRK_INTEGRATOR_YT_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "acados/sim/sim_rk_common_yt.h"
#include "acados/sim/sim_common_yt.h"
#include "acados/utils/types.h"

#include "external/blasfeo/include/blasfeo_target.h"
#include "external/blasfeo/include/blasfeo_common.h"

typedef struct {
    struct d_strmat *JG; // jacobian of G over K (nx*ns, nx*ns)
    struct d_strmat *JGf; // jacobian of G over x and u (nx*ns, nx+nu);
    struct d_strmat *JKf; // jacobian of K over x and u (nx*ns, nx+nu);
    struct d_strmat *JFK; // jacobian of F over K (nx, nx*ns) 
    struct d_strmat *S_forw; // forward sensitivities

    struct d_strvec *rG; // residuals of G (nx*ns)
    struct d_strvec *K; // internal variables (nx*ns)
    struct d_strvec *xt; // temporary x
    struct d_strvec *xn; // x at each integration step
    
    struct d_strvec *lambda; // adjoint seed (nx+nu)
    struct d_strvec *lambdaK; // auxiliary variable (nx*ns)
    
    double *rGt; // temporary residuals of G (nx, 1)
    double *jac_out; // temporary Jacobian of ode (nx, 2*nx+nu)
    double *Jt; // temporary Jacobian of ode (nx, nx)
    double *ode_args; // pointer to ode args
    double *S_adj_w;
    int *ipiv; // index of pivot vector// jacobian of G over x and u (nx*ns, nx+nu);

    struct d_strvec **xn_traj; // xn trajectory
    struct d_strvec **K_traj;  // K trajectory
}sim_irk_memory;

typedef struct {
    double dummy;
} sim_erk_workspace;

int irk_calculate_memory_size(sim_dims *dims, void *opts_);

char *assign_irk_memory(sim_dims *dims, void *opts_, void *raw_memory);

void *sim_irk_create_memory(sim_dims *dims, void *opts_);

int sim_irk_yt(sim_in *in, sim_out *out, void *opts_, void *mem_, void *work_);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // ACADOS_SIM_SIM_IRK_INTEGRATOR_H_