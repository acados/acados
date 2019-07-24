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
#include <stdio.h>
#include <stdlib.h>
// acados
#include "acados/utils/print.h"
#include "acados_c/ocp_nlp_interface.h"
#include "acados_c/external_function_interface.h"
#include "acados/ocp_nlp/ocp_nlp_cost_ls.h"

// blasfeo
#include "blasfeo/include/blasfeo_d_aux.h"
#include "blasfeo/include/blasfeo_d_aux_ext_dep.h"

// example specific
#include "{{ ocp.model.name }}_model/{{ ocp.model.name }}_model.h"
{% if ocp.dims.npd > 0 %}
#include "{{ ocp.con_p.name }}_p_constraint/{{ ocp.con_p.name }}_p_constraint.h"
{% endif %}
{% if ocp.dims.nh > 0 %}
#include "{{ ocp.con_h.name }}_h_constraint/{{ ocp.con_h.name }}_h_constraint.h"
{% endif %}

#include "acados_solver_{{ ocp.model.name }}.h"

{%- for key, value in ocp.constants %}
#define {{ key }} {{ value }}
{%- endfor %}
#define NX_    {{ ocp.dims.nx }}
#define NZ_    {{ ocp.dims.nz }}
#define NU_    {{ ocp.dims.nu }}
#define NP_    {{ ocp.dims.np }}
#define NBX_   {{ ocp.dims.nbx }}
#define NBU_   {{ ocp.dims.nbu }}
#define NSBX_  {{ ocp.dims.nsbx }}
#define NSBU_  {{ ocp.dims.nsbu }}
#define NSH_  {{ ocp.dims.nsh }}
#define NSHN_  {{ ocp.dims.nsh_e }}
#define NSBXN_ {{ ocp.dims.nsbx_e }}
#define NS_    {{ ocp.dims.ns }}
#define NSN_   {{ ocp.dims.ns_e }}
#define NG_    {{ ocp.dims.ng }}
#define NBXN_  {{ ocp.dims.nbx_e }}
#define NGN_   {{ ocp.dims.ng_e }}
#define NY_    {{ ocp.dims.ny }}
#define NYN_   {{ ocp.dims.ny_e }}
#define N_     {{ ocp.dims.N }}
#define NPD_   {{ ocp.dims.npd }}
#define NPDN_  {{ ocp.dims.npd_e }}
#define NH_    {{ ocp.dims.nh }}
#define NHN_   {{ ocp.dims.nh_e }}

#if NX_ < 1
#define NX   1
#else
#define NX   NX_
#endif

#if NZ_ < 1
#define NZ   1
#else
#define NZ   NZ_
#endif

#if NU_ < 1
#define NU   1
#else
#define NU   NU_
#endif

#if NP_ < 1
#define NP   1
#else
#define NP   NP_
#endif

#if NBX_ < 1
#define NBX   1
#else
#define NBX   NBX_
#endif

#if NBU_ < 1
#define NBU   1
#else
#define NBU   NBU_
#endif

#if NSBX_ < 1
#define NSBX   1
#else
#define NSBX   NSBX_
#endif

#if NSBU_ < 1
#define NSBU   1
#else
#define NSBU   NSBU_
#endif

#if NSH_ < 1
#define NSH   1
#else
#define NSH   NSH_
#endif

#if NSHN_ < 1
#define NSHN   1
#else
#define NSHN   NSHN_
#endif

#if NS_ < 1
#define NS   1
#else
#define NS   NS_
#endif

#if NSBXN_ < 1
#define NSBXN   1
#else
#define NSBXN   NSBXN_
#endif

#if NSN_ < 1
#define NSN   1
#else
#define NSN   NSN_
#endif

#if NG_ < 1
#define NG   1
#else
#define NG   NG_
#endif

#if NBXN_ < 1
#define NBXN   1
#else
#define NBXN   NBXN_
#endif

#if NGN_ < 1
#define NGN   1
#else
#define NGN  NGN_
#endif

#if NY_ < 1
#define NY   1
#else
#define NY   NY_
#endif

#if NYN_ < 1
#define NYN   1
#else
#define NYN   NYN_
#endif

#if N_ < 1
#define N   1
#else
#define N   N_
#endif

#if NPD_ < 1
#define NPD   1
#else
#define NPD   NPD_
#endif

#if NPDN_ < 1
#define NPDN   1
#else
#define NPDN   NPDN_
#endif

#if NH_ < 1
#define NH   1
#else
#define NH   NH_
#endif

#if NHN_ < 1
#define NHN   1
#else
#define NHN   NHN_
#endif

int acados_create() {

    int status = 0;

    double Tf = {{ ocp.solver_config.tf }};

    // initial state x0
    int idxbx0[NX];
    {%- for i in range(ocp.dims.nx) %}
    idxbx0[{{ i }}] = {{ i }};
    {%- endfor %}

    double lbx0[NX]; 
    {%- for item in ocp.constraints.x0 %}
    lbx0[{{ loop.index0 }}] = {{ item }};
    {%- endfor %}

    double ubx0[NX];
    {%- for item in ocp.constraints.x0 %}
    ubx0[{{ loop.index0 }}] = {{ item }};
    {%- endfor %}

    // set up bounds for u
    int idxbu[NBU];
    {%- for item in ocp.constraints.idxbu %}
    idxbu[{{ loop.index0 }}] = {{ item }};
    {%- endfor %}

    double lbu[NBU]; 
    {%- for item in ocp.constraints.lbu %}
    lbu[{{ loop.index0 }}] = {{ item }};
    {%- endfor %}

    double ubu[NBU];
    {%- for item in ocp.constraints.ubu %}
    ubu[{{ loop.index0 }}] = {{ item }};
    {%- endfor %}
    
    // set up soft bounds for u
    int idxsbu[NSBU];
    {%- for item in ocp.constraints.idxsbu %}
    idxsbu[{{ loop.index0 }}] = {{ item }};
    {%- endfor %}

    double lsbu[NSBU]; 
    {%- for item in ocp.constraints.lsbu %}
    lsbu[{{ loop.index0 }}] = {{ item }};
    {%- endfor %}

    double usbu[NSBU];
    {%- for item in ocp.constraints.usbu %}
    usbu[{{ loop.index0 }}] = {{ item }};
    {%- endfor %}
    
    // set up soft bounds for nonlinear constraints
    int idxsh[NSH];
    {%- for item in ocp.constraints.idxsh %}
    idxsh[{{ loop.index0 }}] = {{ item }};
    {%- endfor %}

    double lsh[NSH]; 
    {%- for item in ocp.constraints.lsh %}
    lsh[{{ loop.index0 }}] = {{ item }};
    {%- endfor %}

    double ush[NSH];
    {%- for item in ocp.constraints.ush %}
    ush[{{ loop.index0 }}] = {{ item }};
    {%- endfor %}
    
    // bounds on x
    int idxbx[NBX];
    {%- for item in ocp.constraints.idxbx %}
    idxbx[{{ loop.index0 }}] = {{ item }};
    {%- endfor %}

    double lbx[NBX]; 
    {%- for item in ocp.constraints.lbx %}
    lbx[{{ loop.index0 }}] = {{ item }};
    {%- endfor %}

    double ubx[NBX];
    {%- for item in ocp.constraints.ubx %}
    ubx[{{ loop.index0 }}] = {{ item }};
    {%- endfor %}
    
    // soft bounds on x
    int idxsbx[NSBX];
    {%- for item in ocp.constraints.idxsbx %}
    idxsbx[{{ loop.index0 }}] = {{ item }};
    {%- endfor %}

    double lsbx[NSBX]; 
    {%- for item in ocp.constraints.lsbx %}
    lsbx[{{ loop.index0 }}] = {{ item }};
    {%- endfor %}

    double usbx[NSBX];
    {%- for item in ocp.constraints.usbx %}
    usbx[{{ loop.index0 }}] = {{ item }};
    {%- endfor %}

    // set up general constraints for stage 0 to N-1 
    double D[NG*NU];
    double C[NG*NX];
    double lg[NG];
    double ug[NG];

    {%- for item_j in ocp.constraints.D %}
    {% set outer_loop = loop %}
    {%- for item_k in item_j %}
    D[{{ outer_loop.index0 }} + NG * {{ loop.index0 }}] = {{ item_k }}; 
    {%- endfor %}
    {%- endfor %}

    {%- for item_j in ocp.constraints.C %}
    {% set outer_loop = loop %}
    {%- for item_k in item_j %}
    C[{{ outer_loop.index0 }} + NG * {{ loop.index0 }}] = {{ item_k }}; 
    {%- endfor %}
    {%- endfor %}

    {%- for item in ocp.constraints.lg %}
    lg[{{ loop.index0 }}] = {{ item }};
    {%- endfor %}

    {%- for item in ocp.constraints.ug %}
    ug[{{ loop.index0 }}] = {{ item }};
    {%- endfor %}

    // set up nonlinear constraints for stage 0 to N-1 
    double lh[NH];
    double uh[NH];

    {%- for item in ocp.constraints.lh %}
    lh[{{ loop.index0 }}] = {{ item }};
    {%- endfor %}

    {%- for item in ocp.constraints.uh %}
    uh[{{ loop.index0 }}] = {{ item }};
    {%- endfor %}
    
    // set up bounds for last stage
    // x
    int idxbx_e[NBXN];
    {%- for item in ocp.constraints.idxbx_e %}
    idxbx_e[{{ loop.index0 }}] = {{ item }};
    {%- endfor %}

    double lbx_e[NBXN]; 
    {%- for item in ocp.constraints.lbx_e %}
    lbx_e[{{ loop.index0 }}] = {{ item }};
    {%- endfor %}
    
    double ubx_e[NBXN];
    {%- for item in ocp.constraints.ubx_e %}
    ubx_e[{{ loop.index0 }}] = {{ item }};
    {%- endfor %}
    
    // soft bounds on x
    int idxsbx_e[NSBXN];
    {%- for item in ocp.constraints.idxsbx_e %}
    idxsbx_e[{{ loop.index0 }}] = {{ item }};
    {%- endfor %}
    
    double lsbx_e[NSBXN]; 
    {%- for item in ocp.constraints.lsbx_e %}
    lsbx_e[{{ loop.index0 }}] = {{ item }};
    {%- endfor %}
    
    double usbx_e[NBXN];
    {%- for item in ocp.constraints.usbx_e %}
    usbx_e[{{ loop.index0 }}] = {{ item }};
    {%- endfor %}

    // set up soft bounds for nonlinear constraints
    int idxsh_e[NSHN];
    {%- for item in ocp.constraints.idxsh_e %}
    idxsh_e[{{ loop.index0 }}] = {{ item }};
    {%- endfor %}

    double lsh_e[NSHN]; 
    {%- for item in ocp.constraints.lsh_e %}
    lsh[{{ loop.index0 }}] = {{ item }};
    {%- endfor %}

    double ush_e[NSHN];
    {%- for item in ocp.constraints.ush_e %}
    ush_e[{{ loop.index0 }}] = {{ item }};
    {%- endfor %}
    
    // set up general constraints for last stage 
    double C_e[NGN*NX];
    double lg_e[NGN];
    double ug_e[NGN];

    {%- for item_j in ocp.constraints.C_e %}
    {% set outer_loop = loop %}
    {%- for item_k in item_j %}
    C_e[{{ outer_loop.index0 }} + NG * {{ loop.index0 }}] = {{ item_k }}; 
    {%- endfor %}
    {%- endfor %}

    {%- for item in ocp.constraints.lg_e %}
    lg_e[{{ loop.index0 }}] = {{ item }};
    {%- endfor %}

    {%- for item in ocp.constraints.ug_e %}
    ug_e[{{ loop.index0 }}] = {{ item }};
    {%- endfor %}

    // set up nonlinear constraints for last stage 
    double lh_e[NHN];
    double uh_e[NHN];

    {%- for item in ocp.constraints.lh_e %}
    lh_e[{{ loop.index0 }}] = {{ item }};
    {%- endfor %}

    {%- for item in ocp.constraints.uh_e %}
    uh_e[{{ loop.index0}}] = {{ item }};
    {%- endfor %}

    double yref[NY];
    double W[NY*NY];

    double Vx[NY*NX];
    double Vu[NY*NU];
    double Vz[NY*NZ];
    double Zl[NS];
    double Zu[NS];
    double zl[NS];
    double zu[NS];

    double yref_e[NYN];
    double W_e[NYN*NYN];
    double Zl_e[NSN];
    double Zu_e[NSN];
    double zl_e[NSN];
    double zu_e[NSN];

    double Vx_e[NYN*NX];
    
    for (int ii = 0; ii < NU + NX; ii++)
        yref[ii] = 0.0;

    {%- for item_j in ocp.cost.W %}
    {% set outer_loop = loop %}
    {%- for item_k in item_j %}
    W[{{ outer_loop.index0 }} + (NY) * {{ loop.index0 }}] = {{ item_k }}; 
    {%- endfor %}
    {%- endfor %}

    {%- for item_j in ocp.cost.Vx %}
    {% set outer_loop = loop %}
    {%- for item_k in item_j %}
    Vx[{{ outer_loop.index0 }} + (NY) * {{ loop.index0 }}] = {{ item_k }}; 
    {%- endfor %}
    {%- endfor %}

    {%- for item_j in ocp.cost.Vu %}
    {% set outer_loop = loop %}
    {%- for item_k in item_j %}
    Vu[{{ outer_loop.index0 }} + (NY) * {{ loop.index0 }}] = {{ item_k }}; 
    {%- endfor %}
    {%- endfor %}

    {%- for item_j in ocp.cost.Vz %}
    {% set outer_loop = loop %}
    {%- for item_k in item_j %}
    Vz[{{ outer_loop.index0 }} + (NY) * {{ loop.index0 }}] = {{ item_k }}; 
    {%- endfor %}
    {%- endfor %}

    {%- for item_j in ocp.cost.Zl %}
    {% set outer_loop = loop %}
    {%- for item_k in item_j %}
    {% if outer_loop.index0 == loop.index0 %}
    Zl[{{ outer_loop.index0 }}] = {{ item_k }}; 
    {% endif %}
    {%- endfor %}
    {%- endfor %}

    {%- for item_j in ocp.cost.Zu %}
    {% set outer_loop = loop %}
    {%- for item_k in item_j %}
    {% if outer_loop.index0 == loop.index0 %}
    Zu[{{ outer_loop.index0 }}] = {{ item_k }}; 
    {% endif %}
    {%- endfor %}
    {%- endfor %}

    {%- for item  in ocp.cost.zl %}
    zl[{{ loop.index0 }}] = {{ item }}; 
    {%- endfor %}

    {%- for item  in ocp.cost.zu %}
    zu[{{ loop.index0 }}] = {{ item }}; 
    {%- endfor %}

    {%- for item  in ocp.cost.yref %}
    yref[{{ loop.index0 }}] = {{ item }}; 
    {%- endfor %}

    {%- for item_j in ocp.cost.W_e %}
    {% set outer_loop = loop %}
    {%- for item_k in item_j %}
    W_e[{{ outer_loop.index0 }} + (NYN) * {{ loop.index0 }}] = {{ item_k }}; 
    {%- endfor %}
    {%- endfor %}

    {%- for item_j in ocp.cost.Vx_e %}
    {% set outer_loop = loop %}
    {%- for item_k in item_j %}
    Vx_e[{{ outer_loop.index0 }} + (NYN) * {{ loop.index0 }}] = {{ item_k }}; 
    {%- endfor %}
    {%- endfor %}

    {%- for item_j in ocp.cost.Zl_e %}
    {% set outer_loop = loop %}
    {%- for item_k in item_j %}
    {% if outer_loop.index0 == loop.index0 %}
    Zl_e[{{ outer_loop.index0 }}] = {{ item_k }}; 
    {% endif %}
    {%- endfor %}
    {%- endfor %}

    {%- for item_j in ocp.cost.Zu_e %}
    {% set outer_loop = loop %}
    {%- for item_k in item_j %}
    {% if outer_loop.index0 == loop.index0 %}
    Zu_e[{{ outer_loop.index0 }}] = {{ item_k }}; 
    {% endif %}
    {%- endfor %}
    {%- endfor %}

    {%- for item  in ocp.cost.zl_e %}
    zl_e[{{ loop.index0 }}] = {{ item }}; 
    {%- endfor %}

    {%- for item  in ocp.cost.zu_e %}
    zu_e[{{ loop.index0 }}] = {{ item }}; 
    {%- endfor %}

    {%- for item in ocp.cost.yref_e %}
    yref_e[{{ loop.index0 }}] = {{ item }}; 
    {%- endfor %}

    int max_num_sqp_iterations = 100;

    int nx[N+1];
    int nu[N+1];
    int nbx[N+1];
    int nbu[N+1];
    int nsbx[N+1];
    int nsbu[N+1];
    int nsh[N+1];
    int ns[N+1];
    int nb[N+1];
    int ng[N+1];
    int nh[N+1];
    int nz[N+1];
    int nv[N+1];
    int ny[N+1];
    int npd[N+1];
    int npd_e[N+1];

    for(int i = 0; i < N+1; i++) {
        nx[i]  = NX_;
        nu[i]  = NU_;
        nbx[i] = NBX_;
        nbu[i] = NBU_;
        nb[i]  = NBU_ + NBX_;
        nsbx[i]  = NSBX_;
        nsbu[i]  = NSBU_;
        nsh[i]  = NSH_;
        ns[i]  = NS_;
        ng[i]  = NG_;
        nh[i]  = NH_;
        npd[i] = NPD_;
        nz[i]  = NZ_;
        nv[i]  = NX_ + NU_;
        ny[i]  = NY_;
    }

    nbx[0] = NX_;
    nbu[0] = NBU_;
    nb[0]  = NX_ + NBU_;
    nsbx[N]  = NSBXN_;
    nsbu[N]  = 0;
    nsh[N]  = NSHN_;
    ns[N]  = NSN_;

    nu[N]  = 0;
    nx[N]  = NX_;
    nz[N]  = 0;
    nh[N]  = NHN_;
    npd[N]  = NPDN_;
    nv[N]  = NX_; 
    ny[N]  = NYN_;
    nbu[N] = 0;
    nbx[N] = NBXN_;
    ng[N]  = NGN_;
    nb[N]  = NBXN_;

    // Make plan
    nlp_solver_plan = ocp_nlp_plan_create(N);
    {% if ocp.solver_config.nlp_solver_type == "SQP" %}
    nlp_solver_plan->nlp_solver = SQP;
    {% else %}
    nlp_solver_plan->nlp_solver = SQP_RTI;
    {% endif %}
    nlp_solver_plan->ocp_qp_solver_plan.qp_solver = {{ ocp.solver_config.qp_solver }};
    for (int i = 0; i <= N; i++)
        nlp_solver_plan->nlp_cost[i] = LINEAR_LS;
    for (int i = 0; i < N; i++)
    {
        nlp_solver_plan->nlp_dynamics[i] = CONTINUOUS_MODEL;
        nlp_solver_plan->sim_solver_plan[i].sim_solver = {{ ocp.solver_config.integrator_type}};
    }

    for (int i = 0; i < N; i++) {
        {%- if ocp.dims.npd > 0 %}
        nlp_solver_plan->nlp_constraints[i] = BGHP;
        {%- else %}
        nlp_solver_plan->nlp_constraints[i] = BGH;
        {%- endif %}
    }

    {%- if ocp.dims.npd_e > 0 %}
    nlp_solver_plan->nlp_constraints[N] = BGHP;
    {%- else %}
    nlp_solver_plan->nlp_constraints[N] = BGH;
    {%- endif %}

    {% if ocp.solver_config.hessian_approx == "EXACT" %} 
    nlp_solver_plan->regularization = CONVEXIFICATION;
    {% endif %}
    nlp_config = ocp_nlp_config_create(*nlp_solver_plan);

    /* create and set ocp_nlp_dims */
	nlp_dims = ocp_nlp_dims_create(nlp_config);

    ocp_nlp_dims_set_opt_vars(nlp_config, nlp_dims, "nx", nx);
    ocp_nlp_dims_set_opt_vars(nlp_config, nlp_dims, "nu", nu);
    ocp_nlp_dims_set_opt_vars(nlp_config, nlp_dims, "nz", nz);
    ocp_nlp_dims_set_opt_vars(nlp_config, nlp_dims, "ns", ns);

    for (int i = 0; i <= N; i++) {
        ocp_nlp_dims_set_cost(nlp_config, nlp_dims, i, "ny", &ny[i]);
        ocp_nlp_dims_set_constraints(nlp_config, nlp_dims, i, "nbx", &nbx[i]);
        ocp_nlp_dims_set_constraints(nlp_config, nlp_dims, i, "nbu", &nbu[i]);
        ocp_nlp_dims_set_constraints(nlp_config, nlp_dims, i, "nsbx", &nsbx[i]);
        ocp_nlp_dims_set_constraints(nlp_config, nlp_dims, i, "nsbu", &nsbu[i]);
        ocp_nlp_dims_set_constraints(nlp_config, nlp_dims, i, "ng", &ng[i]);
        ocp_nlp_dims_set_constraints(nlp_config, nlp_dims, i, "nh", &nh[i]);
        ocp_nlp_dims_set_constraints(nlp_config, nlp_dims, i, "nsh", &nsh[i]);
    }

    {%- if ocp.dims.npd > 0 %}
    for (int i = 0; i < N; i++) 
        ocp_nlp_dims_set_constraints(nlp_config, nlp_dims, i, "np", &npd[i]);
    {%- endif %}
    {%- if ocp.dims.npd_e > 0 %}
    ocp_nlp_dims_set_constraints(nlp_config, nlp_dims, N, "np", &npd[N]);
    {%- endif %}

    {%- if ocp.dims.npd > 0 %}
    p_constraint = (external_function_casadi *) malloc(sizeof(external_function_casadi)*N);
    for (int i = 0; i < N; ++i) {
        // nonlinear part of convex-composite constraint
        p_constraint[i].casadi_fun = &{{ ocp.con_p.name }}_p_constraint;
        p_constraint[i].casadi_n_in = &{{ ocp.con_p.name }}_p_constraint_n_in;
        p_constraint[i].casadi_n_out = &{{ ocp.con_p.name }}_p_constraint_n_out;
        p_constraint[i].casadi_sparsity_in = &{{ ocp.con_p.name }}_p_constraint_sparsity_in;
        p_constraint[i].casadi_sparsity_out = &{{ ocp.con_p.name }}_p_constraint_sparsity_out;
        p_constraint[i].casadi_work = &{{ ocp.con_p.name }}_p_constraint_work;

        external_function_casadi_create(&p_constraint[i]);
    }
    {%- endif %}

    {%- if ocp.dims.npd_e > 0 %}
	// nonlinear part of convex-composite constraint
	external_function_casadi p_constraint_e;
	p_constraint_e.casadi_fun = &{{ ocp.con_p_e.name }}_p_constraint_e;
	p_constraint_e.casadi_n_in = &{{ ocp.con_p_e.name }}_p_constraint_e_n_in;
	p_constraint_e.casadi_n_out = &{{ ocp.con_p_e.name }}_p_constraint_e_n_out;
	p_constraint_e.casadi_sparsity_in = &{{ ocp.con_p_e.name }}_p_constraint_e_sparsity_in;
	p_constraint_e.casadi_sparsity_out = &{{ ocp.con_p_e.name }}_p_constraint_e_sparsity_out;
	p_constraint_e.casadi_work = &{{ ocp.con_p_e.name }}_p_constraint_e_work;

    external_function_casadi_create(p_constraint_e);
    {%- endif %}

    {%- if ocp.dims.nh > 0 %}
    h_constraint = (external_function_casadi *) malloc(sizeof(external_function_casadi)*N);
    for (int i = 0; i < N; ++i) {
        // nonlinear constraint
        h_constraint[i].casadi_fun = &{{ ocp.con_h.name }}_h_constraint;
        h_constraint[i].casadi_n_in = &{{ ocp.con_h.name }}_h_constraint_n_in;
        h_constraint[i].casadi_n_out = &{{ ocp.con_h.name }}_h_constraint_n_out;
        h_constraint[i].casadi_sparsity_in = &{{ ocp.con_h.name }}_h_constraint_sparsity_in;
        h_constraint[i].casadi_sparsity_out = &{{ ocp.con_h.name }}_h_constraint_sparsity_out;
        h_constraint[i].casadi_work = &{{ ocp.con_h.name }}_h_constraint_work;

        external_function_casadi_create(&h_constraint[i]);
    }
    {%- endif %}

    {%- if ocp.dims.nh_e > 0 %}
	// nonlinear constraint
	external_function_casadi h_constraint_e;
	h_constraint_e.casadi_fun = &{{ ocp.con_h_e.name }}_h_constraint_e;
	h_constraint_e.casadi_n_in = &{{ ocp.con_h_e.name }}_h_constraint_e_n_in;
	h_constraint_e.casadi_n_out = &{{ ocp.con_h_e.name }}_h_constraint_e_n_out;
	h_constraint_e.casadi_sparsity_in = &{{ ocp.con_h_e.name }}_h_constraint_e_sparsity_in;
	h_constraint_e.casadi_sparsity_out = &{{ ocp.con_h_e.name }}_h_constraint_e_sparsity_out;
	p_constraint_e.casadi_work = &{{ ocp.con_h_e.name }}_h_constraint_e_work;

    external_function_casadi_create(h_constraint_e);
    {%- endif %}

    {% if ocp.solver_config.integrator_type == "ERK" %}
    // explicit ode
    {% if ocp.dims.np < 1 %}
    forw_vde_casadi = (external_function_casadi *) malloc(sizeof(external_function_casadi)*N);
    {% else %}
    forw_vde_casadi = (external_function_param_casadi *) malloc(sizeof(external_function_param_casadi)*N);
    {% endif %}

    for (int i = 0; i < N; ++i) {
        forw_vde_casadi[i].casadi_fun = &{{ ocp.model.name }}_expl_vde_forw;
        forw_vde_casadi[i].casadi_n_in = &{{ ocp.model.name }}_expl_vde_forw_n_in;
        forw_vde_casadi[i].casadi_n_out = &{{ ocp.model.name }}_expl_vde_forw_n_out;
        forw_vde_casadi[i].casadi_sparsity_in = &{{ ocp.model.name }}_expl_vde_forw_sparsity_in;
        forw_vde_casadi[i].casadi_sparsity_out = &{{ ocp.model.name }}_expl_vde_forw_sparsity_out;
        forw_vde_casadi[i].casadi_work = &{{ ocp.model.name }}_expl_vde_forw_work;
        external_function_casadi_create(&forw_vde_casadi[i]);
    }

    {% if ocp.solver_config.hessian_approx == "EXACT" %} 
    external_function_casadi * hess_vde_casadi;
    {% if ocp.dims.np < 1 %}
    hess_vde_casadi = (external_function_casadi *) malloc(sizeof(external_function_casadi)*N);
    {% else %}
    hess_vde_casadi = (external_function_param_casadi *) malloc(sizeof(external_function_param_casadi)*N);
    {% endif %}
    for (int i = 0; i < N; ++i) {
        hess_vde_casadi[i].casadi_fun = &{{ ocp.model.name }}_expl_ode_hess;
        hess_vde_casadi[i].casadi_n_in = &{{ ocp.model.name }}_expl_ode_hess_n_in;
        hess_vde_casadi[i].casadi_n_out = &{{ ocp.model.name }}_expl_ode_hess_n_out;
        hess_vde_casadi[i].casadi_sparsity_in = &{{ ocp.model.name }}_expl_ode_hess_sparsity_in;
        hess_vde_casadi[i].casadi_sparsity_out = &{{ ocp.model.name }}_expl_ode_hess_sparsity_out;
        hess_vde_casadi[i].casadi_work = &{{ ocp.model.name }}_expl_ode_hess_work;
        external_function_casadi_create(&hess_vde_casadi[i]);
    }
    {% endif %}
    {% else %}
    {% if ocp.solver_config.integrator_type == "IRK" %}
    // implicit dae
    {% if ocp.dims.np < 1 %}
    impl_dae_fun = (external_function_casadi *) malloc(sizeof(external_function_casadi)*N);
    {% else %}
    impl_dae_fun = (external_function_param_casadi *) malloc(sizeof(external_function_param_casadi)*N);
    {% endif %}
    for (int i = 0; i < N; ++i) {
        impl_dae_fun[i].casadi_fun = &{{ ocp.model.name }}_impl_dae_fun;
        impl_dae_fun[i].casadi_work = &{{ ocp.model.name }}_impl_dae_fun_work;
        impl_dae_fun[i].casadi_sparsity_in = &{{ ocp.model.name }}_impl_dae_fun_sparsity_in;
        impl_dae_fun[i].casadi_sparsity_out = &{{ ocp.model.name }}_impl_dae_fun_sparsity_out;
        impl_dae_fun[i].casadi_n_in = &{{ ocp.model.name }}_impl_dae_fun_n_in;
        impl_dae_fun[i].casadi_n_out = &{{ ocp.model.name }}_impl_dae_fun_n_out;
        // TODO(fix this!!)
        {% if ocp.dims.np < 1 %}
        external_function_casadi_create(&impl_dae_fun[i]);
        {% else %}
        external_function_param_casadi_create(&impl_dae_fun[i], {{ocp.dims.np}});
        {% endif %}
    }

    {% if ocp.dims.np < 1 %}
    impl_dae_fun_jac_x_xdot_z = (external_function_casadi *) malloc(sizeof(external_function_casadi)*N);
    {% else %}
    impl_dae_fun_jac_x_xdot_z = (external_function_param_casadi *) malloc(sizeof(external_function_param_casadi)*N);
    {% endif %}
    for (int i = 0; i < N; ++i) {
        impl_dae_fun_jac_x_xdot_z[i].casadi_fun = &{{ ocp.model.name }}_impl_dae_fun_jac_x_xdot_z;
        impl_dae_fun_jac_x_xdot_z[i].casadi_work = &{{ ocp.model.name }}_impl_dae_fun_jac_x_xdot_z_work;
        impl_dae_fun_jac_x_xdot_z[i].casadi_sparsity_in = &{{ ocp.model.name }}_impl_dae_fun_jac_x_xdot_z_sparsity_in;
        impl_dae_fun_jac_x_xdot_z[i].casadi_sparsity_out = &{{ ocp.model.name }}_impl_dae_fun_jac_x_xdot_z_sparsity_out;
        impl_dae_fun_jac_x_xdot_z[i].casadi_n_in = &{{ ocp.model.name }}_impl_dae_fun_jac_x_xdot_z_n_in;
        impl_dae_fun_jac_x_xdot_z[i].casadi_n_out = &{{ ocp.model.name }}_impl_dae_fun_jac_x_xdot_z_n_out;
        {% if ocp.dims.np < 1 %}
        external_function_casadi_create(&impl_dae_fun_jac_x_xdot_z[i]);
        {% else %}
        external_function_param_casadi_create(&impl_dae_fun_jac_x_xdot_z[i], {{ocp.dims.np}});
        {% endif %}
    }

    {% if ocp.dims.np < 1 %}
    impl_dae_jac_x_xdot_u_z = (external_function_casadi *) malloc(sizeof(external_function_casadi)*N);
    {% else %}
    impl_dae_jac_x_xdot_u_z = (external_function_param_casadi *) malloc(sizeof(external_function_param_casadi)*N);
    {% endif %}
    for (int i = 0; i < N; ++i) {
        impl_dae_jac_x_xdot_u_z[i].casadi_fun = &{{ ocp.model.name }}_impl_dae_jac_x_xdot_u_z;
        impl_dae_jac_x_xdot_u_z[i].casadi_work = &{{ ocp.model.name }}_impl_dae_jac_x_xdot_u_z_work;
        impl_dae_jac_x_xdot_u_z[i].casadi_sparsity_in = &{{ ocp.model.name }}_impl_dae_jac_x_xdot_u_z_sparsity_in;
        impl_dae_jac_x_xdot_u_z[i].casadi_sparsity_out = &{{ ocp.model.name }}_impl_dae_jac_x_xdot_u_z_sparsity_out;
        impl_dae_jac_x_xdot_u_z[i].casadi_n_in = &{{ ocp.model.name }}_impl_dae_jac_x_xdot_u_z_n_in;
        impl_dae_jac_x_xdot_u_z[i].casadi_n_out = &{{ ocp.model.name }}_impl_dae_jac_x_xdot_u_z_n_out;
        {% if ocp.dims.np < 1 %}
        external_function_casadi_create(&impl_dae_jac_x_xdot_u_z[i]);
        {% else %}
        external_function_param_casadi_create(&impl_dae_jac_x_xdot_u_z[i], {{ocp.dims.np}});
        {% endif %}
    }
    {% endif %}
    {% endif %}

    nlp_in = ocp_nlp_in_create(nlp_config, nlp_dims);

    for (int i = 0; i < N; ++i)
        nlp_in->Ts[i] = Tf/N;

    // NLP cost linear least squares
    // C  // TODO(oj) this can be done using
    // // ocp_nlp_cost_set_model(nlp_config, nlp_dims, nlp_in, i, "Cyt", Cyt);
    // ocp_nlp_cost_ls_model **cost_ls = (ocp_nlp_cost_ls_model **) nlp_in->cost;
    // for (int i = 0; i <= N; ++i) {
    //     blasfeo_dgese(nv[i], ny[i], 0.0, &cost_ls[i]->Cyt, 0, 0);
    //     for (int j = 0; j < nu[i]; j++)
    //         BLASFEO_DMATEL(&cost_ls[i]->Cyt, j, nx[i]+j) = 1.0;
    //     for (int j = 0; j < nx[i]; j++)
    //         BLASFEO_DMATEL(&cost_ls[i]->Cyt, nu[i]+j, j) = 1.0;
    // }
    // W
    for (int i = 0; i < N; ++i) {
        ocp_nlp_cost_model_set(nlp_config, nlp_dims, nlp_in, i, "W", W);
    }
    // W_e
    ocp_nlp_cost_model_set(nlp_config, nlp_dims, nlp_in, N, "W", W_e);


	for (int i = 0; i < N; ++i) {
        ocp_nlp_cost_model_set(nlp_config, nlp_dims, nlp_in, i, "Vx", Vx);
        ocp_nlp_cost_model_set(nlp_config, nlp_dims, nlp_in, i, "Vu", Vu);
        ocp_nlp_cost_model_set(nlp_config, nlp_dims, nlp_in, i, "Vz", Vz);
        ocp_nlp_cost_model_set(nlp_config, nlp_dims, nlp_in, i, "yref", yref);
        ocp_nlp_cost_model_set(nlp_config, nlp_dims, nlp_in, i, "Zl", Zl);
        ocp_nlp_cost_model_set(nlp_config, nlp_dims, nlp_in, i, "Zu", Zu);
        ocp_nlp_cost_model_set(nlp_config, nlp_dims, nlp_in, i, "zl", zl);
        ocp_nlp_cost_model_set(nlp_config, nlp_dims, nlp_in, i, "zu", zu);
	}

    ocp_nlp_cost_model_set(nlp_config, nlp_dims, nlp_in, N, "Vx", Vx_e);
    ocp_nlp_cost_model_set(nlp_config, nlp_dims, nlp_in, N, "yref", yref_e);
    ocp_nlp_cost_model_set(nlp_config, nlp_dims, nlp_in, N, "Zl", Zl_e);
    ocp_nlp_cost_model_set(nlp_config, nlp_dims, nlp_in, N, "Zu", Zu_e);
    ocp_nlp_cost_model_set(nlp_config, nlp_dims, nlp_in, N, "zl", zl_e);
    ocp_nlp_cost_model_set(nlp_config, nlp_dims, nlp_in, N, "zu", zu_e);

    // NLP dynamics
    int set_fun_status;
    for (int i = 0; i < N; ++i) {
    {% if ocp.solver_config.integrator_type == "ERK" %} 
        set_fun_status = ocp_nlp_dynamics_model_set(nlp_config, nlp_dims, nlp_in, i, "expl_vde_for", &forw_vde_casadi[i]);
        if (set_fun_status != 0) { printf("Error while setting expl_vde_for[%i]\n", i);  exit(1); }
        {% if ocp.solver_config.hessian_approx == "EXACT" %} 
            set_fun_status = ocp_nlp_dynamics_model_set(nlp_config, nlp_dims, nlp_in, i, "expl_ode_hes", &hess_vde_casadi[i]);
            if (set_fun_status != 0) { printf("Error while setting expl_ode_hes[%i]\n", i);  exit(1); }
        {% endif %}
    {% else %}
    {% if ocp.solver_config.integrator_type == "IRK" %} 
        set_fun_status = ocp_nlp_dynamics_model_set(nlp_config, nlp_dims, nlp_in, i, "impl_ode_fun", &impl_dae_fun[i]);
        if (set_fun_status != 0) { printf("Error while setting impl_dae_fun[%i]\n", i);  exit(1); }
        set_fun_status = ocp_nlp_dynamics_model_set(nlp_config, nlp_dims, nlp_in, i, "impl_ode_fun_jac_x_xdot", &impl_dae_fun_jac_x_xdot_z[i]);
        if (set_fun_status != 0) { printf("Error while setting impl_dae_fun_jac_x_xdot_z[%i]\n", i);  exit(1); }
        set_fun_status = ocp_nlp_dynamics_model_set(nlp_config, nlp_dims, nlp_in, i, "impl_ode_jac_x_xdot_u", &impl_dae_jac_x_xdot_u_z[i]);
        if (set_fun_status != 0) { printf("Error while setting impl_dae_jac_x_xdot_u_z[%i]\n", i);  exit(1); }
    {% endif %}
    {% endif %}
    }

    // NLP constraints
    // TODO(oj) remove this when idxb setter available

    // bounds for stage 0
    ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, 0, "idxbx", idxbx0);
    ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, 0, "lbx", lbx0);
    ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, 0, "ubx", ubx0);

    ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, 0, "idxbu", idxbu);
    ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, 0, "lbu", lbu);
    ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, 0, "ubu", ubu);

    {%- if ocp.dims.nsbx > 0 %} 
    ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, 0, "idxsbx", idxsbx);
    ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, 0, "lsbx", lsbx);
    ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, 0, "usbx", usbx);
    {%- endif %}
    
    {%- if ocp.dims.nsbu > 0 %} 
    ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, 0, "idxsbu", idxsbu);
    ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, 0, "lsbu", lsbu);
    ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, 0, "usbu", usbu);
    {%- endif %}
    
    {%- if ocp.dims.nsh > 0 %} 
    ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, 0, "idxsh", idxsh);
    ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, 0, "lsh", lsh);
    ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, 0, "ush", ush);
    {%- endif %}

    // bounds for intermediate stages
    for (int i = 1; i < N; ++i)
    {
        ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, i, "idxbx", idxbx);
        ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, i, "lbx", lbx);
        ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, i, "ubx", ubx);
        
        ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, i, "idxbu", idxbu);
        ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, i, "lbu", lbu);
        ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, i, "ubu", ubu);

        {%- if ocp.dims.nsbx > 0 %} 
        ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, i, "idxsbx", idxsbx);
        ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, i, "lsbx", lsbx);
        ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, i, "usbx", usbx);
        {%- endif %}
        
        {%- if ocp.dims.nsbu > 0 %} 
        ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, i, "idxsbu", idxsbu);
        ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, i, "lsbu", lsbu);
        ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, i, "usbu", usbu);
        {%- endif %}

        {%- if ocp.dims.nsh > 0 %} 
        ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, i, "idxsh", idxsh);
        ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, i, "lsh", lsh);
        ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, i, "ush", ush);
        {%- endif %}

    }
   
    {%- if ocp.dims.ng > 0 %} 
    // general constraints for stages 0 to N-1
    for (int i = 0; i < N; ++i)
    {
        ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, i, "D", D);
        ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, i, "C", C);
        ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, i, "lg", lg);
        ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, i, "ug", ug);
    }
    {%- endif %}
    
    {%- if ocp.dims.nbx_e > 0 %} 
    // bounds for last
    ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, N, "idxbx", idxbx_e);
    ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, N, "lbx", lbx_e);
    ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, N, "ubx", ubx_e);
    {%- endif %}

    {%- if ocp.dims.nsbx_e > 0 %} 
    ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, N, "idxsbx", idxsbx_e);
    ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, N, "lsbx", lsbx_e);
    ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, N, "usbx", usbx_e);
    {%- endif %}
    
    {%- if ocp.dims.nsh_e > 0 %} 
    ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, N, "idxsh", idxsh_e);
    ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, N, "lsh", lsh_e);
    ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, N, "ush", ush_e);
    {%- endif %}

    {%- if ocp.dims.ng_e > 0 %} 
    // general constraints for last stage
    ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, N, "C", C_e);
    ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, N, "lg", lg_e);
    ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, N, "ug", ug_e);
    {%- endif %}

    
    {%- if ocp.dims.npd > 0 %}
    // convex-composite constraints for stages 0 to N-1
    for (int i = 0; i < N; ++i)
        ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, i, "p", &p_constraint[i]);
    {%- endif %}

    {%- if ocp.dims.npd_e > 0 %}
    // convex-composite constraints for stage N
    ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, N, "p", &p_constraint_e[i]);
    {%- endif %}

    {%- if ocp.dims.nh > 0 %}
    // nonlinear constraints for stages 0 to N-1
    for (int i = 0; i < N; ++i)
    {
        ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, i, "nl_constr_h_fun_jac", &h_constraint[i]);
        ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, i, "lh", lh);
        ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, i, "uh", uh);
    }
    {%- endif %}

    {%- if ocp.dims.nh_e > 0 %}
    // nonlinear constraints for stage N
    ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, N, "nl_constr_h_fun_jac", &h_constraint_e[i]);
    ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, N, "lh", lh_e);
    ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, N, "uh", uh_e);
    {%- endif %}

    nlp_opts = ocp_nlp_opts_create(nlp_config, nlp_dims);
    
    {% if ocp.dims.nz > 0 %}
    bool output_z_val = true; 
    bool sens_algebraic_val = true; 
    int num_steps_val = 1; 
    for (int i = 0; i < N; i++) ocp_nlp_dynamics_opts_set(nlp_config, nlp_opts, i, "output_z", &output_z_val);
    for (int i = 0; i < N; i++) ocp_nlp_dynamics_opts_set(nlp_config, nlp_opts, i, "sens_algebraic", &sens_algebraic_val);
    for (int i = 0; i < N; i++) ocp_nlp_dynamics_opts_set(nlp_config, nlp_opts, i, "num_steps", &num_steps_val);
    {% endif %}
    int ns_val = 1; 
    for (int i = 0; i < N; i++) ocp_nlp_dynamics_opts_set(nlp_config, nlp_opts, i, "ns", &ns_val);
    bool jac_reuse_val = true;
    for (int i = 0; i < N; i++) ocp_nlp_dynamics_opts_set(nlp_config, nlp_opts, i, "jac_reuse", &jac_reuse_val);

    {% if ocp.solver_config.nlp_solver_type == "SQP" %}

    int max_iter = max_num_sqp_iterations;
    double tol_stat = 1e-6;
    double tol_eq   = 1e-6;
    double tol_ineq = 1e-6;
    double tol_comp = 1e-6;

    ocp_nlp_opts_set(nlp_config, nlp_opts, "max_iter", &max_iter);
    ocp_nlp_opts_set(nlp_config, nlp_opts, "tol_stat", &tol_stat);
    ocp_nlp_opts_set(nlp_config, nlp_opts, "tol_eq", &tol_eq);
    ocp_nlp_opts_set(nlp_config, nlp_opts, "tol_ineq", &tol_ineq);
    ocp_nlp_opts_set(nlp_config, nlp_opts, "tol_comp", &tol_comp);


    {% else %}
    // ocp_nlp_sqp_rti_opts *sqp_opts = (ocp_nlp_sqp_rti_opts *) nlp_opts;
    {% endif %}
    {% if ocp.solver_config.hessian_approx == "EXACT" %}
    for (int i = 0; i < N; ++i)
    {
        // TODO(oj) is the following needed, and what does it do? do we
        // ((ocp_nlp_dynamics_cont_opts *) sqp_opts->dynamics[i])->compute_hess = true;
        int num_steps = 5;
        bool sens_hess = true;
        bool sens_adj = true;

        ocp_nlp_dynamics_opts_set(nlp_config, nlp_opts, i, "num_steps", &num_steps);
        ocp_nlp_dynamics_opts_set(nlp_config, nlp_opts, i, "sens_hess", &sens_hess);
        ocp_nlp_dynamics_opts_set(nlp_config, nlp_opts, i, "sens_adj", &sens_adj);
    }
    {% endif %}

    nlp_out = ocp_nlp_out_create(nlp_config, nlp_dims);
    for (int i = 0; i <= N; ++i) {
        blasfeo_dvecse(nu[i]+nx[i], 0.0, nlp_out->ux+i, 0);
    }
    
    nlp_solver = ocp_nlp_solver_create(nlp_config, nlp_dims, nlp_opts);

    // initialize parameters to nominal value
    {% if ocp.dims.np > 0%}
    double p[{{ ocp.dims.np }}];
    {%- for item in ocp.constraints.p %}
    p[{{ loop.index0 }}] = {{ item }};
    {%- endfor %}
    {% if ocp.solver_config.integrator_type == "IRK" %}
    for (int ii = 0; ii < {{ ocp.dims.N }}; ii++) {
    impl_dae_fun[ii].set_param(impl_dae_fun+ii, p);
    impl_dae_fun_jac_x_xdot_z[ii].set_param(impl_dae_fun_jac_x_xdot_z+ii, p);
    impl_dae_jac_x_xdot_u_z[ii].set_param(impl_dae_jac_x_xdot_u_z+ii, p);
    }
    {% else %}
    for (int ii = 0; ii < {{ ocp.dims.N }}; ii++) {
    expl_vde_for[ii].set_param(expl_vde_for+ii, p);
    }
    {% endif %}
    {% endif %}

    return status;
}

int acados_solve() {

    // solve NLP 
    int solver_status = ocp_nlp_solve(nlp_solver, nlp_in, nlp_out);

    return solver_status;
}

int acados_free() {

    // free memory
    ocp_nlp_opts_destroy(nlp_opts);
    ocp_nlp_in_destroy(nlp_in);
    ocp_nlp_out_destroy(nlp_out);
    ocp_nlp_solver_destroy(nlp_solver);
    ocp_nlp_dims_destroy(nlp_dims);
    ocp_nlp_config_destroy(nlp_config);
    ocp_nlp_plan_destroy(nlp_solver_plan);

    // free external function 
    {% if ocp.solver_config.integrator_type == "IRK" %}
    for(int i = 0; i < {{ocp.dims.N}}; i++) {
        {% if ocp.dims.np < 1 %}
        external_function_casadi_free(&impl_dae_fun[i]);
        external_function_casadi_free(&impl_dae_fun_jac_x_xdot_z[i]);
        external_function_casadi_free(&impl_dae_jac_x_xdot_u_z[i]);
        {% else %}
        external_function_param_casadi_free(&impl_dae_fun[i]);
        external_function_param_casadi_free(&impl_dae_fun_jac_x_xdot_z[i]);
        external_function_param_casadi_free(&impl_dae_jac_x_xdot_u_z[i]);
        {% endif %}
    }
    {% else %}
    for(int i = 0; i < {{ocp.dims.N}}; i++) {
        {% if ocp.dims.np < 1 %}
        external_function_casadi_free(&forw_vde_casadi[i]);
        {% else %}
        external_function_param_casadi_free(&forw_vde_casadi[i]);
        {% endif %}
    {% if ocp.solver_config.hessian_approx == "EXACT" %}
        {% if ocp.dims.np < 1 %}
        external_function_casadi_free(&hess_vde_casadi[i]);
        {% else %}
        external_function_param_casadi_free(&hess_vde_casadi[i]);
        {% endif %}
    {% endif %}
    }
    {% endif %}
    
    return 0;
}

ocp_nlp_in * acados_get_nlp_in() { return  nlp_in; }
ocp_nlp_out * acados_get_nlp_out() { return  nlp_out; }
ocp_nlp_solver * acados_get_nlp_solver() { return  nlp_solver; }
ocp_nlp_config * acados_get_nlp_config() { return  nlp_config; }
void * acados_get_nlp_opts() { return  nlp_opts; }
ocp_nlp_dims * acados_get_nlp_dims() { return  nlp_dims; }
