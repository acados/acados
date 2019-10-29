

// standard
#include <stdio.h>
#include <stdlib.h>
// acados
#include "acados/utils/print.h"
#include "acados_c/ocp_nlp_interface.h"
#include "acados_c/external_function_interface.h"

// example specific
#include "{{ model.name }}_model/{{ model.name }}_model.h"
{% if constraints.constr_type == "BGP" %}
#include "{{ con_h.name }}_r_constraint/{{ con_h.name }}_r_constraint.h"
{% endif %}
{% if constraints.constr_type_e == "BGP" %}
#include "{{ con_h_e.name }}_r_e_constraint/{{ con_h_e.name }}_r_e_constraint.h"
{% endif %}
{% if dims.nh > 0 %}
#include "{{ con_h.name }}_h_constraint/{{ con_h.name }}_h_constraint.h"
{% endif %}
{% if dims.nh_e > 0 %}
#include "{{ con_h_e.name }}_h_e_constraint/{{ con_h_e.name }}_h_e_constraint.h"
{% endif %}
{%- if cost.cost_type == "NONLINEAR_LS" %}
#include "{{ cost_phi.name }}_r_cost/{{ cost_r.name }}_r_cost.h"
{% endif %}
{%- if cost.cost_type_e == "NONLINEAR_LS" %}
#include "{{ cost_h_e.name }}_r_e_cost/{{ cost_r_e.name }}_r_e_cost.h"
{% endif %}

#include "acados_solver_{{ model.name }}.h"

#define NX_    {{ dims.nx }}
#define NZ_    {{ dims.nz }}
#define NU_    {{ dims.nu }}
#define NP_    {{ dims.np }}
#define NBX_   {{ dims.nbx }}
#define NBU_   {{ dims.nbu }}
#define NSBX_  {{ dims.nsbx }}
#define NSBU_  {{ dims.nsbu }}
#define NSH_   {{ dims.nsh }}
#define NSHN_  {{ dims.nsh_e }}
#define NSBXN_ {{ dims.nsbx_e }}
#define NS_    {{ dims.ns }}
#define NSN_   {{ dims.ns_e }}
#define NG_    {{ dims.ng }}
#define NBXN_  {{ dims.nbx_e }}
#define NGN_   {{ dims.ng_e }}
#define NY_    {{ dims.ny }}
#define NYN_   {{ dims.ny_e }}
#define N_     {{ dims.N }}
#define NH_    {{ dims.nh }}
#define NHN_   {{ dims.nh_e }}
#define NR_    {{ dims.nr }}
#define NRN_   {{ dims.nr_e }}

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

#if NR_ < 1
#define NR   1
#else
#define NR   NR_
#endif

#if NRN_ < 1
#define NRN   1
#else
#define NRN   NRN_
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

int acados_create()
{

    int status = 0;

    double Tf = {{ solver_config.tf }};

    // set up bounds for stage 0 
    // u
    int idxbu0[NBU];
    {% for i in range(end=dims.nbu) %}
    idxbu0[{{ i }}] = {{ constraints.idxbu[i] }};
    {%- endfor %}
    double lbu0[NBU]; 
    double ubu0[NBU];
    {% for i in range(end=dims.nbu) %}
    lbu0[{{ i }}] = {{ constraints.lbu[i] }};
    ubu0[{{ i }}] = {{ constraints.ubu[i] }};
    {%- endfor %}
    
    // x
    int idxbx0[NX];
    {% for i in range(end=dims.nx) %}
    idxbx0[{{ i }}] = {{ i }};
    {%- endfor %}
    double lbx0[NX]; 
    double ubx0[NX];
    {% for i in range(end=dims.nx) %}
    lbx0[{{ i }}] = {{ constraints.x0[i] }};
    ubx0[{{ i }}] = {{ constraints.x0[i] }};
    {%- endfor %}


    // set up bounds for intermediate stages
    // u
    int idxbu[NBU];
    {% for i in range(end=dims.nbu) %}
    idxbu[{{ i }}] = {{ constraints.idxbu[i] }};
    {%- endfor %}
    double lbu[NBU]; 
    double ubu[NBU];
    {% for i in range(end=dims.nbu) %}
    lbu[{{ i }}] = {{ constraints.lbu[i] }};
    ubu[{{ i }}] = {{ constraints.ubu[i] }};
    {%- endfor %}
    
    // set up soft bounds for u
    int idxsbu[NSBU];
    {% for i in range(end=dims.nsbu) %}
    idsxbu[{{ i }}] = {{ constraints.idxsbu[i] }};
    {%- endfor %}
    double lsbu[NSBU]; 
    double usbu[NSBU];
    {% for i in range(end=dims.nsbu) %}
    lsbu[{{ i }}] = {{ constraints.lsbu[i] }};
    usbu[{{ i }}] = {{ constraints.usbu[i] }};
    {%- endfor %}
    
    // set up soft bounds for nonlinear constraints
    int idxsh[NSBU];
    {% for i in range(end=dims.idxsh) %}
    idsxbu[{{ i }}] = {{ constraints.idxsh[i] }};
    {%- endfor %}
    double lsh[NSH]; 
    double ush[NSH];
    {% for i in range(end=dims.nsh) %}
    lsh[{{ i }}] = {{ constraints.lsh[i] }};
    ush[{{ i }}] = {{ constraints.ush[i] }};
    {%- endfor %}

    // x
    int idxbx[NBX];
    {% for i in range(end=dims.nbx) %}
    idxbx[{{ i }}] = {{ constraints.idxbx[i] }};
    {%- endfor %}
    double lbx[NBX]; 
    double ubx[NBX];
    {% for i in range(end=dims.nbx) %}
    lbx[{{ i }}] = {{ constraints.lbx[i] }};
    ubx[{{ i }}] = {{ constraints.ubx[i] }};
    {%- endfor %}

    // soft bounds on x
    int idxsbx[NSBX];
    {% for i in range(end=dims.nsbx) %}
    idsxbx[{{ i }}] = {{ constraints.idxsbx[i] }};
    {%- endfor %}
    double lsbx[NSBX]; 
    double usbx[NSBX];
    {% for i in range(end=dims.nsbx) %}
    lsbx[{{ i }}] = {{ constraints.lsbx[i] }};
    usbx[{{ i }}] = {{ constraints.usbx[i] }};
    {%- endfor %}
    
    // set up general constraints for stage 0 to N-1 
    double D[NG*NU];
    double C[NG*NX];
    double lg[NG];
    double ug[NG];

    {% for j in range(end=dims.ng) %}
        {%- for k in range(end=dims.nx) %}
    D[{{ j }}+NG * {{ k }}] = {{ constraints.D[j][k] }}; 
        {%- endfor %}
    {%- endfor %}

    {% for j in range(end=dims.ng) %}
        {%- for k in range(end=dims.nu) %}
    C[{{ j }}+NG * {{ k }}] = {{ constraints.C[j][k] }}; 
        {%- endfor %}
    {%- endfor %}

    {% for i in range(end=dims.ng) %}
    lg[{{ i }}] = {{ constraints.lg[i] }};
    {%- endfor %}

    {% for i in range(end=dims.ng) %}
    ug[{{ i }}] = {{ constraints.ug[i] }};
    {%- endfor %}

    // set up nonlinear constraints for stage 0 to N-1 
    double lh[NH];
    double uh[NH];

    {% for i in range(end=dims.nh) %}
    lh[{{ i }}] = {{ constraints.lh[i] }};
    {%- endfor %}

    {% for i in range(end=dims.nh) %}
    uh[{{ i }}] = {{ constraints.uh[i] }};
    {%- endfor %}
    
    // set up bounds for last stage
    // x
    int idxbx_e[NBXN];
    {% for i in range(end=dims.nbx_e) %}
    idxbx_e[{{ i }}] = {{ constraints.idxbx_e[i] }};
    {%- endfor %}
    double lbx_e[NBXN]; 
    double ubx_e[NBXN];
    {% for i in range(end=dims.nbx_e) %}
    lbx_e[{{ i }}] = {{ constraints.lbx_e[i] }};
    ubx_e[{{ i }}] = {{ constraints.ubx_e[i] }};
    {%- endfor %}
    
    // set up general constraints for last stage 
    double C_e[NGN*NX];
    double lg_e[NGN];
    double ug_e[NGN];

    {% for j in range(end=dims.ng_e) %}
        {%- for k in range(end=dims.nu) %}
    C_e[{{ j }}+NG * {{ k }}] = {{ constraints.C_e[j][k] }}; 
        {%- endfor %}
    {%- endfor %}

    {% for i in range(end=dims.ng_e) %}
    lg_e[{{ i }}] = {{ constraints.lg_e[i] }};
    {%- endfor %}

    {% for i in range(end=dims.ng_e) %}
    ug_e[{{ i }}] = {{ constraints.ug_e[i] }};
    {%- endfor %}

    // set up nonlinear constraints for last stage 
    double lh_e[NHN];
    double uh_e[NHN];

    {% for i in range(end=dims.nh_e) %}
    lh_e[{{ i }}] = {{ constraints.lh_e[i] }};
    {%- endfor %}

    {% for i in range(end=dims.nh_e) %}
    uh_e[{{ i }}] = {{ constraints.uh_e[i] }};
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

    {% for j in range(end=dims.ny) %}
        {%- for k in range(end=dims.ny) %}
    W[{{ j }}+(NY) * {{ k }}] = {{ cost.W[j][k] }}; 
        {%- endfor %}
    {%- endfor %}

    {% for j in range(end=dims.ny) %}
        {%- for k in range(end=dims.nx) %}
    Vx[{{ j }}+(NY) * {{ k }}] = {{ cost.Vx[j][k] }}; 
        {%- endfor %}
    {%- endfor %}

    {% for j in range(end=dims.ny) %}
        {%- for k in range(end=dims.nu) %}
    Vu[{{ j }}+(NY) * {{ k }}] = {{ cost.Vu[j][k] }}; 
        {%- endfor %}
    {%- endfor %}

    {% for j in range(end=dims.ny) %}
        {%- for k in range(end=dims.nz) %}
    Vz[{{ j }}+(NY) * {{ k }}] = {{ cost.Vz[j][k] }}; 
        {%- endfor %}
    {%- endfor %}

    {% for j in range(end=dims.ns) %}
    Zl[{{ j }}] = {{ cost.Zl[j] }}; 
    {%- endfor %}

    {% for j in range(end=dims.ns) %}
    Zu[{{ j }}] = {{ cost.Zu[j] }}; 
    {%- endfor %}

    {% for j in range(end=dims.ns) %}
    zl[{{ j }}] = {{ cost.zl[j] }}; 
    {%- endfor %}

    {% for j in range(end=dims.ns) %}
    zu[{{ j }}] = {{ cost.zu[j] }}; 
    {%- endfor %}

    {% for j in range(end=dims.ny) %}
    yref[{{ j }}] = {{ cost.yref[j] }}; 
    {%- endfor %}

    {% for j in range(end=dims.ny_e) %}
        {%- for k in range(end=dims.ny_e) %}
    W_e[{{ j }}+(NYN) * {{ k }}] = {{ cost.W_e[j][k] }}; 
        {%- endfor %}
    {%- endfor %}

    {% for j in range(end=dims.ny_e) %}
        {%- for k in range(end=dims.nx) %}
    Vx_e[{{ j }}+(NYN) * {{ k }}] = {{ cost.Vx_e[j][k] }}; 
        {%- endfor %}
    {%- endfor %}

    {% for j in range(end=dims.ns_e) %}
    Zl_e[{{ j }}] = {{ cost.Zl_e[j] }}; 
    {%- endfor %}

    {% for j in range(end=dims.ns_e) %}
    Zu_e[{{ j }}] = {{ cost.Zu_e[j] }}; 
    {%- endfor %}

    {% for j in range(end=dims.ns_e) %}
    zl_e[{{ j }}] = {{ cost.zl_e[j] }}; 
    {%- endfor %}

    {% for j in range(end=dims.ns_e) %}
    zu_e[{{ j }}] = {{ cost.zu_e[j] }}; 
    {%- endfor %}

    {% for j in range(end=dims.ny_e) %}
    yref_e[{{ j }}] = {{ cost.yref_e[j] }}; 
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
    int nr[N+1];
    int nr_e[N+1];

    for(int i = 0; i < N+1; i++) {
        nx[i]   = NX_;
        nu[i]   = NU_;
        nbx[i]  = NBX_;
        nbu[i]  = NBU_;
        nb[i]   = NBU_ + NBX_;
        nsbx[i] = NSBX_;
        nsbu[i] = NSBU_;
        nsh[i]  = NSH_;
        ns[i]   = NS_;
        ng[i]   = NG_;
        nh[i]   = NH_;
        nr[i]   = NR_;
        nz[i]   = NZ_;
        nv[i]   = NX_ + NU_;
        ny[i]   = NY_;
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
    nr[N]  = NRN_;
    nv[N]  = NX_; 
    ny[N]  = NYN_;
    nbu[N] = 0;

    // Make plan
    nlp_solver_plan = ocp_nlp_plan_create(N);
    {% if solver_config.nlp_solver_type == "SQP" %}
    nlp_solver_plan->nlp_solver = SQP;
    {% else %}
    nlp_solver_plan->nlp_solver = SQP_RTI;
    {% endif %}
    nlp_solver_plan->ocp_qp_solver_plan.qp_solver = {{ solver_config.qp_solver }};
    for (int i = 0; i <= N; i++)
        nlp_solver_plan->nlp_cost[i] = {{ cost.cost_type }};
    for (int i = 0; i < N; i++)
    {
        nlp_solver_plan->nlp_dynamics[i] = CONTINUOUS_MODEL;
        nlp_solver_plan->sim_solver_plan[i].sim_solver = {{ solver_config.integrator_type }};
    }

    for (int i = 0; i < N; i++) {
        {% if constraints.constr_type == "BGP" %}
        nlp_solver_plan->nlp_constraints[i] = BGP;
        {%- else -%}
        nlp_solver_plan->nlp_constraints[i] = BGH;
        {%- endif %}
    }

    {% if constraints.constr_type_e == "BGP" %}
    nlp_solver_plan->nlp_constraints[N] = BGP;
    {% else %}
    nlp_solver_plan->nlp_constraints[N] = BGH;
    {% endif %}

    {% if solver_config.hessian_approx == "EXACT" %} 
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
        ocp_nlp_dims_set_constraints(nlp_config, nlp_dims, i, "ng", &ng[i]);
        ocp_nlp_dims_set_constraints(nlp_config, nlp_dims, i, "nh", &nh[i]);
    }

    {% if constraints.constr_type == "BGP" %}
    for (int i = 0; i < N; i++) 
        ocp_nlp_dims_set_constraints(nlp_config, nlp_dims, i, "nr", &nr[i]);
    {%- endif %}
    {%- if constraints.constr_type_e == "BGP" %}
    ocp_nlp_dims_set_constraints(nlp_config, nlp_dims, N, "nr", &nr[N]);
    {%- endif %}

    {%- if constraints.constr_type == "BGP" %}
    r_constraint = (external_function_param_casadi *) malloc(sizeof(external_function_param_casadi)*N);
    for (int i = 0; i < N; ++i) {
        // nonlinear part of convex-composite constraint
        r_constraint[i].casadi_fun = &{{ con_h.name }}_r_constraint;
        r_constraint[i].casadi_n_in = &{{ con_h.name }}_r_constraint_n_in;
        r_constraint[i].casadi_n_out = &{{ con_h.name }}_r_constraint_n_out;
        r_constraint[i].casadi_sparsity_in = &{{ con_h.name }}_r_constraint_sparsity_in;
        r_constraint[i].casadi_sparsity_out = &{{ con_h.name }}_r_constraint_sparsity_out;
        r_constraint[i].casadi_work = &{{ con_h.name }}_r_constraint_work;

        external_function_param_casadi_create(&r_constraint[i], {{ dims.np }});
    }
    {%- endif %}

    {%- if constraints.constr_type_e == "BGP" %}
    // nonlinear part of convex-composite constraint
    r_e_constraint.casadi_fun = &{{ con_h_e.name }}_r_e_constraint;
    r_e_constraint.casadi_n_in = &{{ con_h_e.name }}_r_e_constraint_n_in;
    r_e_constraint.casadi_n_out = &{{ con_h_e.name }}_r_e_constraint_n_out;
    r_e_constraint.casadi_sparsity_in = &{{ con_h_e.name }}_r_e_constraint_sparsity_in;
    r_e_constraint.casadi_sparsity_out = &{{ con_h_e.name }}_r_e_constraint_sparsity_out;
    r_e_constraint.casadi_work = &{{ con_h_e.name }}_r_e_constraint_work;

    external_function_param_casadi_create(&r_e_constraint, {{ dims.np }});
    {% endif %}

    {% if dims.nh > 0 %}
    h_constraint = (external_function_param_casadi *) malloc(sizeof(external_function_param_casadi)*N);
    for (int i = 0; i < N; ++i) {
        // nonlinear constraint
        h_constraint[i].casadi_fun = &{{ con_h.name }}_h_constraint;
        h_constraint[i].casadi_n_in = &{{ con_h.name }}_h_constraint_n_in;
        h_constraint[i].casadi_n_out = &{{ con_h.name }}_h_constraint_n_out;
        h_constraint[i].casadi_sparsity_in = &{{ con_h.name }}_h_constraint_sparsity_in;
        h_constraint[i].casadi_sparsity_out = &{{ con_h.name }}_h_constraint_sparsity_out;
        h_constraint[i].casadi_work = &{{ con_h.name }}_h_constraint_work;

        external_function_param_casadi_create(&h_constraint[i], {{ dims.np }});
    }
    {% endif %}

    {% if dims.nh_e > 0 %}
	// nonlinear constraint
	h_e_constraint.casadi_fun = &{{ con_h_e.name }}_h_e_constraint;
	h_e_constraint.casadi_n_in = &{{ con_h_e.name }}_h_e_constraint_n_in;
	h_e_constraint.casadi_n_out = &{{ con_h_e.name }}_h_e_constraint_n_out;
	h_e_constraint.casadi_sparsity_in = &{{ con_h_e.name }}_h_e_constraint_sparsity_in;
	h_e_constraint.casadi_sparsity_out = &{{ con_h_e.name }}_h_e_constraint_sparsity_out;
	h_e_constraint.casadi_work = &{{ con_h_e.name }}_h_e_constraint_work;

    external_function_param_casadi_create(&h_e_constraint, {{ dims.np }});
    {% endif %}

    {% if solver_config.integrator_type == "ERK" %}
    // explicit ode
    forw_vde_casadi = (external_function_param_casadi *) malloc(sizeof(external_function_param_casadi)*N);

    for (int i = 0; i < N; ++i) {
        forw_vde_casadi[i].casadi_fun = &{{ model.name }}_expl_vde_forw;
        forw_vde_casadi[i].casadi_n_in = &{{ model.name }}_expl_vde_forw_n_in;
        forw_vde_casadi[i].casadi_n_out = &{{ model.name }}_expl_vde_forw_n_out;
        forw_vde_casadi[i].casadi_sparsity_in = &{{ model.name }}_expl_vde_forw_sparsity_in;
        forw_vde_casadi[i].casadi_sparsity_out = &{{ model.name }}_expl_vde_forw_sparsity_out;
        forw_vde_casadi[i].casadi_work = &{{ model.name }}_expl_vde_forw_work;
        external_function_param_casadi_create(&forw_vde_casadi[i], {{ dims.np }});
    }

    {% if solver_config.hessian_approx == "EXACT" %} 
    external_function_param_casadi * hess_vde_casadi;
    hess_vde_casadi = (external_function_param_casadi *) malloc(sizeof(external_function_param_casadi)*N);
    for (int i = 0; i < N; ++i) {
        hess_vde_casadi[i].casadi_fun = &{{ model.name }}_expl_ode_hess;
        hess_vde_casadi[i].casadi_n_in = &{{ model.name }}_expl_ode_hess_n_in;
        hess_vde_casadi[i].casadi_n_out = &{{ model.name }}_expl_ode_hess_n_out;
        hess_vde_casadi[i].casadi_sparsity_in = &{{ model.name }}_expl_ode_hess_sparsity_in;
        hess_vde_casadi[i].casadi_sparsity_out = &{{ model.name }}_expl_ode_hess_sparsity_out;
        hess_vde_casadi[i].casadi_work = &{{ model.name }}_expl_ode_hess_work;
        external_function_param_casadi_create(&hess_vde_casadi[i], {{ dims.np }});
    }
    {% endif %}

    {% elif solver_config.integrator_type == "IRK" %}
    // implicit dae
    impl_dae_fun = (external_function_param_casadi *) malloc(sizeof(external_function_param_casadi)*N);
    for (int i = 0; i < N; ++i) {
        impl_dae_fun[i].casadi_fun = &{{ model.name }}_impl_dae_fun;
        impl_dae_fun[i].casadi_work = &{{ model.name }}_impl_dae_fun_work;
        impl_dae_fun[i].casadi_sparsity_in = &{{ model.name }}_impl_dae_fun_sparsity_in;
        impl_dae_fun[i].casadi_sparsity_out = &{{ model.name }}_impl_dae_fun_sparsity_out;
        impl_dae_fun[i].casadi_n_in = &{{ model.name }}_impl_dae_fun_n_in;
        impl_dae_fun[i].casadi_n_out = &{{ model.name }}_impl_dae_fun_n_out;
        external_function_param_casadi_create(&impl_dae_fun[i], {{ dims.np }});
    }

    impl_dae_fun_jac_x_xdot_z = (external_function_param_casadi *) malloc(sizeof(external_function_param_casadi)*N);
    for (int i = 0; i < N; ++i) {
        impl_dae_fun_jac_x_xdot_z[i].casadi_fun = &{{ model.name }}_impl_dae_fun_jac_x_xdot_z;
        impl_dae_fun_jac_x_xdot_z[i].casadi_work = &{{ model.name }}_impl_dae_fun_jac_x_xdot_z_work;
        impl_dae_fun_jac_x_xdot_z[i].casadi_sparsity_in = &{{ model.name }}_impl_dae_fun_jac_x_xdot_z_sparsity_in;
        impl_dae_fun_jac_x_xdot_z[i].casadi_sparsity_out = &{{ model.name }}_impl_dae_fun_jac_x_xdot_z_sparsity_out;
        impl_dae_fun_jac_x_xdot_z[i].casadi_n_in = &{{ model.name }}_impl_dae_fun_jac_x_xdot_z_n_in;
        impl_dae_fun_jac_x_xdot_z[i].casadi_n_out = &{{ model.name }}_impl_dae_fun_jac_x_xdot_z_n_out;
        external_function_param_casadi_create(&impl_dae_fun_jac_x_xdot_z[i], {{ dims.np }});
    }

    impl_dae_jac_x_xdot_u_z = (external_function_param_casadi *) malloc(sizeof(external_function_param_casadi)*N);
    for (int i = 0; i < N; ++i) {
        impl_dae_jac_x_xdot_u_z[i].casadi_fun = &{{ model.name }}_impl_dae_jac_x_xdot_u_z;
        impl_dae_jac_x_xdot_u_z[i].casadi_work = &{{ model.name }}_impl_dae_jac_x_xdot_u_z_work;
        impl_dae_jac_x_xdot_u_z[i].casadi_sparsity_in = &{{ model.name }}_impl_dae_jac_x_xdot_u_z_sparsity_in;
        impl_dae_jac_x_xdot_u_z[i].casadi_sparsity_out = &{{ model.name }}_impl_dae_jac_x_xdot_u_z_sparsity_out;
        impl_dae_jac_x_xdot_u_z[i].casadi_n_in = &{{ model.name }}_impl_dae_jac_x_xdot_u_z_n_in;
        impl_dae_jac_x_xdot_u_z[i].casadi_n_out = &{{ model.name }}_impl_dae_jac_x_xdot_u_z_n_out;
        external_function_param_casadi_create(&impl_dae_jac_x_xdot_u_z[i], {{ dims.np }});
    }
    {% endif %}

    nlp_in = ocp_nlp_in_create(nlp_config, nlp_dims);

    for (int i = 0; i < N; ++i)
        nlp_in->Ts[i] = Tf/N;

    // NLP cost linear or nonlinear least squares
    {%- if cost.cost_type == "NONLINEAR_LS" %}
    r_cost = (external_function_param_casadi *) malloc(sizeof(external_function_param_casadi)*N);
    for (int i = 0; i < N; ++i) {
        // residual function
        r_cost[i].casadi_fun = &{{ cost_r.name }}_r_cost;
        r_cost[i].casadi_n_in = &{{ cost_r.name }}_r_cost_n_in;
        r_cost[i].casadi_n_out = &{{ cost_r.name }}_r_cost_n_out;
        r_cost[i].casadi_sparsity_in = &{{ cost_r.name }}_r_cost_sparsity_in;
        r_cost[i].casadi_sparsity_out = &{{ cost_r.name }}_r_cost_sparsity_out;
        r_cost[i].casadi_work = &{{ cost_r.name }}_r_cost_work;

        external_function_param_casadi_create(&r_cost[i], {{ dims.np }});
    }
    {%- endif %}

    {%- if cost.cost_type_e == "NONLINEAR_LS" %}
    // residual function
	r_e_cost.casadi_fun = &{{ cost_r_e.name }}_r_e_cost;
	r_e_cost.casadi_n_in = &{{ cost_r_e.name }}_r_e_cost_n_in;
	r_e_cost.casadi_n_out = &{{ cost_r_e.name }}_r_e_cost_n_out;
	r_e_cost.casadi_sparsity_in = &{{ cost_r_e.name }}_r_e_cost_sparsity_in;
	r_e_cost.casadi_sparsity_out = &{{ cost_r_e.name }}_r_e_cost_sparsity_out;
	r_e_cost.casadi_work = &{{ cost_r_e.name }}_r_e_cost_work;

    external_function_param_casadi_create(&r_e_cost, {{ dims.np }});
    {%- endif %}
	for (int i = 0; i < N; ++i) {
        ocp_nlp_cost_model_set(nlp_config, nlp_dims, nlp_in, i, "W", W);
        ocp_nlp_cost_model_set(nlp_config, nlp_dims, nlp_in, i, "yref", yref);
        ocp_nlp_cost_model_set(nlp_config, nlp_dims, nlp_in, i, "Zl", Zl);
        ocp_nlp_cost_model_set(nlp_config, nlp_dims, nlp_in, i, "Zu", Zu);
        ocp_nlp_cost_model_set(nlp_config, nlp_dims, nlp_in, i, "zl", zl);
        ocp_nlp_cost_model_set(nlp_config, nlp_dims, nlp_in, i, "zu", zu);
        {% if cost.cost_type == "NONLINEAR_LS" %}
        ocp_nlp_cost_model_set(nlp_config, nlp_dims, nlp_in, i, "nls_res_jac", &r_cost[i]);
        {% else %}
        ocp_nlp_cost_model_set(nlp_config, nlp_dims, nlp_in, i, "Vx", Vx);
        ocp_nlp_cost_model_set(nlp_config, nlp_dims, nlp_in, i, "Vu", Vu);
        ocp_nlp_cost_model_set(nlp_config, nlp_dims, nlp_in, i, "Vz", Vz);
        {% endif %}
	}

    ocp_nlp_cost_model_set(nlp_config, nlp_dims, nlp_in, N, "yref", yref_e);
    ocp_nlp_cost_model_set(nlp_config, nlp_dims, nlp_in, N, "Zl", Zl_e);
    ocp_nlp_cost_model_set(nlp_config, nlp_dims, nlp_in, N, "Zu", Zu_e);
    ocp_nlp_cost_model_set(nlp_config, nlp_dims, nlp_in, N, "zl", zl_e);
    ocp_nlp_cost_model_set(nlp_config, nlp_dims, nlp_in, N, "zu", zu_e);
    ocp_nlp_cost_model_set(nlp_config, nlp_dims, nlp_in, N, "W", W_e);
    {% if cost.cost_type_e == "NONLINEAR_LS" %}
    ocp_nlp_cost_model_set(nlp_config, nlp_dims, nlp_in, N, "nls_res_jac", &r_e_cost);
    {% else %}
    ocp_nlp_cost_model_set(nlp_config, nlp_dims, nlp_in, N, "Vx", Vx_e);
    {% endif %}

    // NLP dynamics
    int set_fun_status;
    for (int i = 0; i < N; ++i) {
    {% if solver_config.integrator_type == "ERK" %} 
        set_fun_status = ocp_nlp_dynamics_model_set(nlp_config, nlp_dims, nlp_in, i, "expl_vde_for", &forw_vde_casadi[i]);
        if (set_fun_status != 0) { printf("Error while setting expl_vde_for[%i]\n", i);  exit(1); }
        {% if solver_config.hessian_approx == "EXACT" %} 
            set_fun_status = ocp_nlp_dynamics_model_set(nlp_config, nlp_dims, nlp_in, i, "expl_ode_hes", &hess_vde_casadi[i]);
            if (set_fun_status != 0) { printf("Error while setting expl_ode_hes[%i]\n", i);  exit(1); }
        {% endif %}
    {% elif solver_config.integrator_type == "IRK" %} 
        set_fun_status = ocp_nlp_dynamics_model_set(nlp_config, nlp_dims, nlp_in, i, "impl_ode_fun", &impl_dae_fun[i]);
        if (set_fun_status != 0) { printf("Error while setting impl_dae_fun[%i]\n", i);  exit(1); }
        set_fun_status = ocp_nlp_dynamics_model_set(nlp_config, nlp_dims, nlp_in, i, "impl_ode_fun_jac_x_xdot", &impl_dae_fun_jac_x_xdot_z[i]);
        if (set_fun_status != 0) { printf("Error while setting impl_dae_fun_jac_x_xdot_z[%i]\n", i);  exit(1); }
        set_fun_status = ocp_nlp_dynamics_model_set(nlp_config, nlp_dims, nlp_in, i, "impl_ode_jac_x_xdot_u", &impl_dae_jac_x_xdot_u_z[i]);
        if (set_fun_status != 0) { printf("Error while setting impl_dae_jac_x_xdot_u_z[%i]\n", i);  exit(1); }
    {% endif %}
    }

    // NLP constraints
    // bounds for stage 0
    ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, 0, "idxbx", idxbx0);
    ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, 0, "lbx", lbx0);
    ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, 0, "ubx", ubx0);

    ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, 0, "idxbu", idxbu);
    ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, 0, "lbu", lbu);
    ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, 0, "ubu", ubu);

    {%- if dims.nsbx > 0 %} 
    ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, 0, "idxsbx", idxsbx);
    ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, 0, "lsbx", lsbx);
    ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, 0, "usbx", usbx);
    {%- endif %}
    
    {%- if dims.nsbu > 0 %} 
    ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, 0, "idxsbu", idxsbu);
    ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, 0, "lsbu", lsbu);
    ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, 0, "usbu", usbu);
    {%- endif %}
    
    {%- if dims.nsh > 0 %} 
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

        {%- if dims.nsbx > 0 %} 
        ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, i, "idxsbx", idxsbx);
        ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, i, "lsbx", lsbx);
        ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, i, "usbx", usbx);
        {%- endif %}
        
        {%- if dims.nsbu > 0 %} 
        ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, i, "idxsbu", idxsbu);
        ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, i, "lsbu", lsbu);
        ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, i, "usbu", usbu);
        {%- endif %}

        {%- if dims.nsh > 0 %} 
        ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, i, "idxsh", idxsh);
        ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, i, "lsh", lsh);
        ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, i, "ush", ush);
        {%- endif %}

    }
   
    {%- if dims.ng > 0 %} 
    // general constraints for stages 0 to N-1
    for (int i = 0; i < N; ++i)
    {
        ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, i, "D", D);
        ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, i, "C", C);
        ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, i, "lg", lg);
        ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, i, "ug", ug);
    }
    {%- endif %}
    
    {%- if dims.nbx_e > 0 %} 
    // bounds for last
    ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, N, "idxbx", idxbx_e);
    ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, N, "lbx", lbx_e);
    ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, N, "ubx", ubx_e);
    {%- endif %}

    {%- if dims.nsbx_e > 0 %} 
    ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, N, "idxsbx", idxsbx_e);
    ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, N, "lsbx", lsbx_e);
    ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, N, "usbx", usbx_e);
    {%- endif %}
    
    {%- if dims.nsh_e > 0 %} 
    ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, N, "idxsh", idxsh_e);
    ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, N, "lsh", lsh_e);
    ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, N, "ush", ush_e);
    {%- endif %}

    {%- if dims.ng_e > 0 %} 
    // general constraints for last stage
    ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, N, "C", C_e);
    ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, N, "lg", lg_e);
    ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, N, "ug", ug_e);
    {%- endif %}

    {%- if constraints.constr_type == "BGP" %}
    // convex-composite constraints for stages 0 to N-1
    for (int i = 0; i < N; ++i)
    {
        ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, i, "nl_constr_r_fun_jac", &r_constraint[i]);
    }
    {%- endif %}

    {%- if constraints.constr_type_e == "BGP" %}
    // convex-composite constraints for stage N
    ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, N, "nl_constr_r_fun_jac", &r_e_constraint);
    {%- endif %}

    {% if dims.nh > 0 -%}
    // nonlinear constraints for stages 0 to N-1
    for (int i = 0; i < N; ++i)
    {
        ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, i, "nl_constr_h_fun_jac", &h_constraint[i]);
        ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, i, "lh", lh);
        ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, i, "uh", uh);
    }
    {%- endif %}

    {%- if dims.nh_e > 0 %}
    // nonlinear constraints for stage N
    ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, N, "nl_constr_h_fun_jac", &h_e_constraint);
    ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, N, "lh", lh_e);
    ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, N, "uh", uh_e);
    {% endif %}

    nlp_opts = ocp_nlp_solver_opts_create(nlp_config, nlp_dims);

    {% if dims.nz > 0 -%}
    bool output_z_val = true; 
    bool sens_algebraic_val = true; 

    for (int i = 0; i < N; i++) ocp_nlp_solver_opts_set_at_stage(nlp_config, nlp_opts, i, "dynamics_output_z", &output_z_val);
    for (int i = 0; i < N; i++) ocp_nlp_solver_opts_set_at_stage(nlp_config, nlp_opts, i, "dynamics_sens_algebraic", &sens_algebraic_val);
    {% endif %}

    int num_steps_val = {{ solver_config.sim_method_num_steps }}; 
    for (int i = 0; i < N; i++) ocp_nlp_solver_opts_set_at_stage(nlp_config, nlp_opts, i, "dynamics_num_steps", &num_steps_val);

    int ns_val = {{ solver_config.sim_method_num_stages }}; 
    for (int i = 0; i < N; i++) ocp_nlp_solver_opts_set_at_stage(nlp_config, nlp_opts, i, "dynamics_num_stages", &ns_val);

    bool jac_reuse_val = true;
    for (int i = 0; i < N; i++) ocp_nlp_solver_opts_set_at_stage(nlp_config, nlp_opts, i, "dynamics_jac_reuse", &jac_reuse_val);

    {% if solver_config.nlp_solver_type == "SQP" -%}
    // set SQP specific options
    int max_iter = max_num_sqp_iterations;
    double tol_stat = 1e-6;
    double tol_eq   = 1e-6;
    double tol_ineq = 1e-6;
    double tol_comp = 1e-6;

    ocp_nlp_solver_opts_set(nlp_config, nlp_opts, "max_iter", &max_iter);
    ocp_nlp_solver_opts_set(nlp_config, nlp_opts, "tol_stat", &tol_stat);
    ocp_nlp_solver_opts_set(nlp_config, nlp_opts, "tol_eq", &tol_eq);
    ocp_nlp_solver_opts_set(nlp_config, nlp_opts, "tol_ineq", &tol_ineq);
    ocp_nlp_solver_opts_set(nlp_config, nlp_opts, "tol_comp", &tol_comp);

    {%- else -%}
    // ocp_nlp_sqp_rti_opts *sqp_opts = (ocp_nlp_sqp_rti_opts *) nlp_opts;
    {%- endif %}
    {%- if solver_config.hessian_approx == "EXACT" -%}
    for (int i = 0; i < N; ++i)
    {
        int num_steps = 5;
        bool sens_hess = true;
        bool sens_adj = true;

        ocp_nlp_solver_opts_set_at_stage(nlp_config, nlp_opts, i, "dynamics_num_steps", &num_steps);
        ocp_nlp_solver_opts_set_at_stage(nlp_config, nlp_opts, i, "dynamics_sens_hess", &sens_hess);
        ocp_nlp_solver_opts_set_at_stage(nlp_config, nlp_opts, i, "dynamics_sens_adj", &sens_adj);
    }
    {%- endif %}

    nlp_out = ocp_nlp_out_create(nlp_config, nlp_dims);

    // initialize primal solution
    double x0[{{ dims.nx }}];
    {% for item in constraints.x0 %}
    x0[{{ loop.index0 }}] = {{ item }};
    {%- endfor %}

    double u0[NU];
    {% for i in range(end=dims.nu) %}
    u0[{{ i }}] = 0.0;
    {%- endfor %}

    for (int i = 0; i <= N; ++i)
    {
        // x0
        ocp_nlp_out_set(nlp_config, nlp_dims, nlp_out, i, "x", x0);
        // u0
        ocp_nlp_out_set(nlp_config, nlp_dims, nlp_out, i, "u", u0);
    }
    
    nlp_solver = ocp_nlp_solver_create(nlp_config, nlp_dims, nlp_opts);

    {% if dims.np > 0 %}
    // initialize parameters to nominal value
    double p[{{ dims.np }}];
    {% for i in range(end=dims.np) %}
    p[{{ i }}] = {{ constraints.p[i] }};
    {%- endfor %}
    {% if solver_config.integrator_type == "IRK" %}
    for (int ii = 0; ii < {{ dims.N }}; ii++)
    {
        impl_dae_fun[ii].set_param(impl_dae_fun+ii, p);
        impl_dae_fun_jac_x_xdot_z[ii].set_param(impl_dae_fun_jac_x_xdot_z+ii, p);
        impl_dae_jac_x_xdot_u_z[ii].set_param(impl_dae_jac_x_xdot_u_z+ii, p);
    }
    {% elif solver_config.integrator_type == "ERK" %}
    for (int ii = 0; ii < {{ dims.N }}; ii++)
    {
        forw_vde_casadi[ii].set_param(forw_vde_casadi+ii, p);
    }
    {% endif %}
    for (int ii = 0; ii < {{ dims.N }}; ++ii) {
        {%- if constraints.constr_type == "BGP" %}
        r_constraint[ii].set_param(r_constraint+ii, p);
        h_constraint[ii].set_param(h_constraint+ii, p);
        {% endif %}
        {%- if dims.nh > 0 %}
        h_constraint[ii].set_param(h_constraint+ii, p);
        {% endif %}
    }
    {%- if constraints.constr_type_e == "BGP" %}
    r_e_constraint.set_param(&r_e_constraint, p);
    h_e_constraint.set_param(&h_e_constraint, p);
    {% endif %}
    {%- if dims.nh_e > 0 %}
    h_e_constraint.set_param(&h_e_constraint, p);
    {% endif %}
    {% endif %}

    return status;
}


int acados_update_params(int stage, double *p, int np) {
    int solver_status = 0;
    int casadi_np = 0;
    {% if dims.np > 0 %}
    if (stage < {{ dims.N }}) {
        {% if solver_config.integrator_type == "IRK" %}
        casadi_np = (impl_dae_fun+stage)->np;
        if (casadi_np != np) {
            printf("acados_update_params: trying to set %i parameters "
                "in impl_dae_fun which only has %i. Exiting.\n", np, casadi_np);
            exit(1);
        }
        impl_dae_fun[stage].set_param(impl_dae_fun+stage, p);
        casadi_np = (impl_dae_fun_jac_x_xdot_z+stage)->np;
        if (casadi_np != np) {
            printf("acados_update_params: trying to set %i parameters " 
                "in impl_dae_fun_jac_x_xdot_z which only has %i. Exiting.\n", np, casadi_np);
            exit(1);
        }
        impl_dae_fun_jac_x_xdot_z[stage].set_param(impl_dae_fun_jac_x_xdot_z+stage, p);
        casadi_np = (impl_dae_jac_x_xdot_u_z+stage)->np;
        if (casadi_np != np) {
            printf("acados_update_params: trying to set %i parameters " 
                "in impl_dae_jac_x_xdot_u_z which only has %i. Exiting.\n", np, casadi_np);
            exit(1);
        }
        impl_dae_jac_x_xdot_u_z[stage].set_param(impl_dae_jac_x_xdot_u_z+stage, p);
        {% else %}
        casadi_np = (forw_vde_casadi+stage)->np;
        if (casadi_np != np) {
            printf("acados_update_params: trying to set %i parameters "
                "in forw_vde_casad which only has %i. Exiting.\n", np, casadi_np);
            exit(1);
        }
        forw_vde_casadi[stage].set_param(forw_vde_casadi+stage, p);
        {% endif %}
        {% if constraints.constr_type == "BGP" %}
        casadi_np = (r_constraint+stage)->np;
        if (casadi_np != np) {
            printf("acados_update_params: trying to set %i parameters " 
                "in r_constraint which only has %i. Exiting.\n", np, casadi_np);
            exit(1);
        }
        h_constraint[stage].set_param(h_constraint+stage, p);
        casadi_np = (h_constraint+stage)->np;
        if (casadi_np != np) {
            printf("acados_update_params: trying to set %i parameters " 
                "in h_constraint which only has %i. Exiting.\n", np, casadi_np);
            exit(1);
        }
        h_constraint[stage].set_param(h_constraint+stage, p);
        {% endif %}
        {%- if dims.nh > 0 %}
        casadi_np = (h_constraint+stage)->np;
        if (casadi_np != np) {
            printf("acados_update_params: trying to set %i parameters "
                "in h_constraint which only has %i. Exiting.\n", np, casadi_np);
            exit(1);
        }
        h_constraint[stage].set_param(h_constraint+stage, p);
        {% endif %}
    } else {
    {% if constraints.constr_type_e == "BGP" %}
    casadi_np = (&r_e_constraint)->np;
    if (casadi_np != np) {
        printf("acados_update_params: trying to set %i parameters "
            "in r_e_constraint which only has %i. Exiting.\n", np, casadi_np);
        exit(1);
    }
    r_e_constraint.set_param(&r_e_constraint, p);
    casadi_np = (&h_e_constraint)->np;
    if (casadi_np != np) {
        printf("acados_update_params: trying to set %i parameters " 
            "in h_e_constraint which only has %i. Exiting.\n", np, casadi_np);
        exit(1);
    }
    h_e_constraint.set_param(&h_e_constraint, p);
    {% endif %}
    {%- if dims.nh_e > 0 %}
    casadi_np = (&h_e_constraint)->np;
    if (casadi_np != np) {
        printf("acados_update_params: trying to set %i parameters " 
            "in h_e_constraint which only has %i. Exiting.\n", np, casadi_np);
        exit(1);
    }
    h_e_constraint.set_param(&h_e_constraint, p);
    {% endif %}
    }
    {% endif %}

    return solver_status;
}



int acados_solve() {

    // solve NLP 
    int solver_status = ocp_nlp_solve(nlp_solver, nlp_in, nlp_out);

    return solver_status;
}


int acados_free() {

    // free memory
    ocp_nlp_solver_opts_destroy(nlp_opts);
    ocp_nlp_in_destroy(nlp_in);
    ocp_nlp_out_destroy(nlp_out);
    ocp_nlp_solver_destroy(nlp_solver);
    ocp_nlp_dims_destroy(nlp_dims);
    ocp_nlp_config_destroy(nlp_config);
    ocp_nlp_plan_destroy(nlp_solver_plan);

    // free external function 
    {% if solver_config.integrator_type == "IRK" %}
    for (int i = 0; i < {{ dims.N }}; i++) {
        external_function_param_casadi_free(&impl_dae_fun[i]);
        external_function_param_casadi_free(&impl_dae_fun_jac_x_xdot_z[i]);
        external_function_param_casadi_free(&impl_dae_jac_x_xdot_u_z[i]);
    }
    {% else %}
    for (int i = 0; i < {{ dims.N }}; i++)
    {
        external_function_param_casadi_free(&forw_vde_casadi[i]);
    {% if solver_config.hessian_approx == "EXACT" %}
        external_function_param_casadi_free(&hess_vde_casadi[i]);
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
