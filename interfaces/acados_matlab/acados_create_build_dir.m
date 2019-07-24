function acados_create_build_dir

% clear mex functions (if loaded from previous build)
clear ocp_create
clear ocp_create_ext_fun
clear ocp_destroy
clear ocp_destroy_ext_fun
clear ocp_get
clear ocp_model_ocp_set_ext_fun_cost_0_ext_cost_jac_hes
clear ocp_model_ocp_set_ext_fun_cost_1_ext_cost_jac_hes
clear ocp_model_ocp_set_ext_fun_dyn_0_expl_ode_fun
clear ocp_model_ocp_set_ext_fun_dyn_0_expl_ode_hes
clear ocp_model_ocp_set_ext_fun_dyn_0_expl_vde_adj
clear ocp_model_ocp_set_ext_fun_dyn_0_expl_vde_for
clear ocp_precompute
clear ocp_set
clear ocp_solve

clear sim_create
clear sim_create_ext_fun
clear sim_destroy_ext_fun
clear sim_solve
clear sim_set
clear sim_get
clear sim_precompute

% create build folder and add to path
addpath('build');
rmpath('build')
[~] = rmdir('build', 's');
[~,~] = mkdir('build');
addpath('build');