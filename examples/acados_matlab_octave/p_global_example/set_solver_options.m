%
% Copyright (c) The acados authors.
%
% This file is part of acados.
%
% Licensed under the 2-Clause BSD License.



function ocp = set_solver_options(ocp)
    % set options
    ocp.solver_options.qp_solver = 'PARTIAL_CONDENSING_HPIPM';
    ocp.solver_options.hessian_approx = 'GAUSS_NEWTON'; %'GAUSS_NEWTON'; %'EXACT'; %
    ocp.solver_options.integrator_type = 'ERK';
    ocp.solver_options.print_level = 0;
    ocp.solver_options.nlp_solver_type = 'SQP_RTI';

    % set prediction horizon
    Tf = 1.0;
    N_horizon = 20;
    ocp.solver_options.tf = Tf;
    ocp.solver_options.N_horizon = N_horizon;

    % partial condensing
    ocp.solver_options.qp_solver_cond_N = 5;
    ocp.solver_options.qp_solver_cond_block_size = [3, 3, 3, 3, 7, 1];

    % NOTE: these additional flags are required for code generation of CasADi functions using casadi.blazing_spline
    % These might be different depending on your compiler and operating system.
    flags = ['-I' casadi.GlobalOptions.getCasadiIncludePath ' -O2 -ffast-math -march=native -fno-omit-frame-pointer'];
    ocp.code_gen_options.ext_fun_compile_flags = flags;
    ocp.code_gen_options.ext_fun_expand_dyn = true;
    ocp.code_gen_options.ext_fun_expand_constr = true;
    ocp.code_gen_options.ext_fun_expand_cost = true;
    ocp.code_gen_options.ext_fun_expand_precompute = false;
end
