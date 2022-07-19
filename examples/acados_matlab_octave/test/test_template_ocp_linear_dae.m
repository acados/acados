%
% Copyright 2019 Gianluca Frison, Dimitris Kouzoupis, Robin Verschueren,
% Andrea Zanelli, Niels van Duijkeren, Jonathan Frey, Tommaso Sartor,
% Branimir Novoselnik, Rien Quirynen, Rezart Qelibari, Dang Doan,
% Jonas Koenemann, Yutao Chen, Tobias Schöls, Jonas Schlagenhauf, Moritz Diehl
%
% This file is part of acados.
%
% The 2-Clause BSD License
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
%
% 1. Redistributions of source code must retain the above copyright notice,
% this list of conditions and the following disclaimer.
%
% 2. Redistributions in binary form must reproduce the above copyright notice,
% this list of conditions and the following disclaimer in the documentation
% and/or other materials provided with the distribution.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.;
%

%% test of native matlab interface

addpath('../simple_dae_model/');

for itest = 1:2
    %% options
    compile_interface = 'auto'; % true, false
    codgen_model = 'true'; % true, false

    % ocp
    N = 20;
    nlp_solver = 'sqp'; % sqp, sqp_rti
    nlp_solver_exact_hessian = 'false';
    %nlp_solver_exact_hessian = 'true';
    regularize_method = 'no_regularize';
    %regularize_method = 'project_reduc_hess';
    nlp_solver_max_iter = 100;
    nlp_solver_tol_stat = 1e-12;
    nlp_solver_tol_eq   = 1e-12;
    nlp_solver_tol_ineq = 1e-12;
    nlp_solver_tol_comp = 1e-12;
    nlp_solver_ext_qp_res = 1;
    % qp_solver = 'full_condensing_qpoases'; % partial_condensing_hpipm
    qp_solver = 'partial_condensing_hpipm'; % partial_condensing_hpipm
    qp_solver_cond_N = 5;
    qp_solver_warm_start = 0;
    qp_solver_cond_ric_alg = 1; % 0: dont factorize hessian in the condensing; 1: factorize
    qp_solver_ric_alg = 1; % HPIPM specific
    if itest == 1
        ocp_sim_method = 'irk_gnsf'; % irk, irk_gnsf
        constr_variant = 1; % 0: x bounds; 1: z bounds
    else
        ocp_sim_method = 'irk'; % irk, irk_gnsf
        constr_variant = 1; % 0: x bounds; 1: z bounds
    end
    ocp_sim_method_num_stages = 6; % scalar or vector of size ocp_N;
    ocp_sim_method_num_steps = 4; % scalar or vector of size ocp_N;
    ocp_sim_method_newton_iter = 3; % scalar or vector of size ocp_N;

    % selectors for example variants
    cost_variant = 1; % 0: ls on u,x; 1: ls on u,z; (not implemented yet: 2: nls on u,z)

    name = ['toy_dae' num2str(itest)];

    % get model
    model = simple_dae_model;

    nx = length(model.sym_x);
    nu = length(model.sym_u);
    nz = length(model.sym_z);
    ny = nx+nu;
    ny_e = nx;

    T = 1.0;
    h = T/N;

    Wu = 1e-3*eye(nu);
    Wx = 1e1*eye(nx);
    W = [Wu, zeros(nu,nx); zeros(nx,nu), Wx];
    Vu = [eye(nu); zeros(nx,nu)];
    Vx = [zeros(nu,nx); eye(nx)];
    Vx_e = eye(nx);

    lb = [2; -2];
    ub = [4;  2];

    x0 = [3; -1.8];

    %% acados ocp model
    ocp_model = acados_ocp_model();
    ocp_model.set('T', T);

    % symbolics
    ocp_model.set('sym_x', model.sym_x);
    ocp_model.set('sym_u', model.sym_u);
    ocp_model.set('sym_xdot', model.sym_xdot);
    ocp_model.set('sym_z', model.sym_z);

    % cost
    if cost_variant==0
        ocp_model.set('cost_type', 'linear_ls');
        ocp_model.set('cost_Vu', Vu);
        ocp_model.set('cost_Vx', Vx);
    elseif cost_variant==1
        ocp_model.set('cost_type', 'linear_ls');
        ocp_model.set('cost_Vu', Vu);
        ocp_model.set('cost_Vx', zeros(ny,nx));
        ocp_model.set('cost_Vz', Vx);
    else
        ocp_model.set('cost_type', 'nonlinear_ls');
        ocp_model.set('cost_expr_y', model.expr_y);
    end
    ocp_model.set('cost_type_e', 'linear_ls');
    ocp_model.set('cost_Vx_e', Vx_e);
    ocp_model.set('cost_W', W);
    ocp_model.set('cost_W_e', Wx);
    %ocp_model.set('cost_y_ref', yr);
    %ocp_model.set('cost_y_ref_e', yr_e);

    % dynamics
    ocp_model.set('dyn_type', 'implicit');
    ocp_model.set('dyn_expr_f', model.expr_f_impl);
    ocp_model.set('name', name);

    % constraints
    ocp_model.set('constr_x0', x0);
    if constr_variant==0
        ocp_model.set('constr_Jbx', eye(nx));
        ocp_model.set('constr_lbx', lb);
        ocp_model.set('constr_ubx', ub);
    else
        ocp_model.set('constr_expr_h', model.expr_h);
        ocp_model.set('constr_lh', lb);
        ocp_model.set('constr_uh', ub);
        ocp_model.set('constr_expr_h_e', model.expr_h_e);
        ocp_model.set('constr_lh_e', lb);
        ocp_model.set('constr_uh_e', ub);
    end



    %% acados ocp opts
    ocp_opts = acados_ocp_opts();

    ocp_opts.set('compile_interface', compile_interface);
    ocp_opts.set('codgen_model', codgen_model);
    ocp_opts.set('param_scheme_N', N);
    ocp_opts.set('nlp_solver', nlp_solver);
    ocp_opts.set('nlp_solver_exact_hessian', nlp_solver_exact_hessian);
    ocp_opts.set('regularize_method', regularize_method);
    ocp_opts.set('nlp_solver_ext_qp_res', nlp_solver_ext_qp_res);
    if (strcmp(nlp_solver, 'sqp'))
        ocp_opts.set('nlp_solver_max_iter', nlp_solver_max_iter);
        ocp_opts.set('nlp_solver_tol_stat', nlp_solver_tol_stat);
        ocp_opts.set('nlp_solver_tol_eq', nlp_solver_tol_eq);
        ocp_opts.set('nlp_solver_tol_ineq', nlp_solver_tol_ineq);
        ocp_opts.set('nlp_solver_tol_comp', nlp_solver_tol_comp);
    end
    ocp_opts.set('qp_solver', qp_solver);
    if (strcmp(qp_solver, 'partial_condensing_hpipm'))
        ocp_opts.set('qp_solver_cond_N', qp_solver_cond_N);
        ocp_opts.set('qp_solver_cond_ric_alg', qp_solver_cond_ric_alg);
        ocp_opts.set('qp_solver_ric_alg', qp_solver_ric_alg);
        ocp_opts.set('qp_solver_warm_start', qp_solver_warm_start);
    end
    ocp_opts.set('sim_method', ocp_sim_method);
    ocp_opts.set('sim_method_num_stages', ocp_sim_method_num_stages);
    ocp_opts.set('sim_method_num_steps', ocp_sim_method_num_steps);
    ocp_opts.set('sim_method_newton_iter', ocp_sim_method_newton_iter);
    ocp_opts.set('sim_method_jac_reuse', [0; ones(N-1, 1)]);
    ocp_opts.set('ext_fun_compile_flags', '');


    %% acados ocp
    ocp = acados_ocp(ocp_model, ocp_opts);

    ocp.solve();
    ocp.print('stat')

    status = ocp.get('status');
    sqp_iter = ocp.get('sqp_iter');
    sqp_time = ocp.get('time_tot');


    format short e
    % get solution for initialization of next NLP
    x_traj = ocp.get('x');
    u_traj = ocp.get('u');
    pi_traj = ocp.get('pi');
    z_traj = ocp.get('z');

    diff_x_z = x_traj(:,1:N) - z_traj;
    tol_diff_xz = 1e-15;

    if status ~= 0
        error('test_ocp_linear_dae: ocp_nlp solver exit with nonzero status');
    elseif max(abs(diff_x_z)) > tol_diff_xz
        error(['test_ocp_linear_dae: difference between x and z bigger than',
            num2str(tol_diff_xz, '%e'), ' should be equal'])
    end
    cost_val_ocp = ocp.get_cost();
end % itest

fprintf('\ntest_ocp_linear_dae: success!\n');

% create templated mex solver object: t_ocp
ocp.generate_c_code;
cd c_generated_code/
command = strcat('t_ocp = ', name, '_mex_solver');
eval( command );

t_ocp.set('constr_x0', x0);

t_ocp.solve();
t_ocp.print;
t_x = t_ocp.get('x');
t_u = t_ocp.get('u');
t_z = t_ocp.get('z');

% test setting parameters
t_ocp.set('p',[]);

%
err_x = max(max(abs(x_traj - t_x)))
err_u = max(max(abs(u_traj - t_u)))
err_z = max(max(abs(z_traj - t_z)))

cost_val_t_ocp = t_ocp.get_cost();

if any([err_x, err_u, err_z] > 1e-9)
    error(['test_template_ocp_linear_dae: solution of templated MEX and original MEX',...
         ' differ too much. Should be < 1e-9 ']);
end

if abs(cost_val_ocp - cost_val_t_ocp) > 1e-9
    error(['test_template_ocp_linear_dae: cost function value of templated MEX and original MEX',...
         ' differ too much. Should be < 1e-9, got ', num2str(cost_val_ocp, '%e'), ' template: ', ...
         num2str(cost_val_t_ocp, '%e')]);
else
    disp(['Cost value: ', num2str(cost_val_ocp, '%e'), ' template: ', num2str(cost_val_t_ocp, '%e')]);
end

cd ..
