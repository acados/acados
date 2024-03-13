%
% Copyright (c) The acados authors.
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
clear all

% check that env.sh has been run
env_run = getenv('ENV_RUN');
if (~strcmp(env_run, 'true'))
	error('env.sh has not been sourced! Before executing this example, run: source env.sh');
end
import casadi.*

tic

N = 40;

ncases = 3;
constr_violation = zeros(1, ncases);
constr_vals = zeros(N, ncases);
for i = 1:3
    %% arguments
    compile_interface = 'auto';

    codgen_model = 'true';
    gnsf_detect_struct = 'true';

    % discretization
    h = 0.02;

    nlp_solver = 'sqp';
    %nlp_solver = 'sqp_rti';
    nlp_solver_exact_hessian = 'false';
    regularize_method = 'no_regularize';
    nlp_solver_max_iter = 100;
    tol = 1e-12;
    nlp_solver_tol_stat = tol;
    nlp_solver_tol_eq   = tol;
    nlp_solver_tol_ineq = tol;
    nlp_solver_tol_comp = tol;
    nlp_solver_ext_qp_res = 1;
    qp_solver = 'partial_condensing_hpipm';
%     qp_solver = 'full_condensing_hpipm';
%     qp_solver = 'full_condensing_qpoases';
    qp_solver_cond_N = 5;
    qp_solver_cond_ric_alg = 0;
    qp_solver_ric_alg = 0;
    qp_solver_warm_start = 1;
    sim_method = 'irk';
    sim_method_num_stages = 1;
    sim_method_num_steps = 1;
    sim_method_exact_z_output = 0;

    % cost_type = 'linear_ls';
    cost_type = 'ext_cost';
    model_name = ['ocp_pendulum_' num2str(i)];


    %% create model entries
    switch i
        case 1
            model = pendulum_on_cart_model();
            theta = model.sym_x(2);
            omega = model.sym_x(4);
            model.sym_z = [];
            model.expr_h = cos(theta)*sin(theta)*omega.^2;
            lh = -40;
            uh = 40;
        case {2,3}
            model = pendulum_on_cart_model_dae;
            model.expr_h = model.sym_z;
            lh = -40;
            uh = 40;
            if i == 2
                sim_method_exact_z_output = 1;
            end
    end

    % dims
    T = N*h; % horizon length time
    nx = length(model.sym_x);
    nu = length(model.sym_u);
    if isfield(model, 'sym_z')
        nz = length(model.sym_z);
    else
        nz = 0;
    end

    % constraints
    x0 = [0; pi; 0; 0];
    nbu = nu;
    Jbu = zeros(nbu, nu); for ii=1:nbu; Jbu(ii,ii)=1.0; end
    lbu = -80*ones(nu, 1);
    ubu =  80*ones(nu, 1);


    %% acados ocp model
    ocp_model = acados_ocp_model();
    ocp_model.set('name', model_name);
    ocp_model.set('T', T);
    
    % symbolics
    ocp_model.set('sym_x', model.sym_x);
    if isfield(model, 'sym_u')
        ocp_model.set('sym_u', model.sym_u);
    end
    if isfield(model, 'sym_xdot')
        ocp_model.set('sym_xdot', model.sym_xdot);
    end
    if isfield(model, 'sym_z')
        ocp_model.set('sym_z', model.sym_z);
    end
    % cost
    ocp_model.set('cost_type', cost_type);
    ocp_model.set('cost_type_e', cost_type);

    ocp_model.set('cost_expr_ext_cost', model.expr_ext_cost);
    ocp_model.set('cost_expr_ext_cost_e', model.expr_ext_cost_e);

    % dynamics
    if (strcmp(sim_method, 'erk'))
        ocp_model.set('dyn_type', 'explicit');
        ocp_model.set('dyn_expr_f', model.expr_f_expl);
    else % irk irk_gnsf
        ocp_model.set('dyn_type', 'implicit');
        ocp_model.set('dyn_expr_f', model.expr_f_impl);
    end
  
    % constraints
    ocp_model.set('constr_x0', x0);
    
    nh = length(model.expr_h);
    ocp_model.set('constr_expr_h', model.expr_h);
    ocp_model.set('constr_lh', lh);
    ocp_model.set('constr_uh', uh);

    ocp_model.set('constr_Jbu', Jbu);
    ocp_model.set('constr_lbu', lbu);
    ocp_model.set('constr_ubu', ubu);
%     disp('ocp_model.model_struct')
%     disp(ocp_model.model_struct)


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
        ocp_opts.set('qp_solver_ric_alg', qp_solver_ric_alg);
    end
    ocp_opts.set('qp_solver_cond_ric_alg', qp_solver_cond_ric_alg);
    ocp_opts.set('qp_solver_warm_start', qp_solver_warm_start);
    ocp_opts.set('sim_method', sim_method);
    ocp_opts.set('sim_method_num_stages', sim_method_num_stages);
    ocp_opts.set('sim_method_num_steps', sim_method_num_steps);
    ocp_opts.set('sim_method_exact_z_output', sim_method_exact_z_output);
    
    
    if (strcmp(sim_method, 'irk_gnsf'))
        ocp_opts.set('gnsf_detect_struct', gnsf_detect_struct);
    end

    disp('ocp_opts');
    disp(ocp_opts.opts_struct);


    %% acados ocp
    ocp = acados_ocp(ocp_model, ocp_opts);

    % set trajectory initialization
    x_traj_init = [linspace(0, 0, N+1); linspace(pi, 0, N+1); linspace(0, 0, N+1); linspace(0, 0, N+1)];
    u_traj_init = zeros(nu, N);

    ocp.set('init_x', x_traj_init);
    ocp.set('init_u', u_traj_init);

    ocp.solve();

    % get solution
    u = ocp.get('u');
    x = ocp.get('x');

    %% evaluation
    status = ocp.get('status');
    sqp_iter = ocp.get('sqp_iter');
    time_tot = ocp.get('time_tot');
    time_lin = ocp.get('time_lin');
    time_reg = ocp.get('time_reg');
    time_qp_sol = ocp.get('time_qp_sol');

    fprintf(['\nstatus = %d, sqp_iter = %d, time_int = %f [ms]'...
        ' (time_lin = %f [ms], time_qp_sol = %f [ms], time_reg = %f [ms])\n'],...
        status, sqp_iter, time_tot*1e3, time_lin*1e3, time_qp_sol*1e3, time_reg*1e3 );

    ocp.print('stat');
    
    if i == 1
        % save reference
        xref = x;
        uref = u;
    else
        % compare error w.r.t reference
        err_x(i) = norm(x - xref)
        err_u(i) = norm(u - uref)
    end
    
    % compare z accuracy to respective value obtained from x
    thetas = xref(2,:);
    omegas = xref(4,:);
    z_fromx = cos(thetas).*sin(thetas).*omegas.^2;
    if i > 1
        z = ocp.get('z');
        err_z_zfromx(i) = norm( z - z_fromx(1:end-1) );
    end
    
    % check constraint violation
    theta = model.sym_x(2);
    omega = model.sym_x(4);
    constr_expr = cos(theta)*sin(theta)*omega.^2;
    constr_fun = Function('constr_fun', {model.sym_x, model.sym_u, model.sym_z}, ...
        {constr_expr});

    constr_violation(i) = 0;
    for j=1:N
        valh = full( constr_fun(x(:,j), u(:,j), z_fromx(:,j) ) );
        constr_vals(j,i) = valh;
        violation = max([0, -(valh-lh), valh - uh]);
        constr_violation(i) = max(norm(violation), constr_violation(i));
    end
end

toc

fprintf('\nConstraint values\n');
fprintf('\nODE \t\tDAE-exact_z\tDAE-extrapolation\n');
for i = 1:N
    fprintf('%.4e\t%.4e\t%.4e\n', constr_vals(i,1), constr_vals(i,2),...
        constr_vals(i,3))
end
fprintf('\n\nConstraint violations');
fprintf('\nODE \t\tDAE-exact_z\tDAE-extrapolation\n');
fprintf('%.4e\t%.4e\t%.4e\n', constr_violation(1), constr_violation(2),...
    constr_violation(3))
