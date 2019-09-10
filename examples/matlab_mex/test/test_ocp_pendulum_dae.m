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
clear all

% check that env.sh has been run
env_run = getenv('ENV_RUN');
if (~strcmp(env_run, 'true'))
	disp('ERROR: env.sh has not been sourced! Before executing this example, run:');
	disp('source env.sh');
	return;
end

addpath('../pendulum_dae/');

%% options
compile_mex = 'true'; % true, false
codgen_model = 'true'; % true, false
% compile_mex = 'false'; % true, false
% codgen_model = 'false'; % true, false
% simulation
gnsf_detect_struct = 'true'; % true, false
sim_method = 'irk'; % irk, irk_gnsf, [erk]
sim_sens_forw = 'false'; % true, false
sim_jac_reuse = 'false'; % true, false
sim_num_stages = 3;
sim_num_steps = 3;
sim_newton_iter = 3;
model_name = 'pend_dae';

% ocp
param_scheme = 'multiple_shooting_unif_grid';
ocp_N = 50;
nlp_solver = 'sqp'; % sqp, sqp_rti
nlp_solver_exact_hessian = 'true';
nlp_solver_ext_qp_res = 1;
regularize_method = 'project_reduc_hess'; % no_regularize, project,...
    % project_reduc_hess, mirror, convexify
nlp_solver_max_iter = 100;
qp_solver = 'partial_condensing_hpipm';
        % full_condensing_hpipm, partial_condensing_hpipm
qp_solver_cond_N = 5;
qp_solver_warm_start = 0;
qp_solver_cond_ric_alg = 0; % 0: dont factorize hessian in the condensing; 1: factorize
qp_solver_ric_alg = 0; % HPIPM specific
ocp_sim_method = 'irk'; % irk, irk_gnsf
ocp_sim_method_num_stages = 4;
ocp_sim_method_num_steps = 1;
ocp_sim_method_newton_iter = 3;
cost_type = 'linear_ls'; % linear_ls, ext_cost

%% model
model = pendulum_dae_model;
%  sym_x = [xpos, ypos, alpha, vx, vy, valpha]
% x0 = [1; -5; 1; 0.1; -0.5; 0.1];
length_pendulum = 5;
xsteady = [ 0; -length_pendulum; 0; 0; 0; 0];
xtarget = [ 0; +length_pendulum; pi; 0; 0; 0];
% xtarget = xsteady;
uref = 0;

alpha0 = 0;%.01;
xp0 = length_pendulum * sin(alpha0);
yp0 = - length_pendulum * cos(alpha0);
x0 = [ xp0; yp0; alpha0; 0; 0; 0];
% x0 = xsteady + 1e-4 * ones(nx,1);

h = 0.05;
T = ocp_N*h;

disp('state')
disp(model.sym_x)

nx = length(model.sym_x);
nu = length(model.sym_u);
nz = length(model.sym_z);

ny = nu+nx; % number of outputs in lagrange term
ny_e = nx; % number of outputs in mayer term

ng = 0; % number of general linear constraints intermediate stages
ng_e = 0; % number of general linear constraints final stage
nbx = 0; % number of bounds on state x

nbu = nu; % number of bounds on controls u
nh = 0;
nh_e = 0;

% cost
% linear least square cost: y^T * W * y, where y = Vx * x + Vu * u - y_ref
Vx = eye(ny, nx); % state-to-output matrix in lagrange term
Vu = zeros(ny, nu);
Vu(nx:end, nx:end) = eye(nu); % input-to-output matrix in lagrange term
Vx_e = Vx(1:ny_e,:); % state-to-output matrix in mayer term
% weight matrix in lagrange term
W = diag([1e3 * ones(3,1);... % xpos, ypos, alpha
            ones(3,1); ... %speeds
            1e2 ]); % control
W_e = W(nu+1:nu+nx, nu+1:nu+nx); % weight matrix in mayer term
yr = [xtarget; uref]; % output reference in lagrange term
yr_e = xtarget; % output reference in mayer term

% constraints
Jbu = eye(nbu, nu);
lbu = -80*ones(nu, 1);
ubu =  80*ones(nu, 1);



%% acados ocp model
ocp_model = acados_ocp_model();
ocp_model.set('T', T);

% dims
ocp_model.set('dim_nx', nx);
ocp_model.set('dim_nu', nu);
ocp_model.set('dim_nz', nz);
if (strcmp(cost_type, 'linear_ls'))
	ocp_model.set('dim_ny', ny);
	ocp_model.set('dim_ny_e', ny_e);
end
ocp_model.set('dim_nbx', nbx);
ocp_model.set('dim_nbu', nbu);
ocp_model.set('dim_ng', ng);
ocp_model.set('dim_ng_e', ng_e);
ocp_model.set('dim_nh', nh);
ocp_model.set('dim_nh_e', nh_e);

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
if (strcmp(cost_type, 'linear_ls'))
	ocp_model.set('cost_Vu', Vu);
	ocp_model.set('cost_Vx', Vx);
	ocp_model.set('cost_Vx_e', Vx_e);
	ocp_model.set('cost_W', W);
	ocp_model.set('cost_W_e', W_e);
	ocp_model.set('cost_yr', yr);
	ocp_model.set('cost_yr_e', yr_e);
elseif (strcmp(cost_type, 'ext_cost'))
	ocp_model.set('cost_expr_ext_cost', model.expr_ext_cost);
	ocp_model.set('cost_expr_ext_cost_e', model.expr_ext_cost_e);
end

% dynamics
if (strcmp(ocp_sim_method, 'erk'))
	ocp_model.set('dyn_type', 'explicit');
	ocp_model.set('dyn_expr_f', model.expr_f_expl);
else % irk
	ocp_model.set('dyn_type', 'implicit');
	ocp_model.set('dyn_expr_f', model.expr_f_impl);
end
% constraints
ocp_model.set('constr_x0', x0);
if (ng>0)
	ocp_model.set('constr_C', C);
	ocp_model.set('constr_D', D);
	ocp_model.set('constr_lg', lg);
	ocp_model.set('constr_ug', ug);
	ocp_model.set('constr_C_e', C_e);
	ocp_model.set('constr_lg_e', lg_e);
	ocp_model.set('constr_ug_e', ug_e);
elseif (nh>0)
	ocp_model.set('constr_expr_h', model.expr_h);
	ocp_model.set('constr_lh', lbu);
	ocp_model.set('constr_uh', ubu);
%	ocp_model.set('constr_expr_h_e', model.expr_h_e);
%	ocp_model.set('constr_lh_e', lh_e);
%	ocp_model.set('constr_uh_e', uh_e);
else
%	ocp_model.set('constr_Jbx', Jbx);
%	ocp_model.set('constr_lbx', lbx);
%	ocp_model.set('constr_ubx', ubx);
	ocp_model.set('constr_Jbu', Jbu);
	ocp_model.set('constr_lbu', lbu);
	ocp_model.set('constr_ubu', ubu);
end

ocp_model.model_struct



%% acados ocp opts
ocp_opts = acados_ocp_opts();
ocp_opts.set('compile_mex', compile_mex);
ocp_opts.set('codgen_model', codgen_model);
ocp_opts.set('param_scheme', param_scheme);
ocp_opts.set('param_scheme_N', ocp_N);
ocp_opts.set('nlp_solver', nlp_solver);
ocp_opts.set('nlp_solver_exact_hessian', nlp_solver_exact_hessian);
ocp_opts.set('regularize_method', regularize_method);
if (strcmp(nlp_solver, 'sqp'))
	ocp_opts.set('nlp_solver_max_iter', nlp_solver_max_iter);
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
ocp_opts.set('nlp_solver_ext_qp_res', nlp_solver_ext_qp_res);

ocp_opts.opts_struct



%% acados ocp
% create ocp
ocp = acados_ocp(ocp_model, ocp_opts);


%% acados sim model
sim_model = acados_sim_model();
sim_model.set('name', model_name);
sim_model.set('T', h); % simulation time

sim_model.set('sym_x', model.sym_x);
if isfield(model, 'sym_u')
    sim_model.set('sym_u', model.sym_u);
end
if isfield(model, 'sym_p')
    sim_model.set('sym_p', model.sym_p);
end
sim_model.set('dim_nx', nx);
sim_model.set('dim_nu', nu);


% Note: DAEs can only be used with implicit integrator
sim_model.set('dyn_type', 'implicit');
sim_model.set('dyn_expr_f', model.expr_f_impl);
sim_model.set('sym_xdot', model.sym_xdot);
if isfield(model, 'sym_z')
	sim_model.set('sym_z', model.sym_z);
end
sim_model.set('dim_nz', nz);

%% acados sim opts
sim_opts = acados_sim_opts();
sim_opts.set('compile_mex', compile_mex);
sim_opts.set('codgen_model', codgen_model);
sim_opts.set('num_stages', sim_num_stages);
sim_opts.set('num_steps', sim_num_steps);
sim_opts.set('newton_iter', sim_newton_iter);
sim_opts.set('method', sim_method);
sim_opts.set('sens_forw', sim_sens_forw);
sim_opts.set('sens_adj', 'true');
sim_opts.set('sens_algebraic', 'true');
sim_opts.set('output_z', 'true');
sim_opts.set('sens_hess', 'false');
sim_opts.set('jac_reuse', sim_jac_reuse);
if (strcmp(sim_method, 'irk_gnsf'))
	sim_opts.set('gnsf_detect_struct', gnsf_detect_struct);
end


%% acados sim
% create integrator
sim = acados_sim(sim_model, sim_opts);

%% closed loop simulation
N_sim = 99;

x_sim = zeros(nx, N_sim+1);
x_sim(:,1) = x0; % initial state
u_sim = zeros(nu, N_sim);

% initialization
xdot0 = zeros(nx, 1);
z0 = zeros(nz, 1);

% set trajectory initialization
x_traj_init = repmat(x0, 1, ocp_N + 1);
u_traj_init = zeros(nu, ocp_N);
pi_traj_init = zeros(nx, ocp_N);
z_traj_init = 0*ones(nz, ocp_N);
xdot_traj_init = 0*ones(nx, ocp_N);
if strcmp(ocp_sim_method, 'irk_gnsf')
    gnsf_init_traj = zeros(ocp.model_struct.dim_gnsf_nout, ocp_N);
end

sqp_iter = zeros(N_sim,1);
sqp_time = zeros(N_sim,1);

tic
for ii=1:N_sim

	% set initial state
    ocp.set('constr_x0', x_sim(:,ii));
	% set trajectory initialization (if not, set internally using previous solution)
	ocp.set('init_x', x_traj_init);
	ocp.set('init_u', u_traj_init);
	ocp.set('init_pi', pi_traj_init);
    if ii == 1
        if strcmp(ocp_sim_method, 'irk')
        	ocp.set('init_z', z_traj_init);
        	ocp.set('init_xdot', xdot_traj_init);
        elseif strcmp(ocp_sim_method, 'irk_gnsf')
            ocp.set('init_gnsf_phi', gnsf_init_traj);
        end
    end

	ocp.solve();

    stat = ocp.get('stat');
    fprintf('\niter\tres_g\t\tres_b\t\tres_d\t\tres_m\t\tqp_stat\tqp_iter');
    if size(stat,2)>7
        fprintf('\tqp_res_g\tqp_res_b\tqp_res_d\tqp_res_m');
    end
    fprintf('\n');
    for jj=1:size(stat,1)
        fprintf('%d\t%e\t%e\t%e\t%e\t%d\t%d', stat(jj,1), stat(jj,2), stat(jj,3), stat(jj,4), stat(jj,5), stat(jj,6), stat(jj,7));
        if size(stat,2)>7
            fprintf('\t%e\t%e\t%e\t%e', stat(jj,8), stat(jj,9), stat(jj,10), stat(jj,11));
        end
        fprintf('\n');
    end
        
    status = ocp.get('status');
    sqp_iter(ii) = ocp.get('sqp_iter');
    sqp_time(ii) = ocp.get('time_tot');
    if status ~= 0
        keyboard
    end

	% get solution for initialization of next NLP
	x_traj = ocp.get('x');
	u_traj = ocp.get('u');
	pi_traj = ocp.get('pi');
	z_traj = ocp.get('z');

	% shift trajectory for initialization
	x_traj_init = [x_traj(:,2:end), x_traj(:,end)];
	u_traj_init = [u_traj(:,2:end), u_traj(:,end)];
	pi_traj_init = [pi_traj(:,2:end), pi_traj(:,end)];
	z_traj_init = [z_traj(:,2:end), z_traj(:,end)];

    u_sim(:,ii) = ocp.get('u', 0); % get control input
    % initialize implicit integrator
    if (strcmp(sim_method, 'irk'))
        sim.set('xdot', xdot0);
        sim.set('z', z0);
    elseif (strcmp(sim_method, 'irk_gnsf'))
        import casadi.*
        x01_gnsf = x0(sim.model_struct.dyn_gnsf_idx_perm_x(1:sim.model_struct.dim_gnsf_nx1));
        x01_dot_gnsf = xdot0(sim.model_struct.dyn_gnsf_idx_perm_x(1:sim.model_struct.dim_gnsf_nx1));
        z0_gnsf = z0(sim.model_struct.dyn_gnsf_idx_perm_z( 1:sim.model_struct.dim_gnsf_nz1 ));
        y_in = sim.model_struct.dyn_gnsf_L_x * x01_gnsf ...
                + sim.model_struct.dyn_gnsf_L_xdot * x01_dot_gnsf ...
                + sim.model_struct.dyn_gnsf_L_z * z0_gnsf;
        u_hat = sim.model_struct.dyn_gnsf_L_u * u_sim(:,ii);
        phi_fun = Function([model_name,'_gnsf_phi_fun'],...
                        {sim.model_struct.sym_gnsf_y, sim.model_struct.sym_gnsf_uhat},...
                            {sim.model_struct.dyn_gnsf_expr_phi(:)}); % sim.model_struct.sym_p

        phi_guess = full( phi_fun( y_in, u_hat ) );
        n_out = sim.model_struct.dim_gnsf_nout;
        sim.set('phi_guess', phi_guess(n_out,1));
    end

	sim.set('x', x_sim(:,ii)); 	% set initial state
	sim.set('u', u_sim(:,ii)); 	% set input
	sim.solve();	% simulate state

	% get simulated state
	x_sim(:,ii+1) = sim.get('xn');


end
format short e

toc

if 0
    % trajectory plot
    figure;
    subplot(4, 1, 1);
    plot(1:N_sim+1, x_sim(1:2,:));
    legend('xpos', 'ypos');

    subplot(4, 1, 2);
    plot(1:N_sim+1, x_sim(3,:));
    legend('alpha')

    subplot(4, 1, 3);
    plot(1:N_sim+1, x_sim(4:6,:));
    legend('vx', 'vy', 'valpha');

    subplot(4,1,4);
    plot(1:N_sim+1, [0 u_sim]);
    legend('u')

    if is_octave()
        waitforbuttonpress;
    end


    % iterations, CPU time
    figure();
    subplot(2, 1, 1);
    plot(1:N_sim, sqp_iter)
    ylabel('# iterations')
    xlabel('ocp instance')
    subplot(2, 1, 2);
    plot(1:N_sim, sqp_time)
    ylabel('CPU time')
    xlabel('ocp instance')
    if is_octave()
        waitforbuttonpress;
    end
end

% check consistency
xp = x_sim(1,:);
yp = x_sim(2,:);
check = abs(xp.^2 + yp.^2 - length_pendulum^2);

if any( max(abs(check)) > 1e-11 )
    error('note: check for constant pendulum length failed');
elseif norm( sim.get('xn') - xtarget ) > 1e-13
    error('system should have reached desired state!');
else
    disp('test_ocp_pendulum_dae: SUCCESS');
end
