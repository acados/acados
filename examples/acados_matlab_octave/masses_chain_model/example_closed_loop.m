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

%% example of closed loop simulation
clear all

% check that env.sh has been run
env_run = getenv('ENV_RUN');
if (~strcmp(env_run, 'true'))
	error('env.sh has not been sourced! Before executing this example, run: source env.sh');
end

%% handy arguments
compile_interface = 'auto';
codgen_model = 'true';
% simulation
sim_method = 'irk';
sim_sens_forw = 'false';
sim_num_stages = 4;
sim_num_steps = 4;
% ocp
ocp_N = 40;
%ocp_nlp_solver = 'sqp';
ocp_nlp_solver = 'sqp_rti';
ocp_nlp_solver_exact_hessian = 'false';
%ocp_nlp_solver_exact_hessian = 'true';
regularize_method = 'no_regularize';
%regularize_method = 'project';
%regularize_method = 'project_reduc_hess';
%regularize_method = 'mirror';
%regularize_method = 'convexify';
ocp_nlp_solver_max_iter = 100;
ocp_nlp_solver_ext_qp_res = 1;
ocp_nlp_solver_warm_start_first_qp = 1;
ocp_qp_solver = 'partial_condensing_hpipm';
%ocp_qp_solver = 'full_condensing_hpipm';
%ocp_qp_solver = 'full_condensing_qpoases';
%ocp_qp_solver = 'partial_condensing_osqp';
ocp_qp_solver_cond_N = 5;
%ocp_qp_solver_cond_N = ocp_N;
ocp_qp_solver_cond_ric_alg = 0;
ocp_qp_solver_ric_alg = 0;
ocp_qp_solver_warm_start = 1;
ocp_qp_solver_max_iter = 50;
%ocp_sim_method = 'erk';
ocp_sim_method = 'irk';
ocp_sim_method_num_stages = 4;
ocp_sim_method_num_steps = 2;
ocp_cost_type = 'linear_ls';


%% create model entries
nfm = 4;    % number of free masses
nm = nfm+1; % number of masses
model = masses_chain_model(nfm);
wall = -0.01;


% dims
T = 8.0; % horizon length time
nx = model.nx; % 6*nfm
nu = model.nu; % 3
ny = nu+nx; % number of outputs in lagrange term
ny_e = nx; % number of outputs in mayer term
nbx = nfm;
nbu = nu;
ng = 0;
nh = 0;
nh_e = 0;

% cost
Vx = zeros(ny, nx); for ii=1:nx Vx(ii,ii)=1.0; end % state-to-output matrix in lagrange term
Vu = zeros(ny, nu); for ii=1:nu Vu(nx+ii,ii)=1.0; end % input-to-output matrix in lagrange term
Vx_e = zeros(ny_e, nx); for ii=1:nx Vx_e(ii,ii)=1.0; end % state-to-output matrix in mayer term
W = 10.0*eye(ny); for ii=1:nu W(nx+ii,nx+ii)=1e-2; end % weight matrix in lagrange term
W_e = 10.0*eye(ny_e); % weight matrix in mayer term
yr = [model.x_ref; zeros(nu, 1)]; % output reference in lagrange term
yr_e = model.x_ref; % output reference in mayer term

% constraints
x0 = model.x0;
%x0 = model.x_ref;
Jbx = zeros(nbx, nx); for ii=1:nbx Jbx(ii,2+6*(ii-1))=1.0; end
lbx = wall*ones(nbx, 1);
ubx = 1e+4*ones(nbx, 1);
Jbu = zeros(nbu, nu); for ii=1:nbu Jbu(ii,ii)=1.0; end
lbu = -1.0*ones(nbu, 1);
ubu =  1.0*ones(nbu, 1);


%% acados ocp model
ocp_model = acados_ocp_model();
% dims
ocp_model.set('T', T);

% symbolics
ocp_model.set('sym_x', model.sym_x);
if isfield(model, 'sym_u')
	ocp_model.set('sym_u', model.sym_u);
end
if isfield(model, 'sym_xdot')
	ocp_model.set('sym_xdot', model.sym_xdot);
end
% cost
ocp_model.set('cost_type', ocp_cost_type);
ocp_model.set('cost_type_e', ocp_cost_type);
ocp_model.set('cost_Vu', Vu);
ocp_model.set('cost_Vx', Vx);
ocp_model.set('cost_Vx_e', Vx_e);
ocp_model.set('cost_W', W);
ocp_model.set('cost_W_e', W_e);
ocp_model.set('cost_y_ref', yr);
ocp_model.set('cost_y_ref_e', yr_e);
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
	ocp_model.set('constr_lh', lh);
	ocp_model.set('constr_uh', uh);
	ocp_model.set('constr_expr_h_e', model.expr_h_e);
	ocp_model.set('constr_lh_e', lh_e);
	ocp_model.set('constr_uh_e', uh_e);
else
	ocp_model.set('constr_Jbx', Jbx);
	ocp_model.set('constr_lbx', lbx);
	ocp_model.set('constr_ubx', ubx);
	ocp_model.set('constr_Jbu', Jbu);
	ocp_model.set('constr_lbu', lbu);
	ocp_model.set('constr_ubu', ubu);
end

%ocp_model.model_struct



%% acados ocp opts
ocp_opts = acados_ocp_opts();
ocp_opts.set('compile_interface', compile_interface);
ocp_opts.set('codgen_model', codgen_model);
ocp_opts.set('param_scheme_N', ocp_N);
ocp_opts.set('nlp_solver', ocp_nlp_solver);
ocp_opts.set('nlp_solver_exact_hessian', ocp_nlp_solver_exact_hessian);
ocp_opts.set('regularize_method', regularize_method);
ocp_opts.set('nlp_solver_ext_qp_res', ocp_nlp_solver_ext_qp_res);
ocp_opts.set('nlp_solver_warm_start_first_qp', ocp_nlp_solver_warm_start_first_qp);
if (strcmp(ocp_nlp_solver, 'sqp'))
	ocp_opts.set('nlp_solver_max_iter', ocp_nlp_solver_max_iter);
end
ocp_opts.set('qp_solver', ocp_qp_solver);
ocp_opts.set('qp_solver_iter_max', ocp_qp_solver_max_iter);
ocp_opts.set('qp_solver_warm_start', ocp_qp_solver_warm_start);
ocp_opts.set('qp_solver_cond_ric_alg', ocp_qp_solver_cond_ric_alg);
if (~isempty(strfind(ocp_qp_solver, 'partial_condensing')))
	ocp_opts.set('qp_solver_cond_N', ocp_qp_solver_cond_N);
end
if (strcmp(ocp_qp_solver, 'partial_condensing_hpipm'))
	ocp_opts.set('qp_solver_ric_alg', ocp_qp_solver_ric_alg);
end
ocp_opts.set('sim_method', ocp_sim_method);
ocp_opts.set('sim_method_num_stages', ocp_sim_method_num_stages);
ocp_opts.set('sim_method_num_steps', ocp_sim_method_num_steps);


%% acados ocp
% create ocp
ocp = acados_ocp(ocp_model, ocp_opts);


%% acados sim model
sim_model = acados_sim_model();
% symbolics
sim_model.set('sym_x', model.sym_x);
if isfield(model, 'sym_u')
	sim_model.set('sym_u', model.sym_u);
end
if isfield(model, 'sym_xdot')
	sim_model.set('sym_xdot', model.sym_xdot);
end
% model
sim_model.set('T', T/ocp_N);
if (strcmp(sim_method, 'erk'))
	sim_model.set('dyn_type', 'explicit');
	sim_model.set('dyn_expr_f', model.expr_f_expl);
else % irk
	sim_model.set('dyn_type', 'implicit');
	sim_model.set('dyn_expr_f', model.expr_f_impl);
end



%% acados sim opts
sim_opts = acados_sim_opts();
sim_opts.set('compile_interface', compile_interface);
sim_opts.set('codgen_model', codgen_model);
sim_opts.set('num_stages', sim_num_stages);
sim_opts.set('num_steps', sim_num_steps);
sim_opts.set('method', sim_method);
sim_opts.set('sens_forw', sim_sens_forw);

%% acados sim
% create sim
sim = acados_sim(sim_model, sim_opts);


%% closed loop simulation
n_sim = 50;
x_sim = zeros(nx, n_sim+1);
x_sim(:,1) = x0; % initial state
u_sim = zeros(nu, n_sim);

x_traj_init = repmat(model.x_ref, 1, ocp_N+1);
u_traj_init = zeros(nu, ocp_N);
pi_traj_init = zeros(nx, ocp_N);

%ocp.set('init_x', x_traj_init);
%ocp.set('init_u', u_traj_init);
%ocp.set('init_pi', pi_traj_init);

tic;

for ii=1:n_sim

	% set x0
	ocp.set('constr_x0', x_sim(:,ii));

	% set trajectory initialization (if not, set internally using previous solution)
	ocp.set('init_x', x_traj_init);
	ocp.set('init_u', u_traj_init);
	ocp.set('init_pi', pi_traj_init);

	% solve OCP
	ocp.set('rti_phase', 1);
	ocp.solve();
	ocp.set('rti_phase', 2);
	ocp.solve();

	if 1
		status = ocp.get('status');
		sqp_iter = ocp.get('sqp_iter');
		time_tot = ocp.get('time_tot');
		time_lin = ocp.get('time_lin');
		time_reg = ocp.get('time_reg');
		time_qp_sol = ocp.get('time_qp_sol');
		time_qp_solver_call = ocp.get('time_qp_solver_call');
		qp_iter = ocp.get('qp_iter');

		fprintf('\nstatus = %d, sqp_iter = %d, time_int = %f [ms] (time_lin = %f [ms], time_qp_sol = %f [ms] (time_qp_solver_call = %f [ms]), time_reg = %f [ms])\n', status, sqp_iter, time_tot*1e3, time_lin*1e3, time_qp_sol*1e3, time_qp_solver_call*1e3, time_reg*1e3);
%		fprintf('%e %d\n', time_qp_solver_call, qp_iter);

%		ocp.print('stat');
	end

	% get solution
	x_traj = ocp.get('x');
	u_traj = ocp.get('u');
	pi_traj = ocp.get('pi');

	% shift trajectory for initialization
	x_traj_init = [x_traj(:,2:end), x_traj(:,end)];
	u_traj_init = [u_traj(:,2:end), u_traj(:,end)];
	pi_traj_init = [pi_traj(:,2:end), pi_traj(:,end)];

	% get solution for sim
	u_sim(:,ii) = ocp.get('u', 0);

	% overwrite control to perturb the system
%	if(ii<=5)
%		u_sim(:,ii) = [-1; 1; 1];
%	end

	% set initial state of sim
	sim.set('x', x_sim(:,ii));
	% set input in sim
	sim.set('u', u_sim(:,ii));

	% simulate state
	sim.solve();

	% get new state
	x_sim(:,ii+1) = sim.get('xn');

end

avg_time_solve = toc/n_sim


u_sim;
x_sim;

% print solution
for ii=1:n_sim+1
	cur_pos = x_sim(:,ii);
	visualize;
end



status = ocp.get('status');

if status==0
	fprintf('\nsuccess!\n\n');
else
	fprintf('\nsolution failed!\n\n');
end


if is_octave()
    waitforbuttonpress;
end
