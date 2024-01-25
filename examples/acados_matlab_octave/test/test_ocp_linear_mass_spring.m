
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

addpath('../linear_mass_spring_model/');

%% arguments
compile_interface = 'auto';
codgen_model = 'true';
N = 20;
shooting_nodes = [ linspace(0,1,N/2) linspace(1.1,5,N/2+1) ];

model_name = 'lin_mass';

nlp_solver = 'sqp';
%nlp_solver = 'sqp_rti';
nlp_solver_exact_hessian = 'false';
%nlp_solver_exact_hessian = 'true';
% regularize_method = 'no_regularize';
%regularize_method = 'project';
%regularize_method = 'mirror';
regularize_method = 'convexify';
nlp_solver_max_iter = 100;
nlp_solver_ext_qp_res = 1;
qp_solver = 'partial_condensing_hpipm';
%qp_solver = 'full_condensing_hpipm';
qp_solver_cond_N = 5;
% sim_method = 'irk';
%sim_method = 'irk_gnsf';
sim_method = 'discrete';
sim_method_num_stages = 4 * ones(N,1);
sim_method_num_stages(end) = 2;
sim_method_num_steps = 3;
%cost_type = 'linear_ls';
%cost_type = 'nonlinear_ls';
cost_type = 'ext_cost';

if strcmp(sim_method, 'erk')
    dyn_type = 'explicit';
elseif any(strcmp(sim_method, {'irk','irk_gnsf'}))
	dyn_type = 'implicit';
else
	dyn_type = 'discrete';
end


%% create model entries
model = linear_mass_spring_model();


% dims
T = 10.0; % horizon length time
nx = model.nx;
nu = model.nu;
ny = nu+nx; % number of outputs in lagrange term
ny_e = nx; % number of outputs in mayer term
nbx = 0;
nbu = 0;
ng = 0;
nh = nu+nx;
nh_e = nx;

% cost
Vu = zeros(ny, nu); for ii=1:nu Vu(ii,ii)=1.0; end % input-to-output matrix in lagrange term
Vx = zeros(ny, nx); for ii=1:nx Vx(nu+ii,ii)=1.0; end % state-to-output matrix in lagrange term
Vx_e = zeros(ny_e, nx); for ii=1:nx Vx_e(ii,ii)=1.0; end % state-to-output matrix in mayer term
W = eye(ny); for ii=1:nu W(ii,ii)=2.0; end % weight matrix in lagrange term
W_e = eye(ny_e); % weight matrix in mayer term
yr = zeros(ny, 1); % output reference in lagrange term
yr_e = zeros(ny_e, 1); % output reference in mayer term
% constraints
x0 = zeros(nx, 1); x0(1)=2.5; x0(2)=2.5;
if (ng>0)
	D = zeros(ng, nu); for ii=1:nu D(ii,ii)=1.0; end
	C = zeros(ng, nx); for ii=1:ng-nu C(nu+ii,ii)=1.0; end
	lg = zeros(ng, 1); for ii=1:nu lg(ii)=-0.5; end; for ii=1:ng-nu lg(nu+ii)=-4.0; end
	ug = zeros(ng, 1); for ii=1:nu ug(ii)= 0.5; end; for ii=1:ng-nu ug(nu+ii)= 4.0; end
	C_e = zeros(ng_e, nx); for ii=1:ng_e C_e(ii,ii)=1.0; end
	lg_e = zeros(ng_e, 1); for ii=1:ng_e lg_e(ii)=-4.0; end
	ug_e = zeros(ng_e, 1); for ii=1:ng_e ug_e(ii)= 4.0; end
elseif (nh>0)
	lh = zeros(nh, 1); for ii=1:nu lh(ii)=-0.5; end; for ii=1:nx lh(nu+ii)=-4.0; end
	uh = zeros(nh, 1); for ii=1:nu uh(ii)= 0.5; end; for ii=1:nx uh(nu+ii)= 4.0; end
	lh_e = zeros(nh_e, 1); for ii=1:nh_e lh_e(ii)=-4.0; end
	uh_e = zeros(nh_e, 1); for ii=1:nh_e uh_e(ii)= 4.0; end
else
	Jbx = zeros(nbx, nx); for ii=1:nbx Jbx(ii,ii)=1.0; end
	lbx = -4*ones(nbx, 1);
	ubx =  4*ones(nbx, 1);
	Jbu = zeros(nbu, nu); for ii=1:nbu Jbu(ii,ii)=1.0; end
	lbu = -0.5*ones(nbu, 1);
	ubu =  0.5*ones(nbu, 1);
end


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
% cost
ocp_model.set('cost_type', cost_type);
ocp_model.set('cost_type_e', cost_type);
if (strcmp(cost_type, 'linear_ls'))
	ocp_model.set('cost_Vu', Vu);
	ocp_model.set('cost_Vx', Vx);
	ocp_model.set('cost_Vx_e', Vx_e);
	ocp_model.set('cost_W', W);
	ocp_model.set('cost_W_e', W_e);
	ocp_model.set('cost_y_ref', yr);
	ocp_model.set('cost_y_ref_e', yr_e);
elseif (strcmp(cost_type, 'nonlinear_ls'))
	ocp_model.set('cost_expr_y', model.expr_y);
	ocp_model.set('cost_expr_y_e', model.expr_y_e);
	ocp_model.set('cost_W', W);
	ocp_model.set('cost_W_e', W_e);
	ocp_model.set('cost_y_ref', yr);
	ocp_model.set('cost_y_ref_e', yr_e);
else % if (strcmp(cost_type, 'ext_cost'))
	ocp_model.set('cost_expr_ext_cost', model.expr_ext_cost);
	ocp_model.set('cost_expr_ext_cost_e', model.expr_ext_cost_e);
end
% dynamics
if (strcmp(dyn_type, 'explicit'))
	ocp_model.set('dyn_type', 'explicit');
	ocp_model.set('dyn_expr_f', model.expr_f_expl);
elseif (strcmp(dyn_type, 'implicit'))
	ocp_model.set('dyn_type', 'implicit');
	ocp_model.set('dyn_expr_f', model.expr_f_impl);
else
	ocp_model.set('dyn_type', 'discrete');
	ocp_model.set('dyn_expr_phi', model.expr_phi);
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


%% acados ocp opts
ocp_opts = acados_ocp_opts();
ocp_opts.set('compile_interface', compile_interface);
ocp_opts.set('codgen_model', codgen_model);
ocp_opts.set('param_scheme_N', N);
ocp_opts.set('ext_fun_compile_flags', '');
if (exist('shooting_nodes', 'var'))
	ocp_opts.set('shooting_nodes', shooting_nodes);
end
ocp_opts.set('nlp_solver', nlp_solver);
ocp_opts.set('nlp_solver_exact_hessian', nlp_solver_exact_hessian);
ocp_opts.set('regularize_method', regularize_method);
ocp_opts.set('nlp_solver_ext_qp_res', nlp_solver_ext_qp_res);
ocp_opts.set('nlp_solver_max_iter', nlp_solver_max_iter);
ocp_opts.set('qp_solver', qp_solver);
ocp_opts.set('sim_method', sim_method);
if (strcmp(qp_solver, 'partial_condensing_hpipm'))
	ocp_opts.set('qp_solver_cond_N', qp_solver_cond_N);
end
ocp_opts.set('sim_method_num_stages', sim_method_num_stages);
ocp_opts.set('sim_method_num_steps', sim_method_num_steps);


%% acados ocp
% create ocp
ocp = acados_ocp(ocp_model, ocp_opts);


% set trajectory initialization
x_traj_init = zeros(nx, N+1);
u_traj_init = zeros(nu, N);
ocp.set('init_x', x_traj_init);
ocp.set('init_u', u_traj_init);


% solve
tic;
ocp.solve();
time_ext = toc;

% store and load iterate
filename = 'iterate.json';
ocp.store_iterate(filename, true);
ocp.load_iterate(filename);
delete(filename)

% get solution
u = ocp.get('u');
x = ocp.get('x');

% get info
status = ocp.get('status');
sqp_iter = ocp.get('sqp_iter');
time_tot = ocp.get('time_tot');
time_lin = ocp.get('time_lin');
time_reg = ocp.get('time_reg');
time_qp_sol = ocp.get('time_qp_sol');

fprintf('\nstatus = %d, sqp_iter = %d, time_ext = %f [ms], time_int = %f [ms] (time_lin = %f [ms], time_qp_sol = %f [ms], time_reg = %f [ms])\n', status, sqp_iter, time_ext*1e3, time_tot*1e3, time_lin*1e3, time_qp_sol*1e3, time_reg*1e3);

% print statistics
ocp.print('stat')

if status~=0
    error('ocp_nlp solver returned status nonzero');
elseif sqp_iter > 2
    error('ocp can be solved in 2 iterations!');
else
	fprintf(['\ntest_ocp_linear_mass_spring: success with sim method ', ...
        sim_method, ' !\n']);
end

% plot result
% figure()
% subplot(2, 1, 1)
% plot(0:N, x);
% title('trajectory')
% ylabel('x')
% subplot(2, 1, 2)
% plot(1:N, u);
% ylabel('u')
% xlabel('sample')
% 
% if is_octave()
%     waitforbuttonpress;
% end
