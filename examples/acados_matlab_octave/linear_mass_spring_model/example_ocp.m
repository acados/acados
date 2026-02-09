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

clear all;clc;
check_acados_requirements()

%% acados ocp model
model = get_linear_mass_spring_model();
nx = length(model.x);
nu = length(model.u);
ny = nx + nu;
ny_e = nx;

% constraint formulation
if 1
    % bounds on x and u
	nbx = nx/2;
	nbu = nu;
	ng = 0;
	nh = 0;
elseif 0
    % general linear constraints on x and u
	nbx = 0;
	nbu = 0;
	ng = nu+nx/2;
	ng_e = nx;
	nh = 0;
else
    % external function constraint
	nbx = 0;
	nbu = 0;
	ng = 0;
	nh = nu+nx;
	nh_e = nx;
end
%% set up OCP
ocp = AcadosOcp();
ocp.model = model;

T = 10.0; % horizon length time
N = 20;

ocp.solver_options.tf = T;
ocp.solver_options.N_horizon = N;
ocp.solver_options.nlp_solver_type = 'SQP'; % 'SQP_RTI'
ocp.solver_options.hessian_approx = 'GAUSS_NEWTON'; % 'EXACT', 'GAUSS_NEWTON'
ocp.solver_options.regularize_method = 'NO_REGULARIZE';% NO_REGULARIZE, PROJECT, PROOJECT_REDUC_HESS, MIRROR, CONVEXIFY
ocp.solver_options.nlp_solver_ext_qp_res = 1;
ocp.solver_options.qp_solver = 'PARTIAL_CONDENSING_HPIPM';
ocp.solver_options.qp_solver_cond_N = 5; % for partial condensing
ocp.solver_options.integrator_type = 'ERK'; % 'DISCRETE','ERK','IRK'
ocp.solver_options.sim_method_num_stages = 4;
ocp.solver_options.sim_method_num_steps = 3;
ocp.solver_options.compile_interface = [];

% cost
cost_type = 'AUTO'; % 'LINEAR_LS','NONLINEAR_LS'.'AUTO'
ocp.cost.cost_type_0 = cost_type;
ocp.cost.cost_type = cost_type;
ocp.cost.cost_type_e = cost_type;

Vu = zeros(ny, nu); for ii=1:nu Vu(ii,ii)=1.0; end % input-to-output matrix in lagrange term
Vx = zeros(ny, nx); for ii=1:nx Vx(nu+ii,ii)=1.0; end % state-to-output matrix in lagrange term
Vx_e = zeros(ny_e, nx); for ii=1:nx Vx_e(ii,ii)=1.0; end % state-to-output matrix in mayer term
W = eye(ny); for ii=1:nu W(ii,ii)=2.0; end % weight matrix in lagrange term
W_e = eye(ny_e); % weight matrix in mayer term
yr = zeros(ny, 1); % output reference in lagrange term
yr_e = zeros(ny_e, 1); % output reference in mayer term
    

sym_x = model.x;
sym_u = model.u;
yr_u = zeros(nu, 1);
yr_x = zeros(nx, 1);
dWu = 2*ones(nu, 1);
dWx = ones(nx, 1);
ymyr_0 = sym_u - yr_u;
ymyr = [sym_u; sym_x] - [yr_u; yr_x];
ymyr_e = sym_x - yr_x;
    
if (strcmp(cost_type, 'LINEAR_LS'))
    ocp.cost.Vu_0 = Vu;
    ocp.cost.Vu = Vu;
    ocp.cost.Vx_0 = Vx;
    ocp.cost.Vx = Vx;
    ocp.cost.Vx_e = Vx_e;
    ocp.cost.W_0 = W;
    ocp.cost.W = W;
    ocp.cost.W_e = W_e;
    ocp.cost.yref_0 = yr;
    ocp.cost.yref = yr; 
    ocp.cost.yref_e = yr_e;
elseif (strcmp(cost_type, 'NONLINEAR_LS'))
    cost_expr_y_0 = sym_u;
    cost_expr_y = [sym_u; sym_x];
    cost_expr_y_e = sym_x;

    ocp.model.cost_y_expr_0 = cost_expr_y;
    ocp.model.cost_y_expr = cost_expr_y;
    ocp.model.cost_y_expr_e = cost_expr_y_e;
    ocp.cost.W_0 = W;
    ocp.cost.W = W;
    ocp.cost.W_e = W_e;
    ocp.cost.yref_0 = yr;
    ocp.cost.yref = yr; 
    ocp.cost.yref_e = yr_e;
else
    cost_expr_ext_cost_0 = 0.5 * ymyr_0' * (dWu .* ymyr_0);
    cost_expr_ext_cost = 0.5 * ymyr' * ([dWu; dWx] .* ymyr);
    cost_expr_ext_cost_e = 0.5 * ymyr_e' * (dWx .* ymyr_e);

    ocp.model.cost_expr_ext_cost_0 = cost_expr_ext_cost;
    ocp.model.cost_expr_ext_cost = cost_expr_ext_cost;
    ocp.model.cost_expr_ext_cost_e = cost_expr_ext_cost_e;
end

% constraints
x0 = zeros(nx, 1); x0(1)=2.5; x0(2)=2.5;
ocp.constraints.x0 = x0;
if (ng>0)
	D = zeros(ng, nu); for ii=1:nu D(ii,ii)=1.0; end
	C = zeros(ng, nx); for ii=1:ng-nu C(nu+ii,ii)=1.0; end
	lg = zeros(ng, 1); for ii=1:nu lg(ii)=-0.5; end; for ii=1:ng-nu lg(nu+ii)=-4.0; end
	ug = zeros(ng, 1); for ii=1:nu ug(ii)= 0.5; end; for ii=1:ng-nu ug(nu+ii)= 4.0; end
	C_e = zeros(ng_e, nx); for ii=1:ng_e C_e(ii,ii)=1.0; end
	lg_e = zeros(ng_e, 1); for ii=1:ng_e lg_e(ii)=-4.0; end
	ug_e = zeros(ng_e, 1); for ii=1:ng_e ug_e(ii)= 4.0; end
elseif (nh>0)
    constr_expr_h_0 = sym_u;
    constr_expr_h = [sym_u; sym_x];
    constr_expr_h_e = sym_x;
	lh = zeros(nh, 1); for ii=1:nu lh(ii)=-0.5; end; for ii=1:nx lh(nu+ii)=-4.0; end
	uh = zeros(nh, 1); for ii=1:nu uh(ii)= 0.5; end; for ii=1:nx uh(nu+ii)= 4.0; end
	lh_e = zeros(nh_e, 1); for ii=1:nh_e lh_e(ii)=-4.0; end
	uh_e = zeros(nh_e, 1); for ii=1:nh_e uh_e(ii)= 4.0; end
else
    idxbx = (0:nbx-1)';
	lbx = -4*ones(nbx, 1);
	ubx =  4*ones(nbx, 1);
    idxbu = (0:nbu-1)';
	lbu = -0.5*ones(nbu, 1);
	ubu =  0.5*ones(nbu, 1);
end

if (ng>0)
    ocp.constraints.C = C;
    ocp.constraints.D = D;
    ocp.constraints.lg = lg;
    ocp.constraints.ug = ug;
    ocp.constraints.C = C_e;
    ocp.constraints.D = D_e;
    ocp.constraints.lg = lg_e;
    ocp.constraints.ug = ug_e;
elseif (nh>0)
    ocp.model.con_h_expr_0 = constr_expr_h_0;
    ocp.constraints.lh_0 = lh;
    ocp.constraints.uh_0 = uh;
    ocp.model.con_h_expr = constr_expr_h;
    ocp.constraints.lh = lh;
    ocp.constraints.uh = uh;
    ocp.model.con_h_expr_e = constr_expr_h_e;
    ocp.constraints.lh_e = lh_e;
    ocp.constraints.uh_e = uh_e;
else
    ocp.constraints.idxbx = idxbx;
    ocp.constraints.lbx = lbx;
    ocp.constraints.ubx = ubx;
    ocp.constraints.idxbu = idxbu;
    ocp.constraints.lbu = lbu;
    ocp.constraints.ubu = ubu;
end
%% acados ocp solver
% create ocp
ocp_solver = AcadosOcpSolver(ocp);

% set trajectory initialization
x_traj_init = zeros(nx, N+1);
u_traj_init = zeros(nu, N);
ocp_solver.set('init_x', x_traj_init);
ocp_solver.set('init_u', u_traj_init);

% solve
tic;
ocp_solver.solve();
time_ext = toc;

%x0(1) = 1.5;
%ocp_solver.set('constr_x0', x0);

% if not set, the trajectory is initialized with the previous solution

% get solution
u = ocp_solver.get('u');
x = ocp_solver.get('x');

% get info
status = ocp_solver.get('status');
sqp_iter = ocp_solver.get('sqp_iter');
time_tot = ocp_solver.get('time_tot');
time_lin = ocp_solver.get('time_lin');
time_reg = ocp_solver.get('time_reg');
time_qp_sol = ocp_solver.get('time_qp_sol');

fprintf('\nstatus = %d, sqp_iter = %d, time_ext = %f [ms], time_int = %f [ms] (time_lin = %f [ms], time_qp_sol = %f [ms], time_reg = %f [ms])\n', status, sqp_iter, time_ext*1e3, time_tot*1e3, time_lin*1e3, time_qp_sol*1e3, time_reg*1e3);

% print statistics
ocp_solver.print('stat')

if status==0
	fprintf('\nsuccess!\n\n');
else
	fprintf('\nsolution failed!\n\n');
end

% plot result
figure()
subplot(2, 1, 1)
plot(0:N, x);
title('trajectory')
ylabel('x')
subplot(2, 1, 2)
plot(1:N, u);
ylabel('u')
xlabel('sample')

if is_octave()
    waitforbuttonpress;
end
