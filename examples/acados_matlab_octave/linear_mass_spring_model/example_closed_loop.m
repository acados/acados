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

%% example of closed loop simulation
clear all;clc;
check_acados_requirements()

%% acados ocp model
model = get_linear_mass_spring_model();
% dims
nx = length(model.x);
nu = length(model.u);
ny = nx + nu;
ny_e = nx;
nbx = nx/2;
nbu = nu;
%% set up OCP and OCP solver
ocp = AcadosOcp();
ocp.model = model;

T = 10.0; % horizon length time
ocp_N = 20;

ocp.solver_options.tf = T;
ocp.solver_options.N_horizon = ocp_N;
ocp.solver_options.nlp_solver_type = 'SQP'; % 'SQP_RTI'
ocp.solver_options.regularize_method = 'NO_REGULARIZE';% NO_REGULARIZE, PROJECT, PROOJECT_REDUC_HESS, MIRROR, CONVEXIFY
ocp.solver_options.qp_solver = 'PARTIAL_CONDENSING_HPIPM';
ocp.solver_options.qp_solver_cond_N = 5; % for partial condensing
ocp.solver_options.integrator_type = 'IRK'; % 'DISCRETE','ERK','IRK'
ocp.solver_options.sim_method_num_stages = 2;
ocp.solver_options.sim_method_num_steps = 2;
ocp.solver_options.compile_interface = [];

% cost
cost_type = 'LINEAR_LS'; % 'LINEAR_LS','NONLINEAR_LS'.'EXT_COST'
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
idxbx = (0:nbx-1)';
lbx = -4*ones(nbx, 1);
ubx =  4*ones(nbx, 1);
idxbu = (0:nbu-1)';
lbu = -0.5*ones(nbu, 1);
ubu =  0.5*ones(nbu, 1);
ocp.constraints.idxbx = idxbx;
ocp.constraints.lbx = lbx;
ocp.constraints.ubx = ubx;
ocp.constraints.idxbu = idxbu;
ocp.constraints.lbu = lbu;
ocp.constraints.ubu = ubu;

ocp_solver = AcadosOcpSolver(ocp);

%% set up sim and sim solver
sim = AcadosSim();
sim.model = model;
% simulation
sim.solver_options.integrator_type = 'IRK';
sim.solver_options.Tsim = T/ocp_N;
sim.solver_options.num_stages = 4;
sim.solver_options.num_steps = 4;
sim.solver_options.sens_forw = false;

sim_solver = AcadosSimSolver(sim);
%% closed loop simulation
n_sim = 100;
x_sim = zeros(nx, n_sim+1);
x_sim(:,1) = zeros(nx,1); x_sim(1:2,1) = [3.5; 3.5];
u_sim = zeros(nu, n_sim);

x_traj_init = zeros(nx, ocp_N+1);
u_traj_init = zeros(nu, ocp_N);

tic;

for ii=1:n_sim

	% set x0
	ocp_solver.set('constr_x0', x_sim(:,ii));

	% set trajectory initialization
	ocp_solver.set('init_x', x_traj_init);
	ocp_solver.set('init_u', u_traj_init);

	% solve OCP
	ocp_solver.solve();

	% get solution
	%x_traj = ocp_solver.get('x');
	%u_traj = ocp_solver.get('u');
	u_sim(:,ii) = ocp_solver.get('u', 0);

	% set initial state of sim
	sim_solver.set('x', x_sim(:,ii));
	% set input in sim
	sim_solver.set('u', u_sim(:,ii));

	% simulate state
	sim_solver.solve();

	% get new state
	x_sim(:,ii+1) = sim_solver.get('xn');
end

avg_time_solve = toc/n_sim


% plot result
figure()
subplot(2, 1, 1)
plot(0:n_sim, x_sim);
title('closed loop simulation')
ylabel('x')
subplot(2, 1, 2)
plot(1:n_sim, u_sim);
ylabel('u')
xlabel('sample')


status = ocp_solver.get('status');

if status==0
	fprintf('\nsuccess!\n\n');
else
	fprintf('\nsolution failed!\n\n');
end


if is_octave()
    waitforbuttonpress;
end
