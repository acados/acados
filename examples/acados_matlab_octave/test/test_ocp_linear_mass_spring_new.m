
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

%% 1. Load Model
addpath('.');
model = get_linear_mass_spring_model();

nx = length(model.x);
nu = length(model.u);
ny = nx + nu;
ny_e = nx;

sym_x = model.x;
sym_u = model.u;

%% 2. Set up OCP Object
ocp = AcadosOcp();
ocp.model = model;

% --- Horizon & Shooting Nodes ---
N = 20;
T = 10.0;

shooting_nodes = [ linspace(0,1,N/2) linspace(1.1,5,N/2+1) ];

ocp.solver_options.N_horizon = N;
ocp.solver_options.shooting_nodes = shooting_nodes;
ocp.solver_options.tf = T;

% --- Solver Options ---
ocp.solver_options.nlp_solver_type = 'SQP';
ocp.solver_options.hessian_approx = 'GAUSS_NEWTON';
ocp.solver_options.regularize_method = 'CONVEXIFY';
ocp.solver_options.nlp_solver_ext_qp_res = 1;
ocp.solver_options.nlp_solver_max_iter = 100;

ocp.solver_options.qp_solver = 'PARTIAL_CONDENSING_HPIPM';
ocp.solver_options.qp_solver_cond_N = 5;

% --- Dynamics ---
ocp.solver_options.integrator_type = 'DISCRETE';

%% 3. Cost
%cost_type = 'ext_cost'
ocp.cost.cost_type = 'EXTERNAL';
ocp.cost.cost_type_e = 'EXTERNAL';
ocp.cost.cost_type_0 = 'EXTERNAL';

% y = [u; x]
W_diag = ones(ny, 1);
W_diag(1:nu) = 2.0;
cost_expr = 0.5 * ([sym_u; sym_x] .* W_diag)' * [sym_u; sym_x];
cost_expr_e = 0.5 * sym_x' * eye(nx) * sym_x;

ocp.model.cost_expr_ext_cost = cost_expr;
ocp.model.cost_expr_ext_cost_e = cost_expr_e;
ocp.model.cost_expr_ext_cost_0 = cost_expr;

%% 4. Constraints
nh = nu + nx;
nh_e = nx;

constr_h = [sym_u; sym_x];
constr_h_e = sym_x;

ocp.model.con_h_expr = constr_h;
ocp.model.con_h_expr_e = constr_h_e;

lh = zeros(nh, 1);
for ii=1:nu lh(ii)=-0.5; end;
for ii=1:nx lh(nu+ii)=-4.0; end

uh = zeros(nh, 1);
for ii=1:nu uh(ii)= 0.5; end;
for ii=1:nx uh(nu+ii)= 4.0; end

lh_e = zeros(nh_e, 1); for ii=1:nh_e lh_e(ii)=-4.0; end
uh_e = zeros(nh_e, 1); for ii=1:nh_e uh_e(ii)= 4.0; end

ocp.constraints.lh = lh;
ocp.constraints.uh = uh;
ocp.constraints.lh_e = lh_e;
ocp.constraints.uh_e = uh_e;

x0 = zeros(nx, 1); x0(1)=2.5; x0(2)=2.5;
ocp.constraints.x0 = x0;

%% 5. Acados OCP Solver
ocp_solver = AcadosOcpSolver(ocp);

% Initialization
x_traj_init = zeros(nx, N+1);
u_traj_init = zeros(nu, N);
ocp_solver.set('init_x', x_traj_init);
ocp_solver.set('init_u', u_traj_init);

%% 6. Solve
tic;
ocp_solver.solve();
time_ext = toc;

%% 7. Diagnostics & Strict Checks

% 1. Iterate storing
filename = 'iterate.json';
ocp_solver.store_iterate(filename, true);
ocp_solver.load_iterate(filename);
delete(filename)

% 2. QP Dump
filename = 'qp.json';
try
    ocp_solver.dump_last_qp_to_json(filename);
    delete(filename);
catch
    warning('dump_last_qp_to_json failed or not supported.');
end

% 3. QP Diagnostics
qp_diagnostics_result = ocp_solver.qp_diagnostics();

% 4. Get Info
status = ocp_solver.get('status');
sqp_iter = ocp_solver.get('sqp_iter');
time_tot = ocp_solver.get('time_tot');

fprintf('\nstatus = %d, sqp_iter = %d, time_ext = %f [ms], time_int = %f [ms]\n', ...
    status, sqp_iter, time_ext*1e3, time_tot*1e3);

ocp_solver.print('stat')

% 5. Strict Logic Check
if status ~= 0
    error('ocp_nlp solver returned status nonzero');
elseif sqp_iter > 2
    error('ocp can be solved in 2 iterations! Something is wrong with the setup.');
else
	fprintf(['\ntest_ocp_linear_mass_spring: success with integrator ', ...
        ocp.solver_options.integrator_type, ' !\n']);
end

% Plot (Optional)
% figure()
% subplot(2, 1, 1)
% plot(0:N, x);
% title('trajectory')
% ylabel('x')
% subplot(2, 1, 2)
% plot(1:N, u);
% ylabel('u')
% xlabel('sample')