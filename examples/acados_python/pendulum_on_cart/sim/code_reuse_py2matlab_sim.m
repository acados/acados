%
% Copyright (c) The acados authors.
%
% This file is part of acados.
%
% Licensed under the 2-Clause BSD License.



check_acados_requirements()

json_file = 'codegen_pendulum/pendulum.json';
solver_creation_opts = struct();
solver_creation_opts.json_file = json_file;
solver_creation_opts.generate = false;
solver_creation_opts.build = false;
solver_creation_opts.compile_mex_wrapper = false;

% TODO load sim from json
% sim = [];

% %% create integrator
% sim_solver = AcadosSimSolver(sim, solver_creation_opts);

% % simulation parameters
% N_sim = 100;
% x0 = [0; 1e-1; 0; 0]; % initial state
% u0 = 0; % control input
% nx = length(sim_solver.get('x', 0));
% %% simulate system in loop
% x_sim = zeros(nx, N_sim+1);
% x_sim(:,1) = x0;
% for ii=1:N_sim
%     % set initial state
%     sim_solver.set('x', x_sim(:,ii));
%     sim_solver.set('u', u0);
%     % solve
%     sim_solver.solve();
%     % get simulated state
%     x_sim(:,ii+1) = sim_solver.get('xn');
% end
% disp('simulated state')
% disp(x_sim(:,end))
