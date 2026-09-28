%
% Copyright (c) The acados authors.
%
% This file is part of acados.
%
% Licensed under the 2-Clause BSD License.

import casadi.*

check_acados_requirements()
creation_modes = {'standard', 'precompiled', 'sim_from_json'};

json_file = [];

for i = 1:length(creation_modes)

    % simulation parameters
    N_sim = 100;
    x0 = [0; 1e-1; 0; 0]; % initial state
    u0 = 0; % control input
    nx = length(x0);

    sim_solver = create_sim_solver_code_reuse(creation_modes{i}, json_file);

    % store json file for later run with creation mode sim_from_json
    json_file = sim_solver.sim.code_gen_options.json_file;

    %% simulate system in loop
    x_sim = zeros(nx, N_sim+1);
    x_sim(:,1) = x0;

    for ii=1:N_sim
        x_sim(:,ii+1) = sim_solver.simulate(x_sim(:, ii), u0);
    end

    % forward sensitivities ( dxn_d[x0,u] )
    S_forw = sim_solver.get('S_forw');

    if i == 1
        S_forw_ref = S_forw;
    elseif max(abs(S_forw-S_forw_ref)) > 1e-6
        error('solvers should have the same output independent of compilation options');
    end
    S_forw

    clear sim_solver
end
