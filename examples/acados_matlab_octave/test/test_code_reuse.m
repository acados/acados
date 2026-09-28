%
% Copyright (c) The acados authors.
%
% This file is part of acados.
%
% Licensed under the 2-Clause BSD License.

import casadi.*

check_acados_requirements()
creation_modes = {'standard', 'precompiled', 'force_precompiled', 'ocp_from_json'};

for i = 1:length(creation_modes)
    disp(['testing creation mode ', creation_modes{i}]);
    ocp_solver = create_ocp_solver_code_reuse(creation_modes{i});
    nx = length(ocp_solver.get('x', 0));
    [nu, N] = size(ocp_solver.get('u'));
    T = 1;

    % solver initial guess
    x_traj_init = zeros(nx, N+1);
    u_traj_init = zeros(nu, N);

    %% call ocp solver
    % solve
    ocp_solver.solve();
    % get solution
    sol = ocp_solver.get_iterate();

    status = ocp_solver.get('status'); % 0 - success
    % ocp_solver.print('stat');
    stat = ocp_solver.get('stat');
    if i == 1
        stat_ref = stat;
    elseif max(abs(stat-stat_ref)) > 1e-6
        error('solvers should have the same log independent of compilation options');
    end

    clear ocp_solver
end
