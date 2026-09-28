%
% Copyright (c) The acados authors.
%
% This file is part of acados.
%
% Licensed under the 2-Clause BSD License.

%

check_acados_requirements()

import casadi.*

N = 20; % number of discretization steps
nx = 3;
nu = 3;
np = 10;
[ocp, x0] = create_parametric_ocp_qp(N, np);

%% create ocp solver
ocp_solver = AcadosOcpSolver(ocp);
% NOTE: here we don't perform iterations and just test initialization
% functionality
ocp.solver_options.nlp_solver_max_iter = 0;

%% test parameter setters and getters;
ocp_solver.set('p', zeros(np, 1), 0, N+1);

p_vals_different = reshape(1:(np*(N+1)), [np, N+1]);
ocp_solver.set('p', p_vals_different);

for i_p = 1:N+1;
    p_val = p_vals_different(:, i_p);
    ocp_solver.set('p', p_val, 0);
    p = ocp_solver.get('p', 0);

    if any(p ~= p_val)
        disp('simple parameter setter doesnt work properly');
        exit(1);
    end
end

for stage = 1:N
    p = ocp_solver.get('p', stage);
    if any(p ~= p_vals_different(:, stage+1))
        disp('simple parameter setter doesnt work properly, parameter values should not change after setting at another stage');
        exit(1);
    end
end
disp('simple parameter setter works properly');


%% test sparse parameter update
for stage = [0, 8]
    ocp_solver.set('p', 0 * p_val, stage);

    idx_values = [0, 5, np-1];
    new_p_values = [9, 42, 12];

    ocp_solver.set_params_sparse(idx_values, new_p_values, stage);
    p = ocp_solver.get('p', stage);

    p_val = zeros(np, 1);
    idx_values_matlab = idx_values + 1;
    p_val(idx_values_matlab) = new_p_values;

    if any(p~=p_val)
        disp('sparse parameter update doesnt work properly.');
        exit(1)
    end
end
disp('sparse parameter setter works properly');

clear ocp_solver