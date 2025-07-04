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

clear all;
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
ocp_solver.set('p', zeros(np, 1)); % TODO: this behaviour is only supported for p!

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
