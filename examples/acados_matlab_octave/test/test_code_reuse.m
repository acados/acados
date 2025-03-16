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
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS 'AS IS'
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

import casadi.*

check_acados_requirements()
creation_modes = {'standard', 'precompiled', 'no_ocp'};
for i = 1:length(creation_modes)
    ocp_solver = create_ocp_solver_code_reuse(creation_modes{i});
    nx = length(ocp_solver.get('x', 0));
    [nu, N] = size(ocp_solver.get('u'));
    T = 1;

    % solver initial guess
    x_traj_init = zeros(nx, N+1);
    u_traj_init = zeros(nu, N);

    %% call ocp solver
    % set trajectory initialization
    ocp_solver.set('init_x', x_traj_init); % states
    ocp_solver.set('init_u', u_traj_init); % inputs
    ocp_solver.set('init_pi', zeros(nx, N)); % multipliers for dynamics equality constraints

    % solve
    ocp_solver.solve();
    % get solution
    utraj = ocp_solver.get('u');
    xtraj = ocp_solver.get('x');

    status = ocp_solver.get('status'); % 0 - success
    ocp_solver.print('stat');
    stat = ocp_solver.get('stat');
    if i == 1
        stat_ref = stat;
    elseif max(abs(stat-stat_ref)) > 1e-6
        error('solvers should have the same log independent of compilation options');
    end

    clear ocp_solver
end
