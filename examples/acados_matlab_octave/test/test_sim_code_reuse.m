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
creation_modes = {'standard', 'precompiled', 'no_sim'};
for i = 1:length(creation_modes)

    % simulation parameters
    N_sim = 100;
    x0 = [0; 1e-1; 0; 0]; % initial state
    u0 = 0; % control input

    sim_solver = create_sim_solver_code_reuse(creation_modes{i});
    nx = length(x0);

    %% simulate system in loop
    x_sim = zeros(nx, N_sim+1);
    x_sim(:,1) = x0;

    for ii=1:N_sim

        % set initial state
        sim_solver.set('x', x_sim(:,ii));
        sim_solver.set('u', u0);

        % solve
        sim_solver.solve();

        % get simulated state
        x_sim(:,ii+1) = sim_solver.get('xn');
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
