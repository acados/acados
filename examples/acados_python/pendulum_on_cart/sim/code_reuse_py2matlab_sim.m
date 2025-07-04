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



check_acados_requirements()

json_file = 'acados_sim.json';
solver_creation_opts = struct();
solver_creation_opts.json_file = json_file;
solver_creation_opts.generate = false;
solver_creation_opts.build = false;
solver_creation_opts.compile_mex_wrapper = false;

sim = [];

%% create integrator
sim_solver = AcadosSimSolver(sim, solver_creation_opts);

% simulation parameters
N_sim = 100;
x0 = [0; 1e-1; 0; 0]; % initial state
u0 = 0; % control input
nx = length(sim_solver.get('x', 0));
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
disp('simulated state')
disp(x_sim(:,end))
