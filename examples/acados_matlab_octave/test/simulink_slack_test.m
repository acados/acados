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

%%
N = 20; % number of discretization steps
nx = 3;
nu = 3;
[ocp, simulink_opts] = create_slacked_ocp_qp_solver_formulation(N);
ocp.simulink_opts = simulink_opts;

%% create ocp solver
ocp_solver = AcadosOcpSolver(ocp);

%% test solve
% setup some slack penalty values;
ns = nu;
zl = (1:ns*N);
zu = 2*zl;
Zl = 4*zl;
Zu = 8*zl;
slacks_init = [9e4 * zl(:); 3e2*zl(:)];
% use some initial state such that upper and lower bounds are active
x0 = [1, -3, 2];
ocp_solver.set('constr_x0', x0)
for stage=0:N-1
    ocp_solver.set('cost_zl', zl((1+stage*ns:(stage+1)*ns)), stage);
    ocp_solver.set('cost_zu', zu((1+stage*ns:(stage+1)*ns)), stage);
    ocp_solver.set('cost_Zl', Zl((1+stage*ns:(stage+1)*ns)), stage);
    ocp_solver.set('cost_Zu', Zu((1+stage*ns:(stage+1)*ns)), stage);
    ocp_solver.set('sl', slacks_init((1+2*(stage*ns):2*(stage*ns)+ns)), stage);
    ocp_solver.set('su', slacks_init((1+ns+2*(stage*ns):2*(stage*ns)+2*ns)), stage);
end

ocp_solver.solve();
ocp_solver.print();

xtraj = ocp_solver.get('x');
utraj = ocp_solver.get('u');

%% compile S-function
cd c_generated_code
make_sfun; % ocp solver
cd ..;

%% test
n_sim = 1;
out_sim = sim('block_simulink_slacks', 'SaveOutput', 'on');
element_names = out_sim.logsout.getElementNames();

status_signal = out_sim.logsout.getElement('status');
disp('checking status, should be 0.')
if any(status_signal.Values.Data ~= 0)
    error(['failed. got status values:' mat2str(status_signal.Values.Data)]);
end

utraj_signal = out_sim.logsout.getElement('utraj');
u_simulink = utraj_signal.Values.Data(1, :);
disp('checking u values.')
if any(abs(u_simulink(:) - utraj(:)) > 1e-8)
    error('failed');
end

xtraj_signal = out_sim.logsout.getElement('xtraj');
x_simulink = xtraj_signal.Values.Data(1, :);
disp('checking x values.')
if any(abs(x_simulink(:) - xtraj(:)) > 1e-8)
    error('failed');
end
