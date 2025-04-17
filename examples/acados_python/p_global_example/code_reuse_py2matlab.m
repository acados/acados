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

json_files = {'acados_ocp_pendulum_blazing_True_p_global_True.json', 'mocp.json'};

for i = 1:length(json_files)
    json_file = json_files{i};
    disp('testing solver creation with code reuse with json file: ')
    disp(json_file)
    solver_creation_opts = struct();
    solver_creation_opts.json_file = json_file;
    solver_creation_opts.generate = false;
    solver_creation_opts.build = false;
    solver_creation_opts.compile_mex_wrapper = true;
    ocp = [];

    % create solver
    ocp_solver = AcadosOcpSolver(ocp, solver_creation_opts);

    nx = length(ocp_solver.get('x', 0));
    [nu, N] = size(ocp_solver.get('u'));

    for i = 1:5
        ocp_solver.solve();

        status = ocp_solver.get('status');
        ocp_solver.print('stat');
    end
    clear ocp_solver
end

