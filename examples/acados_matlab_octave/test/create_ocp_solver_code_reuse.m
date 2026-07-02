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


function ocp_solver = create_ocp_solver_code_reuse(creation_mode)

    json_file = 'pendulum_ocp.json';
    solver_creation_opts = struct();
    solver_creation_opts.json_file = json_file;
    if strcmp(creation_mode, 'standard')
        disp('Standard creation mode');
    elseif strcmp(creation_mode, 'ocp_from_json')
        disp('OCP from JSON creation mode');
    elseif strcmp(creation_mode, 'precompiled')
        solver_creation_opts.generate = false;
        solver_creation_opts.build = false;
        solver_creation_opts.compile_mex_wrapper = false;
    elseif strcmp(creation_mode, 'force_precompiled')
        solver_creation_opts.generate = false;
        solver_creation_opts.build = false;
        solver_creation_opts.compile_mex_wrapper = false;
        solver_creation_opts.check_reuse_possible = false;
    else
        error('Invalid creation mode')
    end

    if strcmp(creation_mode, 'ocp_from_json')
        ocp = AcadosOcp.from_json(json_file);
    else
        ocp = create_pendulum_ocp();
    end

    % create solver
    ocp_solver = AcadosOcpSolver(ocp, solver_creation_opts);

    if strcmp(creation_mode, 'precompiled') || strcmp(creation_mode, 'force_precompiled')
        % check if code reuse worked
        if ocp_solver.solver_creation_opts.generate || ocp_solver.solver_creation_opts.build
            error("Code reuse failed, solver was regenerated or rebuilt.");
        end
    end

end
