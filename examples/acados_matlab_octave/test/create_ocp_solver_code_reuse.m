%
% Copyright (c) The acados authors.
%
% This file is part of acados.
%
% Licensed under the 2-Clause BSD License.


function ocp_solver = create_ocp_solver_code_reuse(creation_mode)

    json_file = fullfile('codegen_pendulum_ocp','pendulum_ocp.json');
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
        ocp.name = 'pendulum_ocp';
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
