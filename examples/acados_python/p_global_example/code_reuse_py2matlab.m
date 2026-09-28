%
% Copyright (c) The acados authors.
%
% This file is part of acados.
%
% Licensed under the 2-Clause BSD License.

import casadi.*

check_acados_requirements()

json_files = {'c_generated_code_single_phase/blz_True_pglobal_True.json', 'c_generated_code_multi_phase/mocp_blz_True_pglobal_True_0.json'};

for i = 1:length(json_files)
    json_file = json_files{i};
    disp('testing solver creation with code reuse with json file: ')
    disp(json_file)
    solver_creation_opts = struct();
    solver_creation_opts.json_file = json_file;
    solver_creation_opts.generate = false;
    solver_creation_opts.build = false;
    solver_creation_opts.compile_mex_wrapper = true;

    if ~isempty(strfind(json_file, 'multi'))
        ocp = AcadosMultiphaseOcp.from_json(json_file);
    else
        ocp = AcadosOcp.from_json(json_file);
    end
    % create solver
    ocp_solver = AcadosOcpSolver(ocp, solver_creation_opts);

    % test code reuse was done:
    if ocp_solver.solver_creation_opts.generate || ocp_solver.solver_creation_opts.build
        error('could not load solver from Matlab without rebuilding.');
    end


    nx = length(ocp_solver.get('x', 0));
    [nu, N] = size(ocp_solver.get('u'));

    for i = 1:5
        ocp_solver.solve();

        status = ocp_solver.get('status');
        ocp_solver.print('stat');
    end
    clear ocp_solver
end

