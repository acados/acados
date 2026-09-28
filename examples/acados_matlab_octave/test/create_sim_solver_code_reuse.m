%
% Copyright (c) The acados authors.
%
% This file is part of acados.
%
% Licensed under the 2-Clause BSD License.


function sim_solver = create_sim_solver_code_reuse(creation_mode, json_file)

    addpath('../pendulum_on_cart_model')

    check_acados_requirements()

    solver_creation_opts = struct();
    if strcmp(creation_mode, 'standard')
        disp('Standard creation mode');
    elseif strcmp(creation_mode, 'precompiled') || strcmp(creation_mode, 'sim_from_json')
        solver_creation_opts.generate = false;
        solver_creation_opts.build = false;
        solver_creation_opts.compile_mex_wrapper = false;
    else
        error('Invalid creation mode')
    end

    if strcmp(creation_mode, 'sim_from_json')
        sim = AcadosSim.from_json(json_file);
        [filepath, name, ext] = fileparts(json_file);
        sim.code_gen_options.json_file = json_file;
        sim.code_gen_options.code_export_directory = filepath;
    else
        model = get_pendulum_on_cart_model();
        sim = AcadosSim();
        sim.model = model;
        sim.solver_options.Tsim = 0.1; % simulation time
        sim.solver_options.integrator_type = 'ERK';
    end

    %% create integrator
    sim_solver = AcadosSimSolver(sim, solver_creation_opts);
end
