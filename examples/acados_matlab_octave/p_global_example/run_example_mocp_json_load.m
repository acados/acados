function [state_trajectories, timing] = run_example_mocp_json_load(json_file, code_reuse)
    import casadi.*

    mocp = AcadosMultiphaseOcp.from_json(json_file);

    % MOCP solver
    if code_reuse
        solver_creation_opts.build = false;
        solver_creation_opts.generate = false;
        solver_creation_opts.compile_mex_wrapper = false;
        % TODO: double check if this is what front end should look like.
    else
        solver_creation_opts = struct();
    end
    mocp_solver = AcadosOcpSolver(mocp, solver_creation_opts);

    % check if code reuse worked
    if code_reuse && (mocp_solver.solver_creation_opts.generate || mocp_solver.solver_creation_opts.build)
        error("Code reuse failed, solver was regenerated or rebuilt.");
    else
        disp("Created solver successfully without regenerating.");
        disp("");
    end

    state_trajectories = []; % only for testing purposes

    if ~isempty(mocp.model{1}.p_global)
        disp("Calling precompute.")
        tic
        mocp_solver.set_p_global_and_precompute_dependencies(mocp.p_global_values);
        toc
    end

    timing = 0;
    for i = 1:20
        mocp_solver.solve();
        state_trajectories = [state_trajectories; mocp_solver.get('x')];
        timing = timing + mocp_solver.get('time_lin');
    end
end
