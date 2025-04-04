function sim = create_AcadosSim_from_AcadosOcp(ocp)

    if ~isa(ocp, 'AcadosOcp')
        error('create_AcadosSim_from_AcadosOcp: First argument must be an AcadosOcp object.');
    end
    ocp.make_consistent();
    sim = AcadosSim();
    sim.model = ocp.model;
    % copy all relevant options
    sim.solver_options.integrator_type = ocp.solver_options.integrator_type;
    sim.solver_options.collocation_type = ocp.solver_options.collocation_type;
    sim.solver_options.Tsim = ocp.solver_options.Tsim;
    sim.solver_options.num_stages = ocp.solver_options.sim_method_num_stages(1);
    sim.solver_options.num_steps = ocp.solver_options.sim_method_num_steps(1);
    sim.solver_options.newton_iter = ocp.solver_options.sim_method_newton_iter(1);
    sim.solver_options.newton_tol = ocp.solver_options.sim_method_newton_tol(1);
    sim.solver_options.jac_reuse = ocp.solver_options.sim_method_jac_reuse(1);
    sim.solver_options.ext_fun_compile_flags = ocp.solver_options.ext_fun_compile_flags;
    sim.parameter_values = ocp.parameter_values;
end