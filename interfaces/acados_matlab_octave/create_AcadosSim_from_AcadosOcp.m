%
% Copyright (c) The acados authors.
%
% This file is part of acados.
%
% Licensed under the 2-Clause BSD License.

% creates an AcadosSim object from an AcadosOcp object with the same model and integrator options corresponding to the first stage of the OCP
function sim = create_AcadosSim_from_AcadosOcp(ocp)

    if ~isa(ocp, 'AcadosOcp')
        error('create_AcadosSim_from_AcadosOcp: First argument must be an AcadosOcp object.');
    end
    ocp.make_consistent();
    if strcmp(ocp.solver_options.integrator_type, 'DISCRETE')
        error('create_AcadosSim_from_AcadosOcp: AcadosOcp cannot have integrator_type DISCRETE.');
    end
    if ~ocp.model.p_global.is_empty()
        error('create_AcadosSim_from_AcadosOcp: AcadosOcp cannot have p_global.');
    end
    sim = AcadosSim();
    sim.model = ocp.model.copy();

    if ~sim.model.p_global.is_empty()
        sim.model.p = [sim.model.p; sim.model.p_global];
        sim.model.p_global = [];

        warning('Model contained p_global. Appending p_global to p in the sim model.');
    end

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