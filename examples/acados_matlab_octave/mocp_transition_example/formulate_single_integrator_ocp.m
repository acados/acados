%
% Copyright (c) The acados authors.
%
% This file is part of acados.
%
% Licensed under the 2-Clause BSD License.


function ocp = formulate_single_integrator_ocp(settings, terminal_phase_ocp, cost_type)
    if nargin < 3
        cost_type = 'NONLINEAR_LS'
    end
    if nargin < 2
        terminal_phase_ocp = 0;
    end
    ocp = AcadosOcp();

    ocp.model = get_single_integrator_model();

    ocp.cost.cost_type_e = 'NONLINEAR_LS';
    ocp.cost.yref_e = [0.0];
    ocp.model.cost_y_expr_e = ocp.model.x;
    ocp.cost.W_e = diag([1e1]);
    if strcmp(cost_type, 'NONLINEAR_LS')
        ocp.cost.cost_type = 'NONLINEAR_LS';
        ocp.cost.W = diag([settings.L2_COST_P, settings.L2_COST_V]);
        ocp.cost.yref = [0.0; 0.0];

        ocp.model.cost_y_expr = vertcat(ocp.model.x, ocp.model.u);
    elseif strcmp(cost_type, 'EXTERNAL')
        ocp.cost.cost_type = 'EXTERNAL';

        xu = vertcat(ocp.model.x, ocp.model.u);
        ocp.model.cost_expr_ext_cost = .5 * xu.transpose() * diag([settings.L2_COST_P, settings.L2_COST_V]) * xu;
    else
        error(['cost_type ', cost_type, ' not supported.'])
    end

    u_max = 5.0;
    ocp.constraints.lbu = [-u_max];
    ocp.constraints.ubu = [u_max];
    ocp.constraints.idxbu = [0];

    if settings.WITH_X_BOUNDS
        ocp.constraints.lbx = [-100];
        ocp.constraints.ubx = [100];
        ocp.constraints.idxbx = [0];

        ocp.constraints.lbx_e = [-100];
        ocp.constraints.ubx_e = [100];
        ocp.constraints.idxbx_e = [0];
    end

    % initial state constraint, not needed for terminal phase OCP.
    if ~terminal_phase_ocp
        ocp.constraints.x0 = settings.X0(1);
    end
end



function model = get_single_integrator_model()
    import casadi.*
    model = AcadosModel();
    model.name = 'single_integrator';
    model.x = SX.sym('p');
    model.u = SX.sym('v');
    model.f_expl_expr = model.u;
end

