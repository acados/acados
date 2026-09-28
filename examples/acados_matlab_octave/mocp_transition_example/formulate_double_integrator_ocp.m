%
% Copyright (c) The acados authors.
%
% This file is part of acados.
%
% Licensed under the 2-Clause BSD License.


function ocp = formulate_double_integrator_ocp(settings, first_phase_ocp)
    if nargin < 2
        first_phase_ocp = 0;
    end
    ocp = AcadosOcp();

    ocp.model = get_double_integrator_model();

    ocp.cost.cost_type = 'NONLINEAR_LS';
    ocp.cost.W = diag([settings.L2_COST_P, settings.L2_COST_V, settings.L2_COST_A]);
    ocp.cost.yref = [0.0; 0.0; 0.0];

    ocp.model.cost_y_expr = vertcat(ocp.model.x, ocp.model.u);

    % terminal cost - not needed when formulating first phase of MOCP
    if ~first_phase_ocp
        ocp.cost.cost_type_e = 'NONLINEAR_LS';
        ocp.model.cost_y_expr_e = ocp.model.x;
        ocp.cost.yref_e = [0.0; 0.0];
        ocp.cost.W_e = diag([1e1, 1e1]);
    end

    u_max = 50.0;
    ocp.constraints.lbu = [-u_max];
    ocp.constraints.ubu = [u_max];
    ocp.constraints.idxbu = [0];

    if settings.WITH_X_BOUNDS
        ocp.constraints.lbx = [-100, -10];
        ocp.constraints.ubx = [100, 10];
        ocp.constraints.idxbx = [0, 1];
    end

    ocp.constraints.x0 = settings.X0;

end
