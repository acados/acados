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
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
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


function ocp = formulate_single_integrator_ocp(settings, terminal_phase_ocp)
    if nargin < 2
        terminal_phase_ocp = 0;
    end
    ocp = AcadosOcp();

    ocp.model = get_single_integrator_model();

    ocp.cost.cost_type = 'NONLINEAR_LS';
    ocp.cost.cost_type_e = 'NONLINEAR_LS';
    ocp.cost.W = diag([settings.L2_COST_P, settings.L2_COST_V]);
    ocp.cost.W_e = diag([1e1]);
    ocp.cost.yref = [0.0; 0.0];
    ocp.cost.yref_e = [0.0];

    ocp.model.cost_y_expr = vertcat(ocp.model.x, ocp.model.u);
    ocp.model.cost_y_expr_e = ocp.model.x;

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

