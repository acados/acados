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


function ocp = formulate_double_integrator_ocp(settings)
    ocp = AcadosOcp();

    ocp.model = get_double_integrator_model();

    ocp.dims.N = settings.N_HORIZON;

    ocp.cost.cost_type = 'NONLINEAR_LS';
    ocp.cost.cost_type_e = 'NONLINEAR_LS';
    ocp.cost.W = diag([settings.L2_COST_P, settings.L2_COST_V, settings.L2_COST_A]);
    ocp.cost.W_e = diag([1e1, 1e1]);
    ocp.cost.yref = [0.0; 0.0; 0.0];
    ocp.cost.yref_e = [0.0; 0.0];

    ocp.model.cost_y_expr = vertcat(ocp.model.x, ocp.model.u);
    ocp.model.cost_y_expr_e = ocp.model.x;

    u_max = 50.0;
    ocp.constraints.lbu = [-u_max];
    ocp.constraints.ubu = [u_max];
    ocp.constraints.idxbu = [0];

    ocp.constraints.x0 = settings.X0;

end
