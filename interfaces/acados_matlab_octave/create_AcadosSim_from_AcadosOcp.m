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
%

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