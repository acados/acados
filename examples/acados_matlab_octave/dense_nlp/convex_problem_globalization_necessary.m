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

% Simplest Convex NLP where full-step SQP fails
%
% min log(exp(x) + exp(-x))
%
% Solution is x* = 0.0


import casadi.*
%
check_acados_requirements()

globalization_params = {'FIXED_STEP', 'FUNNEL_L1PEN_LINESEARCH', 'MERIT_BACKTRACKING'}

for num = 1:length(globalization_params)
    globalization = globalization_params{num};

    % set model
    model = AcadosModel();
    x = SX.sym('x');

    % dynamics: identity
    model.x = x;
    model.name = strcat('convex_problem_', globalization);

    %% solver settings
    N_horizon = 0;

    %% OCP formulation object
    ocp = AcadosOcp();
    ocp.model = model;

    % terminal cost term
    ocp.cost.cost_type_e = 'EXTERNAL';
    ocp.model.cost_expr_ext_cost_e = log(exp(model.x) + exp(-model.x));

    % define solver options
    ocp.solver_options.N_horizon = N_horizon;
    ocp.solver_options.qp_solver = 'FULL_CONDENSING_HPIPM';
    ocp.solver_options.hessian_approx = 'EXACT';
    ocp.solver_options.integrator_type = 'DISCRETE';
    ocp.solver_options.print_level = 1;
    ocp.solver_options.nlp_solver_type = 'SQP';
    ocp.solver_options.globalization = globalization;
    ocp.solver_options.qp_solver_iter_max = 400;
    ocp.solver_options.regularize_method = 'MIRROR';
    ocp.solver_options.nlp_solver_max_iter = 100;

    % create solver
    ocp_solver = AcadosOcpSolver(ocp);

    % solver initial guess
    xinit = 1.5;

    % set trajectory initialization
    ocp_solver.set('init_x', xinit, 0);

    % solve
    ocp_solver.solve();
    % get status
    status = ocp_solver.get('status');
    % get solution
    solution = ocp_solver.get('x', 0);

    % compare to analytical solution
    exact_solution = 0;
    sol_err = abs(solution - exact_solution);

    if strcmp(globalization, 'FIXED_STEP')
        assert(status ~= 0, 'FIXED_STEP converged. Theoretically impossible!');
    elseif strcmp(globalization, 'FUNNEL_L1PEN_LINESEARCH')
        assert(status == 0, 'FUNNEL_L1PEN_LINESEARCH did not converge. Algorithm should converge!');
        assert(sol_err <= 1e-5, "numerical solutions do not match to analytical solution with tolerance");
    else
        assert(status == 0, 'MERIT_BACKTRACKING did not converge. Algorithm should converge!');
        assert(sol_err <= 1e-5, "numerical solutions do not match to analytical solution with tolerance");
    end
end

