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
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS 'AS IS'
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


function closed_loop_as_rti()
    algorithms = {'SQP', 'RTI', 'AS-RTI-A', 'AS-RTI-B', 'AS-RTI-C', 'AS-RTI-D'};
    % algorithms = {'AS-RTI-B'};
    for i = 1:length(algorithms)
        main(algorithms{i}, 1);
    end
end

function [ocp_solver, integrator] = setup(x0, Fmax, N_horizon, Tf, algorithm, as_rti_iter)
    if nargin < 6
        as_rti_iter = 1;
    end

    disp(['running with algorithm: ', algorithm, ', as_rti_iter: ', num2str(as_rti_iter)]);

    % create ocp object to formulate the OCP
    ocp = AcadosOcp();

    % set model
    model = get_pendulum_on_cart_model();
    ocp.model = model;

    nx = size(model.x, 1);
    nu = size(model.u, 1);
    ny = nx + nu;
    ny_e = nx;

    ocp.solver_options.N_horizon = N_horizon;

    % set cost module
    Q_mat = 2 * diag([1e3, 1e3, 1e-2, 1e-2]);
    R_mat = 2 * diag([1e-2]);
    if 0
        % NOTE: not yet implemented in MATLAB
        ocp.cost.cost_type = 'NONLINEAR_LS';
        ocp.cost.cost_type_e = 'NONLINEAR_LS';

        ocp.cost.W = blkdiag(Q_mat, R_mat);
        ocp.cost.W_e = Q_mat;

        ocp.model.cost_y_expr = [model.x; model.u];
        ocp.model.cost_y_expr_e = model.x;
        ocp.cost.yref = zeros(ny, 1);
        ocp.cost.yref_e = zeros(ny_e, 1);
        ocp.translate_nls_cost_to_conl();
    else
        ocp.cost.cost_type = 'EXTERNAL';
        ocp.cost.cost_type_e = 'EXTERNAL';

        ocp.model.cost_expr_ext_cost = model.x' * Q_mat * model.x + model.u' * R_mat * model.u;
        ocp.model.cost_expr_ext_cost_e = model.x' * Q_mat * model.x;
    end

    % set constraints
    ocp.constraints.lbu = -Fmax;
    ocp.constraints.ubu = Fmax;

    ocp.constraints.x0 = x0;
    ocp.constraints.idxbu = [0];

    ocp.solver_options.qp_solver = 'PARTIAL_CONDENSING_HPIPM';
    ocp.solver_options.hessian_approx = 'GAUSS_NEWTON';
    ocp.solver_options.integrator_type = 'IRK';
    ocp.solver_options.sim_method_newton_iter = 10;


    if ismember(algorithm, {'RTI', 'AS-RTI-A', 'AS-RTI-B', 'AS-RTI-C', 'AS-RTI-D'})
        ocp.solver_options.nlp_solver_type = 'SQP_RTI';
    elseif strcmp(algorithm, 'SQP')
        ocp.solver_options.nlp_solver_type = 'SQP';
    else
        error(['unknown algorithm: ', algorithm]);
    end

    if strcmp(algorithm, 'AS-RTI-A')
        ocp.solver_options.as_rti_iter = as_rti_iter;
        ocp.solver_options.as_rti_level = 0;
    elseif strcmp(algorithm, 'AS-RTI-B')
        ocp.solver_options.as_rti_iter = as_rti_iter;
        ocp.solver_options.as_rti_level = 1;
    elseif strcmp(algorithm, 'AS-RTI-C')
        ocp.solver_options.as_rti_iter = as_rti_iter;
        ocp.solver_options.as_rti_level = 2;
    elseif strcmp(algorithm, 'AS-RTI-D')
        ocp.solver_options.as_rti_iter = as_rti_iter;
        ocp.solver_options.as_rti_level = 3;
    end

    ocp.solver_options.qp_solver_cond_N = N_horizon;

    % set prediction horizon
    ocp.solver_options.tf = Tf;

    ocp.model.name = strcat(ocp.model.name, '_', strrep(algorithm, '-', '_'));
    ocp_solver = AcadosOcpSolver(ocp);

    % create an integrator with the same settings as used in the OCP solver
    sim = create_AcadosSim_from_AcadosOcp(ocp);
    integrator = AcadosSimSolver(sim);
end

function main(algorithm, as_rti_iter)
    if nargin < 2
        as_rti_iter = 1;
    end

    x0 = [0.0; pi; 0.0; 0.0];
    Fmax = 80;

    Tf = 0.8;
    N_horizon = 40;

    [ocp_solver, integrator] = setup(x0, Fmax, N_horizon, Tf, algorithm, as_rti_iter);

    nx = ocp_solver.ocp.dims.nx;
    nu = ocp_solver.ocp.dims.nu;

    Nsim = 100;
    simX = zeros(Nsim+1, nx);
    simU = zeros(Nsim, nu);

    simX(1, :) = x0';

    if ~strcmp(algorithm, 'SQP')
        t_preparation = zeros(Nsim, 1);
        t_feedback = zeros(Nsim, 1);
    else
        t = zeros(Nsim, 1);
    end

    % closed loop
    for i = 1:Nsim
        if ~strcmp(algorithm, 'SQP')
            % preparation phase
            ocp_solver.set('rti_phase', 1);
            ocp_solver.solve();
            status = ocp_solver.get('status');
            t_preparation(i) = ocp_solver.get('time_tot');

            if ~ismember(status, [0, 2, 5])
                error(['acados returned status ', num2str(status), '. Exiting.']);
            end

            % set initial state
            % NOTE: all bounds can be updated between the phases
            % all other updates, such as parameters are not used in the feedback phase
            ocp_solver.set('constr_lbx', simX(i, :), 0);
            ocp_solver.set('constr_ubx', simX(i, :), 0);

            % feedback phase
            ocp_solver.set('rti_phase', 2);
            ocp_solver.solve();
            t_feedback(i) = ocp_solver.get('time_tot');

            simU(i, :) = ocp_solver.get('u', 0);
        else
            % solve ocp and get next control input
            ocp_solver.set('constr_lbx', simX(i, :), 0);
            ocp_solver.set('constr_ubx', simX(i, :), 0);
            ocp_solver.solve();
            status = ocp_solver.get('status');
            simU(i, :) = ocp_solver.get('u', 0);
            t(i) = ocp_solver.get('time_tot');
        end

        if ~ismember(status, [0, 2, 5])
            error(['acados returned status ', num2str(status), '. Exiting.']);
        end

        % simulate system
        integrator.set('x', simX(i, :));
        integrator.set('u', simU(i, :));
        sim_status = integrator.solve();
        simX(i+1, :) = integrator.get('x');
    end

    % evaluate timings
    if ~strcmp(algorithm, 'SQP')
        % scale to milliseconds
        t_preparation = t_preparation * 1000;
        t_feedback = t_feedback * 1000;
        disp(['Computation time in preparation phase in ms: min ', num2str(min(t_preparation)), ...
              ' median ', num2str(median(t_preparation)), ' max ', num2str(max(t_preparation))]);
        disp(['Computation time in feedback phase in ms: min ', num2str(min(t_feedback)), ...
              ' median ', num2str(median(t_feedback)), ' max ', num2str(max(t_feedback))]);
    else
        % scale to milliseconds
        t = t * 1000;
        disp(['Computation time in ms: min ', num2str(min(t)), ...
              ' median ', num2str(median(t)), ' max ', num2str(max(t))]);
    end

    if 1 % plot
        figure;
        subplot(2,1,1);
        plot(0:Nsim, simX);
        xlim([0 Nsim]);
        legend('p', 'theta', 'v', 'omega');
        subplot(2,1,2);
        plot(0:Nsim-1, simU);
        xlim([0 Nsim]);
        legend('F');
        title(['algorithm: ', algorithm, '-', num2str(as_rti_iter)]);

        % if is_octave()
        %     waitforbuttonpress;
        % end
    end
end
