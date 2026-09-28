%
% Copyright (c) The acados authors.
%
% This file is part of acados.
%
% Licensed under the 2-Clause BSD License.


function [estimator] = setup_estimator(model, h, N)

    T = N * h;

    %% model dynamics
    nx = length(model.x);
    nu = length(model.u);
    ny = 1;
    nw = nu;               % state noise on parameter

    ocp = AcadosOcp();

    ocp.model = model;
    ocp.solver_options.tf = T;

    %% cost
    % weighting matrices
    Q = 1*eye(nw);
    R = 1;
    P0 = 0.1*eye(nx);
    P0(nx, nx) = 0.001;

    W_0 = blkdiag(R, Q, P0);
    W = blkdiag(R, Q);

    ocp.cost.cost_type = 'LINEAR_LS';
    ocp.cost.cost_type_0 = 'LINEAR_LS';
    ocp.cost.cost_type_e = 'LINEAR_LS';

    nout = ny + nu;
    nout_0 = ny + nu + nx;

    Vx = zeros(nout, nx);
    Vx(1, 1) = 1;

    Vu = zeros(nout, nu);
    Vu(ny+1:ny+nu, :) = eye(nu);

    Vx_0 = zeros(nout_0, nx);
    Vx_0(1:ny, :) = eye(ny, nx);
    Vx_0(ny+nu+1:end, :) = eye(nx);

    Vu_0 = zeros(nout_0, nu);
    Vu_0(ny+1:ny+nu, :) = eye(nu);

    yref = zeros(nout, 1);
    yref_0 = zeros(nout_0, 1);

    ocp.cost.Vx = Vx;
    ocp.cost.Vu = Vu;

    ocp.cost.Vx_0 = Vx_0;
    ocp.cost.Vu_0 = Vu_0;

    ocp.cost.yref = yref;
    ocp.cost.yref_0 = yref_0;

    ocp.cost.W = W;
    ocp.cost.W_0 = W_0;

    %% options
    ocp.solver_options.N_horizon = N;
    ocp.solver_options.nlp_solver_type = 'SQP';
    ocp.solver_options.integrator_type = 'ERK';

    ocp.solver_options.sim_method_num_stages = 2;
    ocp.solver_options.sim_method_num_steps = 5;

    ocp.solver_options.qp_solver = 'PARTIAL_CONDENSING_HPIPM';
    ocp.solver_options.qp_solver_cond_N = N;
    ocp.solver_options.print_level = 0;
    ocp.code_gen_opts.ext_fun_compile_flags = ''; % test backward compatibility of deprecated field
    ocp.code_gen_opts.casadi_code_gen_options.static_aux = true; % for testing

    %% create ocp solver
    estimator = AcadosOcpSolver(ocp);
end
