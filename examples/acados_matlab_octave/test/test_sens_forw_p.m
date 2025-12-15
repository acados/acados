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

%% test of native matlab interface (ERK + forward sens + param sens)

addpath('../pendulum_on_cart_model/');

%% integrator / method
method = 'erk';   % only ERK

%% arguments
compile_interface = 'auto';
sens_forw   = 'true';
sens_forw_p = 'true';   % param forward sensitivities
jac_reuse   = 'true';
num_stages  = 3;
num_steps   = 4;
newton_iter = 3;

Ts         = 0.1;
x0         = [1e-1; 1e0; 2e-1; 2e0];   
u          = 0;      
FD_epsilon = 1e-6;

%% model 
model = pendulum_on_cart_model_with_param();
model_name = ['pendulum_sens_p' method];

nx = model.nx;
nu = model.nu;

% detect parameters
if isfield(model, 'sym_p')
    np = length(model.sym_p);
else
    np = 0;
end

if np > 0
    p0 = 1;  
end

%% acados sim model
sim_model = acados_sim_model();
sim_model.set('T', Ts);
sim_model.set('name', model_name);

sim_model.set('sym_x', model.sym_x);
if isfield(model, 'sym_u')
    sim_model.set('sym_u', model.sym_u);
end
if isfield(model, 'sym_p')
    sim_model.set('sym_p', model.sym_p);
end

% ERK: explicit dynamics
sim_model.set('dyn_type', 'explicit');
sim_model.set('dyn_expr_f', model.dyn_expr_f_expl);

%% acados sim opts
sim_opts = acados_sim_opts();
sim_opts.set('compile_interface', compile_interface);
sim_opts.set('num_stages', num_stages);
sim_opts.set('num_steps', num_steps);
sim_opts.set('newton_iter', newton_iter);
sim_opts.set('method', method);
sim_opts.set('sens_forw', sens_forw);
sim_opts.set('sens_forw_p', sens_forw_p);   
sim_opts.set('jac_reuse', jac_reuse);

%% acados sim
sim_solver = acados_sim(sim_model, sim_opts);

% set nominal state, input, parameter
sim_solver.set('x', x0);
sim_solver.set('u', u);
if np > 0
    sim_solver.set('p', p0);
end

% solve once with analytic sensitivities on
sim_solver.solve();

xn         = sim_solver.get('xn');
S_forw_ind = sim_solver.get('S_forw');
if np > 0
    S_p_ind = sim_solver.get('S_p');   
else
    S_p_ind = [];
end

%% --- Forward sensitivities S_forw vs finite differences (x,u) ---

% turn off analytic forward/adjoint for FD runs
sim_opts.set('sens_forw', 'false');
sim_opts.set('sens_adj',  'false');

S_forw_fd = zeros(nx, nx+nu);

% FD w.r.t x
for ii = 1:nx
    dx       = zeros(nx, 1);
    dx(ii)   = 1.0;

    sim_solver.set('x', x0 + FD_epsilon*dx);
    sim_solver.set('u', u);
    if np > 0
        sim_solver.set('p', p0);
    end
    sim_solver.solve();

    xn_tmp = sim_solver.get('xn');
    S_forw_fd(:, ii) = (xn_tmp - xn) / FD_epsilon;
end

% FD w.r.t u
for ii = 1:nu
    du       = zeros(nu, 1);
    du(ii)   = 1.0;

    sim_solver.set('x', x0);
    sim_solver.set('u', u + FD_epsilon*du);
    if np > 0
        sim_solver.set('p', p0);
    end
    sim_solver.solve();

    xn_tmp = sim_solver.get('xn');
    S_forw_fd(:, nx + ii) = (xn_tmp - xn) / FD_epsilon;
end

% error S_forw
error_abs_Sforw = max(max(abs(S_forw_fd - S_forw_ind)));
disp(' ');
disp(['integrator:  ' method]);
disp(['error S_forw (FD vs analytic):   ' num2str(error_abs_Sforw)]);
if error_abs_Sforw > 1e-6
    error(['test_sens_forw FAIL: forward sensitivities error too large for integrator ' method]);
end

%% --- Param sensitivities S_p vs finite differences (p) ---

if np > 0
    % Reset state, input, and parameter to nominal
    sim_solver.set('x', x0);
    sim_solver.set('u', u);
    sim_solver.set('p', p0);
    sim_solver.solve();
    xn_nom = sim_solver.get('xn');

    S_p_fd = zeros(nx, np);

    for jj = 1:np
        dp      = zeros(np, 1);
        dp(jj)  = 1.0;

        sim_solver.set('x', x0);
        sim_solver.set('u', u);
        sim_solver.set('p', p0 + FD_epsilon*dp);
        sim_solver.solve();

        xn_tmp = sim_solver.get('xn');
        S_p_fd(:, jj) = (xn_tmp - xn_nom) / FD_epsilon;
    end

    error_abs_Sp = max(max(abs(S_p_fd - S_p_ind)));
    disp(['error S_p (FD vs analytic):     ' num2str(error_abs_Sp)]);
    if error_abs_Sp > 1e-6
        error(['test_sens_forw_p FAIL: param sensitivities error too large for integrator ' method]);
    end
else
    disp('Model has no parameters (np = 0), skipping S_p test.');
end

fprintf('\nTEST_FORWARD_SENSITIVITIES + PARAM_SENS (ERK): success!\n\n');