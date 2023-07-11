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



function [estimator] = setup_estimator(model)

N = model.N;
T = N * model.h;

nlp_solver = 'sqp'; % sqp, sqp_rti
qp_solver = 'partial_condensing_hpipm';
    % full_condensing_hpipm, partial_condensing_hpipm, full_condensing_qpoases

% integrator type
sim_method = 'erk'; % erk, irk, irk_gnsf

%% model dynamics
nx = model.nx;
nu = model.nu;
ny = model.ny;

ocp_model = acados_ocp_model();
model_name = 'lorentz_model_estimator';

%% acados ocp model
ocp_model.set('name', model_name);
ocp_model.set('T', T);

% symbolics
ocp_model.set('sym_x', model.sym_x);
ocp_model.set('sym_u', model.sym_u);
ocp_model.set('dyn_expr_f', model.f_expl_expr);

% cost
ocp_model.set('cost_type', 'linear_ls');
ocp_model.set('cost_type_0', 'linear_ls');
ocp_model.set('cost_type_e','linear_ls');

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

ocp_model.set('cost_Vx', Vx);
ocp_model.set('cost_Vu', Vu);

ocp_model.set('cost_Vx_0', Vx_0);
ocp_model.set('cost_Vu_0', Vu_0);

ocp_model.set('cost_y_ref', yref);
ocp_model.set('cost_y_ref_0', yref_0);

ocp_model.set('cost_W', model.W);
ocp_model.set('cost_W_0', model.W_0);

% dynamics
ocp_model.set('dyn_type', 'explicit');
ocp_model.set('dyn_expr_f', model.f_expl_expr);


%% acados ocp set opts
ocp_opts = acados_ocp_opts();
ocp_opts.set('param_scheme_N', N);
ocp_opts.set('nlp_solver', nlp_solver);
ocp_opts.set('sim_method', sim_method);

ocp_opts.set('sim_method_num_stages', 2);
ocp_opts.set('sim_method_num_steps', 5);

ocp_opts.set('qp_solver', qp_solver);
ocp_opts.set('qp_solver_cond_N', N);
ocp_opts.set('print_level', 0);
ocp_opts.set('ext_fun_compile_flags', '');

%% create ocp solver
estimator = acados_ocp(ocp_model, ocp_opts);

end

