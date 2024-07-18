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


function generate_c_code_implicit_ode( model, opts, model_dir )

import casadi.*

casadi_opts = struct('mex', false, 'casadi_int', 'int', 'casadi_real', 'double');
check_casadi_version();

%% load model
x = model.x;
u = model.u;
p = model.p;
xdot = model.xdot;
z = model.z;

nx = length(x);
nz = length(z);

% check type
if isa(x(1), 'casadi.SX')
    isSX = true;
else
    isSX = false;
end

model_name = model.name;

if isempty(model.f_impl_expr)
    error("Field `f_impl_expr` is required for integrator type IRK.")
end
f_impl = model.f_impl_expr;

%% generate jacobians
jac_x = jacobian(f_impl, x);
jac_xdot = jacobian(f_impl, xdot);
jac_u = jacobian(f_impl, u);
jac_z = jacobian(f_impl, z);

%% generate hessian
x_xdot_z_u = [x; xdot; z; u];

if isSX
    multiplier  = SX.sym('multiplier', nx + nz);
%    multiply_mat  = SX.sym('multiply_mat', 2*nx+nz+nu, nx + nu);
%    HESS = SX.zeros( length(x_xdot_z_u), length(x_xdot_z_u));
else
    multiplier  = MX.sym('multiplier', nx + nz);
%    multiply_mat  = MX.sym('multiply_mat', 2*nx+nz+nu, nx + nu);
%    HESS = MX.zeros( length(x_xdot_z_u), length(x_xdot_z_u));
end

%% Set up functions
impl_dae_fun = Function([model_name,'_impl_dae_fun'], {x, xdot, u, z, p}, {f_impl});
impl_dae_fun_jac_x_xdot_z = Function([model_name,'_impl_dae_fun_jac_x_xdot_z'], {x, xdot, u, z, p}, {f_impl, jac_x, jac_xdot, jac_z});
impl_dae_jac_x_xdot_u_z = Function([model_name,'_impl_dae_jac_x_xdot_u_z'], {x, xdot, u, z, p}, {jac_x, jac_xdot, jac_u, jac_z});
impl_dae_fun_jac_x_xdot_u = Function([model_name,'_impl_dae_fun_jac_x_xdot_u'], {x, xdot, u, z, p}, {f_impl, jac_x, jac_xdot, jac_u});

%% generate C code in model_dir
return_dir = pwd;
cd(model_dir)

impl_dae_fun.generate([model_name,'_impl_dae_fun'], casadi_opts);
impl_dae_fun_jac_x_xdot_z.generate([model_name,'_impl_dae_fun_jac_x_xdot_z'], casadi_opts);
impl_dae_jac_x_xdot_u_z.generate([model_name,'_impl_dae_jac_x_xdot_u_z'], casadi_opts);
impl_dae_fun_jac_x_xdot_u.generate([model_name,'_impl_dae_fun_jac_x_xdot_u'], casadi_opts);
if opts.generate_hess
    % hessian computed as forward over adjoint !!!
    ADJ = jtimes(f_impl, x_xdot_z_u, multiplier, true);
    HESS = jacobian(ADJ, x_xdot_z_u, struct('symmetric', isSX));

    %HESS_multiplied = multiply_mat' * HESS * multiply_mat;

    %HESS = jtimes(ADJ, x_xdot_z_u, multiply_mat);
    %HESS_multiplied = multiply_mat' * HESS;

    %HESS_multiplied = HESS_multiplied.simplify();
    %HESS_multiplied = HESS; % do the multiplication in BLASFEO !!!
    %    impl_dae_hess = Function([model_name,'_impl_dae_hess'],  {x, xdot, u, z, multiplier, multiply_mat, p}, {HESS_multiplied});

    impl_dae_hess = Function([model_name,'_impl_dae_hess'],...
                             {x, xdot, u, z, multiplier, p}, {HESS});
    impl_dae_hess.generate([model_name,'_impl_dae_hess'], casadi_opts);
end

cd(return_dir);

end
