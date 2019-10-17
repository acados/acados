%
% Copyright 2019 Gianluca Frison, Dimitris Kouzoupis, Robin Verschueren,
% Andrea Zanelli, Niels van Duijkeren, Jonathan Frey, Tommaso Sartor,
% Branimir Novoselnik, Rien Quirynen, Rezart Qelibari, Dang Doan,
% Jonas Koenemann, Yutao Chen, Tobias Schöls, Jonas Schlagenhauf, Moritz Diehl
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

function generate_c_code_implicit_ode( model, opts )

%% import casadi
import casadi.*

casadi_version = CasadiMeta.version();
if strcmp(casadi_version(1:3),'3.4') % require casadi 3.4.x
    casadi_opts = struct('mex', false, 'casadi_int', 'int', 'casadi_real', 'double');
else % old casadi versions
    error('Please download and install CasADi version 3.4.x to ensure compatibility with acados')
end

if isfield(opts, 'generate_hess')
    generate_hess = opts.generate_hess;
else
    generate_hess = 0;
    if opts.print_info
    disp('generate_hess option was not set - default is false')
    end
end

%% load model
x = model.x;
xdot = model.xdot;
u = model.u;
z = model.z;
f_impl = model.f_impl_expr;
model_name = model.name;

%% get model dimensions
nx = length(x);
nu = length(u);
nz = length(z);

%% generate jacobians
jac_x       = jacobian(f_impl, x);
jac_xdot    = jacobian(f_impl, xdot);
jac_u       = jacobian(f_impl, u);
jac_z       = jacobian(f_impl, z);


%% generate hessian
x_xdot_z_u = [x; xdot; z; u];

if class(x(1)) == 'casadi.SX'
    multiplier  = SX.sym('multiplier', length(x) + length(z));
    multiply_mat  = SX.sym('multiply_mat', 2*nx+nz+nu, nx + nu);
    HESS = SX.zeros( length(x_xdot_z_u), length(x_xdot_z_u));
elseif class(x(1)) == 'casadi.MX'
    multiplier  = MX.sym('multiplier', length(x) + length(z));
    multiply_mat  = MX.sym('multiply_mat', 2*nx+nz+nu, nx + nu);
    HESS = MX.zeros( length(x_xdot_z_u), length(x_xdot_z_u));
end

for ii = 1:length(f_impl)
    jac_x_xdot_z = jacobian(f_impl(ii), x_xdot_z_u);
    hess_x_xdot_z = jacobian( jac_x_xdot_z, x_xdot_z_u);
    HESS = HESS + multiplier(ii) * hess_x_xdot_z;
end

HESS = HESS.simplify();
HESS_multiplied = multiply_mat' * HESS * multiply_mat;
HESS_multiplied = HESS_multiplied.simplify();



%% Set up functions
if isfield(model, 'p')
    p = model.p;
    impl_dae_fun = Function([model_name,'_impl_dae_fun'], {x, xdot, u, z, p},...
                             {f_impl});
    impl_dae_fun_jac_x_xdot_z = Function([model_name,'_impl_dae_fun_jac_x_xdot_z'],...
         {x, xdot, u, z, p}, {f_impl, jac_x, jac_xdot, jac_z});
    impl_dae_jac_x_xdot_u_z = Function([model_name,'_impl_dae_jac_x_xdot_u_z'],...
         {x, xdot, u, z, p}, {jac_x, jac_xdot, jac_u, jac_z});
    impl_dae_fun_jac_x_xdot_u_z = Function([model_name,'_impl_dae_fun_jac_x_xdot_u_z'],...
         {x, xdot, u, z, p},...
         {f_impl, jac_x, jac_xdot, jac_u});
    impl_dae_hess = Function([model.name,'_impl_dae_hess'], ...
         {x, xdot, u, z, multiplier, multiply_mat, p}, {HESS_multiplied});
else
    impl_dae_fun = Function([model_name,'_impl_dae_fun'],...
                 {x, xdot, u, z}, {f_impl});
    impl_dae_fun_jac_x_xdot_z = Function([model_name,'_impl_dae_fun_jac_x_xdot_z'],...
         {x, xdot, u, z}, {f_impl, jac_x, jac_xdot, jac_z});
    impl_dae_jac_x_xdot_u_z = Function([model_name,'_impl_dae_jac_x_xdot_u_z'], {x, xdot, u, z},...
             {jac_x, jac_xdot, jac_u, jac_z});
    impl_dae_fun_jac_x_xdot_u_z = Function([model_name,'_impl_dae_fun_jac_x_xdot_u_z'],...
         {x, xdot, u, z}, {f_impl, jac_x, jac_xdot, jac_u});
    impl_dae_hess = Function([model.name,'_impl_dae_hess'], ...
        {x, xdot, u, z, multiplier, multiply_mat}, {HESS_multiplied});
end

%% generate C code
    if ~exist('c_generated_code', 'dir')
        mkdir('c_generated_code');
    end
    cd 'c_generated_code'
    model_dir = [model_name, '_model'];
    if ~exist(model_dir, 'dir')
        mkdir(model_dir);
    end
    model_dir_location = ['./', model_dir];
    cd(model_dir_location);

    fun_name = [model_name, '_impl_dae_fun'];
    impl_dae_fun.generate(fun_name, casadi_opts);

    fun_name = [model_name, '_impl_dae_fun_jac_x_xdot_z'];
    impl_dae_fun_jac_x_xdot_z.generate(fun_name, casadi_opts);
    
    fun_name = [model_name, '_impl_dae_jac_x_xdot_u_z'];
    impl_dae_jac_x_xdot_u_z.generate(fun_name, casadi_opts);

    fun_name = [model_name, '_impl_dae_fun_jac_x_xdot_u_z'];
    impl_dae_fun_jac_x_xdot_u_z.generate(fun_name, casadi_opts);

    if generate_hess
        fun_name = [model_name, '_impl_dae_hess'];
        impl_dae_hess.generate(fun_name, casadi_opts);
    end

    cd '../..'
end
