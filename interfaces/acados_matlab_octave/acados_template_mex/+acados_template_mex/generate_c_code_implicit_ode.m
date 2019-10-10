%   This file is part of acados.
%
%   acados is free software; you can redistribute it and/or
%   modify it under the terms of the GNU Lesser General Public
%   License as published by the Free Software Foundation; either
%   version 3 of the License, or (at your option) any later version.
%
%   acados is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
%   Lesser General Public License for more details.
%
%   You should have received a copy of the GNU Lesser General Public
%   License along with acados; if not, write to the Free Software Foundation,
%   Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
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
