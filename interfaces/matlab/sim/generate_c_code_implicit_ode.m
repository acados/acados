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

if CasadiMeta.version()=='3.4.0'
	% casadi 3.4
	casadi_opts = struct('mex', false, 'casadi_int', 'int', 'casadi_real', 'double');
else
	% old casadi versions
	error('Please download and install Casadi 3.4.0 to ensure compatibility with acados')
end

if nargin > 1
    if isfield(opts, 'generate_hess')
        generate_hess = opts.generate_hess;
    else
        generate_hess = 0;
        if opts.print_info
        disp('generate_hess option was not set - default is false')
        end
    end
else
    generate_hess = 0;
end

%% load model
x = model.x;
if class(x(1)) == 'casadi.SX'
    isSX = true;
else
    isSX = false;
end
xdot = model.xdot;
u = model.u;
if isfield(model, 'z')
    z = model.z;
else
    if class(x(1)) == 'casadi.SX'
        z = SX.sym('z',0, 0);
    else
        z = MX.sym('z',0, 0);
    end
end
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

if isSX
    multiplier  = SX.sym('multiplier', length(x) + length(z));
    multiply_mat  = SX.sym('multiply_mat', 2*nx+nz+nu, nx + nu);
    HESS = SX.zeros( length(x_xdot_z_u), length(x_xdot_z_u));
else
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
% TODO(oj): fix namings such that jac_z is contained!
if isfield(model, 'p')
    p = model.p;
    impl_ode_fun = Function([model_name,'_impl_ode_fun'], {x, xdot, u, z, p},...
                             {f_impl});
    impl_ode_fun_jac_x_xdot = Function([model_name,'_impl_ode_fun_jac_x_xdot'],...
         {x, xdot, u, z, p}, {f_impl, jac_x, jac_xdot, jac_z});
    impl_ode_jac_x_xdot_u = Function([model_name,'_impl_ode_jac_x_xdot_u'],...
         {x, xdot, u, z, p}, {jac_x, jac_xdot, jac_u, jac_z});
    impl_ode_fun_jac_x_xdot_u = Function([model_name,'_impl_ode_fun_jac_x_xdot_u'],...
         {x, xdot, u, z, p},...
         {f_impl, jac_x, jac_xdot, jac_u});
    impl_ode_hess = Function([model.name,'_impl_ode_hess'], ...
         {x, xdot, u, z, multiplier, multiply_mat, p}, {HESS_multiplied});
else
    impl_ode_fun = Function([model_name,'_impl_ode_fun'],...
                 {x, xdot, u, z}, {f_impl});
    impl_ode_fun_jac_x_xdot = Function([model_name,'_impl_ode_fun_jac_x_xdot'],...
         {x, xdot, u, z}, {f_impl, jac_x, jac_xdot, jac_z});
    impl_ode_jac_x_xdot_u = Function([model_name,'_impl_ode_jac_x_xdot_u'], {x, xdot, u, z},...
             {jac_x, jac_xdot, jac_u, jac_z});
    impl_ode_fun_jac_x_xdot_u = Function([model_name,'_impl_ode_fun_jac_x_xdot_u'],...
         {x, xdot, u, z}, {f_impl, jac_x, jac_xdot, jac_u});
    impl_ode_hess = Function([model.name,'_impl_ode_hess'], ...
        {x, xdot, u, z, multiplier, multiply_mat}, {HESS_multiplied});
end

%% generate C code
impl_ode_fun.generate([model_name,'_impl_ode_fun'], casadi_opts);
impl_ode_fun_jac_x_xdot.generate([model_name,'_impl_ode_fun_jac_x_xdot'], casadi_opts);
impl_ode_jac_x_xdot_u.generate([model_name,'_impl_ode_jac_x_xdot_u'], casadi_opts);
impl_ode_fun_jac_x_xdot_u.generate([model_name,'_impl_ode_fun_jac_x_xdot_u'], casadi_opts);
if generate_hess
    impl_ode_hess.generate([model_name,'_impl_ode_hess'], casadi_opts);
end
% keyboard

end
