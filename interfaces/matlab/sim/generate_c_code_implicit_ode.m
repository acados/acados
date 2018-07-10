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

function generate_c_code_implicit_ode(model)

%% import casadi
import casadi.*

if CasadiMeta.version()=='3.4.0'
	% casadi 3.4
	casadi_opts = struct('mex', false, 'casadi_int', 'int', 'casadi_real', 'double');
else
	% old casadi versions
	error('Please download and install Casadi 3.4.0 to ensure compatibility with acados')
end

%% load model
x = model.x;
xdot = model.xdot;
u = model.u;
z = model.z;
f_impl = model.f_impl_expr;
model_name_prefix = model.name;


%% Set up functions
if isfield(model, 'p')
    p = model.p;
    impl_ode_fun = Function([model_name_prefix,'impl_ode_fun'], {x, xdot, u, z, p}, {f_impl});
    impl_ode_fun_jac_x_xdot = Function([model_name_prefix,'impl_ode_fun_jac_x_xdot'], {x, xdot, u, z, p}, {f_impl, jacobian(f_impl, x), jacobian(f_impl, xdot), jacobian(f_impl, z)});
    impl_ode_jac_x_xdot_u = Function([model_name_prefix,'impl_ode_jac_x_xdot_u'], {x, xdot, u, z, p}, {jacobian(f_impl, x), jacobian(f_impl, xdot), jacobian(f_impl, u), jacobian(f_impl, z)});
    impl_ode_fun_jac_x_xdot_u = Function([model_name_prefix,'impl_ode_fun_jac_x_xdot_u'], {x, xdot, u, z, p}, {f_impl, jacobian(f_impl, x), jacobian(f_impl, xdot), jacobian(f_impl, z)});
else
    impl_ode_fun = Function([model_name_prefix,'impl_ode_fun'], {x, xdot, u, z}, {f_impl});
    impl_ode_fun_jac_x_xdot = Function([model_name_prefix,'impl_ode_fun_jac_x_xdot'], {x, xdot, u, z}, {f_impl, jacobian(f_impl, x), jacobian(f_impl, xdot), jacobian(f_impl, z)});
    impl_ode_jac_x_xdot_u = Function([model_name_prefix,'impl_ode_jac_x_xdot_u'], {x, xdot, u, z}, {jacobian(f_impl, x), jacobian(f_impl, xdot), jacobian(f_impl, u), jacobian(f_impl, z)});
    impl_ode_fun_jac_x_xdot_u = Function([model_name_prefix,'impl_ode_fun_jac_x_xdot_u'], {x, xdot, u, z}, {f_impl, jacobian(f_impl, x), jacobian(f_impl, xdot), jacobian(f_impl, z)});
end

%% generate C code
impl_ode_fun.generate([model_name_prefix,'impl_ode_fun'], casadi_opts);
impl_ode_fun_jac_x_xdot.generate([model_name_prefix,'impl_ode_fun_jac_x_xdot'], casadi_opts);
impl_ode_jac_x_xdot_u.generate([model_name_prefix,'impl_ode_jac_x_xdot_u'], casadi_opts);
impl_ode_fun_jac_x_xdot_u.generate([model_name_prefix,'impl_ode_fun_jac_x_xdot_u'], casadi_opts);


end