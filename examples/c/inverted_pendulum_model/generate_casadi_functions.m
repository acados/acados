%
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

clc;
clear all;
close all;

addpath('/home/andrea/casadi-linux-octave-v3.4.0') 
addpath('~/casadi/swig/octave')

% casadi opts for code generation
if CasadiMeta.version()=='3.4.0'
	% casadi 3.4
	opts = struct('mex', false, 'casadi_int', 'int', 'casadi_real', 'double');
else
	% old casadi versions
	error('Please download and install Casadi 3.4.0')
end

NX = 6
NU = 1
NZ = 5

% define model 
dae = export_inverted_pendulum_dae_model();

% Set up DAE
dae.x_dot   = SX.sym('x_dot', NX, 1);
dae.x_dot   = SX.sym('x_dot', NX, 1);
dae.p       = SX.sym('u', NU, 1);

ode_impl    = dae.f_impl_expr; 

jac_x = SX.zeros(NX+NZ, NX) + jacobian(ode_impl, dae.x);
jac_xdot = SX.zeros(NX+NZ, NX) + jacobian(ode_impl, dae.x_dot);
jac_z = SX.zeros(NX+NZ, NZ) + jacobian(ode_impl, dae.z);
jac_u = SX.zeros(NX+NZ, NU) + jacobian(ode_impl, dae.p);

impl_ode_fun = Function('casadi_impl_ode_fun_pendulum_dae', {dae.x, dae.x_dot, dae.z, dae.p}, {ode_impl});
impl_ode_fun_jac_x_xdot = Function('casadi_impl_ode_fun_jac_x_xdot_z_pendulum_dae', {dae.x, dae.x_dot, dae.z, dae.p}, {ode_impl, jac_x, jac_xdot, jac_z});
impl_ode_fun_jac_x_xdot_u = Function('casadi_impl_ode_fun_jac_x_xdot_z_u_pendulum_dae', {dae.x, dae.x_dot, dae.z, dae.p}, {ode_impl, jac_x, jac_xdot, jac_z, jac_u});
impl_ode_jac_x_xdot_u = Function('casadi_impl_ode_jac_x_xdot_z_u_pendulum_dae', {dae.x, dae.x_dot, dae.z, dae.p}, {jac_x, jac_xdot, jac_z, jac_u});

impl_ode_fun.generate('impl_ode_fun_pendulum_dae', opts);
impl_ode_fun_jac_x_xdot.generate('impl_ode_fun_jac_x_xdot_z_pendulum_dae', opts);
impl_ode_fun_jac_x_xdot_u.generate('impl_ode_fun_jac_x_xdot_z_u_pendulum_dae', opts);
impl_ode_jac_x_xdot_u.generate('impl_ode_jac_x_xdot_z_u_pendulum_dae', opts);
