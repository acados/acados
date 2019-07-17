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

function generate_c_code_disc_dyn( model, opts )

%% import casadi
import casadi.*

casadi_version = CasadiMeta.version();
if strcmp(casadi_version(1:3),'3.4') % require casadi 3.4.x
	casadi_opts = struct('mex', false, 'casadi_int', 'int', 'casadi_real', 'double');
else % old casadi versions
	error('Please download and install CasADi version 3.4.x to ensure compatibility with acados')
end

%% load model
% x
x = model.sym_x;
nx = length(x);
% check type
if class(x(1)) == 'casadi.SX'
    isSX = true;
else
    isSX = false;
end
% u
u = model.sym_u;
nu = length(u);
% p
if isfield(model, 'sym_p')
    p = model.sym_p;
	np = length(p);
else
    if isSX
        p = SX.sym('p',0, 0);
    else
        p = MX.sym('p',0, 0);
    end
	np = 0;
end

model_name = model.name;

if isfield(model, 'dyn_expr_phi')
	phi = model.dyn_expr_phi;
	% assume nx1 = nx !!!
	% multipliers for hessian
	if isSX
		lam = SX.sym('lam', nx, 1);
	else
		lam = MX.sym('lam', nx, 1);
	end
	% generate jacobians
	jac_ux = jacobian(phi, [u; x]);
	% generate adjoint
	adj_ux = jtimes(phi, [u; x], lam, true);
	% generate hessian
	hess_ux = jacobian(adj_ux, [u; x]);
	% Set up functions
	phi_fun_jac_ut_xt = Function([model_name,'_dyn_disc_phi_fun_jac'], {x, u, p}, {phi, jac_ux'});
	phi_fun_jac_ut_xt_hess = Function([model_name,'_dyn_disc_phi_fun_jac_hess'], {x, u, lam, p}, {phi, jac_ux', hess_ux});
	% generate C code
	phi_fun_jac_ut_xt.generate([model_name,'_dyn_disc_phi_fun_jac'], casadi_opts);
	phi_fun_jac_ut_xt_hess.generate([model_name,'_dyn_disc_phi_fun_jac_hess'], casadi_opts);
end
