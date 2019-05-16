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

function generate_c_code_nonlinear_constr( model, opts )

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

if isfield(model, 'constr_expr_h')
	h = model.constr_expr_h;
	% multipliers for hessian
	nh = length(h);
	if isSX
		lam_h = SX.sym('lam_h', nh, 1);
	else
		lam_h = MX.sym('lam_h', nh, 1);
	end
	% generate jacobians
	jac_ux = jacobian(h, [u; x]);
	% generate adjoint
	adj_ux = jtimes(h, [u; x], lam_h, true);
	% generate hessian
	hess_ux = jacobian(adj_ux, [u; x]);
	% Set up functions
	if (strcmp(model.constr_param_h, 'true'))
		h_fun_jac_ut_xt = Function([model_name,'_constr_h_fun_jac_ut_xt'], {x, u, p}, {h, jac_ux'});
		h_fun_jac_ut_xt_hess = Function([model_name,'_constr_h_fun_jac_ut_xt_hess'], {x, u, lam_h, p}, {h, jac_ux', hess_ux});
	else
		h_fun_jac_ut_xt = Function([model_name,'_constr_h_fun_jac_ut_xt'], {x, u}, {h, jac_ux'});
		h_fun_jac_ut_xt_hess = Function([model_name,'_constr_h_fun_jac_ut_xt_hess'], {x, u, lam_h}, {h, jac_ux', hess_ux});
	end
	% generate C code
	h_fun_jac_ut_xt.generate([model_name,'_constr_h_fun_jac_ut_xt'], casadi_opts);
	h_fun_jac_ut_xt_hess.generate([model_name,'_constr_h_fun_jac_ut_xt_hess'], casadi_opts);
end

if isfield(model, 'constr_expr_h_e')
	h_e = model.constr_expr_h_e;
	% multipliers for hessian
	nh_e = length(h_e);
	if isSX
		lam_h_e = SX.sym('lam_h', nh_e, 1);
	else
		lam_h_e = MX.sym('lam_h', nh_e, 1);
	end
	% generate jacobians
	jac_x_e     = jacobian(h_e, x);
	% generate adjoint (TODO output also adjoint when hessian is computed ?????)
	adj_ux_e = jtimes(h_e, x, lam_h_e, true);
	% generate hessian
	hess_ux_e = jacobian(adj_ux_e, x);
	% Set up functions
	if (strcmp(model.constr_param_h, 'true'))
		h_e_fun_jac_ut_xt = Function([model_name,'_constr_h_e_fun_jac_ut_xt'], {x, p}, {h_e, jac_x_e'});
		h_e_fun_jac_ut_xt_hess = Function([model_name,'_constr_h_e_fun_jac_ut_xt_hess'], {x, lam_h_e, p}, {h_e, jac_x_e', hess_ux_e});
	else
		h_e_fun_jac_ut_xt = Function([model_name,'_constr_h_e_fun_jac_ut_xt'], {x}, {h_e, jac_x_e'});
		h_e_fun_jac_ut_xt_hess = Function([model_name,'_constr_h_e_fun_jac_ut_xt_hess'], {x, lam_h}, {h_e, jac_x_e', hess_ux_e});
	end
	% generate C code
	h_e_fun_jac_ut_xt.generate([model_name,'_constr_h_e_fun_jac_ut_xt'], casadi_opts);
	h_e_fun_jac_ut_xt_hess.generate([model_name,'_constr_h_e_fun_jac_ut_xt_hess'], casadi_opts);
end
