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

function generate_c_code_ext_cost( model, opts )

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

if isfield(model, 'cost_expr_ext_cost')
	ext_cost = model.cost_expr_ext_cost;
	% generate jacobians
	jac_x = jacobian(ext_cost, x);
	jac_u = jacobian(ext_cost, u);
	% generate hessians
	hes_uu = jacobian(jac_u', u);
	hes_xu = jacobian(jac_u', x);
	hes_ux = jacobian(jac_x', u);
	hes_xx = jacobian(jac_x', x);
	% Set up functions
	if (strcmp(model.cost_param_ext_cost, 'true'))
		ext_cost_jac_hes = Function([model_name,'_cost_ext_cost_jac_hes'], {x, u, p}, {[jac_u'; jac_x'], [hes_uu, hes_xu; hes_ux, hes_xx]});
	else
		ext_cost_jac_hes = Function([model_name,'_cost_ext_cost_jac_hes'], {x, u}, {[jac_u'; jac_x'], [hes_uu, hes_xu; hes_ux, hes_xx]});
	end
	% generate C code
	ext_cost_jac_hes.generate([model_name,'_cost_ext_cost_jac_hes'], casadi_opts);
end

if isfield(model, 'cost_expr_ext_cost_e')
	ext_cost_e = model.cost_expr_ext_cost_e;
	% generate jacobians
	jac_x_e = jacobian(ext_cost_e, x);
	% generate hessians
	hes_xx_e = jacobian(jac_x', x);
	% Set up functions
	if (strcmp(model.cost_param_ext_cost_e, 'true'))
		ext_cost_e_jac_hes = Function([model_name,'_cost_ext_cost_e_jac_hes'], {x, p}, {jac_x_e', hes_xx_e});
	else
		ext_cost_e_jac_hes = Function([model_name,'_cost_ext_cost_e_jac_hes'], {x}, {jac_x_e', hes_xx_e});
	end
	% generate C code
	ext_cost_e_jac_hes.generate([model_name,'_cost_ext_cost_e_jac_hes'], casadi_opts);
end


