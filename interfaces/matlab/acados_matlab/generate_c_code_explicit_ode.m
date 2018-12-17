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

function generate_c_code_explicit_ode( model, opts )

%% import casadi
import casadi.*

if CasadiMeta.version()=='3.4.0'
	% casadi 3.4
	casadi_opts = struct('mex', false, 'casadi_int', 'int', 'casadi_real', 'double');
else
	% old casadi versions
	error('Please download and install Casadi 3.4.0 to ensure compatibility with acados')
end

% TODO check for hessian

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
if isfield(model, 'sym_u')
    u = model.sym_u;
	nu = length(u);
else
    if isSX
        u = SX.sym('u',0, 0);
    else
        u = MX.sym('u',0, 0);
    end
	nu = 0;
end

f_expl = model.dyn_expr;

model_name = model.name;

%% set up functions to be exported
if isSX
    Sx = SX.sym('Sx',nx,nx);
    Sp = SX.sym('Sp',nx,nu);
    lambdaX = SX.sym('lambdaX',nx,1);
    vdeX = SX.zeros(nx,nx);
    vdeP = SX.zeros(nx,nu) + jacobian(f_expl,u);
else
    Sx = MX.sym('Sx',nx,nx);
    Sp = MX.sym('Sp',nx,nu);
    lambdaX = MX.sym('lambdaX',nx,1);
    vdeX = MX.zeros(nx,nx);
    vdeP = MX.zeros(nx,nu) + jacobian(f_expl,u);
end
expl_ode_fun = Function([model_name,'_expl_ode_fun'],{x,u},{f_expl});
% TODO: Polish: get rid of SX.zeros

vdeX = vdeX + jtimes(f_expl,x,Sx);

vdeP = vdeP + jtimes(f_expl,x,Sp);

expl_vde_for = Function([model_name,'_expl_vde_for'],{x,Sx,Sp,u},{f_expl,vdeX,vdeP});

adj = jtimes(f_expl,[x;u],lambdaX,true);

expl_vde_adj = Function([model_name,'_expl_vde_adj'],{x,lambdaX,u},{adj});

S_forw = vertcat(horzcat(Sx, Sp), horzcat(zeros(nu,nx), eye(nu)));
hess = S_forw.'*jtimes(adj,[x;u],S_forw);
hess2 = [];
for j = 1:nx+nu
    for i = j:nx+nu
        hess2 = [hess2; hess(i,j)];
    end
end

expl_ode_hes = Function([model_name,'_expl_ode_hes'],{x,Sx,Sp,lambdaX,u},{adj,hess2});

%% generate C code
expl_ode_fun.generate([model_name,'_expl_ode_fun'], casadi_opts);
expl_vde_for.generate([model_name,'_expl_vde_for'], casadi_opts);
expl_vde_adj.generate([model_name,'_expl_vde_adj'], casadi_opts);
expl_ode_hes.generate([model_name,'_expl_ode_hes'], casadi_opts);

end
