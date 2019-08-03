function [ gnsf ] = determine_trivial_gnsf_transcription(model, print_info)

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
%   Author: Jonathan Frey: jonathanpaulfrey(at)gmail.com

%% Description
% this function takes a model of an implicit ODE/ index-1 DAE and sets up
% an equivalent model in the GNSF structure, with empty linear output
% system and trivial model matrices, i.e. A, B, E, c are zeros, and C is
% eye. - no structure is exploited

% import CasADi
import casadi.*

casadi_version = CasadiMeta.version();
if strcmp(casadi_version(1:3),'3.4') % require casadi 3.4.x
	casadi_opts = struct('mex', false, 'casadi_int', 'int', 'casadi_real', 'double');
else % old casadi versions
	error('Please download and install CasADi version 3.4.x to ensure compatibility with acados')
end

% initial print
disp('*****************************************************************');
disp(' ');
disp(['******      Restructuring ', model.name, ' model    ***********'])
disp(' ');
disp('*****************************************************************');

% load model
f_impl_expr = model.dyn_expr_f;

model_name_prefix = model.name;

% x
x = model.sym_x;
nx = length(x);
% check type
if class(x(1)) == 'casadi.SX'
    isSX = true;
else
    disp('GNSF detection only works for SX CasADi type!!!!!!!!');
	keyboard;
end
% xdot
xdot = model.sym_xdot;
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
% z
if isfield(model, 'sym_z')
    z = model.sym_z;
	nz = length(z);
else
    if isSX
        z = SX.sym('z',0, 0);
    else
        z = MX.sym('z',0, 0);
    end
	nz = 0;
end
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


%% initialize gnsf struct
% dimensions
gnsf = struct('nx', nx, 'nu', nu, 'nz', nz, 'np', np);
gnsf.nx1 = nx;
gnsf.nx2 = 0;
gnsf.nz1 = nz;
gnsf.nz2 = 0;
gnsf.nuhat = nu;
gnsf.ny = 2 * nx + nz;

gnsf.phi_expr = f_impl_expr;
gnsf.A = zeros(nx + nz, nx);
gnsf.B = zeros(nx + nz, nu);
gnsf.E = zeros(nx + nz, nx + nz);
gnsf.c = zeros(nx + nz, 1);
gnsf.C = eye(nx + nz, nx + nz);
gnsf.name = model_name_prefix;

gnsf.x = x;
gnsf.xdot = xdot;
gnsf.z = z;
gnsf.u = u;
gnsf.p = p;

gnsf = determine_input_nonlinearity_function( gnsf );
    
gnsf.A_LO = [];
gnsf.E_LO = [];
gnsf.B_LO = [];
gnsf.c_LO = [];
gnsf.f_lo_expr = [];

% permutation
gnsf.idx_perm_f = 1:nx; % matlab-style
gnsf.ipiv_f = 0:nx-1; % blasfeo-style

check_reformulation(model, gnsf, print_info);
if print_info
    disp(['Success: Set up equivalent GNSF model with trivial matrices']);
    disp(' ');
end

if print_info
    disp('-----------------------------------------------------------------------------------');
    disp(' ');
    disp(['reduced input ny    of phi from  ', num2str(2*nx+nz),'   to  ', num2str( gnsf.ny )]);
    disp(['reduced input nuhat of phi from  ', num2str(nu),'   to  ', num2str( gnsf.nuhat )]);
    disp(' ');
    disp('-----------------------------------------------------------------------------------');
end

end
