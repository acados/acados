%
% Copyright (c) The acados authors.
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

function [ gnsf ] = determine_trivial_gnsf_transcription(model, print_info)

%% Description
% this function takes a model of an implicit ODE/ index-1 DAE and sets up
% an equivalent model in the GNSF structure, with empty linear output
% system and trivial model matrices, i.e. A, B, E, c are zeros, and C is
% eye. - no structure is exploited

% import CasADi
import casadi.*

casadi_opts = struct('mex', false, 'casadi_int', 'int', 'casadi_real', 'double');
check_casadi_version();

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
if isa(x(1), 'casadi.SX')
    isSX = true;
else
    disp('GNSF detection only works for SX CasADi type!!!');
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
        u = SX.sym('u', 0, 0);
    else
        u = MX.sym('u', 0, 0);
    end
    nu = 0;
end
% z
if isfield(model, 'sym_z')
    z = model.sym_z;
    nz = length(z);
else
    if isSX
        z = SX.sym('z', 0, 0);
    else
        z = MX.sym('z', 0, 0);
    end
    nz = 0;
end
% p
if isfield(model, 'sym_p')
    p = model.sym_p;
    np = length(p);
else
    if isSX
        p = SX.sym('p', 0, 0);
    else
        p = MX.sym('p', 0, 0);
    end
    np = 0;
end

% avoid SX of size 0x1
if any(size(u) == 0)
    u = SX.sym('u', 0, 0);
    nu = 0;
end
if any(size(z) == 0)
    z = SX.sym('z', 0, 0);
    nz = 0;
end
if any(size(p) == 0)
    p = SX.sym('p', 0, 0);
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
gnsf.idx_perm_x = 1:nx; % matlab-style
gnsf.ipiv_x = idx_perm_to_ipiv(gnsf.idx_perm_x); % blasfeo-style
gnsf.idx_perm_z = 1:nz;
gnsf.ipiv_z = idx_perm_to_ipiv(gnsf.idx_perm_z);
gnsf.idx_perm_f = 1:(nx+nz);
gnsf.ipiv_f = idx_perm_to_ipiv(gnsf.idx_perm_f);

gnsf.nontrivial_f_LO = 0;

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
