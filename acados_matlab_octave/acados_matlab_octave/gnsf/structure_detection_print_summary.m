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
%   Author: Jonathan Frey: jonathanpaulfrey(at)gmail.com

function structure_detection_print_summary(gnsf, model)

%% Description
% this function prints the most important info after determining a GNSF
% reformulation of the implicit model "initial_model" into "gnsf", which is
% equivalent to the "reordered_model".

% % GNSF
% get dimensions
nx  = gnsf.nx;
nu  = gnsf.nu;
nz  = gnsf.nz;

nx1 = gnsf.nx1;
nx2 = gnsf.nx2;

nz1 = gnsf.nz1;
nz2 = gnsf.nz2;

np = gnsf.np;
n_out = gnsf.n_out;
ny = gnsf.ny;
nuhat = gnsf.nuhat;

%
f_impl_expr = model.dyn_expr_f;
n_nodes_initial = model.dyn_expr_f.n_nodes();
%x_old = model.x;
%f_impl_old = model.f_impl_expr;

x = gnsf.x;
z = gnsf.z;

phi_current = gnsf.phi_expr;

%% PRINT SUMMARY -- STRUCHTRE DETECTION
disp(' ');
disp('*********************************************************************************************');
disp(' ');
disp('******************        SUCCESS: GNSF STRUCTURE DETECTION COMPLETE !!!      ***************');
disp(' ');
disp('*********************************************************************************************');
disp(' ');
disp(['========================= STRUCTURE DETECTION SUMMARY ====================================']);
disp(' ');
disp('-------- Nonlinear Static Feedback type system --------');
disp(' ');
disp(' successfully transcribed dynamic system model into GNSF structure ');
disp(' ');
disp(['reduced dimension of nonlinearity phi from        ', sprintf('%6s', num2str(nx+nz)),      ' to ', sprintf('%6s', num2str(gnsf.n_out))]);
disp(' ');
disp(['reduced input dimension of nonlinearity phi from  ', sprintf('%6s', num2str(2*nx+nz+nu)), ' to ', sprintf('%6s', num2str(gnsf.ny + gnsf.nuhat))]);
disp(' ');
disp(['reduced number of nodes in CasADi expression of']);
disp(['nonlinearity phi from                             '...
    , sprintf('%6s', num2str(n_nodes_initial)), ' to ', sprintf('%6s', num2str(phi_current.n_nodes()))]);
disp(' ');
disp('----------- Linear Output System (LOS) ---------------');
if gnsf.nx2 + gnsf.nz2 >0
    disp(' ');
    disp(['introduced Linear Output System of size           ', sprintf('%6s', num2str(gnsf.nx2 + gnsf.nz2)),'']);
    disp(' ');
end
if gnsf.nx2 >0
    disp('consisting of the states:');
    disp(' ');
    disp(x(gnsf.nx1+1:gnsf.nx));
    disp(' ');
end
if gnsf.nz2 >0
    disp('and algebraic variables:');
    disp(' ');
    disp(z(gnsf.nz1+1:gnsf.nz));
    disp(' ');
end
if gnsf.purely_linear == 1
    disp(' ');
    disp('Model is fully linear!');
    disp(' ');
end

if ~isequal(gnsf.idx_perm_x, [1:nx])
    disp(' ');
    disp('--------------------------------------------------------------------------------------------------');
    disp('NOTE: permuted differential state vector x, such that x_gnsf = x(idx_perm_x) with idx_perm_x =');
    disp(' ');
    disp(gnsf.idx_perm_x)
end

if nz~= 0 && ~isequal(gnsf.idx_perm_z, [1:nz])
    disp(' ');
    disp('--------------------------------------------------------------------------------------------------');
    disp('NOTE: permuted algebraic state vector z, such that z_gnsf = z(idx_perm_z) with idx_perm_z =');
    disp(' ');
    disp(gnsf.idx_perm_z)
end

if ~isequal(gnsf.idx_perm_f, [1:nx+nz])
    disp(' ');
    disp('--------------------------------------------------------------------------------------------------');
    disp('NOTE: permuted rhs expression vector f, such that f_gnsf = f(idx_perm_f) with idx_perm_f =');
    disp(' ');
    disp(gnsf.idx_perm_f)
end

%% print GNSF dimensions
format short
disp('--------------------------------------------------------------------------------------------------------');
disp(' ');
disp('The dimensions of the GNSF reformulated model read as:');
disp(' ');
%T_dim = table(nx, nu, nz, np, nx1, nz1, n_out, ny, nuhat);
%disp( T_dim )
disp(['nx    ' num2str(nx)]);
disp(['nu    ' num2str(nu)]);
disp(['nz    ' num2str(nz)]);
disp(['np    ' num2str(np)]);
disp(['nx1   ' num2str(nx1)]);
disp(['nz1   ' num2str(nz1)]);
disp(['n_out ' num2str(n_out)]);
disp(['ny    ' num2str(ny)]);
disp(['nuhat ' num2str(nuhat)]);

format short e

end
