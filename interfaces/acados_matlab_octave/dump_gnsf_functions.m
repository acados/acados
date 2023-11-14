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


function dump_gnsf_functions(model)

acados_folder = getenv('ACADOS_INSTALL_DIR');
addpath(fullfile(acados_folder, 'external', 'jsonlab'))


%% import casadi
import casadi.*

casadi_version = CasadiMeta.version();
if ~(strcmp(casadi_version(1:3),'3.5') || strcmp(casadi_version(1:3),'3.6'))
    warning('CasADi serialization requires CasADi >= 3.5, you are using: %s.', casadi_version);
end

%% import models
% model matrices
A  = model.dyn_gnsf_A;
B  = model.dyn_gnsf_B;
C  = model.dyn_gnsf_C;
E  = model.dyn_gnsf_E;
c  = model.dyn_gnsf_c;

L_x    = model.dyn_gnsf_L_x;
L_z    = model.dyn_gnsf_L_z;
L_xdot = model.dyn_gnsf_L_xdot;
L_u    = model.dyn_gnsf_L_u;

A_LO = model.dyn_gnsf_A_LO;
E_LO = model.dyn_gnsf_E_LO;
B_LO = model.dyn_gnsf_B_LO;
c_LO = model.dyn_gnsf_c_LO;

% state permutation vector: x_gnsf = dvecpe(x, ipiv)
ipiv_x = model.dyn_gnsf_ipiv_x;
idx_perm_x = model.dyn_gnsf_idx_perm_x;
ipiv_z = model.dyn_gnsf_ipiv_z;
idx_perm_z = model.dyn_gnsf_idx_perm_z;
ipiv_f = model.dyn_gnsf_ipiv_f;
idx_perm_f = model.dyn_gnsf_idx_perm_f;

% CasADi variables and expressions
% x
x = model.sym_x;
x1 = x(model.dyn_gnsf_idx_perm_x(1:model.dim_gnsf_nx1));
% check type
if isa(x(1), 'casadi.SX')
    isSX = true;
else
    isSX = false;
end
% xdot
xdot = model.sym_xdot;
x1dot = xdot(model.dyn_gnsf_idx_perm_x(1:model.dim_gnsf_nx1));
% u
if isfield(model, 'sym_u')
    u = model.sym_u;
else
    if isSX
        u = SX.sym('u',0, 0);
    else
        u = MX.sym('u',0, 0);
    end
end
% z
if isfield(model, 'sym_z')
    z = model.sym_z;
    z1 = model.sym_z(model.dyn_gnsf_idx_perm_z(1:model.dim_gnsf_nz1));
else
    if isSX
        z = SX.sym('z',0, 0);
        z1 = SX.sym('z1',0, 0);
    else
        z = MX.sym('z',0, 0);
        z1 = MX.sym('z1',0, 0);
    end
end
% p
if isfield(model, 'sym_p')
    p = model.sym_p;
else
    if isSX
        p = SX.sym('p',0, 0);
    else
        p = MX.sym('p',0, 0);
    end
end
% y
y = model.sym_gnsf_y;
% uhat
uhat = model.sym_gnsf_uhat;

% expressions
phi = model.dyn_gnsf_expr_phi;
f_lo = model.dyn_gnsf_expr_f_lo;

nontrivial_f_LO = model.dyn_gnsf_nontrivial_f_LO;
purely_linear = model.dyn_gnsf_purely_linear;

% name
model_name = model.name

%% generate functions
if ~purely_linear
    jac_phi_y = jacobian(phi,y);
    jac_phi_uhat = jacobian(phi,uhat);

    phi_fun = Function([model_name,'_gnsf_phi_fun'], {y, uhat, p}, {phi});
    phi_fun_jac_y = Function([model_name,'_gnsf_phi_fun_jac_y'], {y, uhat, p}, {phi, jac_phi_y});
    phi_jac_y_uhat = Function([model_name,'_gnsf_phi_jac_y_uhat'], {y, uhat, p}, {jac_phi_y, jac_phi_uhat});

    if nontrivial_f_LO
        f_lo_fun_jac_x1k1uz = Function([model_name,'_gnsf_f_lo_fun_jac_x1k1uz'], {x1, x1dot, z1, u, p}, ...
            {f_lo, [jacobian(f_lo,x1), jacobian(f_lo,x1dot), jacobian(f_lo,u), jacobian(f_lo,z1)]});
    end
end

% get_matrices function
dummy = x(1);
get_matrices_fun = Function([model_name,'_gnsf_get_matrices_fun'], {dummy},...
     {A, B, C, E, L_x, L_xdot, L_z, L_u, A_LO, c, E_LO, B_LO,...
      nontrivial_f_LO, purely_linear, ipiv_x, ipiv_z, c_LO});

% dump functions to json
out = struct();
out.phi_fun = phi_fun.serialize();
out.phi_fun_jac_y = phi_fun_jac_y.serialize();
out.phi_jac_y_uhat = phi_jac_y_uhat.serialize();
if exist('f_lo_fun_jac_x1k1uz', 'var')
out.f_lo_fun_jac_x1k1uz = f_lo_fun_jac_x1k1uz.serialize();
end
out.get_matrices_fun = get_matrices_fun.serialize();
out.casadi_version = casadi_version;

json_filename = [model.name '_gnsf_functions.json'];
json_string = savejson('', out, 'ForceRootName', 0);

fid = fopen(json_filename, 'w');
if fid == -1, error('Cannot create json file'); end
fwrite(fid, json_string, 'char');
fclose(fid);

disp(['succesfully dumped gnsf model into ', json_filename])
end

