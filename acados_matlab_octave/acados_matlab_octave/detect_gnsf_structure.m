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

function model = detect_gnsf_structure(model, transcribe_opts)

%% Description
% This function takes a CasADi implicit ODE or index-1 DAE model "model" 
% consisting of a CasADi expression f_impl in the symbolic CasADi
% variables x, xdot, u, z, (and possibly parameters p), which are also part
% of the model, as well as a model name.
% It will create a struct "gnsf" containing all information needed to use
% it with the gnsf integrator in acados.
% Additionally it will create the struct "reordered_model" which contains
% the permuted state vector and permuted f_impl, in which additionally some
% functions, which were made part of the linear output system of the gnsf,
% have changed signs.

% Options: transcribe_opts is a Matlab struct consisting of booleans:
%   print_info: if extensive information on how the model is processed
%       is printed to the console.
%   check_E_invertibility: if the transcription method should check if the
%       assumption that the main blocks of the matrix gnsf.E are invertible
%       holds. If not, the method will try to reformulate the gnsf model
%       with a different model, such that the assumption holds.

acados_root_dir = getenv('ACADOS_INSTALL_DIR');
addpath(fullfile(acados_root_dir, 'interfaces', 'acados_matlab_octave', 'gnsf'));

%% load transcribe_opts
if ~exist('transcribe_opts', 'var')
    disp('WARNING: GNSF structure detection called without transcribe_opts');
    disp(' using default settings');
    disp('');
    transcribe_opts = struct;
end

if isfield(transcribe_opts, 'print_info')
    print_info = transcribe_opts.print_info;
else
    print_info = 1;
    disp('print_info option was not set - default is true')
end

if isfield(transcribe_opts, 'detect_LOS')
    detect_LOS = transcribe_opts.detect_LOS;
else
    detect_LOS = 1;
    if print_info
    disp('detect_LOS option was not set - default is true')
    end
end


if isfield(transcribe_opts, 'check_E_invertibility')
    check_E_invertibility = transcribe_opts.check_E_invertibility;
else
    check_E_invertibility = 1;
    if print_info
    disp('check_E_invertibility option was not set - default is true')
    end
end


%% Reformulate implicit index-1 DAE into GNSF form
% (Generalized nonlinear static feedback)
gnsf = determine_trivial_gnsf_transcription( model, print_info );

gnsf = detect_affine_terms_reduce_nonlinearity( gnsf, model, print_info );

if detect_LOS
    gnsf = reformulate_with_LOS( model, gnsf, print_info);
end

if check_E_invertibility
    gnsf = reformulate_with_invertible_E_mat( gnsf, model, print_info);
end

% detect purely linear model
if gnsf.nx1 == 0 && gnsf.nz1 == 0 && gnsf.nontrivial_f_LO == 0
    gnsf.purely_linear = 1;
else
    gnsf.purely_linear = 0;
end

structure_detection_print_summary(gnsf, model);
check_reformulation( model, gnsf, print_info );


%% copy relevant fields from gnsf to model

% dim
model.dim_gnsf_nx1 = gnsf.nx1;
model.dim_gnsf_nx2 = gnsf.nx2;
model.dim_gnsf_nz1 = gnsf.nz1;
model.dim_gnsf_nz2 = gnsf.nz2;
model.dim_gnsf_nuhat = gnsf.nuhat;
model.dim_gnsf_ny = gnsf.ny;
model.dim_gnsf_nout = gnsf.n_out;

% sym
model.sym_gnsf_y = gnsf.y;
model.sym_gnsf_uhat = gnsf.uhat;

% data
model.dyn_gnsf_A = gnsf.A;
model.dyn_gnsf_A_LO = gnsf.A_LO;
model.dyn_gnsf_B = gnsf.B;
model.dyn_gnsf_B_LO = gnsf.B_LO;
model.dyn_gnsf_E = gnsf.E;
model.dyn_gnsf_E_LO = gnsf.E_LO;
model.dyn_gnsf_C = gnsf.C;
model.dyn_gnsf_c = gnsf.c;
model.dyn_gnsf_c_LO = gnsf.c_LO;
model.dyn_gnsf_L_x = gnsf.L_x;
model.dyn_gnsf_L_xdot = gnsf.L_xdot;
model.dyn_gnsf_L_z = gnsf.L_z;
model.dyn_gnsf_L_u = gnsf.L_u;
model.dyn_gnsf_idx_perm_x = gnsf.idx_perm_x;
model.dyn_gnsf_ipiv_x = gnsf.ipiv_x;
model.dyn_gnsf_idx_perm_z = gnsf.idx_perm_z;
model.dyn_gnsf_ipiv_z = gnsf.ipiv_z;
model.dyn_gnsf_idx_perm_f = gnsf.idx_perm_f;
model.dyn_gnsf_ipiv_f = gnsf.ipiv_f;

% flags
model.dyn_gnsf_nontrivial_f_LO = gnsf.nontrivial_f_LO;
model.dyn_gnsf_purely_linear = gnsf.purely_linear;

% casadi expr
model.dyn_gnsf_expr_phi = gnsf.phi_expr;
model.dyn_gnsf_expr_f_lo = gnsf.f_lo_expr;


end

