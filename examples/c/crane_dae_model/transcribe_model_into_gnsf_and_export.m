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


%% Description
% This script analyzes a CasADi expression f_impl in the symbolic CasADi
% variables x, xdot, u, z, which all togehther represent an implicit ODE/
% index-1 DAE.
% The expression and the variables should be provided as in the example
% file, export_pendulum_dae_model;
% It will create a struct "gnsf" containing all information needed to use
% it with the gnsf integrator in acados and can generate the neccessary C
% functions.
% Additionally it will create the struct "reordered_model" which contains
% the permuted state vector and permuted f_impl, in which additionally some
% functions, which were made part of the linear output system of the gnsf,
% have changed signs.
% The C functions to simulate the system as an implicit ODE can also be
% generated

clc;
clear VARIABLES;
close all;

addpath('../../../experimental/interfaces/matlab/external_function_generation/sim/')

%% Set options
print_info = 1;
check_E_invertibility = 1;

generate_reordered_model = 1;
generate_gnsf_model = 1;
generate_hess = 0;

transcribe_opts = struct('print_info', print_info, 'check_E_invertibility',...
    check_E_invertibility, 'generate_reordered_model', generate_reordered_model, ...
    'generate_gnsf_model', generate_gnsf_model);
transcribe_opts.generate_hess = generate_hess;


%% define f_impl
model = export_crane_dae_model();

%% transcribe model into gnsf & export
[ gnsf, reordered_model] = detect_gnsf_structure(model, transcribe_opts);
