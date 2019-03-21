%
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
% This script analyzes a CasADi expression f_impl in the symbolic CasADi
% variables x, xdot, u, z, which all togehther represent an implicit ODE/
% index-1 DAE.
% The expression and the variables should be provided as in the example
% file, export_inverted_pendulum_dae_model;
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
clearvars;
close all;

addpath('../../../interfaces/matlab/external_function_generation/sim/')

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
