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
%% THIS IS A TEMPLATE

clc;
clear all;
close all;

%% Set options
print_info = 0;
generate_c_code = 0;


%% define f_impl
[ model ] = export_inverted_pendulum_dae();
% [ model ] = export_crane_dae_test_problem();
% [ model ] = export_stupid_test_problem();

disp(' ');
disp(['restructuring ', model.name, ' model'])
disp(' ');

initial_model = model;

%% Reformulate implicit index-1 DAE into GNSF form
% (Generalized nonlinear static feedback)
[ gnsf ] = define_equivalent_model_in_gnsf_format(model, print_info);

[ gnsf, reordered_model] = reformulate_with_LOS( model, gnsf, print_info);

[ gnsf ] = reformulate_with_invertible_E_mat(gnsf, reordered_model, print_info);

structure_detection_print_summary(gnsf, initial_model, reordered_model);



%% EXPORT C Code
if generate_c_code
    % generate gnsf model
    gnsf_generate_c_code
    % generate implicit model
    implicit_ode_generate_c_code
end
