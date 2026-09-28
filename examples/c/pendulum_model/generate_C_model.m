%
% Copyright (c) The acados authors.
%
% This file is part of acados.
%
% Licensed under the 2-Clause BSD License.

%


clc;
clear VARIABLES;
close all;

addpath('../../../experimental/interfaces/matlab/external_function_generation/sim/')

model = export_pendulum_ode_model();

%% GNSF Model -- detect structure, reorder model, and generate C Code for
%% GNSF model. --> for more advanded users - uncomment this section
% % Reformulate model as GNSF & Reorder x, xdot, z, f_impl, f_expl
% % accordingly
% transcribe_opts.print_info = 1;
% [ gnsf, reordered_model] = detect_gnsf_structure(model, transcribe_opts);
%     % check output of this function to see if/how the states are reordered
% model = reordered_model;
% generate_c_code_gnsf( gnsf );

%% Implicit Model -- Generate C Code
opts.generate_hess = 1;  % set to 1 if you want to use exact hessian propagation

generate_c_code_implicit_ode( model, opts );

%% Explicit Model -- Generate C Code
generate_c_code_explicit_ode( model );
