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

clc;
clearvars;
close all;

addpath('../../../interfaces/matlab/sim/')
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
opts.generate_hess = 0;  % set to 1 if you want to use exact hessian propagation

generate_c_code_implicit_ode( model, opts );

%% Explicit Model -- Generate C Code
generate_c_code_explicit_ode( model );