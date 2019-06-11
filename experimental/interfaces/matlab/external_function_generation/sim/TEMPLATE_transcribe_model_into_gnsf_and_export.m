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
% This script is a template that shows how the function
% detect_gnsf_structure can be used. For more information look into the
% functions called within this script.

clc;
clearvars;
close all;

% TODO adjust path
% addpath('../../../interfaces/matlab/sim/')

%% Set options
% default is as initialized here
print_info = 1;
check_E_invertibility = 1;
generate_reordered_model = 0;
generate_gnsf_model = 0;
generate_hess = 0;
detect_LOS = 1;

transcribe_opts = struct('print_info', print_info, 'check_E_invertibility',...
    check_E_invertibility, 'generate_reordered_model', generate_reordered_model, ...
    'generate_gnsf_model', generate_gnsf_model);
transcribe_opts.generate_hess = generate_hess;
transcribe_opts.detect_LOS = detect_LOS;
%% define f_impl
% model = TEST_export_inverted_pendulum_dae_model();
% model = TEST_export_inverted_pendulum_dae_model_MX();
% model = TEST_export_crane_dae_test_problem();
% model = TEST_export_stupid_test_problem( );
model = TEST_export_linear_test_problem( );

%% transcribe model into gnsf & export
[ gnsf, reordered_model] = detect_gnsf_structure(model, transcribe_opts);

