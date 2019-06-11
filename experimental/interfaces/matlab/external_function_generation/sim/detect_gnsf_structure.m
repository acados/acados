function [ gnsf, reordered_model] = detect_gnsf_structure(model, transcribe_opts)
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
%   generate_gnsf_model: if the neccessary C functions to simulate the gnsf
%       model with the acados implementation of the GNSF exploiting
%       integrator should be generated.
%   generate_gnsf_model: if the neccessary C functions to simulate the
%       reordered model with the acados implementation of the IRK
%       integrator should be generated.
%   check_E_invertibility: if the transcription method should check if the
%       assumption that the main blocks of the matrix gnsf.E are invertible
%       holds. If not, the method will try to reformulate the gnsf model
%       with a different model, such that the assumption holds.


%% load transcribe_opts
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

if isfield(transcribe_opts, 'generate_gnsf_model')
    generate_gnsf_model = transcribe_opts.generate_gnsf_model;
else
    generate_gnsf_model = 0;
    if print_info
    disp('generate_gnsf_model option was not set - default is false')
    end
end

if isfield(transcribe_opts, 'generate_reordered_model')
    generate_reordered_model = transcribe_opts.generate_reordered_model;
else
    generate_reordered_model = 0;
    if print_info
    disp('generate_reordered_model option was not set - default is false')
    end
end


%% Reformulate implicit index-1 DAE into GNSF form
% (Generalized nonlinear static feedback)
gnsf = determine_trivial_gnsf_transcription( model, print_info );

gnsf = detect_affine_terms_reduce_nonlinearity( gnsf, model, print_info );

if detect_LOS
    [ gnsf, reordered_model] = reformulate_with_LOS( model, gnsf, print_info);
else
    reordered_model = model;
end

if check_E_invertibility
    gnsf = reformulate_with_invertible_E_mat( gnsf, reordered_model, print_info);
end

structure_detection_print_summary(gnsf, model, reordered_model);
check_reformulation( reordered_model, gnsf, print_info );


%% EXPORT C Code
if generate_gnsf_model
    % generate gnsf model
    generate_c_code_gnsf( gnsf );
    disp('Successfully generated C Code to simulate model with acados integrator GNSF');
end

if generate_reordered_model
    % generate implicit model
    generate_c_code_implicit_ode( reordered_model, transcribe_opts );
    disp('Successfully generated C Code to simulate model with acados integrator IRK');
end


end

