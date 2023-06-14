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


clc;
clear VARIABLES;
close all;

addpath('../../../interfaces/acados_matlab_octave/')
addpath('../../matlab_mex/pendulum_dae/')

% define model 
model = pendulum_dae_model();
model.dyn_expr_f = model.expr_f_impl;

%% GNSF Model -- detect structure, reorder model, and generate C Code for
%% GNSF model. --> for more advanded users - uncomment this section
% Reformulate model as GNSF & Reorder x, xdot, z, f_impl, f_expl
% accordingly
transcribe_opts.print_info = 1;
[ gnsf ] = detect_gnsf_structure(model, transcribe_opts);
    % check output of this function to see if/how the states are reordered
generate_c_code_gnsf( gnsf );

%% Implicit Model -- Generate C Code
opts.generate_hess = 1;  % set to 1 if you want to use exact hessian propagation

generate_c_code_implicit_ode( model, opts );

