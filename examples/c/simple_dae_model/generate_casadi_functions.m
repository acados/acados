%
% Copyright (c) The acados authors.
%
% This file is part of acados.
%
% Licensed under the 2-Clause BSD License.

%


clc;
clear all;
close all;

addpath('../../../experimental/interfaces/acados_matlab/') 

import casadi.*

% casadi opts for code generation
if CasadiMeta.version()=='3.4.0'
	% casadi 3.4
	opts = struct('mex', false, 'casadi_int', 'int', 'casadi_real', 'double');
else
	% old casadi versions
	error('Please download and install Casadi 3.4.0')
end

NX = 2;
NU = 2;
NZ = 2;

% define model 
dae     = export_simple_dae_model();
constr  = export_simple_dae_constr();

generate_c_code_implicit_ode( dae )
generate_c_code_nonlinear_constr( constr )

