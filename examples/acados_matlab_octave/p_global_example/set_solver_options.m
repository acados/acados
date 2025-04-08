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



function ocp = set_solver_options(ocp)
    % set options
    ocp.solver_options.qp_solver = 'PARTIAL_CONDENSING_HPIPM';
    ocp.solver_options.hessian_approx = 'GAUSS_NEWTON'; %'GAUSS_NEWTON'; %'EXACT'; %
    ocp.solver_options.integrator_type = 'ERK';
    ocp.solver_options.print_level = 0;
    ocp.solver_options.nlp_solver_type = 'SQP_RTI';

    % set prediction horizon
    Tf = 1.0;
    N_horizon = 20;
    ocp.solver_options.tf = Tf;
    ocp.solver_options.N_horizon = N_horizon;

    % partial condensing
    ocp.solver_options.qp_solver_cond_N = 5;
    ocp.solver_options.qp_solver_cond_block_size = [3, 3, 3, 3, 7, 1];

    % NOTE: these additional flags are required for code generation of CasADi functions using casadi.blazing_spline
    % These might be different depending on your compiler and operating system.
    flags = ['-I' casadi.GlobalOptions.getCasadiIncludePath ' -O2 -ffast-math -march=native -fno-omit-frame-pointer']
    ocp.solver_options.ext_fun_compile_flags = flags;
    ocp.solver_options.ext_fun_expand_dyn = true;
    ocp.solver_options.ext_fun_expand_constr = true;
    ocp.solver_options.ext_fun_expand_cost = true;
    ocp.solver_options.ext_fun_expand_precompute = false;
end
