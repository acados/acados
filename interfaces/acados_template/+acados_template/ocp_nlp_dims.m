%
% Copyright 2019 Gianluca Frison, Dimitris Kouzoupis, Robin Verschueren,
% Andrea Zanelli, Niels van Duijkeren, Jonathan Frey, Tommaso Sartor,
% Branimir Novoselnik, Rien Quirynen, Rezart Qelibari, Dang Doan,
% Jonas Koenemann, Yutao Chen, Tobias Schöls, Jonas Schlagenhauf, Moritz Diehl
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

classdef ocp_nlp_dims < handle
    properties
        nx     % number of states
        nz     % number of algebraic variables
        nu     % number of inputs
        np     % number of parameters
        ny     % number of residuals in Lagrange term
        ny_e   % number of residuals in Mayer term
        npd    % number of positive definite constraints
        npd_e  % number of positive definite constraints at t=T
        nh     % number of nonlinear constraints
        nh_e   % number of nonlinear constraints at t=T
        nbx    % number of state bounds 
        nbx_e  % number of state bounds at t=T 
        nbu    % number of input bounds
        nsbx   % number of soft state bounds 
        nsbu   % number of soft state bounds 
        nsbx_e % number of state bounds at t=T 
        ns     % total number of soft bounds 
        ns_e   % total number of soft bounds at t=T 
        ng     % number of general constraints
        ng_e   % number of general constraints at t=T
        N      % prediction horizon 
    end
    methods 
        function obj = ocp_nlp_dims()
            obj.nx    = []; 
            obj.nz    = 0; 
            obj.nu    = []; 
            obj.np    = 0;
            obj.ny    = []; 
            obj.ny_e   = []; 
            obj.npd   = 0;
            obj.npd_e  = 0; 
            obj.nh    = 0;
            obj.nh_e   = 0;
            obj.nbx   = 0;
            obj.nbu   = 0;
            obj.nbx_e  = 0;
            obj.nsbx  = 0;
            obj.nsbu  = 0;
            obj.nsbx_e = 0;
            obj.ns    = 0;
            obj.ns_e   = 0;
            obj.ng    = 0;
            obj.ng_e   = 0;
            obj.N     = [];
        end
    end
end

