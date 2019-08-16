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

classdef ocp_nlp_constraints < handle
    properties
        % bounds on x and u
        lbx     % lower bounds on x
        lbu     % lower bounds on u
        idxbx   % indexes of bounds on x 
        ubx     % upper bounds on x 
        ubu     % upper bounds on u 
        idxbu   % indexes of bounds on u
        % bounds on x at t=T
        lbx_e    % lower bounds on x at t=T 
        ubx_e    % upper bounds on x at t=T 
        idxbx_e  % indexes for bounds on x at t=T 
        % soft bounds on x and u
        lsbx    % soft lower bounds on x
        lsbu    % soft lower bounds on u
        usbx    % soft upper bounds on x 
        usbu    % soft upper bounds on u 
        idxsbx  % indexes of soft bounds on x 
        idxsbu  % indexes of soft bounds on u
        % soft bounds on x and u at t=T
        lsbx_e   % soft lower bounds on x at t=T
        usbx_e   % soft upper bounds on x at t=T
        idxsbx_e % indexes of soft bounds on x at t=T 
        % polytopic constraints 
        D       % D matrix in lg <= D * u + C * x <= ug
        C       % C matrix in lg <= D * u + C * x <= ug
        lg      % lower bound for general inequalities 
        ug      % upper bound for general inequalities 
        % polytopic constraints at t=T 
        C_e      % C matrix at t=T 
        lg_e     % lower bound on general inequalities at t=T 
        ug_e     % upper bound on general inequalities at t=T 
        % nonlinear constraints
        lh      % lower bound for nonlinear inequalities 
        uh      % upper bound for nonlinear inequalities 
        % nonlinear constraints at t=T
        lh_e     % lower bound on nonlinear inequalities at t=T 
        uh_e     % upper bound on nonlinear inequalities at t=T 
        x0      % initial state 
    end
    methods
        function obj = ocp_nlp_constraints()
            obj.lbx     = [];  
            obj.lbu     = [];
            obj.idxbx   = [];
            obj.ubx     = [];
            obj.ubu     = [];
            obj.idxbu   = [];
            obj.lsbx    = [];  
            obj.lsbu    = [];
            obj.idxsbx  = [];
            obj.usbx    = [];
            obj.usbu    = [];
            obj.idxsbu  = [];
            obj.lsbx_e   = [];  
            obj.idxsbx_e = [];
            obj.lg      = [];
            obj.ug      = [];
            obj.lh      = [];
            obj.uh      = [];
            obj.D       = [];
            obj.C       = [];
            obj.lbx_e    = [];
            obj.ubx_e    = [];
            obj.idxbx_e  = [];
            obj.C_e      = [];
            obj.lg_e     = [];
            obj.ug_e     = [];
            obj.lh_e     = [];
            obj.uh_e     = [];
            obj.x0      = [];
        end
    end
end

