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

classdef ocp_nlp_constraints_json < handle
    properties
        constr_type
        constr_type_e
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
        % bounds on slacks corresponding to softened bounds on x and u
        idxsbx  % indexes of soft bounds on x
        idxsbu  % indexes of soft bounds on u
        lsbx    % lower bounds on slacks corresponding to soft lower bounds on x
        lsbu    % lower bounds on slacks corresponding to soft lower bounds on u
        usbx    % lower bounds on slacks corresponding to soft upper bounds on x
        usbu    % lower bounds on slacks corresponding to soft upper bounds on u
        % bounds on slacks corresponding to softened bounds on t=T
        idxsbx_e % indexes of soft bounds on x at t=T
        lsbx_e   % lower bounds on slacks corresponding to soft lower bounds on x at t=T
        usbx_e   % lower bounds on slacks corresponding to soft upper bounds on x at t=T
        % soft bounds on general linear constraints
        lsg     % lower bounds on slacks corresponding to soft lower bounds for general linear constraints
        usg     % lower bounds on slacks corresponding to soft upper bounds for general linear constraints
        idxsg   % indexes of soft general linear constraints
        % soft bounds on general linear constraints at t=T
        lsg_e     % lower bounds on slacks corresponding to soft lower bounds for general linear constraints
        usg_e     % lower bounds on slacks corresponding to soft upper bounds for general linear constraints
        idxsg_e   % indexes of soft general linear constraints
        % soft bounds on nonlinear constraints
        lsh     % lower bounds on slacks corresponding to soft lower bounds for nonlinear constraints
        ush     % lower bounds on slacks corresponding to soft upper bounds for nonlinear constraints
        idxsh   % indexes of soft nonlinear constraints
        % soft bounds on nonlinear constraints at t=T
        lsh_e    % lower bounds on slacks corresponding to soft lower bounds for nonlinear constraints
        ush_e    % lower bounds on slacks corresponding to soft upper bounds for nonlinear constraints
        idxsh_e  % indexes of soft nonlinear constraints at t=T
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
        % initial state bounds
        lbx_0    % lower bound on initial state
        ubx_0    % upper bound on initial state
        idxbx_0  % bound indices of initial state
        idxbxe_0 %
        % convex over nonlinear constraint (BGP) to work with json
        % TODO: implement in MEX..
        lphi   % lower bound on convex over nonlinear constraint
        uphi   % upper bound on convex over nonlinear constraint
        lphi_e   % lower bound on convex over nonlinear constraint at t=T
        uphi_e   % upper bound on convex over nonlinear constraint at t=T
        lsphi   % lower bounds on slacks corresponding to lower bound on convex over nonlinear constraint
        usphi   % lower bounds on slacks corresponding to upper bound on convex over nonlinear constraint
        lsphi_e   % lower bounds on slacks corresponding to lower bound on convex over nonlinear constraint at t=T
        usphi_e   % lower bounds on slacks corresponding to upper bound on convex over nonlinear constraint at t=T
        idxsphi % indexes of soft convex over nonlinear constraints
        idxsphi_e % indexes of soft convex over nonlinear constraints at t=T

    end
    methods
        function obj = ocp_nlp_constraints_json()
            obj.constr_type     = 'BGH';
            obj.constr_type_e   = 'BGH';
            obj.lbx             = [];
            obj.lbu             = [];
            obj.idxbx           = [];
            obj.ubx             = [];
            obj.ubu             = [];
            obj.idxbu           = [];
            obj.lsbx            = [];
            obj.lsbu            = [];
            obj.idxsbx          = [];
            obj.usbx            = [];
            obj.usbu            = [];
            obj.idxsbu          = [];
            obj.lsbx_e          = [];
            obj.idxsbx_e        = [];
            obj.lg              = [];
            obj.ug              = [];
            obj.lh              = [];
            obj.uh              = [];
            obj.D               = [];
            obj.C               = [];
            obj.lbx_e           = [];
            obj.ubx_e           = [];
            obj.idxbx_e         = [];
            obj.C_e             = [];
            obj.lg_e            = [];
            obj.ug_e            = [];
            obj.lh_e            = [];
            obj.uh_e            = [];
            obj.lbx_0 = [];
            obj.ubx_0 = [];
            obj.idxbx_0 = [];
            obj.lphi            = [];
            obj.uphi            = [];
            obj.lsg = [];
            obj.usg = [];
            obj.idxsg = [];
            obj.lsg_e = [];
            obj.usg_e = [];
            obj.idxsg_e = [];
            obj.lsphi = [];
            obj.usphi = [];
            obj.idxsphi = [];
            obj.lphi_e            = [];
            obj.uphi_e            = [];
            obj.lsphi_e            = [];
            obj.usphi_e            = [];
            obj.idxsphi_e = [];
        end
        function s = struct(self)
            publicProperties = properties(self);
            s = struct();
            for fi = 1:numel(publicProperties)
                s.(publicProperties{fi}) = self.(publicProperties{fi});
            end
        end
    end
end

