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

classdef AcadosOcpConstraints < handle
    properties
        constr_type_0
        constr_type
        constr_type_e

        % bounds on x
        lbx_0   % lower bounds on x
        ubx_0   % lower bounds on u
        idxbx_0
        idxbxe_0
        has_x0

        lbx     % lower bounds on x
        ubx     % upper bounds on x
        idxbx   % indexes of bounds on x

        lbx_e    % lower bounds on x at t=T
        ubx_e    % upper bounds on x at t=T
        idxbx_e  % indexes for bounds on x at t=T

        % bounds on u
        lbu     % lower bounds on u
        ubu     % upper bounds on u
        idxbu   % indexes of bounds on u

        % polytopic constraints
        lg      % lower bound for general inequalities
        ug      % upper bound for general inequalities
        D       % D matrix in lg <= D * u + C * x <= ug
        C       % C matrix in lg <= D * u + C * x <= ug

        % polytopic constraints at t=T
        lg_e     % lower bound on general inequalities at t=T
        ug_e     % upper bound on general inequalities at t=T
        C_e      % C matrix at t=T

        % nonlinear constraints at t=0
        lh_0     % lower bound on nonlinear inequalities at t=0
        uh_0     % upper bound on nonlinear inequalities at t=0
        % nonlinear constraints
        lh      % lower bound for nonlinear inequalities
        uh      % upper bound for nonlinear inequalities
        % nonlinear constraints at t=T
        lh_e     % lower bound on nonlinear inequalities at t=T
        uh_e     % upper bound on nonlinear inequalities at t=T

        % convex-over-nonlinear constraint (BGP) to work with json
        % TODO: implement in MEX..
        lphi_0
        uphi_0
        lphi
        uphi
        lphi_e
        uphi_e

        %  SLACKS
        % bounds on slacks corresponding to softened bounds on x and u
        lsbx    % lower bounds on slacks corresponding to soft lower bounds on x
        usbx    % lower bounds on slacks corresponding to soft upper bounds on x
        idxsbx  % indexes of soft bounds on x

        lsbu    % lower bounds on slacks corresponding to soft lower bounds on u
        usbu    % lower bounds on slacks corresponding to soft upper bounds on u
        idxsbu  % indexes of soft bounds on u

        % bounds on slacks corresponding to softened bounds on t=T
        lsbx_e   % lower bounds on slacks corresponding to soft lower bounds on x at t=T
        usbx_e   % upper bounds on slacks corresponding to soft upper bounds on x at t=T
        idxsbx_e % indexes of soft bounds on x at t=T

        % soft bounds on general linear constraints
        lsg     % lower bounds on slacks corresponding to soft lower bounds for general linear constraints
        usg     % lower bounds on slacks corresponding to soft upper bounds for general linear constraints
        idxsg   % indexes of soft general linear constraints

        % soft bounds on general linear constraints at t=T
        lsg_e     % lower bounds on slacks corresponding to soft lower bounds for general linear constraints
        usg_e     % lower bounds on slacks corresponding to soft upper bounds for general linear constraints
        idxsg_e   % indexes of soft general linear constraints

        % soft bounds on nonlinear constraints at t=0
        lsh_0
        ush_0
        idxsh_0
        % soft bounds on nonlinear constraints
        lsh     % lower bounds on slacks corresponding to soft lower bounds for nonlinear constraints
        ush     % lower bounds on slacks corresponding to soft upper bounds for nonlinear constraints
        idxsh   % indexes of soft nonlinear constraints
        % soft bounds on nonlinear constraints at t=T
        lsh_e
        ush_e
        idxsh_e

        % soft bounds on convex-over-nonlinear constraint (BGP)
        lsphi_0
        usphi_0
        idxsphi_0

        lsphi
        usphi
        idxsphi

        lsphi_e
        usphi_e
        idxsphi_e

    end
    methods
        function obj = AcadosOcpConstraints()
            obj.constr_type_0 = 'BGH';
            obj.constr_type = 'BGH';
            obj.constr_type_e = 'BGH';

            obj.lbx_0 = [];
            obj.ubx_0 = [];
            obj.idxbx_0 = [];
            obj.idxbxe_0 = [];
            obj.has_x0 = false;

            obj.lbx = [];
            obj.ubx = [];
            obj.idxbx = [];

            obj.lbx_e = [];
            obj.ubx_e = [];
            obj.idxbx_e = [];

            obj.lbu = [];
            obj.ubu = [];
            obj.idxbu = [];

            obj.lg = [];
            obj.ug = [];
            obj.D = [];
            obj.C = [];

            obj.lg_e = [];
            obj.ug_e = [];
            obj.C_e = [];

            obj.lh_0 = [];
            obj.uh_0 = [];
            obj.lh = [];
            obj.uh = [];
            obj.lh_e = [];
            obj.uh_e = [];

            obj.lphi_0 = [];
            obj.uphi_0 = [];
            obj.lphi = [];
            obj.uphi = [];
            obj.lphi_e = [];
            obj.uphi_e = [];

            % SLACKS
            obj.lsbx = [];
            obj.usbx = [];
            obj.idxsbx = [];
            obj.lsbu = [];
            obj.usbu = [];
            obj.idxsbu = [];
            obj.lsbx_e = [];
            obj.usbx_e = [];
            obj.idxsbx_e = [];

            obj.lsg = [];
            obj.usg = [];
            obj.idxsg = [];

            obj.lsg_e = [];
            obj.usg_e = [];
            obj.idxsg_e = [];

            obj.lsh_0 = [];
            obj.ush_0 = [];
            obj.idxsh_0 = [];

            obj.lsh = [];
            obj.ush = [];
            obj.idxsh = [];

            obj.lsh_e = [];
            obj.ush_e = [];
            obj.idxsh_e = [];

            obj.lsphi_0 = [];
            obj.usphi_0 = [];
            obj.idxsphi_0 = [];

            obj.lsphi = [];
            obj.usphi = [];
            obj.idxsphi = [];

            obj.lsphi_e = [];
            obj.usphi_e = [];
            obj.idxsphi_e = [];

        end

        function s = struct(self)
            if exist('properties')
                publicProperties = eval('properties(self)');
            else
                publicProperties = fieldnames(self);
            end
            s = struct();
            for fi = 1:numel(publicProperties)
                s.(publicProperties{fi}) = self.(publicProperties{fi});
            end
        end
    end
end

