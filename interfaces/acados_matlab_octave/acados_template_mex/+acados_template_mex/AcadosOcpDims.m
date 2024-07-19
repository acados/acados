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


classdef AcadosOcpDims < handle
    properties
        N      % prediction horizon

        % model
        nx     % number of states
        nu     % number of inputs
        nz     % number of algebraic variables
        np     % number of parameters

        % cost
        ny     % number of residuals in Lagrange term
        ny_0
        ny_e   % number of residuals in Mayer term

        % bounds
        nbu    % number of input bounds
        nbx    % number of state bounds
        nbx_0   % number of state bounds on x0
        nbx_e  % number of state bounds at t=T

        % nonlinear constraints
        nh     % number of nonlinear constraints
        nh_0   % number of nonlinear constraints at t=0
        nh_e   % number of nonlinear constraints at t=T
        nr %
        nr_e %
        nr_0
        nphi
        nphi_e
        nphi_0

        % general linear constraints
        ng     % number of general linear constraints
        ng_e   % number of general linear constraints at t=T

        % soft constraints
        nsbx   % number of soft state bounds
        nsbx_e % number of state bounds at t=T
        nsh    % number of soft bounds on nonlinear constraints
        nsh_0
        nsh_e  % number of soft bounds on nonlinear constraints at t=T
        nsphi
        nsphi_0
        nsphi_e
        nsbu   % number of soft state bounds
        ns     % total number of soft bounds
        ns_0   %
        ns_e   % total number of soft bounds at t=T
        nsg     % number of soft general linear constraints
        nsg_e   % number of soft general linear constraints at t=T
        % equalities within x bounds
        nbxe_0

        % gnsf
        % TODO these dimensions are not part of the corresponding python class (?)
        gnsf_nx1
        gnsf_nz1
        gnsf_nout
        gnsf_ny
        gnsf_nuhat
    end

    methods
        function obj = AcadosOcpDims()
            obj.N     = [];

            obj.nx    = [];
            obj.nu    = [];
            obj.nz    = 0;
            obj.np    = 0;

            obj.ny    = [];
            obj.ny_0 = [];
            obj.ny_e   = [];

            obj.nbu   = 0;
            obj.nbx   = 0;
            obj.nbx_0  = 0;
            obj.nbx_e  = 0;

            obj.nh    = 0;
            obj.nh_0   = 0;
            obj.nh_e   = 0;

            obj.nr = 0;
            obj.nr_0 = 0;
            obj.nr_e = 0;
            obj.nphi = 0;
            obj.nphi_0 = 0;
            obj.nphi_e = 0;
            obj.ng    = 0;
            obj.ng_e   = 0;

            obj.nsbx  = 0;
            obj.nsbx_e = 0;
            obj.nsbu  = 0;
            obj.nsh   = 0;
            obj.nsh_0  = 0;
            obj.nsh_e  = 0;
            obj.nsphi = 0;
            obj.nsphi_0 = 0;
            obj.nsphi_e = 0;
            obj.ns    = 0;
            obj.ns_0 = 0;
            obj.ns_e   = 0;
            obj.nsg   = 0;
            obj.nsg_e  = 0;

            obj.nbxe_0 = 0;

            obj.gnsf_nx1 = 0;
            obj.gnsf_nz1 = 0;
            obj.gnsf_nout = 0;
            obj.gnsf_ny = 0;
            obj.gnsf_nuhat = 0;
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

