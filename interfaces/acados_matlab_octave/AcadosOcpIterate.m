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


classdef AcadosOcpIterate < handle
    properties
        x
        u
        z
        sl
        su
        pi
        lam
    end

    properties (Dependent)
        x_traj
        u_traj
        z_traj
        sl_traj
        su_traj
        pi_traj
        lam_traj
    end

    methods
        function obj = AcadosOcpIterate(x_, u_, z_, sl_, su_, pi_, lam_)
            obj.x = x_;
            obj.u = u_;
            obj.z = z_;
            obj.sl = sl_;
            obj.su = su_;
            obj.pi = pi_;
            obj.lam = lam_;
        end

        function s = to_struct(self)
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

        function val = get.x_traj(obj)
            obj.warn_traj_deprecated('x_traj');
            val = obj.x;
        end

        function val = get.u_traj(obj)
            obj.warn_traj_deprecated('u_traj');
            val = obj.u;
        end

        function val = get.z_traj(obj)
            obj.warn_traj_deprecated('z_traj');
            val = obj.z;
        end

        function val = get.sl_traj(obj)
            obj.warn_traj_deprecated('sl_traj');
            val = obj.sl;
        end

        function val = get.su_traj(obj)
            obj.warn_traj_deprecated('su_traj');
            val = obj.su;
        end

        function val = get.pi_traj(obj)
            obj.warn_traj_deprecated('pi_traj');
            val = obj.pi;
        end

        function val = get.lam_traj(obj)
            obj.warn_traj_deprecated('lam_traj');
            val = obj.lam;
        end

        function set.x_traj(obj, val)
            obj.warn_traj_deprecated('x_traj');
            obj.x = val;
        end

        function set.u_traj(obj, val)
            obj.warn_traj_deprecated('u_traj');
            obj.u = val;
        end

        function set.z_traj(obj, val)
            obj.warn_traj_deprecated('z_traj');
            obj.z = val;
        end

        function set.sl_traj(obj, val)
            obj.warn_traj_deprecated('sl_traj');
            obj.sl = val;
        end

        function set.su_traj(obj, val)
            obj.warn_traj_deprecated('su_traj');
            obj.su = val;
        end

        function set.pi_traj(obj, val)
            obj.warn_traj_deprecated('pi_traj');
            obj.pi = val;
        end

        function set.lam_traj(obj, val)
            obj.warn_traj_deprecated('lam_traj');
            obj.lam = val;
        end
      end
        methods (Access = private)
            function warn_traj_deprecated(~, field_name)
                warning(['The use of the field "' field_name '" is deprecated. Please use the field without "_traj" suffix instead.']);
            end
        end
end

