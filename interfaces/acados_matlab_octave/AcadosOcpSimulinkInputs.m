%
% Copyright (c) The acados authors.
%
% This file is part of acados.
%
% Licensed under the 2-Clause BSD License.

%

classdef AcadosOcpSimulinkInputs < handle

    properties
        lbx_0
        ubx_0
        parameter_traj
        p_global
        y_ref_0
        y_ref
        y_ref_e
        lbx
        ubx
        lbx_e
        ubx_e
        lbu
        ubu
        lg
        ug
        lh
        uh
        lh_0
        uh_0
        lh_e
        uh_e
        cost_W_0
        cost_W
        cost_W_e
        cost_zl
        cost_zu
        cost_Zl
        cost_Zu
        reset_solver
        reset_flags
        ignore_inits
        x_init
        u_init
        pi_init
        slacks_init
        rti_phase
        levenberg_marquardt
        zoRO_payload
    end

    methods
        function obj = AcadosOcpSimulinkInputs()
            obj.lbx_0 = 1;
            obj.ubx_0 = 1;
            obj.parameter_traj = 1;
            obj.p_global = 0;
            obj.y_ref_0 = 1;
            obj.y_ref = 1;
            obj.y_ref_e = 1;
            obj.lbx = 1;
            obj.ubx = 1;
            obj.lbx_e = 1;
            obj.ubx_e = 1;
            obj.lbu = 1;
            obj.ubu = 1;
            obj.lg = 1;
            obj.ug = 1;
            obj.lh = 1;
            obj.uh = 1;
            obj.lh_0 = 1;
            obj.uh_0 = 1;
            obj.lh_e = 1;
            obj.uh_e = 1;
            obj.cost_W_0 = 0;
            obj.cost_W = 0;
            obj.cost_W_e = 0;
            obj.cost_zl = 0;
            obj.cost_zu = 0;
            obj.cost_Zl = 0;
            obj.cost_Zu = 0;
            obj.reset_solver = 0;
            obj.reset_flags = 0;
            obj.ignore_inits = 0;
            obj.x_init = 0;
            obj.u_init = 0;
            obj.pi_init = 0;
            obj.slacks_init = 0;
            obj.rti_phase = 0;
            obj.levenberg_marquardt = 0;
            obj.zoRO_payload = 0;
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

            s = orderfields(s);
        end
        function make_consistent(self)
            s = self.to_struct();
            names = fieldnames(s);
            for i = 1:length(names)
                value = s.(names{i});
                if ~(isequal(value, 0) || isequal(value, 1))
                    error(['AcadosOcpSimulinkInputs.', names{i}, ' must be 0 or 1, got ', num2str(value)]);
                end
            end
        end
    end

    methods (Static)
        function obj = from_struct(data)
            obj = AcadosOcpSimulinkInputs();
            fields = fieldnames(data);
            for i = 1:length(fields)
                f = fields{i};
                try
                    obj.(f) = data.(f);
                catch
                    warning(['Could not assign field ' f ' in AcadosOcpSimulinkInputs.from_struct']);
                end
            end
        end
    end
end