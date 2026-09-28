%
% Copyright (c) The acados authors.
%
% This file is part of acados.
%
% Licensed under the 2-Clause BSD License.

%

classdef AcadosOcpSimulinkOutputs < handle

    properties
        u0
        utraj
        xtraj
        ztraj
        pi_all
        slack_values
        solver_status
        cost_value
        KKT_residual
        KKT_residuals
        x1
        CPU_time
        CPU_time_sim
        CPU_time_qp
        CPU_time_lin
        sqp_iter
        parameter_traj
        zoRO_Pk_matrices
        zoRO_K_matrices
    end

    methods
        function obj = AcadosOcpSimulinkOutputs()
            obj.u0 = 1;
            obj.utraj = 0;
            obj.xtraj = 0;
            obj.ztraj = 0;
            obj.pi_all = 0;
            obj.slack_values = 0;
            obj.solver_status = 1;
            obj.cost_value = 0;
            obj.KKT_residual = 1;
            obj.KKT_residuals = 0;
            obj.x1 = 1;
            obj.CPU_time = 1;
            obj.CPU_time_sim = 0;
            obj.CPU_time_qp = 0;
            obj.CPU_time_lin = 0;
            obj.sqp_iter = 1;
            obj.parameter_traj = 0;
            obj.zoRO_Pk_matrices = 0;
            obj.zoRO_K_matrices = 0;
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
                    error(['AcadosOcpSimulinkOutputs.', names{i}, ' must be 0 or 1, got ', num2str(value)]);
                end
            end
        end
    end

    methods (Static)
        function obj = from_struct(data)
            obj = AcadosOcpSimulinkOutputs();
            fields = fieldnames(data);
            for i = 1:length(fields)
                f = fields{i};
                try
                    obj.(f) = data.(f);
                catch
                    warning(['Could not assign field ' f ' in AcadosOcpSimulinkOutputs.from_struct']);
                end
            end
        end
    end
end