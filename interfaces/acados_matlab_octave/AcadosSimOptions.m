%
% Copyright (c) The acados authors.
%
% This file is part of acados.
%
% Licensed under the 2-Clause BSD License.



classdef AcadosSimOptions < handle
    properties
        integrator_type
        collocation_type
        Tsim
        num_stages
        num_steps
        newton_iter
        newton_tol
        jac_reuse
        sens_forw
        sens_adj
        sens_algebraic
        sens_hess
        output_z
        ext_fun_compile_flags
        ext_fun_expand_dyn
        compile_interface
        with_batch_functionality
        sens_forw_p
    end

    methods
        function obj = AcadosSimOptions()
            obj.integrator_type = 'ERK';
            obj.collocation_type = 'GAUSS_LEGENDRE';
            obj.Tsim = [];
            obj.num_stages = 4;
            obj.num_steps = 1;
            obj.newton_iter = 3;
            obj.newton_tol = 0.;
            obj.sens_forw = true;
            obj.sens_adj = false;
            obj.sens_algebraic = false;
            obj.sens_hess = false;
            obj.output_z = true;
            obj.jac_reuse = 0;
            obj.compile_interface = []; % corresponds to automatic detection, possible values: true, false, []
            obj.with_batch_functionality = false;

            % TODO the options below are deprecated and will be removed
            % check whether flags are provided by environment variable
            env_var = getenv("ACADOS_EXT_FUN_COMPILE_FLAGS");
            if isempty(env_var)
                obj.ext_fun_compile_flags = '-O2';
            else
                obj.ext_fun_compile_flags = env_var;
            end
            obj.ext_fun_expand_dyn = false;
            obj.sens_forw_p = false;     % default: disabled
        end

        function s = to_struct(self)
            if exist('properties')
                publicProperties = eval('properties(self)');
            else
                publicProperties = fieldnames(self);
            end
            s = struct();
            for fi = 1:numel(publicProperties)
                property_name = publicProperties{fi};
                if strcmp(property_name, 'num_stages') || strcmp(property_name, 'num_steps') || strcmp(property_name, 'newton_iter') || ...
                     strcmp(property_name, 'jac_reuse') || strcmp(property_name, 'newton_tol')
                    out_name = strcat('sim_method_', property_name);
                    s.(out_name) = self.(property_name);
                else
                    s.(property_name) = self.(property_name);
                end
            end
        end
    end
    methods (Static)
        function obj = from_struct(s)
            % Create AcadosSimOptions from a struct (e.g. decoded from JSON).
            if ~isstruct(s)
                error('from_struct input must be a struct.');
            end

            obj = AcadosSimOptions();
            sim_method_fields = {'num_stages', 'num_steps', 'newton_iter', 'jac_reuse', 'newton_tol'};
            fields = fieldnames(s);
            for i = 1:length(fields)
                f = fields{i};
                property_name = f;
                for j = 1:length(sim_method_fields)
                    if strcmp(f, ['sim_method_' sim_method_fields{j}])
                        property_name = sim_method_fields{j};
                        break
                    end
                end

                try
                    obj.(property_name) = s.(f);
                catch
                    % ignore unknown fields
                    warning(['Could not assign field ' f ' in AcadosSimOptions.from_struct']);
                end
            end
        end
    end
end
