classdef acados_sim_json < handle

    properties
        % file structure
        acados_include_path
        acados_lib_path
        shared_lib_ext
        json_file
        code_export_directory
        % struct / object
        dims
        model
        sim_options
        % plain data
        parameter_values
        problem_class
    end

    methods
        function obj = acados_sim_json()
            % most fields are initialized as a placeholder
            obj.acados_include_path = [];
            obj.acados_lib_path = [];
            obj.shared_lib_ext = [];
            obj.json_file = [];
            obj.code_export_directory = [];

            obj.dims = struct(
                'nx', [], ...
                'nu', [], ...
                'nz', [], ...
                'np', [] ...
            );
            obj.model = acados_template_mex.acados_model_json();
            obj.sim_options = struct(
                % string
                'integrator_type', [], ...
                'collocation_type', [], ...
                % int
                'sim_method_num_stages', [], ...
                'sim_method_num_steps', [], ...
                'sim_method_newton_iter', [], ...
                % double
                'sim_method_newton_tol', [], ...
                'Tsim', [], ...
                % bool
                'sens_forw', [], ...
                'sens_adj', [], ...
                'sens_algebraic', [], ...
                'sens_hess', [], ...
                'output_z', [], ...
                % extra
                'sim_method_jac_reuse', [], ...
                'ext_fun_compile_flags', [] ...
            );

            obj.parameter_values = [];
            obj.problem_class = 'SIM';
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

end % class
