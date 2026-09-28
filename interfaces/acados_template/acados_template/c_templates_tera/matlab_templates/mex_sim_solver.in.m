%
% Copyright (c) The acados authors.
%
% This file is part of acados.
%
% Licensed under the 2-Clause BSD License.

%

classdef {{ name }}_mex_sim_solver < handle

    properties
        C_sim
        name
        code_gen_dir
    end % properties



    methods

        % constructor
        function obj = {{ name }}_mex_sim_solver()
            make_mex_sim_{{ name }}();
            obj.C_sim = acados_sim_create_{{ name }}();
            % to have path to destructor when changing directory
            addpath('.')
            obj.name = '{{ name }}';
        end

        % set -- borrowed from MEX interface
        function set(obj, field, value)
            if ~isa(field, 'char')
                error('field must be a char vector, use '' ''');
            end
            acados_sim_set_{{ name }}(obj.C_sim, field, value);
        end

        % get -- borrowed from MEX interface
        function value = get(obj, field)
            if ~isa(field, 'char')
                error('field must be a char vector, use '' ''');
            end
            value = sim_get(obj.C_sim, field);
        end

        % solve
        function status = solve(obj)
            status = sim_solve(obj.C_sim);
        end

        % destructor
        function delete(obj)
            disp("delete template...");
            if ~isempty(obj.C_sim)
                acados_sim_free_{{ name }}(obj.C_sim);
            end
            disp("done.");
        end

    end % methods

end % class
