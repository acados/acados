%
% Copyright (c) The acados authors.
%
% This file is part of acados.
%
% Licensed under the 2-Clause BSD License.


classdef AcadosOcpIterates < handle
    properties (Access = public)
        iterates_cell
    end

    properties (Access = private)
        fields = {'x', 'u', 'z', 'sl', 'su', 'lam', 'pi'};
    end

    methods
        function obj = AcadosOcpIterates(iterates_cell_)
            obj.iterates_cell = iterates_cell_;
        end

        function iterate = as_array(obj, field)
            % Return the iterates as matrix of size (nlp_iter + 1, N_horizon (+ 1), n_field)
            % This will fail if the dimension of value `field` is varying stagewise.
            if ~any(strcmp(obj.fields, field))
                error(["Invalid field: got " field]);
            end

            n_iterates = length(obj.iterates_cell); % n_iterates = nlp_iter + 1
            field_iterates_cell = cell(n_iterates, 1);

            iterate = obj.iterates_cell{1};
            traj = iterate.(field);
            num_0 = length(traj);

            try
                % reshape to (num, n_field), num might be either N_horizon or N_horizon + 1
                for i=1:(n_iterates)
                    iterate = obj.iterates_cell{i};
                    traj = iterate.(field);

                    num = length(traj);
                    if num ~= num_0
                        error(['Stage-wise dimensions are not the same for ' field ' trajectory.']);
                    end
                    % NOTE: cannot change reshape order, thus need to transpose afterwards
                    field_iterates_cell{i} = reshape(cell2mat(traj), [], num).';
                end

                iterate = zeros(n_iterates, num_0, size(field_iterates_cell{1}, 2));
                for i=1:n_iterates
                    iterate(i, :, :) = field_iterates_cell{i};
                end
            catch
                error(['Stage-wise dimensions are not the same for ' field ' trajectory.']);
            end
        end

        function s = struct(obj)
            if exist('properties')
                publicProperties = eval('properties(obj)');
            else
                publicProperties = fieldnames(obj);
            end
            s = struct();
            for fi = 1:numel(publicProperties)
                s.(publicProperties{fi}) = obj.(publicProperties{fi});
            end
        end
    end
end