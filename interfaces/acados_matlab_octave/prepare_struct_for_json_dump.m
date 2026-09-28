%
% Copyright (c) The acados authors.
%
% This file is part of acados.
%
% Licensed under the 2-Clause BSD License.

function out = prepare_struct_for_json_dump(out, vector_properties, matrix_properties)
    for i = 1:length(vector_properties)
        prop = vector_properties{i};
        if ~isempty(out.(prop))
            if ~iscell(out.(prop))
                out.(prop) = num2cell(out.(prop));
            end
            out.(prop) = reshape(out.(prop), [1, length(out.(prop))]);
        else
            out.(prop) = [];
        end
    end

    for i = 1:length(matrix_properties)
        prop = matrix_properties{i};
        if ~isempty(out.(prop))
            out.(prop) = num2cell(out.(prop));
            prop_size = size(out.(prop));
            if prop_size(1) == 1
                out.(prop) = {out.(prop)};
            end
        else
            out.(prop) = [];
        end
    end
end
