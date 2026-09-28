%
% Copyright (c) The acados authors.
%
% This file is part of acados.
%
% Licensed under the 2-Clause BSD License.

function out = postprocess_struct_from_json_dump(out, vector_properties, matrix_properties)
    for i = 1:length(vector_properties)
        prop = vector_properties{i};
        if ~isfield(out, prop) || isempty(out.(prop))
            out.(prop) = [];
        else
            if iscell(out.(prop))
                out.(prop) = cell2mat(out.(prop));
            end
            out.(prop) = reshape(out.(prop), [length(out.(prop)), 1]);
        end
    end

    % Handle matrix properties that were transformed by
    % prepare_struct_for_json_dump before JSON serialization.
    for i = 1:length(matrix_properties)
        prop = matrix_properties{i};
        if ~isfield(out, prop) || isempty(out.(prop))
            out.(prop) = [];
        else
            % The prepare function converted matrices using num2cell().
            % Additionally, 1-row matrices were wrapped as a single cell
            % containing the row (i.e. { {a, b, c} }). Undo that here.
            if iscell(out.(prop))
                % unwrap nested single-cell case: { { ... } }
                if numel(out.(prop)) == 1 && iscell(out.(prop){1})
                    out.(prop) = out.(prop){1};
                end
                % convert cell-of-scalars back to numeric matrix
                try
                    out.(prop) = cell2mat(out.(prop));
                catch
                    % If conversion fails, leave as-is to avoid crashing;
                    % caller can handle or report a clearer error.
                end
            end
        end
    end
end
