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
