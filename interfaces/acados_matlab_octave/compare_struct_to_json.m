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

function mismatched_fields = compare_struct_to_json(ocp_struct, json_struct, tol)
    % Compare every entry of an OCP struct to a JSON struct, ignoring certain fields.
    %
    % Args:
    %   ocp_struct: OCP struct (from ocp.to_struct())
    %   json_struct: JSON struct (from loadjson(json_file))
    %   tol: tolerance for floating point comparisons
    % Returns:
    %   mismatched_fields: cell array of field paths that do not match

    % remove keys that should not affect the comparison (same as hash_struct)
    ignored_fields = {'external_function_files_model', 'external_function_files_ocp', 'json_loaded'};
    for i = 1:length(ignored_fields)
        if isfield(ocp_struct, ignored_fields{i})
            ocp_struct = rmfield(ocp_struct, ignored_fields{i});
        end
    end
    % n_global_data is only set during code generation, not in make_consistent
    if isfield(ocp_struct, 'dims') && isfield(ocp_struct.dims, 'n_global_data')
        ocp_struct.dims = rmfield(ocp_struct.dims, 'n_global_data');
        json_struct.dims = rmfield(json_struct.dims, 'n_global_data');
    elseif isfield(ocp_struct, 'phases_dims') && iscell(ocp_struct.phases_dims)
        % for multiphase OCP, phases_dims is a cell of dims, remove in all of them
        for i = 1:length(ocp_struct.phases_dims)
            if isfield(ocp_struct.phases_dims{i}, 'n_global_data')
                ocp_struct.phases_dims{i} = rmfield(ocp_struct.phases_dims{i}, 'n_global_data');
                json_struct.phases_dims{i} = rmfield(json_struct.phases_dims{i}, 'n_global_data');
            end
        end
    end

    mismatched_fields = compare_recursive(ocp_struct, json_struct, '', tol);
end

function mismatched = compare_recursive(ocp_data, json_data, path, tol)
    mismatched = {};
    if isstruct(ocp_data) && isstruct(json_data)
        fields = fieldnames(ocp_data);
        for i = 1:length(fields)
            key = fields{i};
            if isempty(path)
                current_path = key;
            else
                current_path = [path '.' key];
            end
            if ~isfield(json_data, key)
                mismatched{end+1} = current_path;
            else
                sub_mismatched = compare_recursive(ocp_data.(key), json_data.(key), current_path, tol);
                mismatched = [mismatched, sub_mismatched];
            end
        end
    elseif isobject(ocp_data) && isobject(json_data)
        % Handle nested classes
        props = properties(ocp_data);
        props(strcmp(props, 'x0')) = [];
        for i = 1:length(props)
            key = props{i};
            if isempty(path)
                current_path = key;
            else
                current_path = [path '.' key];
            end
            sub_mismatched = compare_recursive(ocp_data.(key), json_data.(key), current_path, tol);
            mismatched = [mismatched, sub_mismatched];
        end
    elseif isempty(ocp_data) && isempty(json_data)
        % treat empty arrays as equal regardless of shape
        return

    elseif iscell(ocp_data) && iscell(json_data)
        if length(ocp_data) ~= length(json_data)
            mismatched{end+1} = path;
        else
            for i = 1:length(ocp_data)
                current_path = sprintf('%s[%d]', path, i);
                sub_mismatched = compare_recursive(ocp_data{i}, json_data{i}, current_path, tol);
                mismatched = [mismatched, sub_mismatched];
            end
        end
    elseif (ischar(ocp_data) || isstring(ocp_data)) && (ischar(json_data) || isstring(json_data))
        if ~strcmp(ocp_data, json_data)
            mismatched{end+1} = path;
        end
    elseif ismatrix(ocp_data) && ismatrix(json_data)
        if ~isequal(size(ocp_data), size(json_data))
            mismatched{end+1} = path;
        elseif ~all(abs(ocp_data - json_data) < tol, 'all')
            % check relative tolerance
            rel_diff = abs(ocp_data - json_data) ./ max(abs(json_data), abs(ocp_data));
            if ~all(rel_diff < tol, 'all')
                mismatched{end+1} = path;
            end
        end
    else
        try
            if ~isequal(ocp_data, json_data)
                mismatched{end+1} = path;
                disp('found mismatch:')
                disp(ocp_data);
                disp(json_data);
            end
        catch
            mismatched{end+1} = path;
        end
    end
end
