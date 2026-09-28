%
% Copyright (c) The acados authors.
%
% This file is part of acados.
%
% Licensed under the 2-Clause BSD License.

function nondefault_fields = find_non_default_fields_of_obj(obj, stage_type)
    if nargin < 2
        stage_type = 'all';
    end

    % Get all public fields of the object
    try
        obj_struct = obj.to_struct();
    catch
        error('The object does not have a to_struct() method. Make sure it is an acados object, got %s.', class(obj));
    end
    all_fields = fieldnames(obj.to_struct());

    % Remove special properties that need translation
    if isa(obj, 'AcadosOcpConstraints')
        all_fields(strcmp(all_fields, 'x0')) = [];
        % all_fields = all_fields(~strncmp(all_fields, 'J', 1));  % Use strncmp for Octave compatibility
    end

    if isa(obj, 'AcadosOcpOptions')
        all_fields(strcmp(all_fields, 'qp_tol')) = [];
        all_fields(strcmp(all_fields, 'tol')) = [];
    end

    % Filter based on stage_type
    switch stage_type
        case 'all'
            % Do nothing
        case 'initial'
            all_fields = all_fields(endsWith_custom(all_fields, '_0'));
        case 'terminal'
            all_fields = all_fields(endsWith_custom(all_fields, '_e'));
        otherwise
            error('stage_type %s not supported.', stage_type);
    end

    % Compare fields with default values
    obj_type = class(obj);
    dummy_obj = feval(obj_type);
    dummy_obj_struct = dummy_obj.to_struct();
    nondefault_fields = {};

    for i = 1:length(all_fields)
        field = all_fields{i};

        % Check if the object has the field (to avoid errors in Octave)
        if isfield(obj_struct, field) && isfield(dummy_obj_struct, field)
            val = obj_struct.(field);
            default_val = dummy_obj_struct.(field);

            if ~isequal(val, default_val)
                nondefault_fields{end+1} = field;
            end
        end
    end
end
