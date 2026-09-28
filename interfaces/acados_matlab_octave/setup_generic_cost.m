%
% Copyright (c) The acados authors.
%
% This file is part of acados.
%
% Licensed under the 2-Clause BSD License.

function setup_generic_cost(context, cost, target_dir, stage_type)
    if strcmp(stage_type, 'initial')
        cost_source_ext_cost = cost.cost_source_ext_cost_0;
    elseif strcmp(stage_type, 'path')
        cost_source_ext_cost = cost.cost_source_ext_cost;
    elseif strcmp(stage_type, 'terminal')
        cost_source_ext_cost = cost.cost_source_ext_cost_e;
    else
        error('Unknown stage_type.')
    end

    check_dir_and_create(target_dir);
    copyfile(fullfile(pwd, cost_source_ext_cost), target_dir);
    context.add_external_function_file(cost_source_ext_cost, target_dir);

end
