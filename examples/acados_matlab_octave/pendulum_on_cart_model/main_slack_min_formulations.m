% Copyright (c) The acados authors.
%
% This file is part of acados.
%
% Licensed under the 2-Clause BSD License.

formulations = {'u_slack', 'u_slack2', 's_slack'};
xtraj_list = {};
for i = 1:length(formulations)
    xtraj = slack_min_formulation(formulations{i});
    if i == 1
        xtraj_ref = xtraj;
    else
        diff_x = max(abs(xtraj(:) - xtraj_ref(:)));
        fprintf('diff xtraj %s vs ref: %g\n', formulations{i}, diff_x);
        if diff_x > 1e-6
            error('xtraj %s differs from reference by %g, expected to be close to zero.', formulations{i}, diff_x);
        end
    end
end
