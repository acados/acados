%
% Copyright (c) The acados authors.
%
% This file is part of acados.
%
% Licensed under the 2-Clause BSD License.


function check_dir_and_create(dir)
    if ~exist(dir, 'dir')
        mkdir(dir);
    end
end
