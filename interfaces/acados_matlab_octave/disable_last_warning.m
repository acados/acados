%
% Copyright (c) The acados authors.
%
% This file is part of acados.
%
% Licensed under the 2-Clause BSD License.

%

function disable_last_warning()
    % to print warning only once:
    w = warning('query','last');
    if ~isempty(w)
        id = w.identifier;
        warning('off',id);
    end
end