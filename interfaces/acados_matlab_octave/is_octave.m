%
% Copyright (c) The acados authors.
%
% This file is part of acados.
%
% Licensed under the 2-Clause BSD License.

%

% function to check if running in octave or matlab
function r = is_octave()
    persistent x;
    if (isempty(x))
        x = exist( 'OCTAVE_VERSION', 'builtin');
    end
    r = x;
end

