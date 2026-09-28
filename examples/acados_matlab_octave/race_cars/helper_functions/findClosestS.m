%
% Copyright (c) The acados authors.
%
% This file is part of acados.
%
% Licensed under the 2-Clause BSD License.

%
% Author: Daniel Kloeser
% Ported by Thomas Jespersen (thomasj@tkjelectronics.dk), TKJ Electronics
%

function idx = findClosestS(si,sref)
    idx = zeros(length(si), 1);
    
    for j = 1:length(si)
        [~,i] = min(abs(sref - si(j)));
        if (i == length(sref))
            i = 1;
        elseif (i == 1)
            i = length(sref);
        end
        idx(j) = i;
    end
end