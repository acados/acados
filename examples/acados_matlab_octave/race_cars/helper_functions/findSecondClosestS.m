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

function idxmindist2 = findSecondClosestS(si,sref,idxmindist)
    sref2 = [sref(end); sref; sref(1)]; % loop/wrap    
    
    d1=abs(si-sref2(idxmindist-1+1));               % distance to node before
    d2=abs(si-sref2(idxmindist+1+1));               % distance to node after

    idx = zeros(length(si), 1);
    idx(d1 > d2) = idxmindist(d1 > d2)+1;
    idx(d1 <= d2) = idxmindist(d1 <= d2)-1;
    
    idx(idx > length(sref)) = idx(idx > length(sref)) - length(sref);
    idx(idx <= 0) = idx(idx <= 0) + length(sref);   
    
    idxmindist2 = idx;
end