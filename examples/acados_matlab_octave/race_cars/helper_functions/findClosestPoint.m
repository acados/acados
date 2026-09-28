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

function idx = findClosestPoint(x,y,xref,yref)
    mindist=1;
    idxmindist=0;
    for i = 1:length(xref)
        dist=dist2D(x,xref(i),y,yref(i));
        if dist<mindist
            mindist=dist;
            idxmindist=i;
        end
    end
    
    idx = idxmindist;
end