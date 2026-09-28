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

function idx = findClosestNeighbour(x,y,xref,yref,idxmindist)
    distBefore=dist2D(x,xref(idxmindist-1),y,yref(idxmindist-1));
    distAfter=dist2D(x,xref(idxmindist+1),y,yref(idxmindist+1));
    if (distBefore<distAfter)
        idxmindist2=idxmindist-1;
    else
        idxmindist2=idxmindist+1;
    end
    if(idxmindist2<0)
        idxmindist2=xref.size-1;
    elseif(idxmindist==length(xref))
        idxmindist2=0;
    end
    
    idx = idxmindist2;
end