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

function [s,n,alpha,v] = transformOrig2Proj(x,y,psi,v,filename)
    [sref,xref,yref,psiref,~]=getTrack(filename);
    idxmindist=findClosestPoint(x,y,xref,yref);
    idxmindist2=findClosestNeighbour(x,y,xref,yref,idxmindist);
    t=findProjection(x,y,xref,yref,sref,idxmindist,idxmindist2);
    s0=(1-t).*sref(idxmindist)+t.*sref(idxmindist2);
    x0=(1-t).*xref(idxmindist)+t.*xref(idxmindist2);
    y0=(1-t).*yref(idxmindist)+t.*yref(idxmindist2);
    psi0=(1-t).*psiref(idxmindist)+t.*psiref(idxmindist2);

    s=s0;
    n=cos(psi0).*(y-y0)-sin(psi0).*(x-x0);
    alpha=psi-psi0;
end