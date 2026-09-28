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

function [x,y,psi,v] = transformProj2Orig(si,ni,alpha,v,filename)
    [sref,xref,yref,psiref,~]=getTrack(filename);
%     tracklength=sref(end);
    %si=tracklength;
    idxmindist=findClosestS(si,sref);
    idxmindist2=findSecondClosestS(si,sref,idxmindist);
    t=(si-sref(idxmindist))./(sref(idxmindist2)-sref(idxmindist));
    x0=(1-t).*xref(idxmindist)+t.*xref(idxmindist2);
    y0=(1-t).*yref(idxmindist)+t.*yref(idxmindist2);
    psi0=(1-t).*psiref(idxmindist)+t.*psiref(idxmindist2);

    x=x0-ni.*sin(psi0);
    y=y0+ni.*cos(psi0);
    psi=psi0+alpha;
end
