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

function t = findProjection(x,y,xref,yref,sref,idxmindist,idxmindist2)
    vabs=abs(sref(idxmindist)-sref(idxmindist2));
    vl=zeros(2,1);
    u=zeros(2,1);
    vl(1)=xref(idxmindist2)-xref(idxmindist);
    vl(2)=yref(idxmindist2)-yref(idxmindist);
    u(1)=x-xref(idxmindist);
    u(2)=y-yref(idxmindist);
    t=(vl(1)*u(1)+vl(2)*u(2))/vabs/vabs;
end