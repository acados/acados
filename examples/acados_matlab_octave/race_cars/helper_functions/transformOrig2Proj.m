%
% Copyright (c) The acados authors.
%
% This file is part of acados.
%
% The 2-Clause BSD License
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
%
% 1. Redistributions of source code must retain the above copyright notice,
% this list of conditions and the following disclaimer.
%
% 2. Redistributions in binary form must reproduce the above copyright notice,
% this list of conditions and the following disclaimer in the documentation
% and/or other materials provided with the distribution.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.;

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