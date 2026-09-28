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

function [sref,xref,yref,psiref,kapparef] = getTrack(filename)
    array=dlmread(filename);
    sref=array(:,1);
    xref=array(:,2);
    yref=array(:,3);
    psiref=array(:,4);
    kapparef=array(:,5);
end