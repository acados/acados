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

function dist = dist2D(x1,x2,y1,y2)
    dist = sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2));
end