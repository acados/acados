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

clear variables
close all
clc

generate_casadi_functions()

M = [   1.05000,   0.00000,   0.00000,  -0.10000; ... 
        0.00000,   1.05000,  -0.10000,   0.00000; ... 
       -0.05000,   0.00000,   1.00000,   0.00000; ... 
        0.00000,  -0.05000,   0.00000,   1.00000 ];

res = [   -1.00000 
          -0.05070; 
           1.00000; 
           0.05070;
      ];

x0 = [   -1.00000; 
         -0.05070 
 ];

xt = [   -1.00000; 
         -0.05070 
 ];

z0 = [0; 0];
u0 = [0; 0];
u0(1) = 0.0; u0(2) = 0.000000;

delta = M\res

K =     0*[    0.94735; 
            -0.04244; 
            -0.95263; 
            -0.05282; 
];

K = K - delta

%% Set up States & Controls
k_sol = K(1:2);
z_sol = K(3:4);

x_sol = x0 + 0.1*0.5*k_sol

f_eval     =        [ k_sol(1)+x_sol(1)-0.1*z_sol(2)-u0(1); ...
                      k_sol(2)+x_sol(2)-0.1*z_sol(1)-u0(2);  ...
                      z_sol(1)-x_sol(1); ...
                      z_sol(2)-x_sol(2)]
