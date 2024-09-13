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
% author: Katrin Baumgaertner





% NOTE: `acados` currently supports both an old MATLAB/Octave interface (< v0.4.0)
% as well as a new interface (>= v0.4.0).

% THIS EXAMPLE still uses the OLD interface. If you are new to `acados` please start
% with the examples that have been ported to the new interface already.
% see https://github.com/acados/acados/issues/1196#issuecomment-2311822122)



function [model] = lorentz_model()

import casadi.*

% system dimensions
nx = 4;               % last state models parameter
nw = 1;               % state noise on parameter
ny = 1;


% dynamics
x_expr = SX.sym('x', nx, 1);
w_expr = SX.sym('w', nw, 1);

f_expl_expr = [10*(x_expr(2) - x_expr(1)); ...
              x_expr(4)*x_expr(1)-x_expr(2)-x_expr(1)*x_expr(3); ...
              -(8/3)*x_expr(3)+x_expr(1)*x_expr(2); ...
              w_expr];

% store eveything in model struct
model = AcadosModel();

model.x = x_expr;
model.u = w_expr;
model.f_expl_expr = f_expl_expr;

model.name = 'lorentz_model_estimator';
end
