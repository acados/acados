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




function generate_c_code_discrete_dynamics( model, opts, model_dir )

import casadi.*

casadi_opts = struct('mex', false, 'casadi_int', 'int', 'casadi_real', 'double');
check_casadi_version();

%% load model
x = model.x;
u = model.u;
p = model.p;
nx = length(x);

if isempty(model.disc_dyn_expr)
    error('Field `disc_dyn_expr` is required for discrete dynamics.')
end
phi = model.disc_dyn_expr;
nx1 = length(phi);

if nx ~= nx1
    disp('Warning: generate_c_code_discrete_dynamics: got nx ~= nx1, this only works for a single shooting interval.');
end

% check type
if isa(x(1), 'casadi.SX')
    isSX = true;
else
    isSX = false;
end

model_name = model.name;

return_dir = pwd;
cd(model_dir)

% multipliers for hessian
if isSX
    lam = SX.sym('lam', nx1, 1);
else
    lam = MX.sym('lam', nx1, 1);
end
% generate jacobians
jac_ux = jacobian(phi, [u; x]);
% generate adjoint
adj_ux = jtimes(phi, [u; x], lam, true);
% generate hessian
hess_ux = jacobian(adj_ux, [u; x], struct('symmetric', isSX));
% Set up functions
phi_fun = Function([model_name,'_dyn_disc_phi_fun'], {x, u, p}, {phi});
phi_fun_jac_ut_xt = Function([model_name,'_dyn_disc_phi_fun_jac'], {x, u, p}, {phi, jac_ux'});
phi_fun_jac_ut_xt_hess = Function([model_name,'_dyn_disc_phi_fun_jac_hess'], {x, u, lam, p}, {phi, jac_ux', hess_ux});

% generate C code
phi_fun.generate([model_name,'_dyn_disc_phi_fun'], casadi_opts);
phi_fun_jac_ut_xt.generate([model_name,'_dyn_disc_phi_fun_jac'], casadi_opts);

if opts.generate_hess
    phi_fun_jac_ut_xt_hess.generate([model_name,'_dyn_disc_phi_fun_jac_hess'], casadi_opts);
end

cd(return_dir);

end