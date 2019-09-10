%
% Copyright 2019 Gianluca Frison, Dimitris Kouzoupis, Robin Verschueren,
% Andrea Zanelli, Niels van Duijkeren, Jonathan Frey, Tommaso Sartor,
% Branimir Novoselnik, Rien Quirynen, Rezart Qelibari, Dang Doan,
% Jonas Koenemann, Yutao Chen, Tobias Schöls, Jonas Schlagenhauf, Moritz Diehl
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


function generate_c_code_implicit_ode( model, opts )

%% import casadi
import casadi.*

casadi_version = CasadiMeta.version();
if strcmp(casadi_version(1:3),'3.4') % require casadi 3.4.x
	casadi_opts = struct('mex', false, 'casadi_int', 'int', 'casadi_real', 'double');
else % old casadi versions
	error('Please download and install CasADi version 3.4.x to ensure compatibility with acados')
end

if nargin > 1
    if isfield(opts, 'sens_hess')
        generate_hess = opts.sens_hess;
    else
        generate_hess = 'false';
%        if opts.print_info
%        disp('generate_hess option was not set - default is false')
%        end
    end
else
    generate_hess = 'false';
end
generate_hess = 'true'; % TODO remove when not needed any more !!!


%% load model
% x
x = model.sym_x;
nx = length(x);
% check type
if class(x(1)) == 'casadi.SX'
    isSX = true;
else
    isSX = false;
end
% xdot
xdot = model.sym_xdot;
% u
if isfield(model, 'sym_u')
    u = model.sym_u;
	nu = length(u);
else
    if isSX
        u = SX.sym('u',0, 0);
    else
        u = MX.sym('u',0, 0);
    end
	nu = 0;
end
% z
if isfield(model, 'sym_z')
    z = model.sym_z;
	nz = length(z);
else
    if isSX
        z = SX.sym('z',0, 0);
    else
        z = MX.sym('z',0, 0);
    end
	nz = 0;
end
% p
if isfield(model, 'sym_p')
    p = model.sym_p;
	np = length(p);
else
    if isSX
        p = SX.sym('p',0, 0);
    else
        p = MX.sym('p',0, 0);
    end
	np = 0;
end


model_name = model.name;

if isfield(model, 'dyn_expr_f')
	f_impl = model.dyn_expr_f;
	model_name = [model_name, '_dyn'];
else
	f_impl = model.expr_f;
end



%% generate jacobians
jac_x       = jacobian(f_impl, x);
jac_xdot    = jacobian(f_impl, xdot);
jac_u       = jacobian(f_impl, u);
jac_z       = jacobian(f_impl, z);


%% generate hessian
x_xdot_z_u = [x; xdot; z; u];

if isSX
    multiplier  = SX.sym('multiplier', length(x) + length(z));
%    multiply_mat  = SX.sym('multiply_mat', 2*nx+nz+nu, nx + nu);
%    HESS = SX.zeros( length(x_xdot_z_u), length(x_xdot_z_u));
else
    multiplier  = MX.sym('multiplier', length(x) + length(z));
%    multiply_mat  = MX.sym('multiply_mat', 2*nx+nz+nu, nx + nu);
%    HESS = MX.zeros( length(x_xdot_z_u), length(x_xdot_z_u));
end

% hessian computed as forward over adjoint !!!
ADJ = jtimes(f_impl, x_xdot_z_u, multiplier, true);
HESS = jacobian(ADJ, x_xdot_z_u);

%HESS_multiplied = multiply_mat' * HESS * multiply_mat;

%HESS = jtimes(ADJ, x_xdot_z_u, multiply_mat);
%HESS_multiplied = multiply_mat' * HESS;

%HESS_multiplied = HESS_multiplied.simplify();
%HESS_multiplied = HESS; % do the multiplication in BLASFEO !!!



%% Set up functions
% TODO(oj): fix namings such that jac_z is contained!
impl_ode_fun = Function([model_name,'_impl_ode_fun'], {x, xdot, u, z, p}, {f_impl});
impl_ode_fun_jac_x_xdot = Function([model_name,'_impl_ode_fun_jac_x_xdot'], {x, xdot, u, z, p}, {f_impl, jac_x, jac_xdot, jac_z});
impl_ode_jac_x_xdot_u = Function([model_name,'_impl_ode_jac_x_xdot_u'], {x, xdot, u, z, p}, {jac_x, jac_xdot, jac_u, jac_z});
impl_ode_fun_jac_x_xdot_u = Function([model_name,'_impl_ode_fun_jac_x_xdot_u'], {x, xdot, u, z, p}, {f_impl, jac_x, jac_xdot, jac_u});
%    impl_ode_hess = Function([model_name,'_impl_ode_hess'],  {x, xdot, u, z, multiplier, multiply_mat, p}, {HESS_multiplied});
impl_ode_hess = Function([model_name,'_impl_ode_hess'],  {x, xdot, u, z, multiplier, p}, {HESS});

%% generate C code
impl_ode_fun.generate([model_name,'_impl_ode_fun'], casadi_opts);
impl_ode_fun_jac_x_xdot.generate([model_name,'_impl_ode_fun_jac_x_xdot'], casadi_opts);
impl_ode_jac_x_xdot_u.generate([model_name,'_impl_ode_jac_x_xdot_u'], casadi_opts);
impl_ode_fun_jac_x_xdot_u.generate([model_name,'_impl_ode_fun_jac_x_xdot_u'], casadi_opts);
if strcmp(generate_hess, 'true')
    impl_ode_hess.generate([model_name,'_impl_ode_hess'], casadi_opts);
end
% keyboard

end
