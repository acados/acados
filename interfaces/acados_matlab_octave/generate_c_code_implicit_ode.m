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


function generate_c_code_implicit_ode( model, opts )

%% import casadi
import casadi.*

casadi_opts = struct('mex', false, 'casadi_int', 'int', 'casadi_real', 'double');
check_casadi_version();

if nargin > 1
    if isfield(opts, 'sens_hess')
        generate_hess = opts.sens_hess;
    elseif isfield(opts, 'nlp_solver_exact_hessian')
        generate_hess = opts.nlp_solver_exact_hessian;
%        if opts.print_info
%        disp('generate_hess option was not set - default is false')
%        end
    else
        generate_hess = 'false';
    end
else
    generate_hess = 'false';
end


%% load model
% x
is_template = false;
if isa(model, 'acados_template_mex.acados_model_json')
    is_template = true;
    % names without sym
    x = model.x;
    nx = length(x);
    % check type
    if isa(x(1), 'casadi.SX')
        isSX = true;
    else
        isSX = false;
    end
    % u
    u = model.u;
    nu = length(u);
    % p
    p = model.p;
    np = length(p);
    % xdot
    xdot = model.xdot;
    % z
    z = model.z;

else
    x = model.sym_x;
    nx = length(x);
    % check type
    if isa(x(1), 'casadi.SX')
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
    else
        if isSX
            z = SX.sym('z',0, 0);
        else
            z = MX.sym('z',0, 0);
        end
    end
    % p
    if isfield(model, 'sym_p')
        p = model.sym_p;
    else
        if isSX
            p = SX.sym('p',0, 0);
        else
            p = MX.sym('p',0, 0);
        end
    end
end
nz = length(z);
np = length(p);

model_name = model.name;

if isfield(model, 'dyn_expr_f')
    f_impl = model.dyn_expr_f;
    model_name = [model_name, '_dyn'];
elseif isfield(model, 'expr_f')
    f_impl = model.expr_f;
else
    f_impl = model.f_impl_expr;
end


%% generate jacobians
jac_x       = jacobian(f_impl, x);
jac_xdot    = jacobian(f_impl, xdot);
jac_u       = jacobian(f_impl, u);
jac_z       = jacobian(f_impl, z);


%% generate hessian
x_xdot_z_u = [x; xdot; z; u];

if isSX
    multiplier  = SX.sym('multiplier', nx + nz);
%    multiply_mat  = SX.sym('multiply_mat', 2*nx+nz+nu, nx + nu);
%    HESS = SX.zeros( length(x_xdot_z_u), length(x_xdot_z_u));
else
    multiplier  = MX.sym('multiplier', nx + nz);
%    multiply_mat  = MX.sym('multiply_mat', 2*nx+nz+nu, nx + nu);
%    HESS = MX.zeros( length(x_xdot_z_u), length(x_xdot_z_u));
end



%% Set up functions
impl_dae_fun = Function([model_name,'_impl_dae_fun'], {x, xdot, u, z, p}, {f_impl});
impl_dae_fun_jac_x_xdot_z = Function([model_name,'_impl_dae_fun_jac_x_xdot_z'], {x, xdot, u, z, p}, {f_impl, jac_x, jac_xdot, jac_z});
impl_dae_jac_x_xdot_u_z = Function([model_name,'_impl_dae_jac_x_xdot_u_z'], {x, xdot, u, z, p}, {jac_x, jac_xdot, jac_u, jac_z});
impl_dae_fun_jac_x_xdot_u = Function([model_name,'_impl_dae_fun_jac_x_xdot_u'], {x, xdot, u, z, p}, {f_impl, jac_x, jac_xdot, jac_u});


if is_template
    if ~exist(fullfile(pwd,'c_generated_code'), 'dir')
        mkdir('c_generated_code');
    end
    cd 'c_generated_code'
    model_dir = [model_name, '_model'];
    if ~exist(fullfile(pwd, model_dir), 'dir')
        mkdir(model_dir);
    end
    cd(model_dir)
end

%% generate C code
impl_dae_fun.generate([model_name,'_impl_dae_fun'], casadi_opts);
impl_dae_fun_jac_x_xdot_z.generate([model_name,'_impl_dae_fun_jac_x_xdot_z'], casadi_opts);
impl_dae_jac_x_xdot_u_z.generate([model_name,'_impl_dae_jac_x_xdot_u_z'], casadi_opts);
impl_dae_fun_jac_x_xdot_u.generate([model_name,'_impl_dae_fun_jac_x_xdot_u'], casadi_opts);
if strcmp(generate_hess, 'true')
    % hessian computed as forward over adjoint !!!
    ADJ = jtimes(f_impl, x_xdot_z_u, multiplier, true);
    HESS = jacobian(ADJ, x_xdot_z_u);

    %HESS_multiplied = multiply_mat' * HESS * multiply_mat;

    %HESS = jtimes(ADJ, x_xdot_z_u, multiply_mat);
    %HESS_multiplied = multiply_mat' * HESS;

    %HESS_multiplied = HESS_multiplied.simplify();
    %HESS_multiplied = HESS; % do the multiplication in BLASFEO !!!
    %    impl_dae_hess = Function([model_name,'_impl_dae_hess'],  {x, xdot, u, z, multiplier, multiply_mat, p}, {HESS_multiplied});

    impl_dae_hess = Function([model_name,'_impl_dae_hess'],...
                             {x, xdot, u, z, multiplier, p}, {HESS});
    impl_dae_hess.generate([model_name,'_impl_dae_hess'], casadi_opts);
end

if is_template
    cd '../..'
end

% keyboard

end
