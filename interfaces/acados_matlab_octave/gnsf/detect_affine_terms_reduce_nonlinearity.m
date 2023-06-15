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

function [ gnsf ] = detect_affine_terms_reduce_nonlinearity( gnsf, model, print_info )

%% Description
% this function takes a gnsf structure with trivial model matrices (A, B,
% E, c are zeros, and C is eye).
% It detects all affine linear terms and sets up an equivalent model in the
% GNSF structure, where all affine linear terms are modeled through the
% matrices A, B, E, c and the linear output system (LOS) is empty.
% NOTE: model is just taken as an argument to check equivalence of the
% models within the function.

import casadi.*

if print_info
disp(' ');
disp('====================================================================');
disp(' ');
disp('============  Detect affine-linear dependencies   ==================');
disp(' ');
disp('====================================================================');
disp(' ');
end

% symbolics
x = gnsf.x;
xdot = gnsf.xdot;
u = gnsf.u;
z = gnsf.z;

% dimensions
nx = gnsf.nx;
nu = gnsf.nu;
nz = gnsf.nz;

ny_old = gnsf.ny;
nuhat_old =gnsf.nuhat;

%% Represent all affine dependencies through the model matrices A, B, E, c
%% determine A
n_nodes_current = gnsf.phi_expr.n_nodes();

for ii = 1:length(gnsf.phi_expr)
    fii = gnsf.phi_expr(ii);
    for ix = 1:nx
        var = x(ix);
        varname = var.name;
        % symbolic jacobian of fii w.r.t. xi;
        jac_fii_xi = jacobian(fii,var);
        if jac_fii_xi.is_constant
            % jacobian value
            jac_fii_xi_fun = Function(['jac_fii_xi_fun'],...
                            {x(1)}, {jac_fii_xi});
            % x(1) as input just to have a scalar input and call the function as follows:
            gnsf.A(ii, ix) = full(jac_fii_xi_fun(0));
        else
            gnsf.A(ii, ix) = 0;
            if print_info
                disp(['phi(' num2str(ii) ') is nonlinear in x(', num2str(ix), ') = ' varname]);
                disp(fii)
                disp('-----------------------------------------------------');
            end
        end
    end
end

f_next = gnsf.phi_expr - gnsf.A * x;
f_next = f_next.simplify();
n_nodes_next = f_next.n_nodes();

if print_info
    fprintf('\n')
    disp(['determined matrix A:']);
    disp(gnsf.A)
    disp(['reduced nonlinearity from  ', num2str(n_nodes_current),...
          ' to ', num2str(n_nodes_next) ' nodes']);
end

% assert(n_nodes_current >= n_nodes_next,'n_nodes_current >= n_nodes_next FAILED')
gnsf.phi_expr = f_next;

check_reformulation(model, gnsf, print_info);


%% determine B
if nu > 0
    n_nodes_current = gnsf.phi_expr.n_nodes();

    for ii = 1:length(gnsf.phi_expr)
        fii = gnsf.phi_expr(ii);
        for iu = 1:nu
            var = u(iu);
            varname = var.name;
            % symbolic jacobian of fii w.r.t. ui;
            jac_fii_ui = jacobian(fii, var);
            if jac_fii_ui.is_constant % i.e. hessian is structural zero
                % jacobian value
                jac_fii_ui_fun = Function(['jac_fii_ui_fun'],...
                                {x(1)}, {jac_fii_ui});
                gnsf.B(ii, iu) = full(jac_fii_ui_fun(0));
            else
                gnsf.B(ii, iu) = 0;
                if print_info
                    disp(['phi(' num2str(ii) ') is nonlinear in u(', num2str(iu), ') = ' varname]);
                    disp(fii)
                    disp('-----------------------------------------------------');
                end
            end
        end
    end
    f_next = gnsf.phi_expr - gnsf.B * u;
    f_next = f_next.simplify();
    n_nodes_next = f_next.n_nodes();

    if print_info
        fprintf('\n')
        disp(['determined matrix B:']);
        disp(gnsf.B)
        disp(['reduced nonlinearity from  ', num2str(n_nodes_current),...
              ' to ', num2str(n_nodes_next) ' nodes']);
    end

    gnsf.phi_expr = f_next;

    check_reformulation(model, gnsf, print_info);
end

%% determine E
n_nodes_current = gnsf.phi_expr.n_nodes();
k = vertcat(xdot, z);

for ii = 1:length(gnsf.phi_expr)
    fii = gnsf.phi_expr(ii);
    for ik = 1:length(k)
        % symbolic jacobian of fii w.r.t. ui;
        var = k(ik);
        varname = var.name;
        jac_fii_ki = jacobian(fii, var);
        if jac_fii_ki.is_constant
            % jacobian value
            jac_fii_ki_fun = Function(['jac_fii_ki_fun'],...
                        {x(1)}, {jac_fii_ki});
            gnsf.E(ii, ik) = - full(jac_fii_ki_fun(0));
        else
            gnsf.E(ii, ik) = 0;
            if print_info
                disp(['phi(' num2str(ii) ') is nonlinear in xdot_z(', num2str(ik), ') = ' varname]);
                disp(fii)
                disp('-----------------------------------------------------');
            end
        end
    end
end

f_next = gnsf.phi_expr + gnsf.E * k;
f_next = f_next.simplify();
n_nodes_next = f_next.n_nodes();

if print_info
    fprintf('\n')
    disp(['determined matrix E:']);
    disp(gnsf.E)
    disp(['reduced nonlinearity from  ', num2str(n_nodes_current),...
          ' to ', num2str(n_nodes_next) ' nodes']);
end

gnsf.phi_expr = f_next;
check_reformulation(model, gnsf, print_info);

%% determine constant term c

n_nodes_current = gnsf.phi_expr.n_nodes();
for ii = 1:length(gnsf.phi_expr)
    fii = gnsf.phi_expr(ii);
    if fii.is_constant
        % function value goes into c
        fii_fun = Function(['fii_fun'],...
                {x(1)}, {fii});
        gnsf.c(ii) = full(fii_fun(0));
    else
        gnsf.c(ii) = 0;
        if print_info
            disp(['phi(',num2str(ii),') is NOT constant']);
            disp(fii)
            disp('-----------------------------------------------------');
        end
    end
end

gnsf.phi_expr = gnsf.phi_expr - gnsf.c;
gnsf.phi_expr = gnsf.phi_expr.simplify();
n_nodes_next = gnsf.phi_expr.n_nodes();

if print_info
    fprintf('\n')
    disp(['determined vector c:']);
    disp(gnsf.c)
    disp(['reduced nonlinearity from  ', num2str(n_nodes_current),...
          ' to ', num2str(n_nodes_next) ' nodes']);
end

check_reformulation(model, gnsf, print_info);


%% determine nonlinearity & corresponding matrix C
%% Reduce dimension of phi
n_nodes_current = gnsf.phi_expr.n_nodes();
ind_non_zero = [];
for ii = 1:length(gnsf.phi_expr)
    fii = gnsf.phi_expr(ii);
    fii = fii.simplify();
    if ~fii.is_zero
        ind_non_zero = union(ind_non_zero, ii);
    end
end

gnsf.phi_expr = gnsf.phi_expr(ind_non_zero);

% C
gnsf.C = zeros(nx+nz, length(ind_non_zero));
for ii = 1:length(ind_non_zero)
    gnsf.C(ind_non_zero(ii), ii) = 1;
end

gnsf = determine_input_nonlinearity_function( gnsf );
n_nodes_next = gnsf.phi_expr.n_nodes();


if print_info
    disp(' ');
    disp('determined matrix C:');
    disp(gnsf.C)
    disp('---------------------------------------------------------------------------------');
    disp('------------- Success: Affine linear terms detected -----------------------------');
    disp('---------------------------------------------------------------------------------');
    disp(['reduced nonlinearity dimension n_out from  ', num2str(nx+nz),'   to  ', num2str(gnsf.n_out)]);
    disp(['reduced nodes in CasADi expr of nonlinearity from  ', num2str(n_nodes_current),...
          ' to ', num2str(n_nodes_next) ' nodes']);
    disp(' ');
    disp('phi now reads as:')
    print_casadi_expression(gnsf.phi_expr);
end

%% determine input of nonlinearity function

check_reformulation(model, gnsf, print_info);

gnsf.ny = length(gnsf.y);
gnsf.nuhat = length(gnsf.uhat);

if print_info
    disp('-----------------------------------------------------------------------------------');
        disp(' ');
    disp(['reduced input ny    of phi from  ', num2str(ny_old),'   to  ', num2str( gnsf.ny )]);
    disp(['reduced input nuhat of phi from  ', num2str(nuhat_old),'   to  ', num2str( gnsf.nuhat )]);
    disp('-----------------------------------------------------------------------------------');
end

end

