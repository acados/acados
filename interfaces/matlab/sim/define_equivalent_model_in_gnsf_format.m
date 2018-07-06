function [ gnsf ] = define_equivalent_model_in_gnsf_format(model, print_info)

%   This file is part of acados.
%
%   acados is free software; you can redistribute it and/or
%   modify it under the terms of the GNU Lesser General Public
%   License as published by the Free Software Foundation; either
%   version 3 of the License, or (at your option) any later version.
%
%   acados is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
%   Lesser General Public License for more details.
%
%   You should have received a copy of the GNU Lesser General Public
%   License along with acados; if not, write to the Free Software Foundation,
%   Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
%
%   Author: Jonathan Frey: jonathanpaulfrey(at)gmail.com

%% Description
% this function takes a model of an implicit ODE/ index-1 DAE and detects
% all linear terms and sets up an equivalent model in the GNSF structure,
% stored in the struct gnsf, where the linear output system (LOS) is empty.

% import CasADi
import casadi.*
if CasadiMeta.version()=='3.4.0'
	% casadi 3.4
	casadi_opts = struct('mex', false, 'casadi_int', 'int', 'casadi_real', 'double');
else
	% old casadi versions
	error('Please download and install Casadi 3.4.0 to ensure compatibility with acados')
end

% load model
f_impl_expr = model.f_impl_expr;

model_name_prefix = model.name;

% get intputs
x = model.x;
xdot = model.xdot;
u = model.u;
z = model.z;

% model dimensions
nx = length(x);
nz  = length(z);
nu  = length(u);

%% initialize gnsf struct
% dimensions
gnsf = struct('nx', nx, 'nu', nu, 'nz', nz);
gnsf.nx1 = nx;
gnsf.nx2 = 0;
gnsf.nuhat = nu;
gnsf.ny = 2 * nx + nz;

gnsf.phi_expr = f_impl_expr;
gnsf.A = zeros(nx + nz, nx);
gnsf.B = zeros(nx + nz, nu);
gnsf.E = zeros(nx + nz, nx + nz);
gnsf.c = zeros(nx + nz, 1);
gnsf.C = eye(nx + nz, nx + nz);
gnsf.name = model_name_prefix;

gnsf.x = x;
gnsf.xdot = xdot;
gnsf.z = z;
gnsf.u = u;

[ gnsf ] = determine_input_nonlinearity_function( gnsf );
    
gnsf.A_LO = [];
gnsf.f_lo_expr = [];
check_reformulation(model, gnsf, print_info);


%% Represent all affine dependencies through the model matrices A, B, E, c
%% determine A
n_nodes_current = gnsf.phi_expr.n_nodes();

for ii = 1:length(gnsf.phi_expr)
    fii = gnsf.phi_expr(ii);
    for ix = 1:nx
        % symbolic jacobian of fii w.r.t. xi;
        jac_fii_xi = jacobian(fii,x(ix));
        if jac_fii_xi.is_constant
            % jacobian value
            jac_fii_xi_fun = Function(['jac_fii_xi_fun'],...
                            {x(1)}, {jac_fii_xi});
            % x(1) as input just to have a scalar input and call the function as follows:
            gnsf.A(ii, ix) = full(jac_fii_xi_fun(0));                
        else
            gnsf.A(ii, ix) = 0;
            if print_info
                disp(['fii is NOT linear in x(ix) ii = ', num2str(ii),',  ix = ', num2str(ix)]);
                disp(fii)
                disp('================================');
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

% keyboard % here it fails for MX

check_reformulation(model, gnsf, print_info);


%% determine B

n_nodes_current = gnsf.phi_expr.n_nodes();

for ii = 1:length(gnsf.phi_expr)
    fii = gnsf.phi_expr(ii);
    for iu = 1:length(u)
        % symbolic jacobian of fii w.r.t. ui;
        jac_fii_ui = jacobian(fii, u(iu));
        if jac_fii_ui.is_constant % i.e. hessian is structural zero
            % jacobian value
            jac_fii_ui_fun = Function(['jac_fii_ui_fun'],...
                            {x(1)}, {jac_fii_ui});
            gnsf.B(ii, iu) = full(jac_fii_ui_fun(0));                
        else
            gnsf.B(ii, iu) = 0;
            if print_info
                disp(['fii is NOT linear in u(iu) ii = ', num2str(ii),',  iu = ', num2str(iu)]);
                disp(fii)
                disp('================================');
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


%% determine E
n_nodes_current = gnsf.phi_expr.n_nodes();
k = vertcat(xdot, z);

for ii = 1:length(gnsf.phi_expr)
    fii = gnsf.phi_expr(ii);
    for ik = 1:length(k)
        % symbolic jacobian of fii w.r.t. ui;
        jac_fii_ki = jacobian(fii, k(ik));
        if jac_fii_ki.is_constant
            % jacobian value
            jac_fii_ki_fun = Function(['jac_fii_ki_fun'],...
                        {x(1)}, {jac_fii_ki});
            gnsf.E(ii, ik) = - full(jac_fii_ki_fun(0));                
        else
            gnsf.E(ii, ik) = 0;
            if print_info
                disp(['fii is NOT linear in k(ik) for ii = ', num2str(ii),',  ik = ', num2str(ik)]);
                disp(fii)
                disp('================================');
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
            disp(['fii is NOT constant for ii = ', num2str(ii)]);
            disp(fii)
            disp('================================');
        end
    end
end

f_next = gnsf.phi_expr - gnsf.c;
% f_next = f_next.simplify();
n_nodes_next = f_next.n_nodes();

if print_info
    fprintf('\n')
    disp(['determined vector c:']);
    disp(gnsf.c)
    disp(['reduced nonlinearity from  ', num2str(n_nodes_current),...
          ' to ', num2str(n_nodes_next) ' nodes']);
end

gnsf.phi_expr = f_next;
check_reformulation(model, gnsf, print_info);


%% determine nonlinearity & corresponding matrix C
%% Reduce dimension of f
n_nodes_current = gnsf.phi_expr.n_nodes();
ind_non_zero = [];
for ii = 1:length(gnsf.phi_expr)
    fii = gnsf.phi_expr(ii);
    fii = fii.simplify();
    if ~fii.is_zero
        ind_non_zero = [ind_non_zero, ii];
    end
end

clear f_next
f_next = gnsf.phi_expr(ind_non_zero);

% C
gnsf.C = zeros(nx+nz, length(ind_non_zero));
for ii = 1:length(ind_non_zero)
    gnsf.C(ind_non_zero(ii), ii) = 1;
end

n_nodes_next = f_next.n_nodes();

% assert(n_nodes_current >= n_nodes_next,'n_nodes_current >= n_nodes_next FAILED')
gnsf.phi_expr = f_next;

check_reformulation(model, gnsf, print_info);

gnsf.n_out = length(gnsf.phi_expr);

if print_info
    fprintf('\n')
    disp(['determined matrix C:']);
    disp(gnsf.C)
    disp(['reduced nonlinearity dimension n_out from  ', num2str(nx+nz),'   to  ', num2str(gnsf.n_out)]);
    disp(['reduced nodes in CasADi expr of nonlinearity from  ', num2str(n_nodes_current),...
          ' to ', num2str(n_nodes_next) ' nodes']);
    disp(gnsf.phi_expr);
end

%% determine input of nonlinearity function
[ gnsf ] = determine_input_nonlinearity_function( gnsf );

check_reformulation(model, gnsf, print_info);

gnsf.ny = length(gnsf.y);
gnsf.nuhat = length(gnsf.uhat);

if print_info
    disp('----------------------------------');
    disp(['reduced input ny    of phi from  ', num2str(2*nx+nz),'   to  ', num2str( gnsf.ny )]);
    disp(['reduced input nuhat of phi from  ', num2str(nu),'   to  ', num2str( gnsf.nuhat )]);
    disp('----------------------------------');
end

