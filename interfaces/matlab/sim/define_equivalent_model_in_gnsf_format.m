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

x = model.x;
xdot = model.xdot;
u = model.u;
z = model.z;
model_name_prefix = model.name;

% set up implicit function
f_impl_fun = Function([model_name_prefix,'f_impl_fun'], {x, xdot, u, z}, {f_impl_expr});

% model dimensions
x1 = x;
x1dot = xdot;
nx = length(x);
nx1 = nx;
nx2 = 0;
nz  = length(z);
nu  = length(u);

nuhat = nu;
ny = 2 * nx + nz;


% set up vectors to eval the functions
x0 = rand(size(x));
x0dot = rand(size(x));
u0 = rand(size(u));
z0 = rand(size(z));


%% initialize gnsf struct
gnsf = struct('nx', nx, 'nu', nu, 'nz', nz, 'nx1', nx1, 'nx2', nx2, ...
    'nuhat', nuhat, 'ny', ny, 'n_out', nx+nz);

x_old = x;
f_impl_old = f_impl_expr;

phi_current = f_impl_expr;
A = zeros(nx + nz, nx);
B = zeros(nx + nz, nu);
E = zeros(nx + nz, nx + nz);
c = zeros(nx + nz, 1);
C = eye(nx + nz, nx + nz);

gnsf.A = A;
gnsf.B = B;
gnsf.C = C;
gnsf.E = E;
gnsf.c = c;


[ gnsf.L_x, gnsf.L_xdot, gnsf.L_z, gnsf.L_u, gnsf.phi_fun, y, uhat ] = ...
        determine_input_nonlinearity_function( x1, x1dot, z, u, phi_current );
gnsf.phi_expr = [];
    
gnsf.A_LO = [];
gnsf.f_lo_fun = Function('f_lo_fun',{x, xdot, z, u}, {[]});

check = check_reformulation(f_impl_fun, gnsf, print_info);


%% Represent all affine dependencies through the model matrices A, B, E, c
%% determine A
n_nodes_current = phi_current.n_nodes();
n_nodes_initial = n_nodes_current;

for ii = 1:length(phi_current)
    fii = phi_current(ii);
    for ix = 1:nx
        % symbolic jacobian of fii w.r.t. xi;
        jac_fii_xi = jacobian(fii,x(ix));
        if jac_fii_xi.is_constant
            % jacobian value
            jac_fii_xi_fun = Function(['jac_fii_xi_fun'],...
                            {x, xdot, u, z}, {jac_fii_xi});
            A(ii, ix) = full(jac_fii_xi_fun(x0,x0dot, u0, z0));                
        else
            A(ii, ix) = 0;
            if print_info
                disp(['fii is NOT linear in x(ix) ii = ', num2str(ii),',  ix = ', num2str(ix)]);
                disp(fii)
                disp('================================');
            end
        end
    end
end

f_next = phi_current - A * x;
f_next = f_next.simplify();
n_nodes_next = f_next.n_nodes();

if print_info
    fprintf('\n')
    disp(['determined matrix A:']);
    disp(A)
    disp(['reduced nonlinearity from  ', num2str(n_nodes_current),...
          ' to ', num2str(n_nodes_next) ' nodes']);
end

% assert(n_nodes_current >= n_nodes_next,'n_nodes_current >= n_nodes_next FAILED')
phi_current = f_next;


gnsf.phi_fun = Function('phi_fun',{y,uhat}, {phi_current});
gnsf.A = A;
check = check_reformulation(f_impl_fun, gnsf, print_info);


%% determine B

n_nodes_current = phi_current.n_nodes();

for ii = 1:length(phi_current)
    fii = phi_current(ii);
    for iu = 1:length(u)
        % symbolic jacobian of fii w.r.t. ui;
        jac_fii_ui = jacobian(fii, u(iu));
        if jac_fii_ui.is_constant % i.e. hessian is structural zero
            % jacobian value
            jac_fii_ui_fun = Function(['jac_fii_ui_fun'],...
                            {x, xdot, u, z}, {jac_fii_ui});
            B(ii, iu) = full(jac_fii_ui_fun(x0,x0dot, u0, z0));                
        else
            B(ii, iu) = 0;
            if print_info
                disp(['fii is NOT linear in u(iu) ii = ', num2str(ii),',  iu = ', num2str(iu)]);
                disp(fii)
                disp('================================');
            end
        end
    end
end

f_next = phi_current - B * u;
f_next = f_next.simplify();
n_nodes_next = f_next.n_nodes();

if print_info
    fprintf('\n')
    disp(['determined matrix B:']);
    disp(B)
    disp(['reduced nonlinearity from  ', num2str(n_nodes_current),...
          ' to ', num2str(n_nodes_next) ' nodes']);
end

phi_current = f_next;

gnsf.phi_fun = Function('phi_fun',{y,uhat}, {phi_current});
gnsf.B = B;
check = check_reformulation(f_impl_fun, gnsf, print_info);


%% determine E
n_nodes_current = phi_current.n_nodes();
k = vertcat(xdot, z);

for ii = 1:length(phi_current)
    fii = phi_current(ii);
    for ik = 1:length(k)
        % symbolic jacobian of fii w.r.t. ui;
        jac_fii_ki = jacobian(fii, k(ik));
        if jac_fii_ki.is_constant
            % jacobian value
            jac_fii_ki_fun = Function(['jac_fii_ki_fun'],...
                        {x, xdot, u, z}, {jac_fii_ki});
            E(ii, ik) = - full(jac_fii_ki_fun(x0,x0dot, u0, z0));                
        else
            E(ii, ik) = 0;
            if print_info
                disp(['fii is NOT linear in k(ik) for ii = ', num2str(ii),',  ik = ', num2str(ik)]);
                disp(fii)
                disp('================================');
            end
        end
    end
end

f_next = phi_current + E * k;
f_next = f_next.simplify();
n_nodes_next = f_next.n_nodes();

if print_info
    fprintf('\n')
    disp(['determined matrix E:']);
    disp(E)
    disp(['reduced nonlinearity from  ', num2str(n_nodes_current),...
          ' to ', num2str(n_nodes_next) ' nodes']);
end

phi_current = f_next;

gnsf.phi_fun = Function('phi_fun',{y,uhat}, {phi_current});
gnsf.E = E;

check = check_reformulation(f_impl_fun, gnsf, print_info);

%% determine constant term c

n_nodes_current = phi_current.n_nodes();
for ii = 1:length(phi_current)
    fii = phi_current(ii);
    if fii.is_constant
        % function value goes into c
        fii_fun = Function(['fii_fun'],...
                {x, xdot, u, z}, {fii});
        c(ii) = full(fii_fun(x0,x0dot, u0, z0));                
    else
        c(ii) = 0;
        if print_info
            disp(['fii is NOT constant for ii = ', num2str(ii)]);
            disp(fii)
            disp('================================');
        end
    end
end

f_next = phi_current - c;
% f_next = f_next.simplify();
n_nodes_next = f_next.n_nodes();

if print_info
    fprintf('\n')
    disp(['determined vector c:']);
    disp(c)
    disp(['reduced nonlinearity from  ', num2str(n_nodes_current),...
          ' to ', num2str(n_nodes_next) ' nodes']);
end

phi_current = f_next;

gnsf.phi_fun = Function('phi_fun',{y,uhat}, {phi_current});
gnsf.c = c;
check = check_reformulation(f_impl_fun, gnsf, print_info);


%% determine nonlinearity & corresponding matrix C
%% Reduce dimension of f
n_nodes_current = phi_current.n_nodes();
ind_non_zero = [];
for ii = 1:length(phi_current)
    fii = phi_current(ii);
    fii = fii.simplify();
    if ~fii.is_zero
        ind_non_zero = [ind_non_zero, ii];
    end
end

clear f_next
f_next = phi_current(ind_non_zero);

% C
C = zeros(nx+nz, length(ind_non_zero));
for ii = 1:length(ind_non_zero)
    C(ind_non_zero(ii), ii) = 1;
end

n_nodes_next = f_next.n_nodes();

% assert(n_nodes_current >= n_nodes_next,'n_nodes_current >= n_nodes_next FAILED')
phi_current = f_next;

gnsf.phi_fun = Function('phi_fun',{y,uhat}, {phi_current});
gnsf.C = C;
check = check_reformulation(f_impl_fun, gnsf, print_info);

gnsf.n_out = length(phi_current);

if print_info
    fprintf('\n')
    disp(['determined matrix C:']);
    disp(C)
    disp(['reduced nonlinearity dimension n_out from  ', num2str(nx+nz),'   to  ', num2str(gnsf.n_out)]);
    disp(['reduced nodes in CasADi expr of nonlinearity from  ', num2str(n_nodes_current),...
          ' to ', num2str(n_nodes_next) ' nodes']);
    disp(phi_current);
end

%% determine input of nonlinearity function
[ gnsf.L_x, gnsf.L_xdot, gnsf.L_z, gnsf.L_u, gnsf.phi_fun, y, uhat ] = ...
       determine_input_nonlinearity_function( x1, x1dot, z, u, phi_current );

check = check_reformulation(f_impl_fun, gnsf, print_info);
ny = length(y);
nuhat = length(uhat);
gnsf.ny = ny;
gnsf.nuhat = nuhat;

gnsf.y = y;
gnsf.uhat = uhat;
gnsf.phi_expr = phi_current;

gnsf.x = x;
gnsf.xdot = xdot;
gnsf.z = z;
gnsf.u = u;

if print_info
    disp('----------------------------------');
    disp(['reduced input ny    of phi from  ', num2str(2*nx+nz),'   to  ', num2str( ny )]);
    disp(['reduced input nuhat of phi from  ', num2str(nu),'   to  ', num2str( nuhat )]);
    disp('----------------------------------');
end

