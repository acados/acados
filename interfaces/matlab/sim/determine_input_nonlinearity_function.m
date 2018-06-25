function [ L_x, L_xdot, L_z, L_u, phi_fun, y, uhat ] = ...
        determine_input_nonlinearity_function( x1, x1dot, z, u, phi )
%
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
% given a CasADi expression phi, which may depend on the variables (x1, 
% x1dot, z, u), this function determines a vector y (uhat) consisting of
% all components of (x, xdot, z) (respectively u) that enter phi.
% Additionally matrices L_x, L_xdot, L_z, L_u are determined such that
%           y    = L_x * x + L_xdot * xdot + L_z * z
%           uhat = L_u * u;
% and a CasADi function phi_fun is created, that maps
%           (y, uhat) -> phi

import casadi.*


nx1 = length(x1);
nz = length(z);
nu = length(u);

%% y
y = [];
% components of x
for ii = 1:nx1
    jac_xi = jacobian(phi, x1(ii));
    if jac_xi.is_zero()  % i.e. f_current does not depend on x(ii);
        % x(ii) is not part of y
    else
        y = vertcat(y, x1(ii));
    end
end
% components of xdot
for ii = 1:nx1
    jac_xidot = jacobian(phi, x1dot(ii));
    if jac_xidot.is_zero()  % i.e. f_current does not depend on xdot(ii);
        % xdot(ii) is not part of y
    else
        y = vertcat(y, x1dot(ii));
    end
end
% components of z
for ii = 1:nz
    jac_zi = jacobian(phi, z(ii));
    if jac_zi.is_zero()  % i.e. f_current does not depend on z(ii);
        % xdot(ii) is not part of y
    else
        y = vertcat(y, z(ii));
    end
end

%% uhat
uhat = [];
% components of u
for ii = 1:nu
    jac_ui = jacobian(phi, u(ii));
    if jac_ui.is_zero()  % i.e. f_current does not depend on u(ii);
        % u(ii) is not part of uhat
    else
        uhat = vertcat(uhat, u(ii));
    end
end

%% generate phi_fun;
phi = phi.simplify();
phi_fun = Function('phi_fun',{y,uhat}, {phi});
% linear input matrices
L_x_fun     = Function('L_x_fun',{x1},{jacobian(y,x1)});
L_xdot_fun  = Function('L_xdot_fun',{x1},{jacobian(y,x1dot)});
L_z_fun     = Function('L_z_fun',{x1},{jacobian(y,z)});
L_u_fun     = Function('L_u_fun',{x1},{jacobian(uhat,u)});

L_x = full(L_x_fun(0));
L_xdot = full(L_xdot_fun(0));
L_u = full(L_u_fun(0));
L_z = full(L_z_fun(0));


end

