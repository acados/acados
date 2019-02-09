function [ gnsf ] = determine_input_nonlinearity_function( gnsf )
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
% this function takes a structure gnsf and updates the matrices L_x,
% L_xdot, L_z, L_u and CasADi vectors y, uhat of this structure as follows:

% given a CasADi expression phi_expr, which may depend on the variables 
% (x1, x1dot, z, u), this function determines a vector y (uhat) consisting 
% of all components of (x1, x1dot, z) (respectively u) that enter phi_expr.
% Additionally matrices L_x, L_xdot, L_z, L_u are determined such that
%           y    = L_x * x + L_xdot * xdot + L_z * z
%           uhat = L_u * u;
% Furthermore the dimensions ny, nuhat, n_out are updated

import casadi.*


%% y
y = [];
% components of x1
for ii = 1:gnsf.nx1
    if gnsf.phi_expr.which_depends(gnsf.x(ii))
        y = vertcat(y, gnsf.x(ii));
    else
        % x(ii) is not part of y
    end
end

% components of x1dot
for ii = 1:gnsf.nx1
    if gnsf.phi_expr.which_depends(gnsf.xdot(ii))
        y = vertcat(y, gnsf.xdot(ii));
    else
        % xdot(ii) is not part of y
    end
end

% components of z
for ii = 1:gnsf.nz1
    if gnsf.phi_expr.which_depends(gnsf.z(ii))
        y = vertcat(y, gnsf.z(ii));
    else
        % z(ii) is not part of y
    end
end

%% uhat
uhat = [];
% components of u
for ii = 1:gnsf.nu
    if gnsf.phi_expr.which_depends(gnsf.u(ii))
        uhat = vertcat(uhat, gnsf.u(ii));
    else
        % u(ii) is not part of uhat
    end
end

%% generate gnsf.phi_expr_fun;
% linear input matrices
if isempty(y)
    gnsf.L_x = [];
    gnsf.L_xdot = [];
    gnsf.L_u = [];
    gnsf.L_z = [];
else
    dummy = SX.sym('dummy_input',0);
    L_x_fun     = Function('L_x_fun', {dummy}, ...
                        {jacobian( y, gnsf.x(1:gnsf.nx1)) });
    L_xdot_fun  = Function('L_xdot_fun', {dummy}, ...
                        {jacobian( y, gnsf.xdot(1:gnsf.nx1) )});
    L_z_fun     = Function('L_z_fun', {dummy},...
                        {jacobian(y, gnsf.z(1:gnsf.nz1) )});
    L_u_fun     = Function('L_u_fun', {dummy},...
                        {jacobian(uhat, gnsf.u)});

    gnsf.L_x = full(L_x_fun(0));
    gnsf.L_xdot = full(L_xdot_fun(0));
    gnsf.L_u = full(L_u_fun(0));
    gnsf.L_z = full(L_z_fun(0));
end

gnsf.y = y;
gnsf.uhat = uhat;

gnsf.ny = length(y);
gnsf.nuhat = length(uhat);
gnsf.n_out = length(gnsf.phi_expr);

end

