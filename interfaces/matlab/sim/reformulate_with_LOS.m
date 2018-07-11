function [ gnsf, reordered_model ] = ...
        reformulate_with_LOS( model, gnsf, print_info)
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

%% Description:
% This function takes an intitial transcription of the implicit ODE model
% "model" into "gnsf" and reformulates "gnsf" with a linear output system
% (LOS), containing as many states of the model as possible.
% Therefore it might be that the state vector and the implicit function
% vector have to be reordered. This reordered model is part of the output,
% namely reordered_model.

%% import CasADi and load models

import casadi.*

% % model
f_impl_expr = model.f_impl_expr;

x = model.x;
xdot = model.xdot;
u = model.u;
z = model.z;

% % GNSF
% get dimensions
nx  = gnsf.nx;
nz  = gnsf.nz;
nx1 = gnsf.nx1;

% get model matrices
A  = gnsf.A;
B  = gnsf.B;
C  = gnsf.C;
E  = gnsf.E;
c  = gnsf.c;

A_LO = gnsf.A_LO;

y = gnsf.y;
x1 = x(1:nx1);
x1dot = xdot(1:nx1);

phi_old = gnsf.phi_expr;

%% build initial I_x1 and I_x2_candidates
% I_x1: all components of x for which either xii or xdot_ii enters y;
% I_x2_candidates: the remaining components

I_x1 = [];
I_x2_candidates = [];

for ii = 1:nx1
    dy_dxii = jacobian(y,x1(ii));
    dy_dxiidot = jacobian(y,x1dot(ii));
    dy_dxii_sparsity = dy_dxii.sparsity;
    dy_dxiidot_sparsity = dy_dxiidot.sparsity;
    if or(dy_dxii_sparsity.nnz, dy_dxiidot_sparsity.nnz)
        % i.e. xii or xiidot are part of y, i.e. enter phi_expr
        if print_info
            disp(['xii is part of x1, ii = ', num2str(ii)]);
        end
        I_x1 = [I_x1, ii];
    else
        % i.e. neither xii nor xiidot are part of y, i.e. enter phi_expr
        I_x2_candidates = [I_x2_candidates, ii];
    end
end

if print_info
disp(' ');
disp('=================================================================');
disp('=== the algorithm will now try to detect linear output system ===');
disp('=================================================================');
disp(' ');
end

%% determine components of Linear Output System
% determine maximal index set I_x2
% such that the components x(I_x2) can be written as a LOS 
%% move successively entries from I_x2_candidates to I_x1, but only if needed..
% idea: assuming x_ii is part of LOS: then xiidot cannot occur in any other
%   equation
% (Simplifying Assumption (SA): for the reformulation we assume that the
%   derivatives of the states that belong to the LOS are given explicitly),
%   i.e. they are of the form:
%           E(i_eq, i_x2) * xdot(i_x2) = ...
%                           A * x + B * u + C * phi + c;
% also: x_ii can not occur in any eq, that defines x1dot, z;
%    thus we will iterate through the candidates for the LOS
%    (I_x2_candidates):
%       - set candidates = [ii];
%       *
%       - check the SA for candidates;
%           - if check failed: move ii to I_x1;
%           - else: check if candidates form a closed LOS, i.e.
%                  A(non_included_equ,candidates) == zero
%               - if closed: add candidates to LOS (I_x2);
%               - else: add new equations (the rows where nonzeros occured)
%                        and go to *
I_x2 = [];
for ii = I_x2_candidates % ensured: xii_dot does not enter phi_expr
    if print_info
    disp(['============= candidate ii = ', num2str(ii), '  ===================']);
    end
    candidates = [ii];
    while 1
        if any(any(candidates' == I_x1)) % some candidate is in I_x1 already
            I_x1 = unique([I_x1, ii]);
            if print_info
                disp(['x1 depends on x_ii, ii = ', num2str(ii)])
            end
            break;
        else
            check = check_SA(E, candidates);
            if check  % check_SA passed
                candidates_dynamics = ...
                    equations_defining_candidates( candidates, E);
                candidates_equs = ...
                    equations_candidates_enter_linearly( candidates, A);
                all_equ = unique( [ candidates_dynamics, candidates_equs]);
                
                if length(candidates_dynamics) == length(all_equ)
                    % candidates form a closed LOS
                    % -> candidates can be made part of LOS
                    I_x2 = unique([I_x2, candidates]);
                    if print_info
                    disp(['the following states form a LOS and are moved to x2']);
                    disp(candidates);
                    end
                    break;
                else
                    % add new candidates
                    candidates = xdot_components_of_equations(all_equ, E);
                end
            else % check_SA failed -> put ii into x1;
                I_x1 = unique([I_x1, ii]);
                if print_info
                    disp(['check_SA failed for x_ii, ii = ', num2str(ii)])
                    disp('move x_ii belongs to x1');
                end
                break;
            end
        end
    end
end

%% permute x, xdot
x1 = x(I_x1);
x2 = x(I_x2);
x1dot = xdot(I_x1);
x2dot = xdot(I_x2);

gnsf.xdot = [x1dot; x2dot];
gnsf.x = [x1; x2];

gnsf.nx1 = length(x1);
gnsf.nx2 = length(x2);

% define reordered_model
reordered_model = model;
reordered_model.x = gnsf.x;
reordered_model.xdot = gnsf.xdot;


%% Set up LOS given the components given in I_x2
% create equivalent equations in LOS
% transcribe the equation
%       E(i_eq, i_x2) * xdot(i_x2) + E(i_eq, I_x1) * xdot( I_x1 ) + ...
%       E(i_eq, I_z) * z = ...
%               A(i_eq,:) * x + B(i_eq,:) * u + C(i_eq,:) * phi + c(i_eq);
% as
%       xdot(i_x2) =
%           (A(i_eq, I_x1) * x1 + B(i_eq,:) * u + C(i_eq,:) * phi + c(i_eq)) ...
%                               / E(i_eq, i_x2) ...
%               + ( A(i_eq, I_x2) / E(i_eq, i_x2) ) * x2;
%
%            =: f_LO(x1, x1dot, u, z) + A_LO( i_LO, :) * x2;
f_LO = [];
I_z = nx1+1:nx1+nz;

equ_changed_sign = [];
for ii_x2 = I_x2
    i_eq = equations_defining_candidates(ii_x2, E);
    f_LO = vertcat(f_LO,...
        ( A(i_eq, I_x1) * x1 + B(i_eq, :) * u + C(i_eq, :) * gnsf.phi_expr + c(i_eq) ...
            -  E(i_eq, I_x1) * x1dot - E(i_eq, I_z) * z) / E(i_eq, ii_x2) );
    i_LO = find( I_x2 == ii_x2 );
    A_LO(i_LO, :) = A(i_eq, I_x2) / E(i_eq, ii_x2) ;
    E_entry = E(i_eq, ii_x2);
    if E_entry < 0
        f_impl_expr( i_eq ) = - f_impl_expr( i_eq );
        equ_changed_sign = [equ_changed_sign, i_eq];
    end
end


f_LO = f_LO.simplify();

% remove corresponding equations from NSF
I_eq_x2 = []; % this must be ordered as I_x2 is; thus for loop
for ii = I_x2
    I_eq_x2 = [I_eq_x2, equations_defining_candidates(ii, E)];
end

I_eq_x1 = 1:nx1+nz;
for i = I_eq_x2
    i_remove = find( I_eq_x1 == i );
    I_eq_x1 = I_eq_x1([1:i_remove-1, i_remove+1:length(I_eq_x1)]);
end

% permute f accordingly
f_permutation = [I_eq_x1, I_eq_x2];
f_impl_expr = f_impl_expr( f_permutation ) ;
f_impl_expr = f_impl_expr.simplify();

reordered_model.f_impl_expr = f_impl_expr;
reordered_model.equ_changed_sign = equ_changed_sign;

gnsf.A_LO = A_LO;
gnsf.f_lo_expr = f_LO;


%% reduce size of first system of GNSF (a)
old_Ix1 = I_x1;

gnsf.A = gnsf.A(I_eq_x1, I_x1);
gnsf.B = gnsf.B(I_eq_x1, :);
gnsf.C = gnsf.C(I_eq_x1, :);
gnsf.E = gnsf.E(I_eq_x1, [I_x1, nx+1 : nx+nz ]);
gnsf.c = gnsf.c(I_eq_x1, :);

C_new = [];
phi_new = [];
for ii = 1:size(gnsf.C, 2) % n_colums of C
    if all(gnsf.C(:,ii) == 0) % if column == 0
    else
        C_new = [C_new, gnsf.C(:,ii)];
        phi_new = [phi_new; gnsf.phi_expr(ii)];
    end
end

gnsf.C = C_new;

gnsf.phi_expr = phi_new;
gnsf.n_out = length(phi_new);

[ gnsf ] = determine_input_nonlinearity_function( gnsf );

check_reformulation(reordered_model, gnsf, print_info);

if print_info
disp('Successfully detected Linear Output System');
disp(['==>>  moved  ', num2str(gnsf.nx2), ' states to the Linear Output System']);
disp(['==>>  recuced output dimension of phi from  ', num2str(length(phi_old)), ' to ', num2str(length(gnsf.phi_expr))]);
end


end

%% auxilary functions

function [ check ] = check_SA(E, candidates)
% (Simplifying Assumption (SA): for the reformulation we assume that the
%   derivatives of the states that belong to the LOS are given as 

% const1 * xdot_ii = vec * xdot_1 + vec * z + A*x + vec * u + vec * f(y,u) 
%           +  const2 + E_vec * x1;

% xdot_ii does only occur in one equation 
% whereby xii and xdot_ii cannot be part of y;

% check is 1, if SA holds for all candidates,
% check is 0, if SA holds,
    for ii = candidates
        eq_xiidot = find(E(:,ii))';
        % contains the indices of equations that xii_dot enters (linearly);
        if length(eq_xiidot) > 1  % put xii into x1, due to (SA)
            check = 0;
            return;
        elseif isempty(eq_xiidot)
            error(['xiidot, ii = ',num2str(ii), ' does not occur in any equation']);
        else
            check = 1;
        end
    end
end

function I_eq = equations_candidates_enter_linearly( candidates, A )
% determine index set of components of x, for which either x_i or xdot_i
% enters at least one of the equations corresponding to I_equ linearly.

I_eq = [];
for ii = candidates
    cand_ii_eq = find(A(:,ii)');
    I_eq = [I_eq, cand_ii_eq];
end
I_eq = unique(I_eq); 

end

function I_equ = equations_defining_candidates( candidates, E)
% determine the index set of equations that the candidates or its
% derivatives enter linearly, i.e. through A, E (as f is excluded)
    I_equ = [];
    for ii = candidates
        % find indices of equations that xii enters
        eq_xiidot = find(E(:,ii))';
        I_equ = [I_equ, eq_xiidot];
    end
    I_equ = unique(I_equ);
end

function I = xdot_components_of_equations(I_eq, E);

components = [];
for ii = I_eq
    components = [components, find(E(ii,:))];
end
I = unique(components);
    

end

