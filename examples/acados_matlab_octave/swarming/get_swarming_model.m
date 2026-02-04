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

% Author: Enrica


% NOTE: `acados` currently supports both an old MATLAB/Octave interface (< v0.4.0)
% as well as a new interface (>= v0.4.0).

% THIS EXAMPLE still uses the OLD interface. If you are new to `acados` please start
% with the examples that have been ported to the new interface already.
% see https://github.com/acados/acados/issues/1196#issuecomment-2311822122)



function model = get_swarming_model(S)

% SWARMING_MODEL - Function that describes the dynamics of the swarm and
% the cost function for controlling it.
%
% Swarming  is the behavior of collective and ordered motion. It can be
% obtained through the combination of the following rules:
%
% Separation rule: drives the agents to a reference inter-agent ...
%       distance (d_ref)
% Direction rule: make the agents' velocity converge to a ...
%       reference direction (u_ref)
% Navigation rule: make the agents' speed converge to a reference ...
%       value (v_ref)
%

import casadi.*

%% Rename swarming parameters

N = S.N; % number of agents in the swarm
d_ref = S.d_ref; % reference distance among every couple of neighboring agents
u_ref = S.u_ref; % reference direction of velocity for all agents
v_ref = S.v_ref; % reference speed for all agents

%% System dimensions

nx = 6 * N; % nb of state variables
nu = 3 * N; % nb of control inputs

%% Named symbolic variables

p = SX.sym('p', 3*N); % 3D positions of the agents [m]
v = SX.sym('v', 3*N); % 3D velocities of the agents [m/s]
u = SX.sym('a', 3*N); % 3D acceleration to apply to agents [m/s^2]

%% Unnamed symbolic variables

sym_x = vertcat(p, v);
sym_xdot = SX.sym('xdot', nx, 1);
sym_u = u;

%% Dynamics

expr_f_expl = vertcat(v, u);
expr_f_impl = expr_f_expl - sym_xdot;

%% Nonlinear least squares

% Weights
W_sep = 1;
W_dir = 1;
W_nav = 2;
W_u = 1e-1; % Penalization of high values of the control input variables

sym_sep = SX.zeros(N*(N-1),1);
sym_dir = SX.zeros(N,1);
sym_nav = SX.zeros(N,1);

%ny = N*(N+1);

% Neighborhood matrix
M = ones(N,N) - eye(N,N);

% For every agent define the nonlinear_ls terms
for agent = 1:N

    % Get the index triplet related to the current agent
    agent_idx = [1,2,3]' + 3*(agent-1)*ones(3,1);

    % For every neighbor, compute the distance to the current agent
    for neig = 1:(N-1)
        if neig < agent
            neig_idx = [1,2,3]' + 3*(neig-1)*ones(3,1);
        else
            neig_idx = [1,2,3]' + 3*(neig)*ones(3,1);
        end
        % Separation term
        p_rel = p(neig_idx)-p(agent_idx);
        sym_sep((agent-1)*(N-1)+neig) = 1/(N-1)*(p_rel'*p_rel - d_ref^2);
    end
    v_agent = v(agent_idx);
    % Direction term
    sym_dir(agent) = 1 - (v_agent'*u_ref)^2/(v_agent'*v_agent);
    % sym_dir(agent) = (v_agent - v_ref*u_ref)'*(v_agent - v_ref*u_ref);
    % Navigation term
    sym_nav(agent) = v_agent'*v_agent - v_ref^2;
end

sym_sep = W_sep * sym_sep;
sym_dir = W_dir * sym_dir;
sym_nav = W_nav * sym_nav;

% Assemble expr_y
expr_y = vertcat(sym_sep, sym_dir, sym_nav, W_u*sym_u);
expr_y_e = vertcat(sym_sep, sym_dir, sym_nav);
%% Populate structure
model = AcadosModel();
model.x = sym_x;
model.xdot = sym_xdot;
model.u = sym_u;
model.f_expl_expr = expr_f_expl;
model.f_impl_expr = expr_f_impl;
model.cost_y_expr = expr_y;
model.cost_y_expr_e = expr_y_e;
end
