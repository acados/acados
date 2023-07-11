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

clc;
close all;

%addpath('../../../external/casadi-octave-v3.2.3')
import casadi.*


% casadi opts for code generation
if CasadiMeta.version()=='3.4.0'
	% casadi 3.4
	opts = struct('mex', false, 'casadi_int', 'int', 'casadi_real', 'double');
else
	% old casadi versions
	error('Please download and install Casadi 3.4.0')
end


for Nm = 2:10
    
    disp(['---- Nm value = ' num2str(Nm) '----']);
    
    % Environment
    g = 9.81;     % [N/kg]
    L = 0.033;
    D = 1.0;
    m = 0.03;
    x0 = zeros(3,1);
    xN = [1 0 0].';
    
    wall_pos = -0.01;
    
    T = 5.0;
    N = 20;
    
    % Number of variables
    nx = (Nm-1)*2*3;
    nu = 3;
    
    % State variables
    u = SX.sym('u',3);
    dae.p = u;
    
    dae.x = [];
    states = [];
    for i = 1:Nm-1
        p = SX.sym(['p' num2str(i)],3);
        v = SX.sym(['v' num2str(i)],3);
        
        x_struct = struct('p',p,'v',v);
        states = [states; x_struct];
        dae.x = [dae.x; casadi_struct2vec(x_struct)];
    end
    
    % Compute forces
    F = {};
    for i = 1:Nm-1
        if i == 1
            dist = states(1).p-x0;
        else
            dist = states(i).p-states(i-1).p;
        end
        tmp = D*(1 - L/sqrt(dist.'*dist));
        F = {F{:}, tmp*dist};
    end
    
    % Set up ODE
    dae.ode = [];
    for i = 1:Nm-2
        f = 1/m*(F{i+1} - F{i}) - [0;0;g];
        dae.ode = [dae.ode; casadi_vec(x_struct,'p',states(i).v,'v',f)];
    end
    dae.ode = [dae.ode; casadi_vec(x_struct,'p',states(end).v,'v',u)];
        
    dae.x_dot = SX.sym('x_dot',nx,1);
    
    ode_impl = dae.x_dot - dae.ode;
    jac_x = SX.zeros(nx,nx) + jacobian(ode_impl, dae.x);
    jac_xdot = SX.zeros(nx,nx) + jacobian(ode_impl, dae.x_dot);
    jac_u = SX.zeros(nx,nu) + jacobian(ode_impl, dae.p);
    
    impl_ode_fun = Function(['casadi_impl_ode_fun_chain_nm' num2str(Nm)], {dae.x, dae.x_dot, dae.p}, {ode_impl});
	impl_ode_fun_jac_x_xdot = Function(['casadi_impl_ode_fun_jac_x_xdot_chain_nm',num2str(Nm)], {dae.x, dae.x_dot, dae.p}, {ode_impl, jac_x, jac_xdot});
    impl_ode_fun_jac_x_xdot_u = Function(['casadi_impl_ode_fun_jac_x_xdot_u_chain_nm',num2str(Nm)], {dae.x, dae.x_dot, dae.p}, {ode_impl, jac_x, jac_xdot, jac_u});
	impl_ode_jac_x_xdot_u = Function(['casadi_impl_ode_jac_x_xdot_u_chain_nm',num2str(Nm)], {dae.x, dae.x_dot, dae.p}, {jac_x, jac_xdot, jac_u});
    
    impl_ode_fun.generate(['impl_ode_fun_chain_nm' num2str(Nm)], opts);
	impl_ode_fun_jac_x_xdot.generate(['impl_ode_fun_jac_x_xdot_chain_nm',num2str(Nm)], opts);
    impl_ode_fun_jac_x_xdot_u.generate(['impl_ode_fun_jac_x_xdot_u_chain_nm',num2str(Nm)], opts);
	impl_ode_jac_x_xdot_u.generate(['impl_ode_jac_x_xdot_u_chain_nm',num2str(Nm)], opts);

end
