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


function generate_c_code_explicit_ode( model )

%% import casadi
import casadi.*

    casadi_version = CasadiMeta.version();
    if strcmp(casadi_version(1:3),'3.4') % require casadi 3.4.x
        casadi_opts = struct('mex', false, 'casadi_int', 'int', 'casadi_real', 'double');
    else % old casadi versions
        error('Please download and install CasADi version 3.4.x to ensure compatibility with acados')
    end

%% load model
x = model.x;
u = model.u;
if class(x(1)) == 'casadi.SX'
    isSX = true;
else
    isSX = false;
end

f_expl = model.f_expl_expr;
model_name = model.name;

%% get model dimensions
nx = length(x);
nu = length(u);

%% set up functions to be exported
if isSX
    Sx = SX.sym('Sx',nx,nx);
    Sp = SX.sym('Sp',nx,nu);
    lambdaX = SX.sym('lambdaX',nx,1);
    vdeX = SX.zeros(nx,nx);
    vdeP = SX.zeros(nx,nu) + jacobian(f_expl,u);
else
    Sx = MX.sym('Sx',nx,nx);
    Sp = MX.sym('Sp',nx,nu);
    lambdaX = MX.sym('lambdaX',nx,1);
    vdeX = MX.zeros(nx,nx);
    vdeP = MX.zeros(nx,nu) + jacobian(f_expl,u);
end
expl_ode_fun = Function([model_name,'_expl_ode_fun'],{x,u},{f_expl});
% TODO: Polish: get rid of SX.zeros

vdeX = vdeX + jtimes(f_expl,x,Sx);

vdeP = vdeP + jtimes(f_expl,x,Sp);

expl_vde_forw = Function([model_name,'_expl_vde_forw'],{x,Sx,Sp,u},{f_expl,vdeX,vdeP});

adj = jtimes(f_expl,[x;u],lambdaX,true);

expl_vde_adj = Function([model_name,'_expl_vde_adj'],{x,lambdaX,u},{adj});

S_forw = vertcat(horzcat(Sx, Sp), horzcat(zeros(nu,nx), eye(nu)));
hess = mtimes(S_forw.',jtimes(adj,[x;u],S_forw));
hess2 = [];
for j = 1:nx+nu
    for i = j:nx+nu
        hess2 = [hess2; hess(i,j)];
    end
end

expl_ode_hess = Function([model_name,'_expl_ode_hess'],{x,Sx,Sp,lambdaX,u},{adj,hess2});
%% generate C code
    if ~isfolder('c_generated_code')
        mkdir('c_generated_code');
    end

    cd('c_generated_code');
    model_dir = [model_name, '_model'];
    
    if ~isfolder(model_dir)
        mkdir(model_dir);
    end
    
    model_dir_location = ['./', model_dir];
    cd(model_dir_location);
    fun_name = [model_name, '_expl_ode_fun'];
    expl_ode_fun.generate(fun_name, casadi_opts);

    fun_name = [model_name, '_expl_vde_forw'];
    expl_vde_forw.generate(fun_name, casadi_opts);
    
    fun_name = [model_name, '_expl_vde_adj'];
    expl_vde_adj.generate(fun_name, casadi_opts);
    
    fun_name = [model_name, '_expl_ode_hess'];
    expl_ode_hess.generate(fun_name, casadi_opts);
    cd('../..')

end
