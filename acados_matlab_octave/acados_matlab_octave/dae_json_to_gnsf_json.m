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

function dae_json_to_gnsf_json(json_filename)

import casadi.*

casadi_version = CasadiMeta.version();
acados_folder = getenv('ACADOS_INSTALL_DIR');
addpath(fullfile(acados_folder, 'external', 'jsonlab'))

loaded_struct = loadjson(fileread(json_filename));

if ~strcmp(casadi_version, loaded_struct.casadi_version)
    error(['Current casadi_version ', casadi_version, ' does not match' ...
        ' casadi_version of dumped dae model, ', loaded_struct.casadi_version]);
end

%% set up model
impl_dae_fun = Function.deserialize(loaded_struct.str_impl_dae_fun);
% created as:
% impl_dae_fun = Function(fun_name, [x, xdot, u, z, p], [f_impl])

model_name = strrep(impl_dae_fun.name, '_impl_dae_fun', '');

size_x = impl_dae_fun.size_in(0);
nx = size_x(1);

size_xdot = impl_dae_fun.size_in(1);
nxdot = size_xdot(1);

if nx ~= nxdot
    error('nx != nxdot loading impl_dae_fun')
end

size_u = impl_dae_fun.size_in(2);
nu = size_u(1);

size_z = impl_dae_fun.size_in(3);
nz = size_z(1);

size_p = impl_dae_fun.size_in(4);
np = size_p(1);

model.name = model_name;
x = SX.sym('x', nx, 1);
xdot = SX.sym('xdot', nx, 1);
u = SX.sym('u', nu, 1);
z = SX.sym('z', nz, 1);
p = SX.sym('p', np, 1);

model.sym_x = x;
model.sym_xdot = xdot;
model.sym_u = u;
model.sym_z = z;
model.sym_p = p;
model.dyn_expr_f = impl_dae_fun(model.sym_x, model.sym_xdot,...
                                model.sym_u, model.sym_z, model.sym_p);

%%
model = detect_gnsf_structure(model);

%%
dump_gnsf_functions(model)

