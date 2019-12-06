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

function check = sim_check_dims(model)

check = 1;

if isfield(model, 'sym_x')
    if all(size(model.sym_x))
        nx = length(model.sym_x);
    else
        nx = 0;
    end
    if nx ~= model.dim_nx
        check = 0;
        fail = 'x';
    end
end

if isfield(model, 'sym_u')
    if all(size(model.sym_u))
        nu = length(model.sym_u);
    else
        nu = 0;
    end
    if nu ~= model.dim_nu
        check = 0;
        fail = 'u';
    end
end


if isfield(model, 'sym_p')
    if all(size(model.sym_p))
        np = length(model.sym_p);
    else
        np = 0;
    end
    if np ~= model.dim_np
        check = 0;
        fail = 'p';
    end
end


if isfield(model, 'sym_z')
    if all(size(model.sym_z))
        nz = length(model.sym_z);
    else
        nz = 0;
    end
    if nz ~= model.dim_nz
        check = 0;
        fail = 'z';
    end
end


if check == 0
    message = strcat('\nSIM_DIM_CHECK FAIL: check consistency of dim_',...
        fail, ' with CasADi symbolic sym_', fail, '!\n\n');
    fprintf(message);
    error(message);
end
