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

function model  = detect_dims_sim(model, opts)

    %% general
    model.dim_nx = length(model.sym_x);

    if isfield(model, 'sym_u')
        model.dim_nu = length(model.sym_u);
    else
        model.dim_nu = 0;
    end

    if isfield(model, 'sym_z')
        model.dim_nz = length(model.sym_z);
    else
        model.dim_nz = 0;
    end

    if isfield(model, 'sym_p')
        model.dim_np = length(model.sym_p);
    else
        model.dim_np = 0;
    end

    if isfield(model, 'sym_xdot')
        if numel(model.sym_xdot) ~= model.dim_nx
            warning('sym_xdot is not of shape nx');
        end
    end

    if ~isempty(opts.num_stages)
        if(strcmp(opts.method,"erk"))
            if(opts.num_stages == 1 || opts.num_stages == 2 || ...
                opts.num_stages == 3 || opts.num_stages == 4)
            else
                error(['ERK: num_stages = ', num2str(opts.num_stages) ' not available. Only number of stages = {1,2,3,4} implemented!']);
            end
        end
    end
end
