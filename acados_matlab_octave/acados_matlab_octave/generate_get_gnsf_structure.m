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

function generate_get_gnsf_structure(model, opts)

output_dir = opts.output_dir;

fileID = fopen(fullfile(output_dir, 'get_gnsf_structure.m'), 'w');

fprintf(fileID, 'function model = get_gnsf_structure(model)\n');
fprintf(fileID, 'model.dim_gnsf_nx1 = %d;\n', model.dim_gnsf_nx1);
%fprintf(fileID, 'model.dim_gnsf_nx2 = %d;\n', model.dim_gnsf_nx2);
fprintf(fileID, 'model.dim_gnsf_nz1 = %d;\n', model.dim_gnsf_nz1);
%fprintf(fileID, 'model.dim_gnsf_nz2 = %d;\n', model.dim_gnsf_nz2);
fprintf(fileID, 'model.dim_gnsf_nuhat = %d;\n', model.dim_gnsf_nuhat);
fprintf(fileID, 'model.dim_gnsf_ny = %d;\n', model.dim_gnsf_ny);
fprintf(fileID, 'model.dim_gnsf_nout = %d;\n', model.dim_gnsf_nout);

fclose(fileID);
