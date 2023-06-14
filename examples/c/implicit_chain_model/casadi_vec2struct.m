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

function [out] = casadi_vec2struct(s,vec)
  import casadi.*
  assert(isvector(vec))
  try
      vec.sparsity();
  catch
      vec = DM(vec);
  end
  flat = {};
  if isstruct(s)
    out = struct;
    sizes = {0};
    for f=fieldnames(s)'
      dim = size(casadi_struct2vec(s.(f{1})));
      sizes = {sizes{:} sizes{end}+dim(1)};
    end
    comps = vertsplit(vec,sizes);
    i = 1;
    for f=fieldnames(s)'
      out.(f{1}) = casadi_vec2struct(s.(f{1}),comps{i});
      i = i+1;
    end
  else if iscell(s)
    out = cell(size(s));
    sizes = {0};
    for i=1:length(s)
      n = size(casadi_struct2vec(s{i}),1);
      sizes = {sizes{:} sizes{end}+n};
    end
    comps = vertsplit(vec,sizes);
    for i=1:length(s)
      out{i} = casadi_vec2struct(s{i},comps{i});
    end
      else
    out = reshape(vec,size(s));
  end
end