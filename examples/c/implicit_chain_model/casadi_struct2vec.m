%
% Copyright (c) The acados authors.
%
% This file is part of acados.
%
% Licensed under the 2-Clause BSD License.

%

function [out] = casadi_struct2vec(s)
  flat = {};
  if isstruct(s)
    for f=fieldnames(s)'
      flat = {flat{:} casadi_struct2vec(s.(f{1}))};
    end
    out = vertcat(flat{:});
  else if iscell(s)
    for i=1:length(s)
       flat = {flat{:} casadi_struct2vec(s{i})};
    end
    out = vertcat(flat{:});
  else
    try
      out = vec(s);
    catch
      import casadi.*
      out = vec(DM(s));
    end
  end
end