%
% Copyright (c) The acados authors.
%
% This file is part of acados.
%
% Licensed under the 2-Clause BSD License.

%

function [out] = casadi_vec(varargin)
  out = casadi_struct2vec(casadi_struct(varargin{:}));
end
