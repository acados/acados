function [out] = casadi_vec(varargin)
  out = casadi_struct2vec(casadi_struct(varargin{:}));
end
