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