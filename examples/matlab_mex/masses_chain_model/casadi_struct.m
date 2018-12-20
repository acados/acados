function [out] = casadi_struct(s,varargin)
  if ischar(varargin{1})
      default = 0;
  else
      default = varargin{1};
      varargin = varargin(2:end);
  end
  import casadi.*
  out = struct;
  for k=fieldnames(s)'
    k = k{1};
    found = -1;
    for l=1:length(varargin)/2
        if strcmp(varargin{2*l-1},k)
           found = l;
           break;
        end
    end
    if found>0
      e = varargin{2*found};
      if isscalar(e)
        dims = size(s.(k));
        e = repmat(e,dims(1),dims(2));
      end
      dims = size(s.(k));
      assert(size(e,1)==dims(1))
      assert(size(e,2)==dims(2))
    else
      dims = size(s.(k));
      e = default*DM.ones(dims(1),dims(2));
    end
    out.(k) = e;
  end
end