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