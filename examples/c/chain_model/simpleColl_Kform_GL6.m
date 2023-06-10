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

function [G] = simpleColl_Kform_GL6(dae,tau_root,h)
  import casadi.*
  daefun = Function('fun',dae,char('x','p'),char('ode','quad'));
  % Degree of interpolating polynomial
  d = 3;

  AA(1,1) = 5.0/36.0;
  AA(1,2) = 2.0/9.0-1.0/15.0*sqrt(15.0);
  AA(1,3) = 5.0/36.0-1.0/30.0*sqrt(15.0);
  AA(2,1) = 5.0/36.0+1.0/24.0*sqrt(15.0);
  AA(2,2) = 2.0/9.0;
  AA(2,3) = 5.0/36.0-1.0/24.0*sqrt(15.0);
  AA(3,1) = 5.0/36.0+1.0/30.0*sqrt(15.0);
  AA(3,2) = 2.0/9.0+1.0/15.0*sqrt(15.0);
  AA(3,3) = 5.0/36.0;
  
  bb(1) = 5.0/18.0;
  bb(2) = 4.0/9.0;
  bb(3) = 5.0/18.0;
  
  % State variable
  X  = MX.sym('x',dae.x.size1(),1);
  
  % Helper state variables
  K = MX.sym('k',dae.x.size1(),d);
  
  % Fixed parameters (controls)
  CVp  = MX.sym('p',dae.p.size1());

  g = {};

  % For all collocation points
  for j=1:d
        
    % Get an expression for the state at the collocation point
    xp_jk = X;
    for r=1:d
      xp_jk = xp_jk + h*AA(j,r)*K(:,r);
    end
    % Add collocation equations to the NLP
    out = daefun({xp_jk,CVp});
    ode = out{1};
    g = {g{:} ode - K(:,j)};
  end
  % Get an expression for the state at the end of the finite element
  xf_k = X;
  for r=1:d
    xf_k = xf_k + h*bb(r)*K(:,r);
  end
  G = Function('G',{X,K,CVp},{xf_k,vertcat(g{:})});
  
end
