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

function [G] = simpleColl(dae,tau_root,h)
  import casadi.*
  daefun = Function('fun',dae,char('x','p'),char('ode','quad'));
  % Degree of interpolating polynomial
  d = length(tau_root)-1;

  % Coefficients of the collocation equation
  C = zeros(d+1,d+1);

  % Coefficients of the continuity equation
  D = zeros(d+1,1);

  % Dimensionless time inside one control interval
  tau = SX.sym('tau');

  % For all collocation points
  for j=1:d+1
    % Construct Lagrange polynomials to get the polynomial basis at the collocation point
    L = 1;
    for r=1:d+1
      if r ~= j
        L = L * (tau-tau_root{r})/(tau_root{j}-tau_root{r});
      end
    end
    lfcn = Function('lfcn', {tau},{L});
    out = lfcn({1.0});
    
    % Evaluate the polynomial at the final time to get the coefficients of the continuity equation
    D(j) = full(out{1});

    % Evaluate the time derivative of the polynomial at all collocation points to get the coefficients of the continuity equation
    tfcn = lfcn.tangent();
    for r=1:d+1
      out = tfcn({tau_root{r}});
      C(j,r) = full(out{1});
    end
  end
  
  % State variable
  CVx  = MX.sym('x',dae.x.size1(),1);
  
  % Helper state variables
  CVCx = MX.sym('x',dae.x.size1(),d);
  
  % Fixed parameters (controls)
  CVp  = MX.sym('p',dae.p.size1());

  X = [CVx CVCx];
  g = {};

  % For all collocation points
  for j=2:d+1
        
    % Get an expression for the state derivative at the collocation point
    xp_jk = 0;
    for r=1:d+1
      xp_jk = xp_jk + C(r,j)*X(:,r);
    end
    % Add collocation equations to the NLP
    out = daefun({CVCx(:,j-1),CVp});
    ode = out{1};
    g = {g{:} h*ode - xp_jk};
  end
  % Get an expression for the state at the end of the finite element
  xf_k = 0;
  for r=1:d+1
    xf_k = xf_k + D(r)*X(:,r);
  end
  G = Function('G',{CVx,CVCx,CVp},{xf_k,vertcat(g{:})});
  
end
