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

function [G] = simpleColl_Kform_GL8(dae,tau_root,h)
  addpath('../../external/casadi-octave-v3.2.2')
  import casadi.*
  daefun = Function('fun',dae,char('x','p'),char('ode','quad'));
  % Degree of interpolating polynomial
  d = 4;
  
  AA(1,1) = (1/144)*sqrt(30)+1/8;
  AA(1,2) = -(1/840)*sqrt(525-70*sqrt(30))*sqrt(30)+(1/144)*sqrt(30)-(1/105)*sqrt(525-70*sqrt(30))+1/8;
  AA(1,3) = (1/2352)*sqrt(525+70*sqrt(30))*sqrt(30)-(1/144)*sqrt(30)+(1/1680)*sqrt(525-70*sqrt(30))*sqrt(30)+1/8+(1/1470)*sqrt(525+70*sqrt(30))-(1/420)*sqrt(525-70*sqrt(30));
  AA(1,4) = -(1/2352)*sqrt(525+70*sqrt(30))*sqrt(30)-(1/144)*sqrt(30)+(1/1680)*sqrt(525-70*sqrt(30))*sqrt(30)+1/8-(1/1470)*sqrt(525+70*sqrt(30))-(1/420)*sqrt(525-70*sqrt(30));
  AA(2,1) = (1/840)*sqrt(525-70*sqrt(30))*sqrt(30)+(1/144)*sqrt(30)+(1/105)*sqrt(525-70*sqrt(30))+1/8;
  AA(2,2) = (1/144)*sqrt(30)+1/8;
  AA(2,3) = (1/2352)*sqrt(525+70*sqrt(30))*sqrt(30)-(1/144)*sqrt(30)-(1/1680)*sqrt(525-70*sqrt(30))*sqrt(30)+1/8+(1/1470)*sqrt(525+70*sqrt(30))+(1/420)*sqrt(525-70*sqrt(30));
  AA(2,4) = -(1/2352)*sqrt(525+70*sqrt(30))*sqrt(30)-(1/144)*sqrt(30)-(1/1680)*sqrt(525-70*sqrt(30))*sqrt(30)+1/8-(1/1470)*sqrt(525+70*sqrt(30))+(1/420)*sqrt(525-70*sqrt(30));
  AA(3,1) = -(1/2352)*sqrt(525-70*sqrt(30))*sqrt(30)+(1/144)*sqrt(30)-(1/1680)*sqrt(525+70*sqrt(30))*sqrt(30)+(1/1470)*sqrt(525-70*sqrt(30))+1/8-(1/420)*sqrt(525+70*sqrt(30));
  AA(3,2) = (1/2352)*sqrt(525-70*sqrt(30))*sqrt(30)+(1/144)*sqrt(30)-(1/1680)*sqrt(525+70*sqrt(30))*sqrt(30)-(1/1470)*sqrt(525-70*sqrt(30))+1/8-(1/420)*sqrt(525+70*sqrt(30));
  AA(3,3) = -(1/144)*sqrt(30)+1/8;
  AA(3,4) = (1/840)*sqrt(525+70*sqrt(30))*sqrt(30)-(1/144)*sqrt(30)-(1/105)*sqrt(525+70*sqrt(30))+1/8;
  AA(4,1) = -(1/2352)*sqrt(525-70*sqrt(30))*sqrt(30)+(1/144)*sqrt(30)+(1/1680)*sqrt(525+70*sqrt(30))*sqrt(30)+(1/1470)*sqrt(525-70*sqrt(30))+1/8+(1/420)*sqrt(525+70*sqrt(30));
  AA(4,2) = (1/2352)*sqrt(525-70*sqrt(30))*sqrt(30)+(1/144)*sqrt(30)+(1/1680)*sqrt(525+70*sqrt(30))*sqrt(30)-(1/1470)*sqrt(525-70*sqrt(30))+1/8+(1/420)*sqrt(525+70*sqrt(30));
  AA(4,3) = -(1/840)*sqrt(525+70*sqrt(30))*sqrt(30)-(1/144)*sqrt(30)+(1/105)*sqrt(525+70*sqrt(30))+1/8;
  AA(4,4) = -(1/144)*sqrt(30)+1/8;
  
  bb(1) = (1/72)*sqrt(30)+1/4;
  bb(2) = (1/72)*sqrt(30)+1/4;
  bb(3) = -(1/72)*sqrt(30)+1/4;
  bb(4) = -(1/72)*sqrt(30)+1/4;
  
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
    ode = daefun(xp_jk,CVp);
    g = {g{:} ode - K(:,j)};
  end
  % Get an expression for the state at the end of the finite element
  xf_k = X;
  for r=1:d
    xf_k = xf_k + h*bb(r)*K(:,r);
  end
  G = Function('G',{X,K,CVp},{xf_k,vertcat(g{:})});
  
end
