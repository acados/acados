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
