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
