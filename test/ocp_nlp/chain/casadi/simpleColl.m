function [G] = simpleColl(dae,tau_root,h)
  import casadi.*
  daefun = Function('fun',dae,char('x','p'),char('ode','quad'));
  % Degree of interpolating polynomial
  tau_root2 = [0 tau_root];
  d = length(tau_root2)-1;

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
        L = L * (tau-tau_root2(r))/(tau_root2(j)-tau_root2(r));
      end
    end
    lfcn = Function('lfcn', {tau},{L});
    out = lfcn(1.0);
    
    % Evaluate the polynomial at the final time to get the coefficients of the continuity equation
    D(j) = full(out);

    % Evaluate the time derivative of the polynomial at all collocation points to get the coefficients of the continuity equation
    tfcn = Function('tfcn', {tau}, {tangent(L,tau)});
    for r=1:d+1
      out = tfcn(tau_root2(r));
      C(j,r) = full(out);
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
    ode = daefun(CVCx(:,j-1),CVp);
    g = {g{:} h*ode - xp_jk};
  end
  % Get an expression for the state at the end of the finite element
  xf_k = 0;
  for r=1:d+1
    xf_k = xf_k + D(r)*X(:,r);
  end
  G = Function('G',{CVx,CVCx,CVp},{xf_k,vertcat(g{:})});
  
end
