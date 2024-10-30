try
    import casadi.*
    x = SX.sym('x')
    disp('casadi is available');
    jac_xsquare = jacobian(x^2, x)
catch e
    disp('casadi is not available');
end
