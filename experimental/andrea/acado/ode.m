function dx = ode(t,x,u)

    x1 = x(1);
    x2 = x(2);
    u = u(1);
    
    dx =   [    x2+u*(0.5+0.5*x1); ...
                x1+u*(0.5-2*x2)];

end

