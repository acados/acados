function [ x_next ] = dynamics( x,u,h )
% Discrete time dynamics obtained by applying the RK4 scheme
k1      = c_dynamics(x,u);
k2      = c_dynamics(x + h/2*k1,u);
k3      = c_dynamics(x + h/2*k2,u);
k4      = c_dynamics(x + h*k3,u);
x_next  = x + h/6*(k1 + 2*k2 + 2*k3 + k4);
end

