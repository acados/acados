function dx = dyn( x, u )

% constants
M = 1;
m = 0.1;
g = 9.81;
l = 0.8;

x1 = x(1);
theta = x(2);
v1 = x(3);
dtheta = x(4);

F = u(1);

dx = [  v1; ...
        dtheta; ...
        (- l*m*sin(theta)*dtheta^2 + F + g*m*cos(theta)*sin(theta))/(M + m - m*cos(theta)^2); ...
        (- l*m*cos(theta)*sin(theta)*dtheta^2 + F*cos(theta) + g*m*sin(theta) + M*g*sin(theta))/(l*(M + m - m*cos(theta)^2)) ];

end

