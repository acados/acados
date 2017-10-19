function dx = dyn( x,u )
rho = 1.23;
A = 0.1;
Cl = 0.25;
Cd = 0.3*Cl;
m = 10;
g = 9.81;
L = 0.5;
L2 = 1; 
J1 = 0.25;
J2 = 0.25;
J3 = 1;

q1 = x(1);
q2 = x(2);
q3 = x(3);
q4 = x(4);

Omega1 = x(5);
Omega2 = x(6);
Omega3 = x(7);

W1 = x(8);
W2 = x(9);
W3 = x(10);
W4 = x(11);

alpha = 0.0;
rW1 = u(1);
rW2 = u(2);
rW3 = u(3);
rW4 = u(4);

dx = ...
    [  
    - (Omega1*q2)/2 - (Omega2*q3)/2 - (Omega3*q4)/2 - (alpha*q1*(q1*q1 + q2*q2 + q3*q3 + q4*q4 - 1))/(q1*q1 + q2*q2 + q3*q3 + q4*q4);
    (Omega1*q1)/2 - (Omega3*q3)/2 + (Omega2*q4)/2 - (alpha*q2*(q1*q1 + q2*q2 + q3*q3 + q4*q4 - 1))/(q1*q1 + q2*q2 + q3*q3 + q4*q4);
    (Omega2*q1)/2 + (Omega3*q2)/2 - (Omega1*q4)/2 - (alpha*q3*(q1*q1 + q2*q2 + q3*q3 + q4*q4 - 1))/(q1*q1 + q2*q2 + q3*q3 + q4*q4);
    (Omega3*q1)/2 - (Omega2*q2)/2 + (Omega1*q3)/2 - (alpha*q4*(q1*q1 + q2*q2 + q3*q3 + q4*q4 - 1))/(q1*q1 + q2*q2 + q3*q3 + q4*q4);
    (J3*Omega2*Omega3 - J2*Omega2*Omega3 + (A*Cl*L*rho*(W2*W2 - W4*W4))/2)/J1;
    -(J3*Omega1*Omega3 - J1*Omega1*Omega3 + (A*Cl*L*rho*(W1*W1 - W3*W3))/2)/J2;
    (J2*Omega1*Omega2 - J1*Omega1*Omega2 + (A*Cd*L2*rho*(W1*W1 - W2*W2 + W3*W3 - W4*W4))/2)/J3;
    rW1;
    rW2;
    rW3;
    rW4];
end

