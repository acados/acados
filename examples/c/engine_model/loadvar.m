% Simulink initialization

x0 = [50, 50, 1.14275, 1.53787];

load reference.mat;
reference = [(0:0.05:30).', reference];
