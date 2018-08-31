clear all;

close all;

clc

% restoredefaultpath

% perpare the MATLAB path to include CASADI

% tmp = path();

% if(isempty(strfind(strrep(tmp,'\','/'),'D:\GitLab\casadi')))
% 
%     addpath(path,'D:\GitLab\casadi')
% 
% end;

import casadi.* 

x1 = MX.sym('x1');
x2 = MX.sym('x2');


x = [x1; x2];

f = x(1)^2;
 

jac_phi_x1 = jacobian(f,x1) % funktioniert

jac_phi_x = jacobian(f,x) % funktioniert