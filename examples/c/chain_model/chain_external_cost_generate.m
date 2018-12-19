%% create external cost function for nonlinear chain model
% import casadi
clearvars

import casadi.*
addpath('../../../interfaces/matlab/external_function_generation');

for nm = 2:6
    source = strcat('xN_nm', num2str(nm), '.txt');
    xref = load(source);
    nx = 6 * (nm-1);
    nu = 3;
    uxref = zeros(nu + nx, 1);
    for i = 1:nx
        uxref(i+nu) = xref(i);
    end

    model.x = SX.sym('x',nx,1);
    model.u = SX.sym('u',nu,1);
    cost.name = strcat('chain_nm_', num2str(nm));
    
    ux = [model.u; model.x];

    % define least square cost as general expression
    cost.general_expr = SX.zeros(1,1);
    for j = 1:nu
        cost.general_expr = cost.general_expr + 0.5 * (ux(j) - uxref(j))^2;
    end
    for j = 1:nx
        cost.general_expr = cost.general_expr + 0.5 * 1e-2 * (ux(j+nu) - uxref(j+nu))^2;
    end
    
    generate_c_code_external_cost(model, cost)
    
end
    