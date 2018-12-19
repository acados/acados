%% create external cost function for nonlinear chain model
% import casadi
clearvars

import casadi.*
addpath('../../../interfaces/matlab/external_function_generation');

for nm = 2:6
    source = strcat('xN_nm', num2str(nm), '.txt');
    ref = load(source);
    nx = 6 * (nm-1);
    nu = 3;
    for i = length(ref)+1:nx+nu
        ref(i) = 0;
    end

    model.x = SX.sym('x',nx,1);
    model.u = SX.sym('u',nu,1);
    cost.name = strcat('chain_nm_', num2str(nm));
    
    ux = [model.x; model.u];

    % define least square cost as general expression
    cost.general_expr = SX.zeros(1,1);
    for j = 1:nu
        cost.general_expr = cost.general_expr + 0.5 * (ux(j) - ref(j))^2;
    end
    for j = 1:nx
        cost.general_expr = cost.general_expr + 0.5 * 1e-2 * (ux(j+nu) - ref(j+nu))^2;
    end
    
    generate_c_code_external_cost(model, cost)
    
end
    