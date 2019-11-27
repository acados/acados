function generate_c_code_external_cost(  model, cost, opts )
    %% import casadi
    import casadi.*

    casadi_version = CasadiMeta.version();
    if ( strcmp(casadi_version(1:3),'3.4') || strcmp(casadi_version(1:3),'3.5')) % require casadi 3.4.x
        casadi_opts = struct('mex', false, 'casadi_int', 'int', 'casadi_real', 'double');
    else % old casadi versions
        error('Please download and install CasADi version 3.4.x to ensure compatibility with acados')
    end


    % load cost, model
    if nargin > 2 && opts.is_terminal == 1 %% TERMINAL stage
        x = model.x;
        u = [];
        name = [cost.name, '_final'];
    else  %% INTERMEDIATE stage
        x = model.x;
        u = model.u;
        name = cost.name;
    end

    external_cost = cost.general_expr;

    cost_jac = jacobian(external_cost, [u; x]);
    cost_hes = jacobian( cost_jac, [u;x]);
    
    external_cost_fun = Function( [name, '_external_cost'], {x, u}, ...
        { external_cost, cost_jac', cost_hes });
    external_cost_fun.generate( [name, '_external_cost'], casadi_opts );

    
%     keyboard
%     external_cost_fun(zeros(size([u;x])))
end
