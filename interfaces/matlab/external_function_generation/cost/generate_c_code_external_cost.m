function generate_c_code_external_cost(  model, cost, opts )
    %% import casadi
    import casadi.*

    if CasadiMeta.version()=='3.4.0'
        % casadi 3.4
        casadi_opts = struct('mex', false, 'casadi_int', 'int', 'casadi_real', 'double');
    else
        % old casadi versions
        error('Please download and install CasADi 3.4.0 to ensure compatibility with acados')
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
        { cost_jac', cost_hes });
    external_cost_fun.generate( [name, '_external_cost'], casadi_opts );

    
%     keyboard
%     external_cost_fun(zeros(size([u;x])))
end