% generate_c_code_nls
function generate_c_code_nls_cost( model, cost, opts )

    %% import casadi
    import casadi.*

    if CasadiMeta.version()=='3.4.0'
        % casadi 3.4
        casadi_opts = struct('mex', false, 'casadi_int', 'int', 'casadi_real', 'double');
    else
        % old casadi versions
        error('Please download and install Casadi 3.4.0 to ensure compatibility with acados')
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

    nls_cost = cost.nls_expr;


    nls_cost_jac = jacobian(nls_cost, [u; x]);

    
    nls_cost_fun = Function( [name, '_nls_cost_fun_jac'], {x, u}, ...
        { nls_cost, nls_cost_jac });
    nls_cost_fun.generate( [name, '_nls_cost_fun_jac'], casadi_opts );

    

    % % Trivial NLS (=LS) cost function
    % cost_x = SX.sym('x',nx);
    % cost_u = SX.sym('u',nu);
    % cost_ux = [cost_u; cost_x];
    % cost_y = [cost_x; cost_u];
    % cost_jac_y = jacobian(cost_y, cost_ux);
    % 
    % ls_cost = Function(['ls_cost_nm' num2str(Nm)], {[cost_u; cost_x]}, {cost_y, cost_jac_y'});
    % ls_cost.generate(['ls_cost_nm' num2str(Nm)], opts);

end