% generate_c_code_nls
function generate_c_code_nonlinear_constraints( model, constraints, opts )
%% TODO: test this

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
        name = [constraints.name, '_final'];
    else  %% INTERMEDIATE stage
        x = model.x;
        u = model.u;
        name = constraints.name;
    end

    nonl_constr = constraints.constraint_expr;


    nonl_constr_jac_tran = jacobian(nonl_constr, [u; x])';

    
    nonl_constr_fun_jac_tran = Function( [name, '_nonl_constr_fun_jac_tran'], {x, u}, ...
        { nonl_constr, nonl_constr_jac_tran });
    nonl_constr_fun_jac_tran.generate( [name, '_nonl_constr_fun_jac_tran'], casadi_opts );

end