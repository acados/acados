%
% Copyright (c) The acados authors.
%
% This file is part of acados.
%
% Licensed under the 2-Clause BSD License.

%


function generate_c_code_explicit_ode(context, model, model_dir)

    if isempty(model.f_expl_expr)
        error("Field `f_expl_expr` is required for integrator type ERK.")
    end

    add_explicit_ode_function_definitions(context, model, model_dir, model.x, model.f_expl_expr);
end

