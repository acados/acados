%
% Copyright (c) The acados authors.
%
% This file is part of acados.
%
% Licensed under the 2-Clause BSD License.


function model = get_transition_model()
    import casadi.*
    model = AcadosModel();
    model.name = 'transition_model';
    % set up states
    p = SX.sym('p');
    v = SX.sym('v');
    model.x = vertcat(p, v);
    model.disc_dyn_expr = p;
end
