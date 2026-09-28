%
% Copyright (c) The acados authors.
%
% This file is part of acados.
%
% Licensed under the 2-Clause BSD License.


function solver = acados_sim(model, opts)

    sim = setup_AcadosSim_from_legacy_sim_description(model, opts);
    solver = AcadosSimSolver(sim, struct('output_dir', opts.opts_struct.output_dir));
    % warning('In acados v0.4.0, many changes to the MATLAB/Octave interface of acados have been introduced.', ...
    % 'We recommend directly using the new AcadosSimSolver and to check the examples for the intended use.')
end