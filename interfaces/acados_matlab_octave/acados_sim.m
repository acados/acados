%
% Copyright (c) The acados authors.
%
% This file is part of acados.
%
% Licensed under the 2-Clause BSD License.


function solver = acados_sim(model, opts)

    warning('acados_sim will be deprecated in the future and is not tested anymore. Use AcadosSimSolver instead. For more information on the major acados MATLAB interface overhaul, see https://github.com/acados/acados/releases/tag/v0.4.0');

    sim = setup_AcadosSim_from_legacy_sim_description(model, opts);
    solver = AcadosSimSolver(sim, struct('output_dir', opts.opts_struct.output_dir));
end