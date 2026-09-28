%
% Copyright (c) The acados authors.
%
% This file is part of acados.
%
% Licensed under the 2-Clause BSD License.


function solver = acados_ocp(model, opts, simulink_opts)

    warning('acados_ocp will be deprecated in the future. Use AcadosOcpSolver instead. For more information on the major acados MATLAB interface overhaul, see https://github.com/acados/acados/releases/tag/v0.4.0');

    if nargin < 3
        simulink_opts = AcadosOcpSimulinkOptions();
    end

    ocp = setup_AcadosOcp_from_legacy_ocp_description(model, opts, simulink_opts);
    solver = AcadosOcpSolver(ocp, struct('output_dir', opts.opts_struct.output_dir));

end
