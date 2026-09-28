%
% Copyright (c) The acados authors.
%
% This file is part of acados.
%
% Licensed under the 2-Clause BSD License.

%

function simulink_opts = get_acados_simulink_opts()
    warning("Function get_acados_simulink_opts() is deprecated in acados v0.5.6, please instead use AcadosOcpSimulinkOptions().");
    simulink_opts = AcadosOcpSimulinkOptions();
end
