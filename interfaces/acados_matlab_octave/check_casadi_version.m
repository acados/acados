%
% Copyright (c) The acados authors.
%
% This file is part of acados.
%
% Licensed under the 2-Clause BSD License.

%
function check_casadi_version()
    import casadi.*
    casadi_version = CasadiMeta.version();
    if ~(strcmp(casadi_version(1:3),'3.7'))
        if ~(strcmp(casadi_version(1:3),'3.4') || strcmp(casadi_version(1:3),'3.5') || strcmp(casadi_version(1:3),'3.6'))
            warning('Tested CasADi versions are 3.4, 3.5, 3.6, 3.7 you are using: %s.', casadi_version);
        else
            warning('Recommended CasADi version to be used in acados is 3.7 for a full featured experience. Versions 3.4 to 3.6 work for most problem formulations. You are using %s.', casadi_version);
        end
    end
end
