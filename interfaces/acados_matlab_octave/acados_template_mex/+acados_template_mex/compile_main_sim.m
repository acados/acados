%
% Copyright (c) The acados authors.
%
% This file is part of acados.
%
% The 2-Clause BSD License
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
%
% 1. Redistributions of source code must retain the above copyright notice,
% this list of conditions and the following disclaimer.
%
% 2. Redistributions in binary form must reproduce the above copyright notice,
% this list of conditions and the following disclaimer in the documentation
% and/or other materials provided with the distribution.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.;

%

function compile_main_sim(export_dir)
    return_dir = pwd;
    cd(export_dir);
    %% build main file
    if isunix
        [ status, result ] = system('make sim_shared_lib');
        if status
            cd(return_dir);
            error('building templated code as shared library failed.\nGot status %d, result: %s',...
                  status, result);
        end
        fprintf('Successfully built main file!\n');
    else
        % compile on Windows platform
        [ status, result ] = system('cmake -G "MinGW Makefiles" -DCMAKE_BUILD_TYPE=Release -DBUILD_ACADOS_SIM_SOLVER_LIB=ON -S . -B .');
        if status
            cd(return_dir);
            error('Building templated code failed.\nGot status %d, result: %s',...
                  status, result);
        end
        [ status, result ] = system('cmake --build . --config Release');
        if status
            cd(return_dir);
            error('Building templated code as shared library failed.\nGot status %d, result: %s',...
                  status, result);
        end
        fprintf('Successfully built main file!\n');
    end

    cd(return_dir);
end
