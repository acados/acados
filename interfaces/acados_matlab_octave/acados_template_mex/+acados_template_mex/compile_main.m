%
% Copyright 2019 Gianluca Frison, Dimitris Kouzoupis, Robin Verschueren,
% Andrea Zanelli, Niels van Duijkeren, Jonathan Frey, Tommaso Sartor,
% Branimir Novoselnik, Rien Quirynen, Rezart Qelibari, Dang Doan,
% Jonas Koenemann, Yutao Chen, Tobias Schöls, Jonas Schlagenhauf, Moritz Diehl
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

function compile_main()
    return_dir = pwd;
    cd c_generated_code
    %% build main file
    if isunix || ismac 
        % compile if on Mac or Unix platform
        [ status, result ] = system('make');
        if status
            cd(return_dir);
            error('building templated code failed.\nGot status %d, result: %s',...
                  status, result);
        end
        [ status, result ] = system('make shared_lib');
        if status
            cd(return_dir);
            error('building templated code as shared library failed.\nGot status %d, result: %s',...
                  status, result);
        end
        fprintf('Successfully built main file!\n');
    else
        % disp(['Commandline compilation of generated C code not yet supported under Windows.', ...
        %     'Please consider building the code in the c_generated_code folder from Windows Subsystem for Linux.'])
        % compile if on Windows platform
        [ status, result ] = system('mingw32-make.exe');
        if status
            cd(return_dir);
            error('building templated code failed.\nGot status %d, result: %s',...
                  status, result);
        end
        [ status, result ] = system('mingw32-make.exe shared_lib');
        if status
            cd(return_dir);
            error('building templated code as shared library failed.\nGot status %d, result: %s',...
                  status, result);
        end
        fprintf('Successfully built main file!\n');
    end
    % cd .. % fails if directory junctions are used (Matlab bug)
    cd(return_dir);
end