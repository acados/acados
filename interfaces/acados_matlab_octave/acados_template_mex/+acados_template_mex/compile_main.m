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
        % compile if on Windows platform
        disp(['Compilation of generated C code main file not thoroughly tested under Windows. Attempting to continue.'])

        % check compiler
        use_msvc = false;
        if ~is_octave()
            mexOpts = mex.getCompilerConfigurations('C', 'Selected');
            if contains(mexOpts.ShortName, 'MSVC')
                use_msvc = true;
            end
        end

        % get compiler
        if use_msvc
            % get env vars for MSVC
            msvc_env = fullfile(mexOpts.Location, 'VC\Auxiliary\Build\vcvars64.bat');
            assert(isfile(msvc_env), 'Cannot find definition of MSVC env vars.');

            make_cmd = sprintf('"%s" & nmake', msvc_env);

            % TODO
            warning('Templated Makefile not (yet) implemented for MSVC compiler.')
            cd(return_dir);
            return;
        else
            % using MinGW
            make_cmd = 'mingw32-make.exe';
        end

        % compile
        [ status, result ] = system(make_cmd);
        if status
            cd(return_dir);
            error('Building templated code failed.\nGot status %d, result: %s',...
                  status, result);
        end
        [ status, result ] = system(sprintf('%s shared_lib', make_cmd));
        if status
            cd(return_dir);
            error('Building templated code as shared library failed.\nGot status %d, result: %s',...
                  status, result);
        end
        fprintf('Successfully built main file!\n');
    end
    cd(return_dir);
end