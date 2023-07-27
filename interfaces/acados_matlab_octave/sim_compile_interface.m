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

function sim_compile_interface(opts)

% get acados folder
acados_folder = getenv('ACADOS_INSTALL_DIR');
mex_flags = getenv('ACADOS_MEX_FLAGS');

% set paths
acados_mex_folder = fullfile(acados_folder, 'interfaces', 'acados_matlab_octave');
acados_include = ['-I' acados_folder];
acados_interfaces_include = ['-I' fullfile(acados_folder, 'interfaces')];
acados_lib_path = ['-L' fullfile(acados_folder, 'lib')];

mex_names = { ...
    'sim_solve', ...
    'sim_get', ...
};

mex_files = cell(length(mex_names), 1);
for k=1:length(mex_names)
    mex_files{k} = fullfile(acados_mex_folder, [mex_names{k}, '.c']);
end


%% compile mex
if is_octave()
    if ~exist(fullfile(opts.output_dir, 'cflags_octave.txt'), 'file')
        diary(fullfile(opts.output_dir, 'cflags_octave.txt'));
        diary on
        mkoctfile -p CFLAGS
        diary off
        input_file = fopen(fullfile(opts.output_dir, 'cflags_octave.txt'), 'r');
        cflags_tmp = fscanf(input_file, '%[^\n]s');
        fclose(input_file);
        cflags_tmp = [cflags_tmp, ' -std=c99'];
        input_file = fopen(fullfile(opts.output_dir, 'cflags_octave.txt'), 'w');
        fprintf(input_file, '%s', cflags_tmp);
        fclose(input_file);
    end
    % read cflags from file
    input_file = fopen(fullfile(opts.output_dir, 'cflags_octave.txt'), 'r');
    cflags_tmp = fscanf(input_file, '%[^\n]s');
    fclose(input_file);
    setenv('CFLAGS', cflags_tmp);
end

for ii=1:length(mex_files)
    disp(['compiling ', mex_files{ii}])
    if is_octave()
%        mkoctfile -p CFLAGS
        mex(acados_include, acados_interfaces_include, acados_lib_path,...
            '-lacados', '-lhpipm', '-lblasfeo', mex_files{ii})
    else
        FLAGS = 'CFLAGS=$CFLAGS -std=c99';
        mex(mex_flags, FLAGS, acados_include, acados_interfaces_include, acados_lib_path, ...
             '-lacados', '-lhpipm', '-lblasfeo', mex_files{ii}, '-outdir', opts.output_dir)
    end
end


if is_octave()
    octave_version = OCTAVE_VERSION();
    if octave_version < 5
        movefile('*.o', opts.output_dir);
    end

    %system(['mv -f *.mexa64 ', opts.output_dir])
    for k=1:length(mex_names)
        clear(mex_names{k})
    %    [status, message] = movefile([mex_names{k}, '.', mexext], opts.output_dir)
        [status, message] = copyfile([mex_names{k}, '.', mexext], opts.output_dir);
        delete([mex_names{k}, '.', mexext]);
    end
end

end
