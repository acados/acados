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

classdef AcadosCodeGenOpts < handle

    properties
        cython_include_dirs
        shared_lib_ext
        acados_include_path
        acados_lib_path
        os
        acados_link_libs
        json_file
        code_export_directory
        acados_version
    end

    methods
        function obj = AcadosCodeGenOpts()
            obj.cython_include_dirs = []; % just for python compatibility

            acados_folder = getenv('ACADOS_INSTALL_DIR');
            obj.acados_include_path = [acados_folder, '/include'];
            obj.acados_link_libs = struct();

            if ismac
                obj.shared_lib_ext = '.dylib';
            else
                obj.shared_lib_ext = '.so';
            end

            if ismac
                obj.os = 'mac';
            elseif isunix
                obj.os = 'unix';
            else
                obj.os = 'pc';
            end

            % public
            obj.acados_lib_path = [acados_folder, '/lib'];
            obj.json_file = '';
            obj.code_export_directory = '';
            obj.acados_version = '';

        end

        function make_consistent(obj)
            acados_folder = getenv('ACADOS_INSTALL_DIR');
            addpath(fullfile(acados_folder, 'external', 'jsonlab'));
            libs = loadjson(fileread(fullfile(obj.acados_lib_path, 'link_libs.json')));
            obj.acados_link_libs = orderfields(libs);
            obj.json_file = absolute_path(obj.json_file);

            obj.acados_version = fileread(fullfile(obj.acados_lib_path, 'git_commit_hash'));
            if isempty(obj.code_export_directory)
                obj.code_export_directory = 'c_generated_code';
            end
            obj.code_export_directory = absolute_path(obj.code_export_directory);
        end
        function s = to_struct(self)
            if exist('properties')
                publicProperties = eval('properties(self)');
            else
                publicProperties = fieldnames(self);
            end
            s = struct();
            for fi = 1:numel(publicProperties)
                s.(publicProperties{fi}) = self.(publicProperties{fi});
            end
        end
    end
    methods (Static)
        function obj = from_struct(s)
            % Create AcadosCodeGenOpts from a struct (e.g. decoded from JSON).
            obj = AcadosCodeGenOpts();
            fields = fieldnames(s);
            for i = 1:length(fields)
                f = fields{i};
                % direct assignment for simple fields
                try
                    obj.(f) = s.(f);
                catch
                    % ignore unknown fields
                    warning(['Could not assign field ' f ' in AcadosCodeGenOpts.from_struct']);
                end
            end
        end
    end
end
