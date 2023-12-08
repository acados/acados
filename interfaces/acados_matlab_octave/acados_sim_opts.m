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

classdef acados_sim_opts < handle



    properties
        opts_struct
    end % properties



    methods

        function obj = acados_sim_opts()
            obj.opts_struct = struct;
            obj.opts_struct.compile_interface = 'auto'; % auto, true, false
            obj.opts_struct.codgen_model = 'true';
            obj.opts_struct.compile_model = 'true';
            obj.opts_struct.method = 'irk';
            obj.opts_struct.collocation_type = 'gauss_legendre';
            obj.opts_struct.num_stages = 4;
            obj.opts_struct.num_steps = 1;
            obj.opts_struct.newton_iter = 3;
            obj.opts_struct.newton_tol = 0.0;
            obj.opts_struct.sens_forw = 'false';
            obj.opts_struct.sens_adj = 'false';
            obj.opts_struct.sens_hess = 'false';
            obj.opts_struct.sens_algebraic = 'false';
            obj.opts_struct.output_z = 'false';
            obj.opts_struct.jac_reuse = 'false';
            obj.opts_struct.gnsf_detect_struct = 'true';
            obj.opts_struct.output_dir = fullfile(pwd, 'build');
            obj.opts_struct.ext_fun_compile_flags = '-O2';
            obj.opts_struct.parameter_values = [];
        end


        function obj = set(obj, field, value)
            % convert Matlab strings to char arrays
            if isstring(value)
                value = char(value);
            end

            if (strcmp(field, 'compile_interface'))
                obj.opts_struct.compile_interface = value;
            elseif (strcmp(field, 'ext_fun_compile_flags'))
                obj.opts_struct.ext_fun_compile_flags = value;
            elseif (strcmp(field, 'codgen_model'))
                obj.opts_struct.codgen_model = value;
            elseif (strcmp(field, 'compile_model'))
                obj.opts_struct.compile_model = value;
            elseif (strcmp(field, 'num_stages'))
                obj.opts_struct.num_stages = value;
            elseif (strcmp(field, 'num_steps'))
                obj.opts_struct.num_steps = value;
            elseif (strcmp(field, 'newton_iter'))
                obj.opts_struct.newton_iter = value;
            elseif (strcmp(field, 'newton_tol'))
                obj.opts_struct.newton_iter = value;
            elseif (strcmp(field, 'method'))
                obj.opts_struct.method = value;
            elseif (strcmp(field, 'sens_forw'))
                obj.opts_struct.sens_forw = value;
            elseif (strcmp(field, 'sens_adj'))
                obj.opts_struct.sens_adj = value;
            elseif (strcmp(field, 'sens_hess'))
                obj.opts_struct.sens_hess = value;
            elseif (strcmp(field, 'sens_algebraic'))
                obj.opts_struct.sens_algebraic = value;
            elseif (strcmp(field, 'output_z'))
                obj.opts_struct.output_z = value;
            elseif (strcmp(field, 'jac_reuse'))
                obj.opts_struct.jac_reuse = value;
            elseif (strcmp(field, 'gnsf_detect_struct'))
                obj.opts_struct.gnsf_detect_struct = value;
            elseif (strcmp(field, 'output_dir'))
                obj.opts_struct.output_dir = value;
            elseif (strcmp(field, 'collocation_type'))
                obj.opts_struct.collocation_type = value;
            elseif (strcmp(field, 'compile_mex'))
                disp(['Option compile_mex is not supported anymore,'...
                    'please use compile_interface instead or dont set the option.', ...
                    'options are: true, false, auto.']);
                keyboard
            elseif (strcmp(field, 'parameter_values'))
                obj.opts_struct.parameter_values = value;
            else
                disp(['acados_sim_opts: set: wrong field: ', field]);
            end
        end

    end % methods



end % class
