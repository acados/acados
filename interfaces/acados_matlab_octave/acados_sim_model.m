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

classdef acados_sim_model < handle
    


    properties
        model_struct
    end %properties



    methods
        

        function obj = acados_sim_model()
            obj.model_struct = struct;
            obj.model_struct.name = 'sim_model';
            obj.model_struct.ext_fun_type = 'casadi'; % generic
        end


        function obj = set(obj, field, value)
            % misc
            if (strcmp(field, 'name'))
                obj.model_struct.name = value;
            elseif (strcmp(field, 'ext_fun_type'))
                obj.model_struct.ext_fun_type = value;
            % dims
            elseif (strcmp(field, 'dim_nx'))
                obj.model_struct.dim_nx = value;
            elseif (strcmp(field, 'dim_nu'))
                obj.model_struct.dim_nu = value;
            elseif (strcmp(field, 'dim_nz'))
                obj.model_struct.dim_nz = value;
            elseif (strcmp(field, 'dim_np'))
                obj.model_struct.dim_np = value;
            % symbolics
            elseif (strcmp(field, 'sym_x'))
                obj.model_struct.sym_x = value;
            elseif (strcmp(field, 'sym_xdot'))
                obj.model_struct.sym_xdot = value;
            elseif (strcmp(field, 'sym_u'))
                obj.model_struct.sym_u = value;
            elseif (strcmp(field, 'sym_z'))
                obj.model_struct.sym_z = value;
            elseif (strcmp(field, 'sym_p'))
                obj.model_struct.sym_p = value;
            % dynamics
            elseif (strcmp(field, 'dyn_type'))
                obj.model_struct.dyn_type = value;
            elseif (strcmp(field, 'dyn_expr_f'))
                obj.model_struct.dyn_expr_f = value;
            elseif (strcmp(field, 'T'))
                obj.model_struct.T = value;
            elseif (strcmp(field, 'seed_adj'))
                obj.model_struct.seed_adj = value;
            else
                disp(['acados_sim_model: set: wrong field: ', field]);
				keyboard;
            end
        end

    end % methods



end % class

