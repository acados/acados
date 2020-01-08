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

SOURCES = [ 'acados_solver_sfunction_{{ ocp.model_name }}.c ', ...
            'acados_solver_{{ ocp.model_name }}.c ', ...
            {%- if  ocp.solver_config.integrator_type == 'ERK': %}
            '{{ ocp.model_name }}_model/{{ ocp.model_name }}_expl_ode_fun.c ', ...
            '{{ ocp.model_name }}_model/{{ ocp.model_name }}_expl_vde_forw.c ',...
            {% if ocp.solver_config.hessian_approx == 'EXACT': -%} 
            {% endif -%}
            {% else: %}
            '{{ ocp.model_name }}_model/{{ ocp.model_name }}_impl_dae_fun.c ', ...
            '{{ ocp.model_name }}_model/{{ ocp.model_name }}_impl_dae_fun_jac_x_xdot_z.c ', ...
            '{{ ocp.model_name }}_model/{{ ocp.model_name }}_impl_dae_jac_x_xdot_u_z.c ', ...
            {% endif -%}
            {% if ocp.dims.npd > 0: -%}
            '{{ ocp.con_p_name }}_p_constraint/{{ ocp.con_p_name }}_p_constraint.c ', ...
            {% endif -%}
            {% if ocp.dims.nh > 0: -%}
            '{{ ocp.con_h_name }}_h_constraint/{{ ocp.con_h_name }}_h_constraint.c ', ...
            {% endif -%}
          ];

INC_PATH = '{{ ocp.acados_include_path }}';

INCS = [ ' -I', INC_PATH, '/blasfeo/include/ ', ...
          '-I', INC_PATH, ' -I', INC_PATH, '/acados/ ', ...
          '-I', INC_PATH, '/qpOASES_e/' ];

CFLAGS  = ' -O';

{% if  'QPOASES' in ocp.solver_config.qp_solver: %}
CFLAGS = [ CFLAGS, ' -DACADOS_WITH_QPOASES ' ];
{% endif %}

LIB_PATH = '{{ ocp.acados_lib_path }}';

LIBS = '-lacados -lhpipm -lblasfeo -lqpOASES_e -lm'; 
    
eval( [ 'mex -v -output  acados_solver_sfunction_{{ ocp.model_name }} ', ...
    CFLAGS, INCS, ' ', SOURCES, ' -L', LIB_PATH, ' ', LIBS ]);

disp( [ 'acados_solver_sfunction_{{ ocp.model_name }}', '.', ...
    eval('mexext'), ' successfully created!'] );
