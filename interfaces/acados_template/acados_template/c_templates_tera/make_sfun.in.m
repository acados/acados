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

SOURCES = [ ...
        {%- if solver_options.integrator_type == 'ERK' %}
            '{{ model.name }}_model/{{ model.name }}_expl_ode_fun.c ', ...
            '{{ model.name }}_model/{{ model.name }}_expl_vde_forw.c ',...
            {%- if solver_options.hessian_approx == 'EXACT' %}
            '{{ model.name }}_model/{{ model.name }}_expl_ode_hess.c ',...
            {%- endif %}
        {%- elif solver_options.integrator_type == "IRK" %}
            '{{ model.name }}_model/{{ model.name }}_impl_dae_fun.c ', ...
            '{{ model.name }}_model/{{ model.name }}_impl_dae_fun_jac_x_xdot_z.c ', ...
            '{{ model.name }}_model/{{ model.name }}_impl_dae_jac_x_xdot_u_z.c ', ...
            {%- if solver_options.hessian_approx == 'EXACT' %}
            '{{ model.name }}_model/{{ model.name }}_impl_dae_hess.c ',...
            {%- endif %}
        {%- elif solver_options.integrator_type == "GNSF" %}
            '{{ model.name }}_model/{{ model.name }}_gnsf_phi_fun.c ',...
            '{{ model.name }}_model/{{ model.name }}_gnsf_phi_fun_jac_y.c ',...
            '{{ model.name }}_model/{{ model.name }}_gnsf_phi_jac_y_uhat.c ',...
            '{{ model.name }}_model/{{ model.name }}_gnsf_f_lo_fun_jac_x1k1uz.c ',...
            '{{ model.name }}_model/{{ model.name }}_gnsf_get_matrices_fun.c ',...
        {%- endif %}
        {%- if cost.cost_type == "NONLINEAR_LS" %}
            '{{ model.name }}_cost/{{ model.name }}_cost_y_fun.c ',...
            '{{ model.name }}_cost/{{ model.name }}_cost_y_fun_jac_ut_xt.c ',...
            '{{ model.name }}_cost/{{ model.name }}_cost_y_hess.c ',...
        {%- elif cost.cost_type == "EXTERNAL" %}
            '{{ model.name }}_cost/{{ model.name }}_cost_ext_cost_fun.c ',...
            '{{ model.name }}_cost/{{ model.name }}_cost_ext_cost_fun_jac.c ',...
            '{{ model.name }}_cost/{{ model.name }}_cost_ext_cost_fun_jac_hess.c ',...
        {%- endif %}
        {%- if cost.cost_type_e == "NONLINEAR_LS" %}
            '{{ model.name }}_cost/{{ model.name }}_cost_y_e_fun.c ',...
            '{{ model.name }}_cost/{{ model.name }}_cost_y_e_fun_jac_ut_xt.c ',...
            '{{ model.name }}_cost/{{ model.name }}_cost_y_e_hess.c ',...
        {%- elif cost.cost_type_e == "EXTERNAL" %}
            '{{ model.name }}_cost/{{ model.name }}_cost_ext_cost_e_fun.c ',...
            '{{ model.name }}_cost/{{ model.name }}_cost_ext_cost_e_fun_jac.c ',...
            '{{ model.name }}_cost/{{ model.name }}_cost_ext_cost_e_fun_jac_hess.c ',...
        {%- endif %}
        {%- if constraints.constr_type == "BGH"  and dims.nh > 0 %}
            '{{ model.name }}_constraints/{{ model.name }}_constr_h_fun.c ', ...
            '{{ model.name }}_constraints/{{ model.name }}_constr_h_fun_jac_uxt_hess.c ', ...
            '{{ model.name }}_constraints/{{ model.name }}_constr_h_fun_jac_uxt_zt.c ', ...
        {%- elif constraints.constr_type == "BGP" and dims.nphi > 0 %}
            '{{ model.name }}_constraints/{{ model.name }}_phi_constraint.c ', ...
        {%- endif %}
        {%- if constraints.constr_type_e == "BGH"  and dims.nh_e > 0 %}
            '{{ model.name }}_constraints/{{ model.name }}_constr_h_e_fun.c ', ...
            '{{ model.name }}_constraints/{{ model.name }}_constr_h_e_fun_jac_uxt_zt_hess.c ', ...
            '{{ model.name }}_constraints/{{ model.name }}_constr_h_e_fun_jac_uxt_zt.c ', ...
        {%- elif constraints.constr_type_e == "BGP" and dims.nphi_e > 0 %}
            '{{ model.name }}_constraints/{{ model.name }}_phi_e_constraint.c ', ...
        {%- endif %}
            'acados_solver_sfunction_{{ model.name }}.c ', ...
            'acados_solver_{{ model.name }}.c '
          ];

INC_PATH = '{{ acados_include_path }}';

INCS = [ ' -I', fullfile(INC_PATH, 'blasfeo', 'include'), ...
         ' -I', fullfile(INC_PATH, 'hpipm', 'include'), ...
        ' -I', INC_PATH, ' -I', fullfile(INC_PATH, 'acados'), ' '];

{% if  solver_options.qp_solver == "FULL_CONDENSING_QPOASES" %}
INCS = strcat(INCS, '-I', fullfile(INC_PATH, 'qpOASES_e') )
{% endif %}

CFLAGS  = ' -O';

{% if  solver_options.qp_solver == "FULL_CONDENSING_QPOASES" %}
CFLAGS = [ CFLAGS, ' -DACADOS_WITH_QPOASES ' ];
{% endif %}

LIB_PATH = '{{ acados_lib_path }}';

LIBS = '-lacados -lhpipm -lblasfeo';

{% if  solver_options.qp_solver == "FULL_CONDENSING_QPOASES" %}
LIBS = strcat(LIBS, ' -lqpOASES_e');
{% endif %}

eval( [ 'mex -v -output  acados_solver_sfunction_{{ model.name }} ', ...
    CFLAGS, INCS, ' ', SOURCES, ' -L', LIB_PATH, ' ', LIBS ]);

fprintf( [ '\n\nSuccessfully created sfunction:\nacados_solver_sfunction_{{ model.name }}', '.', ...
    eval('mexext')] );


%% print note on usage of s-function
fprintf('\n\nNote: Usage of Sfunction is as follows:\n')
input_note = 'Inputs are:\n1) x0, initial state, size [{{ dims.nx }}]\n ';
i_in = 2;
{%- if dims.ny > 0 %}
input_note = strcat(input_note, num2str(i_in), ') y_ref - concatenated for intermediate stages,',...
                    ' size [{{ dims.N * dims.ny }}]\n ');
i_in = i_in + 1;
{%- endif %}

{%- if dims.ny > 0 %}
input_note = strcat(input_note, num2str(i_in), ') y_ref_e, size [{{ dims.ny_e }}]\n ');
i_in = i_in + 1;
{%- endif %}

{%- if dims.np > 0 %}
input_note = strcat(input_note, num2str(i_in), ') parameters - concatenated for all stages,',...
                    ' size [{{ (dims.N+1)*dims.np }}]\n ');
i_in = i_in + 1;
{%- endif %}

{%- if dims.nbx > 0 %}
input_note = strcat(input_note, num2str(i_in), ') lbx, size [{{ dims.nbx }}]\n ');
i_in = i_in + 1;
input_note = strcat(input_note, num2str(i_in), ') ubx, size [{{ dims.nbx }}]\n ');
i_in = i_in + 1;
{%- endif %}

{%- if dims.nbu > 0 %}
input_note = strcat(input_note, num2str(i_in), ') lbu, size [{{ dims.nbu }}]\n ');
i_in = i_in + 1;
input_note = strcat(input_note, num2str(i_in), ') ubu, size [{{ dims.nbu }}]\n ');
i_in = i_in + 1;
{%- endif %}

{%- if dims.ng > 0 %}
input_note = strcat(input_note, num2str(i_in), ') lg, size [{{ dims.ng }}]\n ');
i_in = i_in + 1;
input_note = strcat(input_note, num2str(i_in), ') ug, size [{{ dims.ng }}]\n ');
i_in = i_in + 1;
{%- endif %}

{%- if dims.nh > 0 %}
input_note = strcat(input_note, num2str(i_in), ') lh, size [{{ dims.nh }}]\n ');
i_in = i_in + 1;
input_note = strcat(input_note, num2str(i_in), ') uh, size [{{ dims.nh }}]\n ');
i_in = i_in + 1;
{%- endif %}

fprintf(input_note)

disp(' ')

output_note = strcat('Outputs are:\n', ...
                '1) u0 - optimal input, size [{{ dims.nu }}]\n',...
                '2) acados solver status (0 = SUCCESS)\n3) KKT residual\n',...
                '4) first state \n5) CPU time\n6) sqp iter\n');

fprintf(output_note)
