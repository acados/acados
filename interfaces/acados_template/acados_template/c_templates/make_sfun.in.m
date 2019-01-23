SOURCES = [ 'acados_solver_sfunction_{{ ra.model_name }}.c ', ...
            'acados_solver_{{ ra.model_name }}.c ', ...
            {%- if  ra.solver_config.integrator_type == 'ERK': %}
            '{{ ra.model_name }}_model/{{ ra.model_name }}_expl_ode_fun.c ', ...
            '{{ ra.model_name }}_model/{{ ra.model_name }}_expl_vde_forw.c ',...
            {% if ra.solver_config.hessian_approx == 'EXACT': -%} 
            {% endif -%}
            {% else: %}
            '{{ ra.model_name }}_model/{{ ra.model_name }}_impl_dae_fun.c ', ...
            '{{ ra.model_name }}_model/{{ ra.model_name }}_impl_dae_fun_jac_x_xdot_z.c ', ...
            '{{ ra.model_name }}_model/{{ ra.model_name }}_impl_dae_jac_x_xdot_u_z.c '
            {% endif -%}
          ];

INC_PATH = '{{ ra.acados_include_path }}';

INCS = [ ' -I', INC_PATH, '/blasfeo/include/ ', ...
          '-I', INC_PATH, ' -I', INC_PATH, '/acados/ ', ...
          '-I', INC_PATH, '/qpOASES_e/' ];

CFLAGS  = ' -O';

{% if  'QPOASES' in ra.solver_config.qp_solver: %}
CFLAGS = [ CFLAGS, ' -DACADOS_WITH_QPOASES ' ];
{% endif %}

LIB_PATH = '{{ ra.acados_lib_path }}';

LIBS = '-lacados -lhpipm -lblasfeo -lqpOASES_e -lm'; 
    
eval( [ 'mex -v -output  acados_solver_sfunction_{{ ra.model_name }} ', ...
    CFLAGS, INCS, ' ', SOURCES, ' -L', LIB_PATH, ' ', LIBS ]);

disp( [ 'acados_solver_sfunction_{{ ra.model_name }}', '.', ...
    eval('mexext'), ' successfully created!'] );
