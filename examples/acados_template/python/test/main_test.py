import itertools as it
import os

FORMULATION_values = ['LS', 'NLS']
QP_SOLVER_values = ['PARTIAL_CONDENSING_HPIPM', 'FULL_CONDENSING_HPIPM', 'FULL_CONDENSING_QPOASES']
INTEGRATOR_TYPE_values = ['ERK', 'IRK']
SOLVER_TYPE_values = ['SQP', 'SQP_RTI']

test_parameters = { 'FORMULATION_values': FORMULATION_values, 
                    'QP_SOLVER_values': QP_SOLVER_values, 
                    'INTEGRATOR_TYPE_values': INTEGRATOR_TYPE_values,
                    'SOLVER_TYPE_values': SOLVER_TYPE_values} 

all_test_parameters = sorted(test_parameters)
combinations = list(it.product(*(test_parameters[Name] for Name in all_test_parameters)))

# combinations = [('LS', 'IRK', 'PARTIAL_CONDENSING_HPIPM', 'SQP')]
for parameters in combinations:
    os_cmd = ("python generate_c_code.py" + 
        " --FORMULATION {}".format(parameters[0]) + 
        " --INTEGRATOR_TYPE {}".format(parameters[1]) +   
        " --QP_SOLVER {}".format(parameters[2]) +
        " --SOLVER_TYPE {}".format(parameters[3])) 
    status = os.system(os_cmd)
    if status != 0:
        raise Exception("acados status  = {} on test {}. Exiting\n".format(status, parameters))

