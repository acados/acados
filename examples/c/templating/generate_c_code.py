from jinja2 import Environment, FileSystemLoader
from ocp_nlp_dims import *

# setting up loader and environment
file_loader = FileSystemLoader('c_templates')
env = Environment(loader = file_loader)
template = env.get_template('template_example.in.c')

# create render arguments
ra = ocp_nlp_render_arguments()

# set model_name 
ra.model_name = 'pendulum'

# set ocp_nlp_dimensions
nlp_dims = ra.dims
nlp_dims.nx = 4
nlp_dims.nu = 1
nlp_dims.N  = 100

# set weighting matrices

# set constants
const1 = ocp_nlp_constant()
const1.name  = 'PI'
const1.value = 3.1415926535897932
ra.constants = [const1]
# set QP solver

ra.solver_config.qp_solver = 'PARTIAL_CONDENSING_HPIPM'
# ra.solver_config.qp_solver = 'FULL_CONDENSING_QPOASES'

# render template
output = template.render(ra=ra)

# output file
out_file = open('./c_generated_code/template_example.c', 'w+')
out_file.write(output)
