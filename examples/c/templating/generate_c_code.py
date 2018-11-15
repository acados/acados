from jinja2 import Environment, FileSystemLoader
from ocp_nlp_dims import *

# setting up loader and environment
file_loader = FileSystemLoader('c_templates')
env = Environment(loader = file_loader)
template = env.get_template('template_example.in.c')

# set ocp_nlp_dimensions
nlp_dims = ocp_nlp_dims()
nlp_dims.nx = 4
nlp_dims.nu = 1
nlp_dims.N  = 100

# render template
output = template.render(ocp_nlp_dims=nlp_dims)

# output file
out_file = open('./c_generated_code/template_example.c', 'w+')
out_file.write(output)
