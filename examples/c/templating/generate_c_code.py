from jinja2 import Environment, FileSystemLoader

file_loader = FileSystemLoader('c_templates')

env = Environment(loader = file_loader)

template = env.get_template('template_example.in.c')

output = template.render(nx=4)

out_file = open('./c_generated_code/template_example.c', 'w+')
out_file.write(output)
