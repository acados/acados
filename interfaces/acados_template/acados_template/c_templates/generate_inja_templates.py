import glob, os
for template_name in glob.glob("*.in.*"):
    with open(template_name, 'U') as f:
        newText=f.read()
    
        while '{% set outer_loop = loop %}' in newText:
            newText=newText.replace('{% set outer_loop = loop %}', '')    


        while 'outer_loop.index0' in newText:
            newText=newText.replace('outer_loop.index0', 'loop.parent.index0')    

    with open('../c_templates_inja/' + template_name, "w") as f:
        f.write(newText)
