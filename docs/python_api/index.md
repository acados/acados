# Python API

<!-- ``` eval_rst
.. automodule:: acados_template.
    :members:
    :private-members:
    :undoc-members:
``` -->

## Overview
The following image shows an overview of the available classes in acados and their dependencies.
``` eval_rst
.. graphviz:: py_acados_classes.dot
    :name: sphinx.ext.graphviz.pyclasses
    :caption: Python API classes overview
    :alt: Overview of acados Python classes
    :align: center
```

## acados OCP solver
``` eval_rst
.. automodule:: acados_template.acados_ocp_solver
    :members:
    :private-members:
    :exclude-members: make_ocp_dims_consistent, get_ocp_nlp_layout
```


<!-- ## acados OCP -->
``` eval_rst
.. automodule:: acados_template.acados_ocp
    :members:
    :private-members:
    :exclude-members:
```


<!-- ## acados model -->
``` eval_rst
.. automodule:: acados_template.acados_model
    :members:
    :private-members:
    :exclude-members:
```



## acados integrator interface
The class `AcadosSim` can be used to formulate a simulation problem, for which an acados integrator (`AcadosSimSolver`) can be created.
``` eval_rst
.. automodule:: acados_template.acados_sim
    :members:
    :private-members:
    :exclude-members:
```

<!-- ## acados sim solver -->
``` eval_rst
.. automodule:: acados_template.acados_sim_solver
    :members:
    :private-members:
    :exclude-members: make_ocp_dims_consistent, get_ocp_nlp_layout
```
