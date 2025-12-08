import json, os
import warnings
import numpy as np
from acados_template import AcadosOcp, AcadosOcpSolver


def create_solver_from_json(json_file):

    # load json, store options in object
    with open(json_file, 'r') as f:
        acados_ocp_json = json.load(f)
    acados_ocp_json['json_file'] = os.path.abspath(json_file)

    warnings.filterwarnings("ignore", message=".*not in dictionary.*", category=UserWarning)
    ocp = AcadosOcp().from_dict(acados_ocp_json)

    ocp.render_templates()

    solver = AcadosOcpSolver(ocp)
    solver.solve()
    solver.print_statistics()

    print("create_solver_from_json passed.")


if __name__ == "__main__":
    json_file = "acados_ocp.json"
    create_solver_from_json(json_file)
