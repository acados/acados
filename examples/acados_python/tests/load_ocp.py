import json, os
import warnings
import numpy as np
from acados_template import AcadosOcp, AcadosOcpSolver


def create_solver_from_json(json_file):

    warnings.filterwarnings("ignore", message=".*not in dictionary.*", category=UserWarning)
    ocp = AcadosOcp.from_json(json_file)

    solver = AcadosOcpSolver(ocp)
    solver.solve()
    solver.print_statistics()

    print("create_solver_from_json passed.")


if __name__ == "__main__":
    json_file = "acados_ocp.json"
    create_solver_from_json(json_file)
