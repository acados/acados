import warnings
import numpy as np
from acados_template import AcadosOcpSolver, AcadosMultiphaseOcp
# Test assuming main_multiphase_ocp.m was run just before.

def main_mocp_json_load(json_file: str):
    warnings.filterwarnings("ignore", message=".*not in dictionary.*", category=UserWarning)

    for code_gen in [False, True]:
        mocp = AcadosMultiphaseOcp.from_json(json_file)
        ocp_solver = AcadosOcpSolver(mocp, generate=code_gen, build=code_gen)

        ocp_solver.set_p_global_and_precompute_dependencies(mocp.p_global_values)

        status = ocp_solver.solve()
        assert status == 0
        ocp_solver.print_statistics()

        print(f"test with {code_gen=} worked.")

if __name__ == "__main__":
    main_mocp_json_load('ocp_mocp.json')
