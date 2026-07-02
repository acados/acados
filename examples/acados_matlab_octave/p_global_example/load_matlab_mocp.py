import warnings
import numpy as np
from acados_template import AcadosOcpSolver, AcadosMultiphaseOcp

def main_mocp_json_load(json_file: str):
    warnings.filterwarnings("ignore", message=".*not in dictionary.*", category=UserWarning)

    for code_gen in [False, True]:
        mocp = AcadosMultiphaseOcp.from_json(json_file)
        ocp_solver = AcadosOcpSolver(mocp, generate=code_gen, build=code_gen)

        ocp_solver.set_p_global_and_precompute_dependencies(mocp.p_global_values)

        for iter in range(30):
            status = ocp_solver.solve()
            # ocp_solver.print_statistics()
            assert status == 0
        ocp_solver.print_statistics()
        qp_iter = ocp_solver.get_stats('qp_iter')
        assert qp_iter[-1] == 3
        print(f"test with {code_gen=} worked.")


if __name__ == "__main__":
    main_mocp_json_load('mocp_blz_true_pglbl_true_lut_true.json')
