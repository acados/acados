import warnings

from acados_template import AcadosOcpSolver, AcadosOcp

def main_mocp_json_load(json_file: str):
    warnings.filterwarnings("ignore", message=".*not in dictionary.*", category=UserWarning)

    mocp = AcadosOcp.from_json(json_file)
    ocp_solver = AcadosOcpSolver(mocp, generate=False, build=False)


    status = ocp_solver.solve()
    ocp_solver.print_statistics()


if __name__ == "__main__":
    main_mocp_json_load('pendulum_on_cart_ocp.json')
