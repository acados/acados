from acados_template import AcadosOcpQp, AcadosOcpQpSolver
import time
import numpy as np

from acados_template.acados_ocp_iterate import AcadosOcpIterate


sqp_qp_sol_pairs = [
        ('last_qp_nonuniform_pendulum.json', 'sqp_sol_nonuniform_pendulum.json'),
                    ('last_qp_one_sided_test.json', 'sqp_sol_one_sided_test.json')
                    ]

# NOTE: create test cases in SQP example using:
    # id = 'xxx'
    # ocp_solver.dump_last_qp_to_json(filename=f'last_qp_{id}.json', overwrite=True)
    # ocp_solver.store_iterate(filename=f'sqp_sol_{id}.json', overwrite=True)


for qp_json_file, sqp_sol_file in sqp_qp_sol_pairs:

    sqp_sol = AcadosOcpIterate.from_json(sqp_sol_file)

    start_time = time.time()
    qp = AcadosOcpQp.from_json(qp_json_file)
    end_time = time.time()
    print(f"Loading took {end_time - start_time:.4f} seconds")


    solver = AcadosOcpQpSolver(qp)

    solver.solve()
    sol = solver.get_iterate()

    sqp_sol = AcadosOcpIterate.from_json(sqp_sol_file)
    
    # print(f"sol x: {sol.x}")
    # print(f"sol u: {sol.u}")
    # breakpoint()

    # # assert duals (lam and pi) are the same
    # print(f"sol lam: {sol.lam}")
    # print(f"sqp_sol lam: {sqp_sol.lam}")
    # print(f"sol pi: {sol.pi}")
    # print(f"sqp_sol pi: {sqp_sol.pi}")
    tol = 1e-5
    for stage in range(1,len(sol.lam)):
        assert np.allclose(sol.lam[stage], sqp_sol.lam[stage], atol=tol), \
            f"lam mismatch at stage {stage}: max diff = {np.max(np.abs(sol.lam[stage] - sqp_sol.lam[stage]))}"
    
    for stage in range(1,len(sol.pi)):
        assert np.allclose(sol.pi[stage], sqp_sol.pi[stage], atol=tol), \
            f"pi mismatch at stage {stage}: max diff = {np.max(np.abs(sol.pi[stage] - sqp_sol.pi[stage]))}"
    print(f"QP solution matches SQP solution for {qp_json_file} and {sqp_sol_file}.")

# solver.print_statistics()



# breakpoint()
