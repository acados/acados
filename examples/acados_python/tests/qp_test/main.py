from acados_template import AcadosOcpQp, AcadosOcpQpSolver, AcadosOcpIterate, AcadosOcpQpOptions
import time
import numpy as np

sqp_qp_sol_pairs = [
        ('last_qp_nonuniform_pendulum.json', 'sqp_sol_nonuniform_pendulum.json'),
        ('last_qp_one_sided_test.json', 'sqp_sol_one_sided_test.json')
                    ]

# NOTE: create test cases in SQP example using:
    # id = 'xxx'
    # ocp_solver.dump_last_qp_to_json(filename=f'last_qp_{id}.json', overwrite=True)
    # ocp_solver.store_iterate(filename=f'sqp_sol_{id}.json', overwrite=True)

for qp_solver in ['PARTIAL_CONDENSING_HPIPM', 'FULL_CONDENSING_HPIPM', 'FULL_CONDENSING_DAQP']:
    n_iter_list = []
    timings_list = []
    for qp_json_file, sqp_sol_file in sqp_qp_sol_pairs:

        sqp_sol = AcadosOcpIterate.from_json(sqp_sol_file)

        # start_time = time.time()
        qp = AcadosOcpQp.from_json(qp_json_file)
        # end_time = time.time()
        # print(f"Loading took {end_time - start_time:.4f} seconds")


        opts = AcadosOcpQpOptions()
        opts.iter_max = 500
        # opts.print_level = 1
        # opts.hpipm_mode = "ROBUST"
        opts.qp_solver = qp_solver
        solver = AcadosOcpQpSolver(qp, opts=opts)

        status = solver.solve()

        assert status == 0, f"QP solver returned status {status} != 0"
        sol = solver.get_iterate()
        n_iter_list.append(solver.get_stats("iter"))
        timings_list.append(solver.get_stats("time_tot"))

        sqp_sol = AcadosOcpIterate.from_json(sqp_sol_file)

        tol = 1e-5
        for stage in range(1,len(sol.lam)):
            assert np.allclose(sol.lam[stage], sqp_sol.lam[stage], atol=tol), \
                f"lam mismatch at stage {stage}: max diff = {np.max(np.abs(sol.lam[stage] - sqp_sol.lam[stage]))}"

        for stage in range(1,len(sol.pi)):
            assert np.allclose(sol.pi[stage], sqp_sol.pi[stage], atol=tol), \
                f"pi mismatch at stage {stage}: max diff = {np.max(np.abs(sol.pi[stage] - sqp_sol.pi[stage]))}"
        print(f"QP solution matches SQP solution for {qp_json_file} and {sqp_sol_file}.")


    print(f"QP solver used: {qp_solver}")
    print(f"Number of iterations for each QP: {n_iter_list}")
    print(f"Timings for each QP: {np.array(timings_list)*1e3} ms")
# solver.print_statistics()



# breakpoint()
