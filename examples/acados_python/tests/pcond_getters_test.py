import os
import numpy as np
import casadi as ca
from acados_template import AcadosModel, AcadosOcp, AcadosOcpSolver



if __name__ == "__main__":
    N = 40
    Tf = 10.0
    cond_N = 10

    # Problem data
    nx, nu = 2, 1
    Q = np.diag([5.0, 1.0])
    R = np.diag([0.2])
    x0 = np.zeros(nx)
    x_ref = np.array([10.0, 0.0])
    lbu = np.array([-5.5])
    ubu = np.array([5.5])

    # Build OCP
    ocp = AcadosOcp()
    ocp.solver_options.tf = Tf
    ocp.solver_options.N_horizon = N

    ocp.model = AcadosModel()
    ocp.model.name = "double_integrator"

    s, v = ca.SX.sym("s"), ca.SX.sym("v")
    a = ca.SX.sym("a")
    x = ca.vertcat(s, v)
    u = ca.vertcat(a)
    xdot = ca.vertcat(v, a)

    ocp.model.x = x
    ocp.model.u = u
    ocp.model.xdot = xdot
    ocp.model.f_expl_expr = xdot
    ocp.model.f_impl_expr = xdot

    # Bounds
    ocp.constraints.x0 = x0
    ocp.constraints.idxbu = np.array([0], dtype=int)
    ocp.constraints.lbu = lbu
    ocp.constraints.ubu = ubu

    # Linear least-squares cost
    ny = nx + nu
    W = np.zeros((ny, ny))
    W[:nx, :nx] = Q
    W[nx:, nx:] = R

    ocp.cost.cost_type = "LINEAR_LS"
    ocp.cost.cost_type_0 = "LINEAR_LS"
    ocp.cost.cost_type_e = "LINEAR_LS"
    ocp.cost.W = W
    ocp.cost.W_0 = W
    ocp.cost.W_e = Q

    ocp.cost.Vx = np.vstack((np.eye(nx), np.zeros((nu, nx))))
    ocp.cost.Vx_0 = np.vstack((np.eye(nx), np.zeros((nu, nx))))
    ocp.cost.Vx_e = np.eye(nx)
    ocp.cost.Vu = np.vstack((np.zeros((nx, nu)), np.eye(nu)))
    ocp.cost.Vu_0 = np.vstack((np.zeros((nx, nu)), np.eye(nu)))

    ocp.cost.yref = np.hstack((x_ref, np.zeros(nu)))
    ocp.cost.yref_0 = np.hstack((x_ref, np.zeros(nu)))
    ocp.cost.yref_e = x_ref

    # Solver options
    ocp.solver_options.qp_solver = "PARTIAL_CONDENSING_HPIPM"
    ocp.solver_options.qp_solver_cond_N = cond_N
    ocp.solver_options.hessian_approx = "EXACT"
    ocp.solver_options.reg_epsilon = 0.0
    ocp.solver_options.nlp_solver_type = "SQP_RTI"

    ocp.solver_options.qp_solver_cond_N = cond_N

    ocp.solver_options.qp_solver_cond_block_size = (cond_N) * [N // ((cond_N))] + [N - (cond_N * (N // cond_N))]
    # ocp.solver_options.qp_solver_cond_block_size = (cond_N) * [1] + [N-((cond_N))]

    uniq = f"double_integrator_{os.getpid()}"
    ocp.code_export_directory = f"c_generated_code/{uniq}"
    json_path = f".solver_files/ACADOS_SOLVER_{uniq}.json"

    os.makedirs(os.path.dirname(ocp.code_export_directory), exist_ok=True)
    os.makedirs(os.path.dirname(json_path), exist_ok=True)

    solver = AcadosOcpSolver(ocp, json_file=json_path)

    solver.solve()
    cond_N = ocp.solver_options.qp_solver_cond_N

    FIELDS_PCOND = [
        "pcond_A", "pcond_B", "pcond_b",
        "pcond_Q", "pcond_R", "pcond_S", "pcond_q", "pcond_r",
        "pcond_C", "pcond_D", "pcond_lg", "pcond_ug",
        "pcond_lbx", "pcond_ubx", "pcond_lbu", "pcond_ubu"
            ]
    qp = {"N": cond_N, "stages": []}
    for k in range(cond_N + 1):
        stage = {}
        for field in FIELDS_PCOND:
            print(f"trying to get {field} at stage {k}")
            try:
                arr = solver.get_from_qp_in(k, field)
                stage[field.split("_", 1)[1]] = np.array(arr).tolist()
            except Exception as e:
                print(f"[warn] stage {k}: couldn't get {field}: {e}")
        qp["stages"].append(stage)


    from pprint import pprint
    pprint(qp["stages"])
