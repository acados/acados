import sys
sys.path.insert(0, '../getting_started')
import numpy as np
import casadi as ca
from typing import Union

from acados_template import AcadosOcp, AcadosOcpSolver, AcadosCasadiOcpSolver
from pendulum_model import export_pendulum_ode_model

def formulate_ocp(bu_cons=True)-> AcadosOcp:
    # create ocp object to formulate the OCP
    ocp = AcadosOcp()
    Tf: float = 1.0
    N: int = 20

    # set model
    model = export_pendulum_ode_model()
    # add dummy control
    model.u = ca.vertcat(model.u, ca.SX.sym('dummy_u'))
    ocp.model = model

    nx = model.x.rows()
    nu = model.u.rows()
    
    # set prediction horizon
    ocp.solver_options.N_horizon = N
    ocp.solver_options.tf = Tf

    # cost matrices
    Q_mat = 2*np.diag([1e3, 1e3, 1e-2, 1e-2])
    R_mat = 2*np.diag([1e-2, 1e-2])

    # path cost
    ocp.cost.cost_type = 'NONLINEAR_LS'
    ocp.model.cost_y_expr = ca.vertcat(model.x, model.u)
    ocp.cost.yref = np.zeros((nx+nu,))
    ocp.cost.W = ca.diagcat(Q_mat, R_mat).full()

    # terminal cost
    ocp.cost.cost_type_e = 'NONLINEAR_LS'
    ocp.cost.yref_e = np.zeros((nx,))
    ocp.model.cost_y_expr_e = model.x
    ocp.cost.W_e = Q_mat

    # set constraints
    ocp.constraints.x0 = np.array([0, np.pi, 0, 0])  # initial state

    # path constraints
    Fmax = 80
    if bu_cons:
        ocp.constraints.lbu = np.array([-Fmax])
        ocp.constraints.ubu = np.array([+Fmax])
        ocp.constraints.idxbu = np.array([0])
    else:
        ocp.model.con_h_expr_0 = model.u[0]
        ocp.constraints.lh_0 = np.array([-Fmax])
        ocp.constraints.uh_0 = np.array([Fmax])
        ocp.constraints.idxh = np.array([0])

        ocp.model.con_h_expr = model.u[0]
        ocp.constraints.lh = np.array([-Fmax])
        ocp.constraints.uh = np.array([Fmax])
        ocp.constraints.idxh = np.array([0])

    # set partial bounds for state
    ocp.constraints.lbx = np.array([-10,-10])
    ocp.constraints.ubx = np.array([10, 10])
    ocp.constraints.idxbx = np.array([0,3])

    # set options
    ocp.solver_options.qp_solver = 'PARTIAL_CONDENSING_HPIPM' # FULL_CONDENSING_QPOASES
    ocp.solver_options.hessian_approx = 'GAUSS_NEWTON' # 'GAUSS_NEWTON', 'EXACT'
    ocp.solver_options.integrator_type = 'ERK'
    ocp.solver_options.nlp_solver_type = 'SQP' # SQP_RTI, SQP
    ocp.solver_options.globalization = 'MERIT_BACKTRACKING' # turns on globalization
    
    return ocp

def main(bu_cons: bool = True):
    N_horizon = 20
    F_max_new = 50.0
    x_max_new = 5.0

    ocp = formulate_ocp(bu_cons=True) # need to reformulate to set lh and uh
    initial_iterate = ocp.create_default_initial_iterate()
    ## solve using acados
    # create acados solver
    ocp_solver = AcadosOcpSolver(ocp,verbose=False)
    ocp_solver.set_iterate(initial_iterate)
    for i in range(1, N_horizon):
        ocp_solver.set(i, "lbu", np.array([-F_max_new]))
        ocp_solver.set(i, "ubu", np.array([F_max_new]))
        ocp_solver.set(i, "lbx", np.array([-x_max_new, -x_max_new]))
        ocp_solver.set(i, "ubx", np.array([x_max_new, x_max_new]))

    # solve with acados
    status = ocp_solver.solve()
    acados_x = np.array([ocp_solver.get(i, "x") for i in range(N_horizon+1)])
    acados_u = np.array([ocp_solver.get(i, "u") for i in range(N_horizon)])
    lam = np.concatenate([ocp_solver.get(i, "lam") for i in range(N_horizon+1)])
    pi = np.concatenate([ocp_solver.get(i, "pi") for i in range(N_horizon)])
    result = ocp_solver.get_iterate()

    # ## solve using casadi
    if bu_cons:
        pass
    else:
        ocp = formulate_ocp(bu_cons=False) # need to reformulate to set lh and uh

    casadi_ocp_solver = AcadosCasadiOcpSolver(ocp=ocp,solver="ipopt",verbose=False)
    # casadi_ocp_solver = AcadosOcpSolver(ocp,verbose=False)
    casadi_ocp_solver.set_iterate(result)

    for i in range(1, N_horizon):
        if bu_cons:
            casadi_ocp_solver.set(i, "lbu", np.array([-F_max_new]))
            casadi_ocp_solver.set(i, "ubu", np.array([F_max_new]))
        else:
            casadi_ocp_solver.set(i, "lh", np.array([-F_max_new]))
            casadi_ocp_solver.set(i, "uh", np.array([F_max_new]))
        casadi_ocp_solver.set(i, "lbx", np.array([-x_max_new, -x_max_new]))
        casadi_ocp_solver.set(i, "ubx", np.array([x_max_new, x_max_new]))

    casadi_ocp_solver.solve()
    casadi_x = np.array([casadi_ocp_solver.get(i, "x") for i in range(N_horizon+1)])
    casadi_u = np.array([casadi_ocp_solver.get(i, "u") for i in range(N_horizon)])
    lam_casadi = np.concatenate([casadi_ocp_solver.get(i, "lam") for i in range(N_horizon+1)])
    pi_casadi = np.concatenate([casadi_ocp_solver.get(i, "pi") for i in range(N_horizon)])
    result_casadi = casadi_ocp_solver.get_iterate()

    # evaluate difference
    diff_x = np.linalg.norm(casadi_x - acados_x)
    assert np.allclose(casadi_x, acados_x, atol=1e-5, rtol=1e-5), f"x mismatch with error {diff_x}"
    print(f"Difference between casadi and acados solution in x: {diff_x}")
    diff_u = np.linalg.norm(casadi_u - acados_u)
    assert np.allclose(casadi_u, acados_u, atol=1e-5, rtol=1e-5), f"u mismatch with error {diff_u}"
    print(f"Difference between casadi and acados solution in u: {diff_u}")
    if bu_cons:
        diff_lam = np.linalg.norm(lam_casadi - lam)
        assert np.allclose(lam_casadi, lam, atol=1e-5, rtol=1e-5), f"lam mismatch with error {diff_lam}"
        print(f"Difference between casadi and acados solution in lam: {diff_lam}")
    else:
        pass # skip lam comparison since lam corresponds to different constraints in the two formulations
    diff_pi = np.linalg.norm(pi_casadi - pi)
    assert np.allclose(pi_casadi, pi, atol=1e-5, rtol=1e-5), f"pi mismatch with error {diff_pi}"
    print(f"Difference between casadi and acados solution in pi: {diff_pi}")

    print("Test passed for bu_cons=", bu_cons)

if __name__ == "__main__":
    main(bu_cons=True)
    main(bu_cons=False)