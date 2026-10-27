import sys
sys.path.insert(0, '../getting_started')
import numpy as np
import casadi as ca
from typing import Union

from acados_template import AcadosOcp, AcadosModel, AcadosOcpSolver, AcadosCasadiOcpSolver
from pendulum_model import export_pendulum_ode_model

Fmax = 80

def formulate_ocp(Tf: float = 1.0, N: int = 20)-> AcadosOcp:
    ocp = AcadosOcp()

    # set model
    model = export_pendulum_ode_model()
    ocp.model = model

    nx = model.x.rows()
    nu = model.u.rows()
    
    # set prediction horizon
    ocp.solver_options.N_horizon = N
    ocp.solver_options.tf = Tf

    # cost matrices
    Q_mat = 2*np.diag([1e3, 1e3, 1e-2, 1e-2])
    R_mat = 2*np.diag([1e-2])

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
    ocp.constraints.lbu = np.array([-Fmax])
    ocp.constraints.ubu = np.array([+Fmax])
    ocp.constraints.idxbu = np.array([0])

    # set initial bounds and state
    ocp.constraints.lbx_0 = np.array([0, np.pi, -0.2, 0])
    ocp.constraints.ubx_0 = np.array([0, np.pi, 0.2, 0])
    ocp.constraints.idxbx_0 = np.array([0, 1, 2, 3])

    # set linear constraints 
    ocp.constraints.C = np.array([[0, 0, 0, 0]])
    ocp.constraints.D = np.array([[0.1]])
    ocp.constraints.lg = np.array([-8])
    ocp.constraints.ug = np.array([8])

    # set x_1 at the end of the horizon
    ocp.constraints.C_e = np.array([[1, 0, 0, 0]])
    ocp.constraints.lg_e = np.array([0.3])
    ocp.constraints.ug_e = np.array([0.5])

    # set options
    ocp.solver_options.qp_solver = 'PARTIAL_CONDENSING_HPIPM' # FULL_CONDENSING_QPOASES
    ocp.solver_options.hessian_approx = 'GAUSS_NEWTON' # 'GAUSS_NEWTON', 'EXACT'
    ocp.solver_options.integrator_type = 'ERK'
    ocp.solver_options.nlp_solver_type = 'SQP' # SQP_RTI, SQP
    ocp.solver_options.globalization = 'MERIT_BACKTRACKING' # turns on globalization
    
    return ocp

def formulate_scalar_ocp(x0: float, W_x: float, W_e: float) -> AcadosOcp:
    model = AcadosModel()
    model.name = 'scalar_ocp'
    x = ca.SX.sym('x')
    u = ca.SX.sym('u')
    p = ca.SX.sym('p')
    model.x = x
    model.u = u
    model.p = p
    model.disc_dyn_expr = x + p * u

    ocp = AcadosOcp()
    ocp.model = model
    ocp.parameter_values = np.array([1.0])
    ocp.solver_options.N_horizon = 2
    ocp.solver_options.tf = 2.0
    ocp.solver_options.integrator_type = 'DISCRETE'

    ocp.cost.cost_type = 'NONLINEAR_LS'
    ocp.model.cost_y_expr = ca.vertcat(x, u)
    ocp.cost.W = np.diag([W_x, 1e-3])
    ocp.cost.yref = np.array([10.0, 0.0])
    ocp.cost.cost_type_e = 'NONLINEAR_LS'
    ocp.model.cost_y_expr_e = x
    ocp.cost.W_e = np.array([[W_e]])
    ocp.cost.yref_e = np.array([-10.0])

    ocp.constraints.x0 = np.array([x0])
    ocp.constraints.idxbx = np.array([0])
    ocp.constraints.lbx = np.array([-1.0])
    ocp.constraints.ubx = np.array([1.0])
    ocp.constraints.idxbu = np.array([0])
    ocp.constraints.lbu = np.array([-1.0])
    ocp.constraints.ubu = np.array([1.0])
    return ocp

def test_LICQ_across_stages():
    #Constraint x1 = x0 + u0

    # Case 1: x0 = 0.5. At optimum, x1 = 1 (upper bound) and u1 = -1 (lower bound). LICQ holds.
    ocp=formulate_scalar_ocp(x0=0.5, W_x=100., W_e=1.)
    casadi_ocp_solver = AcadosCasadiOcpSolver(ocp, verbose=False)
    casadi_ocp_solver.solve()
    assert np.allclose([casadi_ocp_solver.get(1, 'x')[0], casadi_ocp_solver.get(1, 'u')[0]], [1.0, -1.0], atol=1e-6)
    if not casadi_ocp_solver.satisfies_LICQ():
        raise ValueError("LICQ should hold if state and control bounds are active at the same stage.")

    # Case 2: x0 = 0, We = 0. At optimum, x1 = 1 and u0 = 1, both at upper bounds. 
    # Jacobian of the active constraints goes:
    #                      u0  x1  u1  x2
    #   u0 <= 1          [ 1   0   0   0 ]
    #   x1 <= 1          [ 0   1   0   0 ]
    #   x1 - x0 - u0 = 0 [-1   1   0   0 ]
    #   x2 - x1 - u1 = 0 [ 0  -1  -1   1 ]
    # LICQ fails.
    ocp=formulate_scalar_ocp(x0=0.0, W_x=1., W_e=0.)
    casadi_ocp_solver = AcadosCasadiOcpSolver(ocp, verbose=False)
    casadi_ocp_solver.solve()
    assert np.allclose([casadi_ocp_solver.get(0, 'u')[0], casadi_ocp_solver.get(1, 'x')[0]], [1.0, 1.0], atol=1e-6)
    if casadi_ocp_solver.satisfies_LICQ():
        raise ValueError("LICQ should be violated by active constraints coupled across stages.")

def main():
    N_horizon = 3
    Tf = 1.0
    ocp = formulate_ocp(Tf, N_horizon)

    initial_iterate = ocp.create_default_initial_iterate()

    ## solve using acados
    # create acados solver
    ocp_solver = AcadosOcpSolver(ocp,verbose=False)
    ocp_solver.set_iterate(initial_iterate)
    # solve with acados
    status = ocp_solver.solve()
    # get solution
    result = ocp_solver.get_iterate()

    # ## solve using casadi
    casadi_ocp_solver = AcadosCasadiOcpSolver(ocp=ocp,solver="ipopt",verbose=False)
    casadi_ocp_solver.set_iterate(result)
    casadi_ocp_solver.solve()

    assert np.isclose(casadi_ocp_solver.get(2, 'u')[0], Fmax, atol=1e-6)
    if casadi_ocp_solver.satisfies_LICQ():
        raise ValueError("LICQ should be violated by the redundant bound and linear constraint on u.")

if __name__ == "__main__":
    main()
    test_LICQ_across_stages()