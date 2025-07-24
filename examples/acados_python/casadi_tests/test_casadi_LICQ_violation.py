import sys
sys.path.insert(0, '../getting_started')
import numpy as np
import casadi as ca
from typing import Union

from acados_template import AcadosOcp, AcadosOcpSolver, AcadosCasadiOcpSolver
from pendulum_model import export_pendulum_ode_model

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
    Fmax = 80
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

def main():
    N_horizon = 3
    Tf = 1.0
    ocp = formulate_ocp(Tf, N_horizon)

    initial_iterate = ocp.create_default_initial_iterate()

    ## solve using acados
    # create acados solver
    ocp_solver = AcadosOcpSolver(ocp,verbose=False)
    ocp_solver.load_iterate_from_obj(initial_iterate)
    # solve with acados
    status = ocp_solver.solve()
    # get solution
    result = ocp_solver.store_iterate_to_obj()

    # ## solve using casadi
    casadi_ocp_solver = AcadosCasadiOcpSolver(ocp=ocp,solver="ipopt",verbose=False)
    casadi_ocp_solver.load_iterate_from_obj(result)
    casadi_ocp_solver.solve()
    licq = casadi_ocp_solver.satisfies_LICQ()

    # Check for violation of specific stage, at least one stage should violate LICQ
    if licq:
        raise ValueError("LICQ condition is not violated in the solution.")

if __name__ == "__main__":
    main()