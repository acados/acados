import sys
import os

import numpy as np
import scipy.linalg
from acados_template import AcadosOcp, AcadosOcpSolver

# import zoro
local_path = os.path.dirname(os.path.abspath(__file__))
zoro_source_dir = os.path.join(local_path, '..')
sys.path.append(zoro_source_dir)
from zoro_description import ZoroDescription, process_zoro_description
from zoro_utils import samplesFromEllipsoid

# same as in normal pendulum model
pendulum_source_dir = os.path.join(local_path, '..', '..', 'pendulum_on_cart', 'common')
sys.path.append(pendulum_source_dir)
from pendulum_model import export_pendulum_ode_model
from utils import plot_pendulum


def main():
    # create ocp object to formulate the OCP
    ocp = AcadosOcp()

    # set model
    model = export_pendulum_ode_model()
    ocp.model = model

    Tf = 1.0
    nx = model.x.size()[0]
    nu = model.u.size()[0]
    ny = nx + nu
    ny_e = nx
    N = 20

    # set dimensions
    ocp.dims.N = N

    # set cost
    Q = 2 * np.diag([1e3, 1e3, 1e-2, 1e-2])
    R = 2 * np.diag([1e-2])

    ocp.cost.W = scipy.linalg.block_diag(Q, R)
    ocp.cost.W_e = Q

    ocp.cost.cost_type = 'LINEAR_LS'
    ocp.cost.cost_type_e = 'LINEAR_LS'

    Vu = np.zeros((ny, nu))
    Vu[4,0] = 1.0
    ocp.cost.Vu = Vu
    ocp.cost.Vx = np.zeros((ny, nx))
    ocp.cost.Vx[:nx,:nx] = np.eye(nx)

    ocp.cost.Vx_e = np.eye(nx)

    ocp.cost.yref  = np.zeros((ny, ))
    ocp.cost.yref_e = np.zeros((ny_e, ))

    # set constraints
    # bound on u
    Fmax = 40
    ocp.constraints.lbu = np.array([-Fmax])
    ocp.constraints.ubu = np.array([+Fmax])
    ocp.constraints.idxbu = np.array([0])

    # bound on x
    theta_min = -np.pi * 0.1
    theta_max = np.pi * 0.3
    ocp.constraints.lbx = np.array([theta_min])
    ocp.constraints.ubx = np.array([theta_max])
    ocp.constraints.idxbx = np.array([1])

    # bound on the terminal state
    ocp.constraints.lbx_e = np.array([theta_min])
    ocp.constraints.ubx_e = np.array([theta_max])
    ocp.constraints.idxbx_e = np.array([1])

    # initial state
    x_init = np.array([0.0, 0.15*np.pi, 0.0, 0.0])
    ocp.constraints.x0 = x_init

    # set options
    ocp.solver_options.qp_solver = 'PARTIAL_CONDENSING_HPIPM' # FULL_CONDENSING_QPOASES
    ocp.solver_options.hessian_approx = 'GAUSS_NEWTON'
    ocp.solver_options.integrator_type = 'ERK'
    ocp.solver_options.print_level = 0
    ocp.solver_options.nlp_solver_type = 'SQP_RTI' # SQP_RTI, SQP

    # set prediction horizon
    ocp.solver_options.tf = Tf

    # custom update: disturbance propagation
    ocp.solver_options.custom_update_filename = 'custom_update_function.c'
    ocp.solver_options.custom_update_header_filename = 'custom_update_function.h'

    ocp.solver_options.custom_update_copy = False
    ocp.solver_options.custom_templates = [
        ('custom_update_function_zoro_template.in.c', 'custom_update_function.c'),
        ('custom_update_function_zoro_template.in.h', 'custom_update_function.h'),
    ]

    # zoro stuff
    zoro_description = ZoroDescription()
    zoro_description.P0_mat = np.zeros((nx, nx))
    zoro_description.fdbk_K_mat = np.array([[0.0, 0.0, 10.0, 10.0]])
    # zoro_description.fdbk_K_mat = np.zeros((nu, nx))
    zoro_description.W_mat = np.diag([5*1e-4, 5*1e-4, 5*1e-3, 5*1e-3])
    zoro_description.idx_lbu_t = [0]
    zoro_description.idx_ubu_t = [0]
    zoro_description.idx_lbx_t = [0]
    zoro_description.idx_lbx_e_t = [0]
    ocp.zoro_description = process_zoro_description(zoro_description)

    ocp_solver = AcadosOcpSolver(ocp, json_file = 'acados_ocp.json')

    Nsim = 100
    simX = np.ndarray((Nsim+1, nx))
    simU = np.ndarray((Nsim, nu))
    simX[0,:] = x_init

    # zoro parameters
    max_zoro_iter = 100
    zoro_tol = 1e-5

    # sample disturbances
    np.random.seed(1)
    dist_samples = samplesFromEllipsoid(Nsim, np.zeros((nx,)), zoro_description.W_mat)

    for idx_sim in range(Nsim):
        residuals = np.inf
        # set initial state
        ocp_solver.set(0, "lbx", simX[idx_sim,:])
        ocp_solver.set(0, "ubx", simX[idx_sim,:])

        for idx_iter in range(max_zoro_iter):
            # constraint tightening
            ocp_solver.custom_update([])
            # call SQP_RTI solver
            status = ocp_solver.solve()

            if status != 0:
                print(f"{simU[idx_sim, :]=}")
                raise Exception(f'acados returned status {status} at idx_sim={idx_sim} for initial state {simX[idx_sim,:]}.')

            residuals = ocp_solver.get_residuals()
            if max(residuals) <= zoro_tol:
                break

        if max(residuals) > zoro_tol:
            print("zoro does not converge")

        # get solution
        simU[idx_sim,:] = ocp_solver.get(0, "u")
        simX[idx_sim+1,:] = ocp_solver.get(1, "x")
        simX[idx_sim+1,:] += dist_samples[idx_sim]

        if idx_sim == 0:
            # print backoffs
            print("backoff in the terminal state constraint=", \
                (theta_max - theta_min) \
                    - (ocp_solver.get_from_qp_in(N, 'ubx') - ocp_solver.get_from_qp_in(N, 'lbx')))
            print("backoff in the N-1 input constraint=", \
                2 * Fmax \
                    - (ocp_solver.get_from_qp_in(N-1, 'ubu') - ocp_solver.get_from_qp_in(N-1, 'lbu')))
            # print costs
            cost = ocp_solver.get_cost()
            print("cost function value of solution = ", cost)

    plot_pendulum(np.linspace(0, Tf, Nsim+1), Fmax, simU, simX, latexify=False)


if __name__ == "__main__":
    main()