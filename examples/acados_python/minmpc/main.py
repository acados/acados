#
# Copyright (c) The acados authors.
#
# This file is part of acados.
#
# Licensed under the 2-Clause BSD License.


"""
Reference tracking of  x_dot = x^3 - u,  u in {0,1}

    min   sum_{k=0}^{N} (x_k - x_ref)^2
    s.t.  x_0 = 0.9
          x_{k+1} = F_RK4(x_k, u_k),   k = 0..N-1
          0 <= u_k <= 1

N = 30, t_s = 0.05 s  ->  T = 1.5 s

Example taken from Ghezzi et al.:
A sequential Benders-based mixed-integer quadratic programming algorithm and its implementation in the CAMINO toolbox

inspired by a textbook example of the MPC book by Rawlings, Mayne, Diehl
"""

import numpy as np
import casadi as ca
import matplotlib.pyplot as plt
from acados_template import AcadosOcp, AcadosOcpSolver, AcadosModel, plot_trajectories, ACADOS_INFTY, AcadosSim, AcadosSimSolver
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator


N = 30
TS = 0.05
T = N * TS
X0 = 0.9
X_REF = 0.7

N_SIM = 51

def create_model() -> AcadosModel:
    x = ca.SX.sym("x")
    u = ca.SX.sym("u")

    f_expl = x**3 - u

    model = AcadosModel()
    model.name = "unstable_tracking_relaxed"
    model.x = x
    model.u = u
    model.f_expl_expr = f_expl

    model.cost_y_expr = x
    model.cost_y_expr_e = x
    return model

def create_plant() -> ca.Function:
    sim = AcadosSim()
    sim.model = create_model()
    sim.solver_options.T = TS
    sim.name = sim.model.name + "_plant"

    return AcadosSimSolver(sim, verbose=False)

def create_solver() -> AcadosOcpSolver:
    ocp = AcadosOcp()
    ocp.name = "nonlinear_unstable_miocp"
    ocp.solver_options.N_horizon = N
    ocp.solver_options.tf = T

    ocp.model = create_model()
    ocp.model.reformulate_with_discrete_dynamics(dt=TS)
    ocp.model.f_expl_expr = []

    # state augmentation for minimum-uptime constraints
    u_prev = ca.SX.sym("u_prev", 1)
    u_prev_prev = ca.SX.sym("u_prev_prev", 1)
    ocp.model.x = ca.vertcat(ocp.model.x, u_prev_prev, u_prev)
    ocp.model.disc_dyn_expr = ca.vertcat(ocp.model.disc_dyn_expr, u_prev, ocp.model.u)

    ocp.model.con_h_expr_0 = ocp.model.u - u_prev + u_prev_prev
    ocp.model.con_h_expr = ocp.model.u - u_prev + u_prev_prev

    ocp.constraints.lh_0 = np.zeros((1,))
    ocp.constraints.lh = np.zeros((1,))
    ocp.constraints.uh_0 = ACADOS_INFTY * np.ones((1,))
    ocp.constraints.uh = ACADOS_INFTY * np.ones((1,))

    ocp.constraints.idxbinu = np.array([0])

    # cost: 0.5 * W * (x - x_ref)^2 with W = 2  ->  (x - x_ref)^2
    ocp.cost.cost_type = "NONLINEAR_LS"
    ocp.cost.cost_type_e = "NONLINEAR_LS"
    ocp.cost.W = np.array([[2.0]])
    ocp.cost.W_e = np.array([[2.0]])
    ocp.cost.yref = np.array([X_REF])
    ocp.cost.yref_e = np.array([X_REF])

    # constraints: initial state plus augmented states
    ocp.constraints.x0 = np.array([X0, 0, 0])
    ocp.constraints.idxbu = np.array([0])
    ocp.constraints.lbu = np.array([0.0])
    ocp.constraints.ubu = np.array([1.0])

    # solver options
    ocp.solver_options.integrator_type = "DISCRETE"

    ocp.solver_options.nlp_solver_type = "SQP_RTI"
    ocp.solver_options.hessian_approx = "GAUSS_NEWTON"

    # u is not penalized in the cost, so the Gauss-Newton Hessian has a zero R-block
    ocp.solver_options.levenberg_marquardt = 1e-6
    ocp.solver_options.qp_solver = "PARTIAL_CONDENSING_DAQP"

    ocp.solver_options.qp_solver_cond_N = N
    ocp.solver_options.qp_solver_iter_max = 100
    ocp.solver_options.cost_scaling = np.ones((N+1,))

    return AcadosOcpSolver(ocp, verbose=False)


def set_uptime_constraint(solver: AcadosOcpSolver, active: bool, N_bin: int = N):
    """Switch the minimum-uptime constraint on / off."""
    lh = np.zeros((1,)) if active else -ACADOS_INFTY * np.ones((1,))
    for n in range(N_bin):
        solver.constraints_set(n, "lh", lh)

    for n in range(N_bin, N):
        solver.constraints_set(n, "lh", -ACADOS_INFTY * np.ones((1,)))


def plot(x_traj_list, u_traj_list, labels_list, t_traj_list, x_ref_traj):
    fig = plot_trajectories(x_traj_list,
                            u_traj_list,
                            labels_list,
                            t_traj_list,
                            X_ref=x_ref_traj,
                            show_legend=False,
                            show_plot=False
                            )
    fig.axes[0].legend(loc="upper left")
    fig.axes[0].xaxis.set_minor_locator(MultipleLocator(TS))
    fig.axes[0].grid(which='minor', linestyle=':', linewidth=0.5, alpha=0.5)
    fig.axes[1].grid(which='minor', linestyle=':', linewidth=0.5, alpha=0.5)
    fig.tight_layout()


def solve_open_loop(solver: AcadosOcpSolver):

    # solve relaxed problem without minimum uptime constraints
    solver.options_set("qp_relax_binaries", True)
    set_uptime_constraint(solver, False)

    num_iter_relaxed = 10
    for _ in range(num_iter_relaxed):
        solver.solve()

    iterate_relaxed = solver.get_iterate()

    # solve MIQP
    solver.options_set("qp_relax_binaries", False)
    set_uptime_constraint(solver, True)

    solver.solve()

    iterate = solver.get_iterate()

    print(f"open-loop cost: {solver.get_cost():.6f} (globally optimal: 0.1765)")

    x_traj = np.array(iterate.x)[:, :1]
    u_traj = np.array(iterate.u)

    x_traj_relaxed = np.array(iterate_relaxed.x)[:, :1]
    u_traj_relaxed = np.array(iterate_relaxed.u)

    plot([x_traj_relaxed, x_traj],
         [u_traj_relaxed, u_traj],
         ["relaxed solution", "solution"],
         2*[np.arange(N + 1) * TS],
         X_REF*np.ones((N+1, 1)),
        )


def run_closed_loop(solver: AcadosOcpSolver):

    plant = create_plant()

    x_ref_traj = np.ones((N_SIM+N+1, 1))
    x_ref_traj[:N+1] = X_REF
    x_ref_traj[N+1:] = X_REF + 0.3

    x_traj_all = []
    u_traj_all = []
    labels_all = []

    for N_bin in [N, N // 2]:
        x_cl = np.zeros((N_SIM + 1, 1))
        u_cl = np.zeros((N_SIM, 1))
        x_cl[0] = X0

        u_prev_prev, u_prev = 0.0, 0.0

        solver.reset()

        # warm start
        x0_aug = np.array([x_cl[0, 0], u_prev_prev, u_prev])

        for n in range(N+1):
            solver.set(n, "yref", x_ref_traj[n])

        solver.options_set("qp_relax_binaries", True)

        set_uptime_constraint(solver, False)

        num_iter_relaxed = 3
        for _ in range(num_iter_relaxed):
            solver.solve_for_x0(x0_aug)

        solver.options_set("qp_binary_horizon", N_bin)

        # closed-loop simulation: solve one MIQP per time step
        solver.options_set("qp_relax_binaries", False)
        set_uptime_constraint(solver, True, N_bin=N_bin)

        for k in range(N_SIM):
            x0_aug = np.array([x_cl[k, 0], u_prev_prev, u_prev])

            for n in range(N+1):
                solver.set(n, "yref", x_ref_traj[k+n])

            u0 = solver.solve_for_x0(x0_aug)

            if not (np.allclose(u0, np.ones((1,))) or np.allclose(u0, np.zeros((1,)))):
                print(f"u not integer at iteration {k}, {u0=}")
                solver.print_statistics()
            u_cl[k] = u0

            # simulate the plant
            x_cl[k + 1] = plant.simulate(x_cl[k], u0)

            u_prev_prev, u_prev = u_prev, u0.item()

        cl_cost = np.sum((x_cl - x_ref_traj[:N_SIM+1]) ** 2)
        print(f"{N_bin =}: closed-loop cost over {N_SIM} steps: {cl_cost:.6f}")

        x_traj_all.append(x_cl)
        u_traj_all.append(u_cl)
        labels_all.append(r"$N_{\mathrm{bin}} =" + f"{N_bin}$")

    plot(x_traj_all, u_traj_all, labels_all,
          [np.arange(N_SIM + 1) * TS]*2,
          x_ref_traj[:N_SIM+1]
    )


def main():
    solver = create_solver()

    solve_open_loop(solver)
    run_closed_loop(solver)

    plt.show()


if __name__ == "__main__":
    main()