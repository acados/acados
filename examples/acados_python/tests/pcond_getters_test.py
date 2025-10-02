import os
import numpy as np
import casadi as ca
from acados_template import AcadosModel, AcadosOcp, AcadosOcpSolver, latexify_plot
import matplotlib.pyplot as plt


def plot_qp_sparsity(qp, fig_filename=None, title=None, with_legend=True):
    latexify_plot()

    plt.figure(figsize=(4.0, 4.0))
    nx_total = 0
    nu_total = 0
    for stage in qp["stages"]:
        # cost Hessian
        nx = stage["Q"].shape[0]
        nu = stage["R"].shape[0]
        nx_total += nx
        nu_total += nu
    nv_total = nx_total + nu_total

    H_x = np.zeros((nv_total, nv_total))
    H_u = np.zeros((nv_total, nv_total))
    H_xu = np.zeros((nv_total, nv_total))

    offset = 0
    for stage in qp["stages"]:
        nx = stage["Q"].shape[0]
        nu = stage["R"].shape[0]

        H_x[offset:offset + nx, offset:offset + nx] = stage["Q"]
        H_u[offset + nx:offset + nx + nu, offset + nx:offset + nx + nu] = stage["R"]
        if nx > 0 and nu > 0:
            H_xu[offset+nx:offset + nx + nu, offset:offset + nx] = stage["S"]
            H_xu[offset:offset + nx, offset+nx:offset + nx + nu] = stage["S"].T
        offset += nx + nu
    markersize = 5 * (40 / nv_total)
    plt.spy(H_xu, markersize=markersize, label='$S$', color='C2')
    plt.spy(H_x, markersize=markersize, label='$Q$', color='C0')
    plt.spy(H_u, markersize=markersize, label='$R$', color='C1')
    # remove xticks
    plt.xticks([])
    if with_legend:
        plt.legend(handletextpad=0.3)
    if title is not None:
        plt.title(title)
    if fig_filename is not None:
        plt.savefig(fig_filename, bbox_inches='tight')
    plt.show()

def main(N=20, cond_N=10, qp_solver_cond_block_size=None, x0_elimination=True, fig_title=None, with_legend=True, fig_filename=None):
    N = 20
    Tf = 10.0

    # Problem data
    nx, nu = 5, 2
    Q = np.diag([float(i+1) for i in range(nx)])
    R = np.diag([float(i+1) for i in range(nu)])
    x0 = np.zeros(nx)
    x_ref = np.ones(nx)
    lbu = -np.ones(nu)
    ubu = -lbu

    # Build OCP
    ocp = AcadosOcp()
    ocp.solver_options.tf = Tf
    ocp.solver_options.N_horizon = N

    ocp.model = AcadosModel()

    A = np.eye(nx) + 0.1 * np.ones((nx, nx))
    B = np.ones((nx, nu))
    u = ca.SX.sym('u', nu, 1)
    x = ca.SX.sym('x', nx, 1)
    f_expl_expr = A @ x + B @ u

    ocp.model.x = x
    ocp.model.u = u
    ocp.model.name = 'model'
    ocp.model.f_expl_expr = f_expl_expr

    # Bounds
    if x0_elimination:
        ocp.constraints.x0 = x0
    else:
        ocp.constraints.lbx_0 = x0
        ocp.constraints.ubx_0 = x0
        ocp.constraints.idxbx_0 = np.arange(nx)
    ocp.constraints.idxbu = np.arange(nu)
    ocp.constraints.lbu = lbu
    ocp.constraints.ubu = ubu

    # Linear least-squares cost
    ny = nx + nu
    W = np.zeros((ny, ny))
    W[:nx, :nx] = Q
    W[nx:, nx:] = R
    S = 1e-3 * np.ones((nx, nu))
    W[:nx, nx:] = S
    W[nx:, :nx] = S.T

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

    if qp_solver_cond_block_size is not None:
        ocp.solver_options.qp_solver_cond_block_size = qp_solver_cond_block_size

    uniq = f"double_integrator_{os.getpid()}"
    ocp.code_export_directory = f"c_generated_code/{uniq}"
    json_path = f".solver_files/ACADOS_SOLVER_{uniq}.json"

    os.makedirs(os.path.dirname(ocp.code_export_directory), exist_ok=True)
    os.makedirs(os.path.dirname(json_path), exist_ok=True)

    solver = AcadosOcpSolver(ocp, json_file=json_path)

    solver.solve()
    cond_N = ocp.solver_options.qp_solver_cond_N

    dyn_fields = ["A", "B", "b"]
    cost_constr_fields = ["Q", "R", "S", "q", "r", "C", "D", "lg", "ug", "lbx", "ubx", "lbu", "ubu"]
    fields = dyn_fields + cost_constr_fields
    qp = {"N": cond_N, "stages": []}
    for k in range(cond_N + 1):
        stage = {}
        for field in fields:
            if not (k == cond_N and field in dyn_fields):
                arr = solver.get_from_qp_in(k, "pcond_" + field)
                stage[field] = np.array(arr)
        qp["stages"].append(stage)

    plot_qp_sparsity(qp, fig_filename=fig_filename, title=fig_title, with_legend=with_legend)

if __name__ == "__main__":
    # main(N=20, cond_N=20, x0_elimination=False, fig_title="Original QP $N=20$ -- no condensing", fig_filename="sparsity_no_condensing.pdf")
    # main(N=20, cond_N=1, with_legend=False, fig_title="Condensing, $N_{\mathrm{cond}}=1$ (fully condensed)", fig_filename="sparsity_fcond.pdf")
    # main(N=20, cond_N=5, with_legend=False, fig_title="Partial condensing, $N_{\mathrm{cond}}=5$", fig_filename="sparsity_pcond.pdf")
    # main(N=20, cond_N=5, qp_solver_cond_block_size=[6, 5, 4, 3, 2, 0], with_legend=False, fig_title="Condensing, $N_{\mathrm{cond}}=5$, custom block sizes", fig_filename="sparsity_pcond_custom_block_sizes.pdf")

    main(N=20, cond_N=5, qp_solver_cond_block_size=[6, 5, 4, 2, 2, 1])
