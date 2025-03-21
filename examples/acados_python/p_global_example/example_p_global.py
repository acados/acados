#
# Copyright (c) The acados authors.
#
# This file is part of acados.
#
# The 2-Clause BSD License
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice,
# this list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright notice,
# this list of conditions and the following disclaimer in the documentation
# and/or other materials provided with the distribution.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.;
#

from acados_template import AcadosOcp, AcadosOcpSolver, AcadosModel, AcadosMultiphaseOcp, get_simulink_default_opts
import numpy as np
import scipy.linalg
from utils import plot_pendulum

from casadi import MX, vertcat, sin, cos
import casadi as ca
import time

# NOTE: This example requires CasADi version nightly-se2 or later.
# Furthermore, this example uses additional flags for the CasADi code generation,
# cf. the solver option ext_fun_compile_flags, which you might need to adapt based
# on your compiler and operating system.

LARGE_SCALE = False
PLOT = False
np.random.seed(1)

if LARGE_SCALE:
    knots = [np.arange(200),np.arange(200)]
    data = np.random.random((38416,)).ravel(order='F')
else:
    knots = [np.arange(20),np.arange(20)]
    data = 0.1 + 0.*np.random.random((256,)).ravel(order='F')

def create_p_global(lut=True):
    m = MX.sym("m")
    l = MX.sym("l")
    p_global = [m, l]
    p_global_values = np.array([0.1, 0.8])

    if lut:
        # Coefficient of B-spline
        C = MX.sym("C", data.shape[0], 1)
        p_global += [C]
        p_global_values = np.concatenate([p_global_values, data])
    else:
        C = None
    p_global = ca.vcat(p_global)

    return p_global, m, l, C, p_global_values


def export_pendulum_ode_model(p_global, m, l, C, lut=True, blazing=True) -> AcadosModel:
    model_name = f'pendulum_blazing_{blazing}'

    # constants
    m_cart = 1. # mass of the cart [kg]

    # parameters
    g = MX.sym("g")
    p = g

    # set up states & controls
    x1      = MX.sym('x1')
    theta   = MX.sym('theta')
    v1      = MX.sym('v1')
    dtheta  = MX.sym('dtheta')

    x = vertcat(x1, theta, v1, dtheta)

    F = MX.sym('F')
    u = vertcat(F)

    # xdot
    nx = x.shape[0]
    xdot = MX.sym('xdot', nx)

    # dynamics
    cos_theta = cos(theta)
    sin_theta = sin(theta)
    denominator = m_cart + m - m*cos_theta**2
    f_expl = vertcat(v1,
                     dtheta,
                     (-m*l*sin_theta*dtheta**2 + m*g*cos_theta*sin_theta+F)/denominator,
                     (-m*l*cos_theta*sin_theta*dtheta**2 + F*cos_theta+(m_cart+m)*g*sin_theta)/(l*denominator)
                     )

    if lut:
        x_in = ca.vertcat(u/100+0.5,theta/np.pi+0.5)

        if blazing:
            # Disturb the dynamics by a sprinkle of bspline
            spline_fun = ca.blazing_spline('blazing_spline', knots)
            f_expl[3] += 0.01*spline_fun(x_in, C)
        else:
            f_expl[3] += 0.01*ca.bspline(x_in, C, knots, [3, 3], 1)

    f_impl = xdot - f_expl

    model = AcadosModel()

    model.f_impl_expr = f_impl
    model.f_expl_expr = f_expl
    model.x = x
    model.xdot = xdot
    model.u = u
    # model.z = z
    model.p = p
    model.p_global = p_global
    model.name = model_name

    # store meta information
    model.x_labels = ['$x$ [m]', r'$\theta$ [rad]', '$v$ [m]', r'$\dot{\theta}$ [rad/s]']
    model.u_labels = ['$F$']
    model.t_label = '$t$ [s]'

    return model


def create_ocp_formulation_without_opts(p_global, m, l, C, lut=True, use_p_global=True, blazing=True) -> AcadosOcp:

    # create ocp object to formulate the OCP
    ocp = AcadosOcp()

    # set model
    model = export_pendulum_ode_model(p_global, m, l, C, lut=lut, blazing=blazing)
    model.p_global = p_global
    model.name += f'_p_global_{use_p_global}'
    ocp.model = model

    # dimensions
    nx = model.x.rows()
    nu = model.u.rows()
    ny = nx + nu
    ny_e = nx

    # set cost
    Q = np.diag([1e3, 1e3, 1e-2, 1e-2])
    R = np.diag([1e-2])

    ocp.cost.W_e = Q
    ocp.cost.W = scipy.linalg.block_diag(Q, R)

    ocp.cost.cost_type = 'LINEAR_LS'
    ocp.cost.cost_type_e = 'LINEAR_LS'

    ocp.cost.Vx = np.zeros((ny, nx))
    ocp.cost.Vx[:nx,:nx] = np.eye(nx)

    Vu = np.zeros((ny, nu))
    Vu[4,0] = 1.0
    ocp.cost.Vu = Vu

    ocp.cost.Vx_e = np.eye(nx)

    ocp.cost.yref  = np.zeros((ny, ))
    ocp.cost.yref_e = np.zeros((ny_e, ))

    # set constraints
    Fmax = 80
    ocp.constraints.lbu = np.array([-Fmax])
    ocp.constraints.ubu = np.array([+Fmax])
    ocp.constraints.idxbu = np.array([0])

    ocp.constraints.x0 = np.array([0.0, np.pi, 0.0, 0.0])

    ocp.parameter_values = np.array([9.81])

    if not use_p_global:
        model.p = ca.vertcat(model.p, p_global)
        model.p_global = MX.sym('p_global', 0, 1)

    return ocp


def main(use_cython=False, lut=True, use_p_global=True, blazing=True, with_matlab_templates=False, code_export_directory=None):

    print(f"\n\nRunning example with lut={lut}, use_p_global={use_p_global}, {blazing=}")
    p_global, m, l, C, p_global_values = create_p_global(lut=lut)

    # create ocp
    ocp = create_ocp_formulation_without_opts(p_global, m, l, C, lut=lut, use_p_global=use_p_global, blazing=blazing)

    if not use_p_global:
        ocp.parameter_values = np.concatenate([ocp.parameter_values, p_global_values])
    else:
        ocp.p_global_values = p_global_values

    Tf = 1.0
    N_horizon = 20

    # set options
    ocp.solver_options.qp_solver = 'PARTIAL_CONDENSING_HPIPM' # FULL_CONDENSING_QPOASES
    ocp.solver_options.hessian_approx = 'GAUSS_NEWTON'
    ocp.solver_options.integrator_type = 'ERK'
    ocp.solver_options.print_level = 0
    ocp.solver_options.nlp_solver_type = 'SQP_RTI' # SQP_RTI, SQP
    ocp.solver_options.ext_fun_compile_flags += ' -I' + ca.GlobalOptions.getCasadiIncludePath() + ' -ffast-math -march=native'
    if code_export_directory is not None:
        ocp.code_export_directory = code_export_directory

    # set prediction horizon
    ocp.solver_options.tf = Tf
    ocp.solver_options.N_horizon = N_horizon

    # create ocp solver
    print(f"Creating ocp solver with p_global = {ocp.model.p_global}, p = {ocp.model.p}")
    if with_matlab_templates:
        ocp.simulink_opts = get_simulink_default_opts()
    solver_json = 'acados_ocp_' + ocp.model.name + '.json'
    if use_cython:
        AcadosOcpSolver.generate(ocp, json_file=solver_json)
        AcadosOcpSolver.build(ocp.code_export_directory, with_cython=True)
        ocp_solver = AcadosOcpSolver.create_cython_solver(solver_json)
    else:
        ocp_solver = AcadosOcpSolver(ocp, json_file = solver_json, generate=True, build=True)

    # call SQP_RTI solver in the loop:
    residuals = []
    if use_p_global:
        t_start = time.time()
        ocp_solver.set_p_global_and_precompute_dependencies(p_global_values)
        t_elapsed = time.time() - t_start

        print(f"Precompute {t_elapsed}.")

    timing = 0
    for i in range(20):
        status = ocp_solver.solve()
        # ocp_solver.print_statistics() # encapsulates: stat = ocp_solver.get_stats("statistics")
        residuals+= list(ocp_solver.get_residuals())
        timing += ocp_solver.get_stats("time_lin")

    # plot results
    if PLOT:
        u_traj = np.array([ocp_solver.get(i, "u") for i in range(N_horizon)])
        x_traj = np.array([ocp_solver.get(i, "x") for i in range(N_horizon+1)])
        plot_pendulum(ocp.solver_options.shooting_nodes, ocp.constraints.ubu[0], u_traj, x_traj, x_labels=ocp.model.x_labels, u_labels=ocp.model.u_labels)

    return residuals, timing


def main_mocp(lut=True, use_p_global=True, with_matlab_templates=False):
    print(f"\n\nRunning multi-phase example with lut={lut}, use_p_global={use_p_global}")
    p_global, m, l, C, p_global_values = create_p_global(lut=lut)

    Tf = 1.0
    N_horizon = 20

    # create ocp
    n_phases = 2
    mocp = AcadosMultiphaseOcp(N_list=[10, 10])

    ocp_phase_1 = create_ocp_formulation_without_opts(p_global, m, l, C, lut=lut, use_p_global=use_p_global)
    ocp_phase_2 = create_ocp_formulation_without_opts(p_global, m, l, C, lut=lut, use_p_global=use_p_global)

    mocp.set_phase(ocp_phase_1, 0)
    mocp.set_phase(ocp_phase_2, 1)

    if not use_p_global:
        for ip in range(n_phases):
            mocp.parameter_values[ip] = np.concatenate([mocp.parameter_values[ip], p_global_values])
    else:
        mocp.p_global_values = p_global_values

    # set options
    mocp.solver_options.qp_solver = 'PARTIAL_CONDENSING_HPIPM' # FULL_CONDENSING_QPOASES
    mocp.solver_options.hessian_approx = 'GAUSS_NEWTON'
    mocp.solver_options.integrator_type = 'ERK'
    mocp.solver_options.print_level = 0
    mocp.solver_options.nlp_solver_type = 'SQP_RTI' # SQP_RTI, SQP

    if lut:
        # NOTE: these additional flags are required for code generation of CasADi functions using ca.blazing_spline
        mocp.solver_options.ext_fun_compile_flags = '-I' + ca.GlobalOptions.getCasadiIncludePath() + ' -ffast-math -march=native'

    # set prediction horizon
    mocp.solver_options.tf = Tf
    mocp.solver_options.N_horizon = N_horizon

    # create ocp solver
    print(f"Creating ocp solver with p_global = {mocp.model[0].p_global}, p_phase_1 = {mocp.model[0].p}, p_phase_2 = {mocp.model[1].p}")

    if with_matlab_templates:
        mocp.simulink_opts = get_simulink_default_opts()
    ocp_solver = AcadosOcpSolver(mocp, generate=True, build=True)

    # call SQP_RTI solver in the loop:
    residuals = []

    if use_p_global:
        ocp_solver.set_p_global_and_precompute_dependencies(p_global_values)

    timing = 0
    for i in range(20):
        status = ocp_solver.solve()
        # ocp_solver.print_statistics() # encapsulates: stat = ocp_solver.get_stats("statistics")
        residuals+= list(ocp_solver.get_residuals())
        timing += ocp_solver.get_stats('time_lin')

    return residuals, timing


if __name__ == "__main__":

    # OCP with lookuptable, comparing blazing, bspline, p_global
    ref_lut_no_blazing, t_lin_lut_no_blazing_ref = main(use_cython=False, use_p_global=False, lut=True, blazing=False)
    res_lut_no_blazing, t_lin_lut_no_blazing = main(use_cython=False, use_p_global=True, lut=True, blazing=False)

    ref_lut, t_lin_lut_ref = main(use_cython=False, use_p_global=False, lut=True)
    res_lut, t_lin_lut = main(use_cython=False, use_p_global=True, lut=True)

    print(f"\t\t bspline \t blazing")
    print(f"ref\t\t {t_lin_lut_no_blazing_ref:.5f} \t {t_lin_lut_ref:.5f}")
    print(f"p_global\t {t_lin_lut_no_blazing:.5f} \t {t_lin_lut:.5f}")

    np.testing.assert_almost_equal(ref_lut, res_lut)
    np.testing.assert_almost_equal(ref_lut_no_blazing, res_lut_no_blazing)

    np.testing.assert_almost_equal(ref_lut, ref_lut_no_blazing)
    np.testing.assert_almost_equal(res_lut, res_lut_no_blazing)
    np.testing.assert_almost_equal(ref_lut, res_lut_no_blazing)

    ref_nolut, _ = main(use_cython=False, use_p_global=False, lut=False)
    res_nolut, _ = main(use_cython=False, use_p_global=True, lut=False)
    np.testing.assert_almost_equal(ref_nolut, res_nolut)

    # MOCP tests
    res_mocp_nolut_p, _ = main_mocp(use_p_global=False, lut=False)
    res_mocp_nolut_p_global, _ = main_mocp(use_p_global=True, lut=False)
    np.testing.assert_almost_equal(ref_nolut, res_mocp_nolut_p)
    np.testing.assert_almost_equal(ref_nolut, res_mocp_nolut_p_global)

    res_mocp_lut_p, _ = main_mocp(use_p_global=False, lut=True)
    res_mocp_lut_p_global, _ = main_mocp(use_p_global=True, lut=True, with_matlab_templates=True)
    np.testing.assert_almost_equal(ref_lut, res_mocp_lut_p)
    np.testing.assert_almost_equal(ref_lut, res_mocp_lut_p_global)

    with np.testing.assert_raises(Exception):
        np.testing.assert_almost_equal(ref_lut, ref_nolut)

    # to test transfer to MATLAB/Octave
    res_lut, t_lin_lut = main(use_cython=False, use_p_global=True, lut=True, with_matlab_templates=True, code_export_directory='c_generated_code_single_phase')
