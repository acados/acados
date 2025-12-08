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

import numpy as np
from acados_template import AcadosSim, AcadosSimSolver, AcadosModel
from casadi import SX, vertcat, sin, cos

def export_pendulum_model_with_M_param() -> AcadosModel:
    model_name = 'pendulum_with_param_M'

    # Parameters
    M = SX.sym('M')  # mass of the cart [kg]
    m = 0.1          # mass of the ball [kg]
    l = 0.8          # length of the rod [m]
    g = 9.81         # gravity constant [m/s^2]

    # States
    p_pos   = SX.sym('p_pos')
    theta   = SX.sym('theta')
    v       = SX.sym('v')
    dtheta  = SX.sym('dtheta')

    x = vertcat(p_pos, theta, v, dtheta)

    # Controls
    F = SX.sym('F')
    u = vertcat(F)

    # Dynamics
    sin_theta = sin(theta)
    cos_theta = cos(theta)
    denominator = M + m - m * cos_theta**2
    
    f_expl = vertcat(
        v,
        dtheta,
        (-l*m*sin_theta*dtheta**2 + F + g*m*cos_theta*sin_theta) / denominator,
        (-l*m*cos_theta*sin_theta*dtheta**2 + F*cos_theta + g*m*sin_theta + M*g*sin_theta) / (l*denominator)
    )

    xdot = SX.sym('xdot', 4, 1)
    f_impl = xdot - f_expl

    model = AcadosModel()
    model.f_impl_expr = f_impl
    model.f_expl_expr = f_expl
    model.x = x
    model.xdot = xdot
    model.u = u
    model.p = vertcat(M)
    model.name = model_name

    return model


def main():
    # Setup
    model = export_pendulum_model_with_M_param()
    
    T_s = 0.1
    nx = model.x.size()[0]
    np_param = model.p.size()[0]

    # Inputs
    x0 = np.array([1e-1, 1e0, 2e-1, 2e0])
    u0 = np.array([0.0])
    p0 = np.array([1.0])

    # Simulation config
    sim = AcadosSim()
    sim.model = model
    sim.solver_options.T = T_s
    sim.solver_options.num_stages = 3
    sim.solver_options.integrator_type = 'ERK'
    sim.solver_options.collocation_type = 'GAUSS_RADAU_IIA'
    
    # Enable sensitivities
    sim.solver_options.sens_forw = True
    sim.solver_options.sens_forw_p = True

    sim.parameter_values = np.array([1.0])

    # Initialize solver
    acados_integrator = AcadosSimSolver(sim)

    # Solve nominal
    acados_integrator.set('x', x0)
    acados_integrator.set('u', u0)
    acados_integrator.set('p', p0)

    if acados_integrator.solve() != 0:
        raise Exception('acados returned error status.')

    S_p_solver = acados_integrator.get('S_p')

    # Finite Differences Verification
    FD_epsilon = 1e-6
    S_p_fd = np.zeros((nx, np_param))

    for jj in range(np_param):
        p_pert = p0.copy()
        p_pert[jj] += FD_epsilon

        acados_integrator.set('x', x0)
        acados_integrator.set('u', u0)
        acados_integrator.set('p', p_pert)

        acados_integrator.solve()
        xn_tmp = acados_integrator.get('x')

        xn_nom = S_p_solver  # Not used in FD, just placeholder
        # Re-solve nominal for diff (optional if we stored xn_nom earlier, 
        # but cleaner to just use the solver loop logic if strictly FD)
        # To strictly compute FD, we need xn_nom.
        
    # Re-compute nominal x for FD check
    acados_integrator.set('x', x0)
    acados_integrator.set('u', u0)
    acados_integrator.set('p', p0)
    acados_integrator.solve()
    xn_nom = acados_integrator.get('x')
    
    # Compute FD
    for jj in range(np_param):
        p_pert = p0.copy()
        p_pert[jj] += FD_epsilon
        acados_integrator.set('p', p_pert)
        acados_integrator.solve()
        xn_pert = acados_integrator.get('x')
        S_p_fd[:, jj] = (xn_pert - xn_nom) / FD_epsilon

    # Validation
    error_abs_Sp = np.max(np.abs(S_p_fd - S_p_solver))
    print(f"Max Error S_p (FD vs Analytic): {error_abs_Sp:.2e}")

    if error_abs_Sp > 1e-5:
        raise Exception("Failure: parameter sensitivity error too large.")
    
    print("Success: sens_forw_p matches finite differences.")

if __name__ == "__main__":
    main()