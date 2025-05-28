# -*- coding: future_fstrings -*-
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

from utils import plot_furuta_pendulum, plot_time_per_solve
from furuta_common import get_furuta_model, setup_ocp_solver
from integrator_experiment import setup_acados_integrator, IntegratorSetting
import numpy as np

def get_plant_integrator_settings():
    integrator_settings = IntegratorSetting(integrator_type="IRK",
                                            num_stages=2,
                                            num_steps=2,
                                            newton_iter=20,
                                            newton_tol=1e-10,
                                            )
    return integrator_settings

def main(use_RTI=False, timeout_max_time=0., heuristic="ZERO"):

    x0 = np.array([0.0, np.pi, 0.0, 0.0])
    umax = .45

    Tf = .350       # total prediction time
    N_horizon = 8   # number of shooting intervals
    dt_0 = 0.025    # sampling time = length of first shooting interval

    ocp_solver = setup_ocp_solver(x0, umax, dt_0, N_horizon, Tf, use_RTI, timeout_max_time, heuristic)
    # setup plant simulator
    integrator_settings = get_plant_integrator_settings()
    model = get_furuta_model()
    integrator = setup_acados_integrator(model, dt_0, integrator_settings)

    nx = ocp_solver.acados_ocp.dims.nx
    nu = ocp_solver.acados_ocp.dims.nu

    Nsim = 50
    simX = np.zeros((Nsim+1, nx))
    simU = np.zeros((Nsim, nu))

    simX[0,:] = x0

    if use_RTI:
        t_preparation = np.zeros((Nsim))
        t_feedback = np.zeros((Nsim))

    else:
        t = np.zeros((Nsim))

    # do some initial iterations to start with a good initial guess
    num_iter_initial = 10
    for _ in range(num_iter_initial):
        ocp_solver.solve_for_x0(x0_bar = x0, fail_on_nonzero_status=False)

    # closed loop
    for i in range(Nsim):

        if use_RTI:
            # preparation phase
            ocp_solver.options_set('rti_phase', 1)
            status = ocp_solver.solve()
            t_preparation[i] = ocp_solver.get_stats('time_tot')

            # set initial state
            ocp_solver.set(0, "lbx", simX[i, :])
            ocp_solver.set(0, "ubx", simX[i, :])

            # feedback phase
            ocp_solver.options_set('rti_phase', 2)
            status = ocp_solver.solve()
            t_feedback[i] = ocp_solver.get_stats('time_tot')

            simU[i, :] = ocp_solver.get(0, "u")

        else:
            # solve ocp and get next control input
            simU[i,:] = ocp_solver.solve_for_x0(x0_bar = simX[i, :], fail_on_nonzero_status=False)
            t[i] = ocp_solver.get_stats('time_tot')

        # simulate system
        simX[i+1, :] = integrator.simulate(x=simX[i, :], u=simU[i,:])

    # evaluate timings
    if use_RTI:
        # scale to milliseconds
        t_preparation *= 1000
        t_feedback *= 1000
        print(f'Computation time in preparation phase in ms: \
                min {np.min(t_preparation):.3f} median {np.median(t_preparation):.3f} max {np.max(t_preparation):.3f}')
        print(f'Computation time in feedback phase in ms:    \
                min {np.min(t_feedback):.3f} median {np.median(t_feedback):.3f} max {np.max(t_feedback):.3f}')
    else:
        # scale to milliseconds
        t *= 1000
        plot_time_per_solve(t, timeout_max_time*1000, heuristic, plt_show=False, store_figure=False)
        print(f'Computation time in ms: min {np.min(t):.3f} median {np.median(t):.3f} max {np.max(t):.3f}')

    # plot results
    plot_furuta_pendulum(np.linspace(0, (Tf/N_horizon)*Nsim, Nsim+1),simX, simU, umax, plt_show=True)

    ocp_solver = None


if __name__ == '__main__':
    main(use_RTI=False, timeout_max_time=0.)
    # main(use_RTI=False, timeout_max_time=1*1e-3, heuristic="ZERO")
    # main(use_RTI=False, timeout_max_time=1*1e-3, heuristic="LAST")
    # main(use_RTI=False, timeout_max_time=1*1e-3, heuristic="MAX_CALL")
    # main(use_RTI=False, timeout_max_time=1*1e-3, heuristic="MAX_OVERALL")

    # main(use_RTI=True) # timeout not implemented for RTI

