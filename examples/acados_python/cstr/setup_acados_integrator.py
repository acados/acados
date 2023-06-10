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

# authors: Katrin Baumgaertner, Jonathan Frey

from acados_template import AcadosSim, AcadosSimSolver
import numpy as np


def setup_acados_integrator(
    model,
    dt,
    cstr_param,
    sensitivity_propagation=False,
    num_stages=4,
    num_steps=100,
    integrator_type="ERK",
):

    sim = AcadosSim()

    # set model
    sim.model = model

    # set simulation time
    sim.solver_options.T = dt

    ## set options
    sim.solver_options.integrator_type = integrator_type
    sim.solver_options.num_stages = num_stages
    sim.solver_options.num_steps = num_steps
    # for implicit integrator
    sim.solver_options.newton_iter = 10
    sim.solver_options.newton_tol = 1e-8
    sim.solver_options.collocation_type = "GAUSS_LEGENDRE"
    # sensitivity_propagation
    sim.solver_options.sens_adj = sensitivity_propagation
    sim.solver_options.sens_forw = sensitivity_propagation
    sim.solver_options.sens_hess = sensitivity_propagation

    # nominal parameter values
    sim.parameter_values = np.array([cstr_param.F0])

    # create
    acados_integrator = AcadosSimSolver(sim)

    return acados_integrator


if __name__ == "__main__":
    setup_acados_integrator()
