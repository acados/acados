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

from dataclasses import dataclass, field, fields, is_dataclass, InitVar
from typing import Any, Dict
import warnings

NONSUPPORTED_MOCP_INPUTS = ('y_ref', 'lg', 'ug', 'cost_W_0', 'cost_W', 'cost_W_e')
DEFAULT_OFF_MOCP_INPUTS = ('lbx', 'ubx', 'lbx_e', 'ubx_e', 'lh', 'uh', 'y_ref_0', 'y_ref_e')

@dataclass
class AcadosOcpSimulinkInputs:
    """
    Class to toggle which optional input ports are generated on the acados Simulink S-function block.
    All fields must be set to `0` or `1`: a value of `1` enables the corresponding input port, a value of `0` disables it.

    This is the Python equivalent of the MATLAB class ``AcadosOcpSimulinkInputs``.

    Attributes:
        lbx_0: Lower bound on x for stage 0.
            Shape: `nbx_0`. MOCP support: Yes.
        ubx_0: Upper bound on x for stage 0.
            Shape: `nbx_0`. MOCP support: Yes.
        parameter_traj: Parameters - concatenated for all stages 0 to N.
            Shape: `sum(np_i), i=0,..., N`. MOCP support: Yes.
        p_global: Global parameters - first value indicates if update should be performed (0 = no update), followed by new numerical values of `p_global`.
            Shape: `1 + np_global`. MOCP support: Yes.
        y_ref_0: Reference `y_ref` at stage 0.
            Shape: `ny_0`. MOCP support: Yes.
        y_ref: `y_ref` concatenated for stages 1 to N-1.
            Shape: `(N-1) * ny`. MOCP support: No.
        y_ref_e: Reference `y_ref` at stage N.
            Shape: `ny_e`. MOCP support: Yes.
        lbx: Lower bound x values concatenated for stages 1 to N-1.
            Shape: `sum(nbx_i), for i = 1,..., N-1`. MOCP support: Yes.
        ubx: Upper bound x values concatenated for stages 1 to N-1.
            Shape: `sum(nbx_i), for i = 1,..., N-1`. MOCP support: Yes.
        lbx_e: Lower bound x at shooting node N.
            Shape: `nbx_e`. MOCP support: Yes.
        ubx_e: Upper bound x at shooting node N.
            Shape: `nbx_e`. MOCP support: Yes.
        lbu: Lower bound u values concatenated for stages 0 to N-1.
            Shape: `sum(nbu_i), for i = 0,..., N-1`. MOCP support: Yes.
        ubu: Upper bound u values concatenated for stages 0 to N-1.
            Shape: `sum(nbu_i), for i = 0,..., N-1`. MOCP support: Yes.
        lg: Lower bound g values concatenated for stages 0 to N-1.
            Shape: `sum(ng_i), for i = 0,..., N-1`. MOCP support: No.
        ug: Upper bound g values concatenated for stages 0 to N-1.
            Shape: `sum(ng_i), for i = 0,..., N-1`. MOCP support: No.
        lh: Lower bound h values concatenated for stages 1 to N-1.
            Shape: `sum(nh_i), for i = 1,..., N-1`. MOCP support: Yes.
        uh: Upper bound h values concatenated for stages 1 to N-1.
            Shape: `sum(nh_i), for i = 1,..., N-1`. MOCP support: Yes.
        lh_0: Lower bound h at stage 0.
            Shape: `nh_0`. MOCP support: Yes.
        uh_0: Upper bound h at stage 0.
            Shape: `nh_0`. MOCP support: Yes.
        lh_e: Lower bound h at stage N.
            Shape: `nh_e`. MOCP support: Yes.
        uh_e: Upper bound h at stage N.
            Shape: `nh_e`. MOCP support: Yes.
        cost_W_0: Cost matrix `W_0` in column-major format.
            Shape: `ny_0 * ny_0`. MOCP support: No.
        cost_W: Cost matrix `W` in column-major format.
            Shape: `ny * ny`. MOCP support: No.
        cost_W_e: Cost matrix `W_e` in column-major format.
            Shape: `ny_e * ny_e`. MOCP support: No.
        cost_zl: Cost `zl` for all nodes 0 to N.
            Shape: `sum(ns_i) for i = 0,..., N`. MOCP support: Yes.
        cost_zu: Cost `zu` for all nodes 0 to N.
            Shape: `sum(ns_i) for i = 0,..., N`. MOCP support: Yes.
        cost_Zl: Cost `Zl` for all nodes 0 to N.
            Shape: `sum(ns_i) for i = 0,..., N`. MOCP support: Yes.
        cost_Zu: Cost `Zu` for all nodes 0 to N.
            Shape: `sum(ns_i) for i = 0,..., N`. MOCP support: Yes.
        reset_solver: Determines if the solver's iterate is set to all zeros before other initializations.
            Shape: `1`. MOCP support: Yes.
        reset_flags: Additional flags for the solver reset, `[reset_qp_solver, reset_numerical_values, reset_solver_options, reset_x_to_x0_bar]`.
            Shape: `4`. MOCP support: Yes.
        ignore_inits: Determines if initialization (`x_init`, `u_init`, `pi_init`, `slacks_init`) is set (0) or ignored (1).
            Shape: `1`. MOCP support: Yes.
        x_init: Initialization of x for all stages.
            Shape: `sum(nx_i), i=0,..., N`. MOCP support: Yes.
        u_init: Initialization of u for stages 0 to N-1.
            Shape: `sum(nu_i), i=0,..., N-1`. MOCP support: Yes.
        pi_init: Initialization of pi for stages 0 to N-1.
            Shape: `sum(npi_i), i=0,..., N-1`. MOCP support: Yes.
        slacks_init: Initialization of slack values for all stages (0 to N).
            Shape: `2 * ns_total`. MOCP support: Yes.
        rti_phase: Real-time iteration phase.
            Shape: `1`. MOCP support: Yes.
        levenberg_marquardt: Factor for LM regularization.
            Shape: `1`. MOCP support: Yes.
        zoRO_payload: vector forwarded to custom_update_function. Only supported for zoRO. Shape depends on zoRO configuration.
    """
    lbx_0: int = 1
    ubx_0: int = 1
    parameter_traj: int = 1
    p_global: int = 0
    y_ref_0: int = 1
    y_ref: int = 1
    y_ref_e: int = 1
    lbx: int = 1
    ubx: int = 1
    lbx_e: int = 1
    ubx_e: int = 1
    lbu: int = 1
    ubu: int = 1
    lg: int = 1
    ug: int = 1
    lh: int = 1
    uh: int = 1
    lh_0: int = 1
    uh_0: int = 1
    lh_e: int = 1
    uh_e: int = 1
    cost_W_0: int = 0
    cost_W: int = 0
    cost_W_e: int = 0
    cost_zl: int = 0
    cost_zu: int = 0
    cost_Zl: int = 0
    cost_Zu: int = 0
    reset_solver: int = 0
    reset_flags: int = 0
    ignore_inits: int = 0
    x_init: int = 0
    u_init: int = 0
    pi_init: int = 0
    slacks_init: int = 0
    rti_phase: int = 0
    levenberg_marquardt: int = 0
    zoRO_payload: int = 0

    def to_dict(self) -> Dict[str, Any]:
        return {f.name: getattr(self, f.name) for f in fields(self)}

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "AcadosOcpSimulinkInputs":
        obj = cls()
        for key, value in data.items():
            if hasattr(obj, key):
                setattr(obj, key, value)
            else:
                warnings.warn(
                    f"Could not assign field {key} in AcadosOcpSimulinkInputs.from_dict()",
                    UserWarning,
                    stacklevel=2,
                )
        return obj


@dataclass
class AcadosOcpSimulinkOutputs:
    """
    Class to toggle which optional output ports are generated on the acados Simulink S-function block.
    All fields must be set to `0` or `1`: a value of `1` enables the corresponding output port, a value of `0` disables it.

    This is the Python equivalent of the MATLAB class ``AcadosOcpSimulinkOutputs``.

    Attributes:
        u0: Control input at node 0.
            Shape: `nu_0`. MOCP support: Yes.
        utraj: Control input concatenated for nodes 0 to N-1.
            Shape: `sum(nu_i), i=0,..., N-1`. MOCP support: Yes.
        xtraj: State concatenated for nodes 0 to N.
            Shape: `sum(nx_i), i=0,..., N`. MOCP support: Yes.
        ztraj: Algebraic states concatenated for nodes 0 to N-1.
            Shape: `sum(nz_i), i=0,..., N-1`. MOCP support: Yes.
        pi_all: Equality Lagrange multipliers concatenated for nodes 0 to N-1.
            Shape: `sum(npi_i), i=0,..., N-1`. MOCP support: Yes.
        slack_values: Slack values concatenated in order [sl_0, su_0, ..., sl_N, su_N].
            Shape: `2 * ns_total`. MOCP support: Yes.
        solver_status: Acados solver status (0 = SUCCESS).
            Shape: `1`. MOCP support: Yes.
        cost_value: Cost function value.
            Shape: `1`. MOCP support: Yes.
        KKT_residual: KKT residual.
            Shape: `1`. MOCP support: Yes.
        KKT_residuals: KKT residuals, size [4] (stat, eq, ineq, comp).
            Shape: `4`. MOCP support: Yes.
        x1: State at node 1.
            Shape: `nx_1`. MOCP support: Yes.
        CPU_time: CPU time.
            Shape: `1`. MOCP support: Yes.
        CPU_time_sim: CPU time for integrator.
            Shape: `1`. MOCP support: Yes.
        CPU_time_qp: CPU time for QP solution.
            Shape: `1`. MOCP support: Yes.
        CPU_time_lin: CPU time for linearization (including integrator).
            Shape: `1`. MOCP support: Yes.
        sqp_iter: NLP solver iterations.
            Shape: `1`. MOCP support: Yes.
        parameter_traj: Parameter trajectory.
            Shape: `sum(np_i), i=0,..., N`. MOCP support: Yes.
        zoRO_Pk_matrices: P_k matrices within zoRO, only supported if zoRO custom update is used.
            Shape: `sum(nx_i * nx_i), i=0,...,N`. MOCP support: No.
        zoRO_K_matrices: K matrices within zoRO, only supported if zoRO custom update is used.
            Shape: `sum(nx_i * nu_i), i=0,...,N-1`. MOCP support: No.
    """
    u0: int = 1
    utraj: int = 0
    xtraj: int = 0
    ztraj: int = 0
    pi_all: int = 0
    slack_values: int = 0
    solver_status: int = 1
    cost_value: int = 0
    KKT_residual: int = 1
    KKT_residuals: int = 0
    x1: int = 1
    CPU_time: int = 1
    CPU_time_sim: int = 0
    CPU_time_qp: int = 0
    CPU_time_lin: int = 0
    sqp_iter: int = 1
    parameter_traj: int = 0
    zoRO_Pk_matrices: int = 0
    zoRO_K_matrices: int = 0

    def to_dict(self) -> Dict[str, Any]:
        return {f.name: getattr(self, f.name) for f in fields(self)}

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "AcadosOcpSimulinkOutputs":
        obj = cls()
        for key, value in data.items():
            if hasattr(obj, key):
                setattr(obj, key, value)
            else:
                warnings.warn(
                    f"Could not assign field {key} in AcadosOcpSimulinkOutputs.from_dict()",
                    UserWarning,
                    stacklevel=2,
                )
        return obj


@dataclass
class AcadosOcpSimulinkOptions:
    """
    Class containing the options that configure the acados Simulink block,
    i.e. which inputs/outputs are exposed as ports, as well as block-level
    settings such as the sampling time.
    """
    problem_class: InitVar[str] = 'OCP'
    outputs: AcadosOcpSimulinkOutputs = field(default_factory=AcadosOcpSimulinkOutputs)
    inputs: AcadosOcpSimulinkInputs = field(default_factory=AcadosOcpSimulinkInputs)
    samplingtime: str = 't0'
    show_port_info: int = 1
    zoro_iterations: int = 1
    generate_simulink_block: int = 1
    customizable_inputs: Dict[str, Any] = field(default_factory=dict)

    def __post_init__(self, problem_class: str):
        if problem_class not in ('OCP', 'MOCP'):
            raise ValueError(f"problem_class must be 'OCP' or 'MOCP', got '{problem_class}'")

        if problem_class == 'MOCP':
            for input_name in NONSUPPORTED_MOCP_INPUTS + DEFAULT_OFF_MOCP_INPUTS:
                setattr(self.inputs, input_name, 0)

    def _add_customizable_input(self, input_name: str, input_spec: Dict[str, Any]) -> None:
        """Register an additional customizable (e.g. sparse parameter) input port."""
        self.customizable_inputs[input_name] = input_spec

    def to_dict(self) -> Dict[str, Any]:
        d: Dict[str, Any] = {}
        for f in fields(self):
            value = getattr(self, f.name)
            if is_dataclass(value):
                d[f.name] = value.to_dict()
            else:
                d[f.name] = value

        if not d.get('customizable_inputs'):
            d.pop('customizable_inputs', None)

        return d

    def make_consistent(self, solver_options, problem_class: str) -> None:
        if self.inputs.rti_phase and solver_options.nlp_solver_type != 'SQP_RTI':
            raise Exception('rti_phase is only supported for SQP_RTI')

        if self.outputs.KKT_residuals and solver_options.nlp_solver_type == 'SQP_RTI':
            warnings.warn(
                "KKT_residuals now computes the residuals of the output iterate in SQP_RTI, "
                "this leads to increased computation time, turn off this port if it is not needed. "
                "See https://github.com/acados/acados/pull/1346."
            )

        if problem_class == 'MOCP':
            for input_name in NONSUPPORTED_MOCP_INPUTS:
                if getattr(self.inputs, input_name):
                    warnings.warn(
                        f"Simulink input {input_name} is not supported for MOCP, turning it off."
                    )
                    setattr(self.inputs, input_name, 0)

        # validate that all inputs/outputs are 0 or 1
        for group_name, group in (('inputs', self.inputs), ('outputs', self.outputs)):
            for f in fields(group):
                value = getattr(group, f.name)
                if value not in (0, 1):
                    raise ValueError(
                        f"AcadosOcpSimulinkOptions.{group_name}.{f.name} must be 0 or 1, got {value}"
                    )

    def add_sparse_param_port(self, idx_p, port_name: str, stage_idx_0: int, stage_idx_e: int) -> None:
        """
        Allows one to specify information for an input port of the simulink block corresponding
        to an acados OCP solver.

        :param idx_p: 0-based list/array of parameter indices to be updated by the port.
        :param port_name: used to identify the port and print information.
        :param stage_idx_0: first stage for which the parameters should be updated by the port.
        :param stage_idx_e: last stage for which the parameters should be updated by the port.
        """
        if stage_idx_0 > stage_idx_e:
            raise ValueError("stage_idx_0 > stage_idx_e")

        input_name = f'sparse_parameter_{port_name}'

        input_spec = {
            'parameter_indices': list(idx_p),
            'stage_idx_0': stage_idx_0,
            'stage_idx_e': stage_idx_e,
        }

        self._add_customizable_input(input_name, input_spec)

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "AcadosOcpSimulinkOptions":
        obj = cls()
        for key, value in data.items():
            if key == 'outputs':
                obj.outputs = AcadosOcpSimulinkOutputs.from_dict(value)
            elif key == 'inputs':
                obj.inputs = AcadosOcpSimulinkInputs.from_dict(value)
            elif hasattr(obj, key):
                setattr(obj, key, value)
            else:
                warnings.warn(
                    f"Could not assign field {key} in AcadosOcpSimulinkOptions.from_dict()",
                    UserWarning,
                    stacklevel=2,
                )

        if obj.customizable_inputs is None:
            obj.customizable_inputs = {}

        return obj


def get_simulink_default_opts() -> AcadosOcpSimulinkOptions:
    warnings.warn(
        "The function get_simulink_default_opts has been changed in acados v0.5.6. It returns an AcadosOcpSimulinkOptions object instead of a dict. Creating acados Simulink options should now be done using AcadosOcpSimulinkOptions().",
        UserWarning,
        stacklevel=2,
    )
    return AcadosOcpSimulinkOptions()
